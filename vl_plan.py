"""Voronoi Lattice - porosity <-> cell size, and what a run will cost. STDLIB ONLY.

Imported by vl_core (inside Fusion, no numpy) for the live numbers in the dialog,
and by vl_engine (py -3.10) for its first guess, so both start from the same model.

Units are millimetres throughout.

The model: a lattice is a set of struts of diameter d = 2r with total centreline
length L, meeting at J junctions. Its solid volume is

    V = pi r^2 L - BETA r^3 J

the second term being the volume the struts share where they meet. For a
Voronoi lattice of cell size s both L and J follow from the cell count:

    volume lattice    L = C_VOL  V_env / s^2        J = J_VOL  V_env / s^3
    surface lattice   L = C_SURF A / s              J = J_SURF A / s^2

with s the cube (or square) root of the volume (area) per cell. The constants
were measured on Lloyd-relaxed seeds (tests/calibrate.py) and are only a first
guess: the engine measures the real L and J of the graph it built and corrects s
until the predicted porosity matches, then reports the porosity measured off the
final mesh.

Porosity means, for a VOLUME lattice, the void fraction of the lattice envelope;
for a SURFACE lattice, the open fraction of the surface - the part a light shone
along the normal would get through.
"""

import math

# measured by tests/calibrate.py (Lloyd regularity 5, struts 0.3-0.6 mm, cells
# 2-4.5 mm): 30 mm box, 60 mm sheet, 25 mm sphere. C_VOL falls from 5.25 to 4.4
# as cells grow against a 30 mm box - clipped boundary cells - so it is a first
# guess only. BETA is the node overlap: ~5.1 where four struts meet in 3D, ~2.0
# where three meet in a surface, flat across strut and cell size.
C_VOL = 5.0
J_VOL = 5.3
C_SURF = 1.87
J_SURF = 1.9
BETA = {'volume': 5.1, 'surface': 2.0}

PLACEMENTS = ('flush', 'middle', 'external')
RESOLUTION = {'coarse': 3.0, 'medium': 4.0, 'fine': 6.0}  # voxels per diameter
MAX_VOXELS_SLAB = 24_000_000


def mesh_organic(level, d):
    """ORGANIC 1..10 for the mesh output: the smooth-union radius of the strut
    distance field (0 = struts simply intersect) and Taubin smoothing steps."""
    o = min(max(float(level), 1.0), 10.0)
    return {'blend': 0.06 * (o - 1.0) * d, 'smooth': int(round(10 + 3 * (o - 1.0)))}


class PlanError(ValueError):
    pass


def shift_for(placement, d):
    """Offset of the strut CENTRELINES from the selected surface, along its
    outward normal: flush keeps the struts inside, tangent to the surface."""
    r = 0.5 * d
    return {'flush': -r, 'middle': 0.0, 'external': r}[placement]


def envelope(kind, placement, d, volume, area):
    """What porosity is measured against: the volume the lattice's outer
    envelope encloses (first order in the offset), or the surface area."""
    if kind == 'volume':
        return volume + area * (shift_for(placement, d) + 0.5 * d)
    return area


def solid_fraction(kind, s, d):
    """Strut volume over envelope volume (volume lattice), or covered area over
    surface area (surface lattice), for cell size s."""
    r = 0.5 * d
    if kind == 'volume':
        return math.pi * r * r * C_VOL / s ** 2 - BETA['volume'] * r ** 3 * J_VOL / s ** 3
    vol_per_area = math.pi * r * r * C_SURF / s - BETA['surface'] * r ** 3 * J_SURF / s ** 2
    return 4.0 * vol_per_area / (math.pi * d)


def porosity_of(kind, s, d):
    return 1.0 - solid_fraction(kind, s, d)


def cell_for(kind, porosity, d):
    """Cell size giving this porosity. The overlap term makes the model's solid
    fraction peak and then fall as cells shrink below a few strut widths, where
    the model stops meaning anything; above that peak it falls monotonically
    with s, so bisect from the peak up."""
    if not 0.0 < porosity < 1.0:
        raise PlanError('porosity must be between 0 and 100 %')
    grid = [d * 1.02 ** k for k in range(400)]
    lo = max(grid, key=lambda s: solid_fraction(kind, s, d))
    hi = 1000.0 * d
    if porosity_of(kind, lo, d) > porosity:
        raise PlanError('%.0f %% porosity is too low for %.3g mm struts - the '
                        'cells would close up' % (porosity * 100.0, d))
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if porosity_of(kind, mid, d) < porosity:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def estimate(kind, s, d, volume, area, placement, resolution='medium'):
    """Rough cost of a run, for the dialog: cells, strut length, triangles."""
    env = envelope(kind, placement, d, volume, area)
    if kind == 'volume':
        cells = env / s ** 3
        length = C_VOL * env / s ** 2
    else:
        cells = area / s ** 2
        length = C_SURF * area / s
    h = d / RESOLUTION[resolution]
    skin = 2.0 * math.pi * 0.5 * d * length
    tris = 2.2 * skin / (h * h)
    return {'cells': cells, 'length_mm': length, 'voxel_mm': h,
            'raw_triangles': tris, 'pore_mm': max(s - d, 0.0),
            'envelope': env}
