"""Lattice - porosity <-> cell size, and what a run will cost. STDLIB ONLY.

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

CUBIC = {'volume': (3.0, 1.0, 8.0 * math.sqrt(2.0)), 'surface': (2.0, 1.0, 16.0 / 3.0)}
CELL_TYPES = ('voronoi', 'cubic', 'bcc', 'fcc', 'octet', 'diamond', 'kelvin')
# every cell but Voronoi tiles a cube of side s: one unit cell, built once
REGULAR = CELL_TYPES[1:]
# on a surface only Voronoi and the cubic grid's trace mean anything; the other
# cells are three-dimensional and fill a volume
SURFACE_TYPES = ('voronoi', 'cubic')

# strut centreline length per cubic cell, in units of s
LENGTH_PER_CELL = {'cubic': 3.0, 'bcc': 4.0 * math.sqrt(3.0), 'fcc': 6.0 * math.sqrt(2.0),
                   'octet': 12.0 * math.sqrt(2.0), 'diamond': 4.0 * math.sqrt(3.0),
                   'kelvin': 6.0 * math.sqrt(2.0)}

# Solid fraction of the periodic unit cell (union of round-ended struts of
# diameter d) at d/s = FRACTION_T0 + FRACTION_DT * i. No closed form survives the
# node overlaps of these cells, so it was measured once: 4 M jittered points per
# cell against the struts sampled every 0.002 s (the same measurement gives the
# cubic cell's exact value to 5e-5). BCC and diamond agree to 1e-4 everywhere:
# both are <111> rods at the same density. Past the last entry the cell is solid.
FRACTION_T0 = 0.02
FRACTION_DT = 0.02
FRACTION = {
    'bcc': [0.0021, 0.0084, 0.0186, 0.0323, 0.0495, 0.0699, 0.0932, 0.1192, 0.1477, 0.1784,
            0.2112, 0.2457, 0.2817, 0.3190, 0.3575, 0.3967, 0.4364, 0.4766, 0.5170, 0.5570,
            0.5969, 0.6361, 0.6746, 0.7119, 0.7480, 0.7825, 0.8153, 0.8461, 0.8746, 0.9007,
            0.9241, 0.9446, 0.9619, 0.9757, 0.9859, 0.9925, 0.9964, 0.9986, 0.9996, 1.0000],
    'fcc': [0.0026, 0.0102, 0.0225, 0.0392, 0.0598, 0.0842, 0.1119, 0.1427, 0.1760, 0.2119,
            0.2498, 0.2895, 0.3306, 0.3726, 0.4154, 0.4587, 0.5020, 0.5451, 0.5876, 0.6293,
            0.6697, 0.7085, 0.7455, 0.7803, 0.8126, 0.8419, 0.8681, 0.8907, 0.9095, 0.9249,
            0.9376, 0.9483, 0.9573, 0.9649, 0.9713, 0.9768, 0.9815, 0.9856, 0.9889, 0.9917,
            0.9940, 0.9958, 0.9972, 0.9982, 0.9990, 0.9995, 0.9998, 0.9999],
    'octet': [0.0052, 0.0201, 0.0438, 0.0754, 0.1140, 0.1586, 0.2083, 0.2621, 0.3192, 0.3786,
            0.4394, 0.5007, 0.5616, 0.6210, 0.6781, 0.7320, 0.7817, 0.8263, 0.8648, 0.8965,
            0.9206, 0.9389, 0.9532, 0.9645, 0.9734, 0.9805, 0.9862, 0.9907, 0.9940, 0.9964,
            0.9981, 0.9991, 0.9997, 0.9999],
    'diamond': [0.0021, 0.0084, 0.0185, 0.0323, 0.0495, 0.0699, 0.0932, 0.1192, 0.1477, 0.1784,
            0.2112, 0.2457, 0.2818, 0.3191, 0.3574, 0.3966, 0.4365, 0.4767, 0.5169, 0.5571,
            0.5969, 0.6361, 0.6746, 0.7119, 0.7480, 0.7825, 0.8153, 0.8461, 0.8747, 0.9007,
            0.9241, 0.9445, 0.9619, 0.9757, 0.9859, 0.9925, 0.9964, 0.9986, 0.9996, 1.0000],
    'kelvin': [0.0026, 0.0102, 0.0224, 0.0388, 0.0591, 0.0829, 0.1099, 0.1397, 0.1719, 0.2062,
            0.2423, 0.2796, 0.3179, 0.3569, 0.3961, 0.4353, 0.4739, 0.5117, 0.5488, 0.5849,
            0.6200, 0.6542, 0.6870, 0.7186, 0.7487, 0.7771, 0.8037, 0.8283, 0.8507, 0.8708,
            0.8886, 0.9041, 0.9180, 0.9302, 0.9410, 0.9506, 0.9591, 0.9664, 0.9728, 0.9783,
            0.9830, 0.9870, 0.9902, 0.9929, 0.9950, 0.9967, 0.9979, 0.9988, 0.9994, 0.9997,
            0.9999],
}

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


def solid_fraction(kind, s, d, cell_type='voronoi'):
    """Strut volume over envelope volume (volume lattice), or covered area over
    surface area (surface lattice), for cell size s."""
    r = 0.5 * d
    if cell_type in FRACTION:
        if kind != 'volume':
            raise PlanError('a %s cell fills a volume; on a surface use Voronoi or '
                            'Cubic' % cell_type)
        return _table_fraction(cell_type, d / s)
    if cell_type == 'cubic':
        c, j, beta = CUBIC[kind]
        if kind == 'volume':
            return (math.pi * r * r * c * s - beta * r ** 3 * j) / s ** 3
        return 4.0 * (math.pi * r * r * c * s - beta * r ** 3 * j) / (math.pi * d * s * s)
    if kind == 'volume':
        return math.pi * r * r * C_VOL / s ** 2 - BETA['volume'] * r ** 3 * J_VOL / s ** 3
    vol_per_area = math.pi * r * r * C_SURF / s - BETA['surface'] * r ** 3 * J_SURF / s ** 2
    return 4.0 * vol_per_area / (math.pi * d)


def _table_fraction(cell_type, t):
    tab = FRACTION[cell_type]
    x = (t - FRACTION_T0) / FRACTION_DT
    if x <= 0.0:
        # below the first entry the struts barely meet: pi r^2 L / s^3
        return math.pi * 0.25 * t * t * LENGTH_PER_CELL[cell_type]
    i = int(x)
    if i >= len(tab) - 1:
        return tab[-1]
    return tab[i] + (x - i) * (tab[i + 1] - tab[i])


def porosity_of(kind, s, d, cell_type='voronoi'):
    return 1.0 - solid_fraction(kind, s, d, cell_type)


def cell_for(kind, porosity, d, cell_type='voronoi'):
    """Cell size giving this porosity. The overlap term makes the model's solid
    fraction peak and then fall as cells shrink below a few strut widths, where
    the model stops meaning anything; above that peak it falls monotonically
    with s, so bisect from the peak up."""
    if not 0.0 < porosity < 1.0:
        raise PlanError('porosity must be between 0 and 100 %')
    if cell_type in FRACTION:
        if kind != 'volume':
            raise PlanError('a %s cell fills a volume; on a surface use Voronoi or '
                            'Cubic' % cell_type)
        # the measured fraction rises monotonically with d/s: bisect on it
        t_hi = FRACTION_T0 + FRACTION_DT * (len(FRACTION[cell_type]) - 1)
        if 1.0 - _table_fraction(cell_type, t_hi) > porosity:
            raise PlanError('%.0f %% porosity is too low for a %s cell - it would '
                            'close up' % (porosity * 100.0, cell_type))
        lo, hi = 1e-4, t_hi
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if 1.0 - _table_fraction(cell_type, mid) > porosity:
                lo = mid
            else:
                hi = mid
        return d / (0.5 * (lo + hi))
    grid = [d * 1.02 ** k for k in range(400)]
    if cell_type == 'cubic':
        grid = [g for g in grid if g > d]      # the model holds while d < s
    lo = max(grid, key=lambda s: solid_fraction(kind, s, d, cell_type))
    hi = 1000.0 * d
    if porosity_of(kind, lo, d, cell_type) > porosity:
        raise PlanError('%.0f %% porosity is too low for %.3g mm struts - the '
                        'cells would close up' % (porosity * 100.0, d))
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if porosity_of(kind, mid, d, cell_type) < porosity:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def estimate(kind, s, d, volume, area, placement, resolution='medium',
             cell_type='voronoi'):
    """Rough cost of a run, for the dialog: cells, strut length, triangles."""
    env = envelope(kind, placement, d, volume, area)
    if kind == 'volume':
        cells = env / s ** 3
        length = LENGTH_PER_CELL.get(cell_type, C_VOL) * env / s ** 2
    else:
        cells = area / s ** 2
        length = (CUBIC['surface'][0] if cell_type in REGULAR else C_SURF) * area / s
    h = d / RESOLUTION[resolution]
    skin = 2.0 * math.pi * 0.5 * d * length
    tris = 2.2 * skin / (h * h)
    return {'cells': cells, 'length_mm': length, 'voxel_mm': h,
            'raw_triangles': tris, 'pore_mm': max(s - d, 0.0),
            'envelope': env}
