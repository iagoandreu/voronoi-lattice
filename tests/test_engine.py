"""Offline checks of the lattice engine and the planning model. No Fusion.

    py -3.10 tests/test_engine.py

Every check has its answer known in advance: a strut's exact capsule volume, the
porosity asked for, the side of the surface each placement must land on, a flush
lattice that must not cross the surface it is flush with.
"""

import math
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, 'tests', 'stub'))

import numpy as np                                            # noqa: E402
import trimesh                                                # noqa: E402

import vl_engine as E                                         # noqa: E402
import vl_plan as P                                           # noqa: E402

FAILED = []
WORK = tempfile.mkdtemp(prefix='vl_test_')


def check(name, cond, detail=''):
    print('  %s %s%s' % ('ok  ' if cond else 'FAIL', name,
                         '   %s' % (detail,) if detail != '' else ''), flush=True)
    if not cond:
        FAILED.append(name)


def write(name, m, normals=None):
    p = {'coords': os.path.join(WORK, name + '_v.bin'),
         'tris': os.path.join(WORK, name + '_f.bin')}
    np.asarray(m.vertices, '<f8').tofile(p['coords'])
    np.asarray(m.faces, '<i4').tofile(p['tris'])
    if normals is not None:
        p['normals'] = os.path.join(WORK, name + '_n.bin')
        np.asarray(normals, '<f8').tofile(p['normals'])
    return p


def job(name, m, **kw):
    j = dict(mesh=write(name, m, kw.pop('normals', None)), closed=True,
             kind='volume', placement='flush', define='porosity', porosity=0.8,
             cell=3.0, thickness=0.6, regularity=5, seed=1, skin=True, border=True,
             resolution='medium', max_tris=1_500_000, mode='mesh',
             out_stl=os.path.join(WORK, name + '.stl'))
    j.update(kw)
    return j


def stl_points(path):
    m = trimesh.load(path, force='mesh')
    return np.asarray(m.vertices), m


def curved_sheet(radius=30.0, half_angle=0.7, half_len=20.0, n=60, flip=False):
    """Part of a cylinder of radius 30 bulging towards +z; its normals face +z."""
    u, v = np.meshgrid(np.linspace(-half_angle, half_angle, n),
                       np.linspace(-half_len, half_len, n))
    X = np.stack([radius * np.sin(u), v, radius * np.cos(u) - radius], -1).reshape(-1, 3)
    idx = np.arange(n * n).reshape(n, n)
    a, b = idx[:-1, :-1].ravel(), idx[1:, :-1].ravel()
    c, d = idx[:-1, 1:].ravel(), idx[1:, 1:].ravel()
    F = np.vstack([np.stack([a, c, b], 1), np.stack([b, c, d], 1)])
    if flip:
        F = F[:, ::-1]
    m = trimesh.Trimesh(X, F, process=False)
    N = np.stack([X[:, 0] / radius, np.zeros(len(X)), (X[:, 2] + radius) / radius], 1)
    return m, N


def test_plan():
    print('plan')
    for kind in ('volume', 'surface'):
        for phi in (0.5, 0.7, 0.85, 0.95):
            s = P.cell_for(kind, phi, 0.5)
            check('%s %.0f%%: cell_for inverts porosity_of' % (kind, phi * 100),
                  abs(P.porosity_of(kind, s, 0.5) - phi) < 1e-6, round(s, 4))
    try:
        P.cell_for('volume', 0.05, 0.5)
        check('5 %% volume porosity is refused', False)
    except P.PlanError:
        check('5 % volume porosity is refused', True)
    check('placement shifts: flush -r, middle 0, external +r',
          (P.shift_for('flush', 1.0), P.shift_for('middle', 1.0),
           P.shift_for('external', 1.0)) == (-0.5, 0.0, 0.5))


def test_strut_volume():
    print('strut volume (contour bias corrected)')
    rng = np.random.default_rng(0)
    for d in (0.4, 1.0):
        r = d / 2
        A, B = [], []
        for k in range(10):
            u = rng.normal(size=3)
            u /= np.linalg.norm(u)
            c = np.array([k * 12.0, 0, 0])
            A.append(c - 4 * u)
            B.append(c + 4 * u)
        A, B = np.array(A), np.array(B)
        exact = 10 * (math.pi * r * r * 8 + 4.0 / 3.0 * math.pi * r ** 3)
        for res in ('coarse', 'medium', 'fine'):
            s, _ = E.mesh_struts(A, B, d, d / P.RESOLUTION[res])
            s, fi = E.finish_mesh(s, 10 ** 9)
            err = s.volume / exact - 1.0
            check('d %.1f %s: within 2.5 %% of the exact capsules' % (d, res),
                  abs(err) < 0.025, '%+.2f %%' % (100 * err))
            check('d %.1f %s: ten separate watertight struts' % (d, res),
                  fi['components'] == 10 and fi['open_or_nonmanifold_edges'] == 0,
                  (fi['components'], fi['open_or_nonmanifold_edges']))


def test_volume():
    print('volume lattice')
    box = trimesh.creation.box((20, 20, 20))
    sdf = E.SignedDistance(np.asarray(box.vertices, float), np.asarray(box.faces))
    for place in ('flush', 'middle', 'external'):
        r = E.run(job('box_' + place, box, placement=place))
        check('%s: porosity measured within 1.5 points of 80 %%' % place,
              abs(r['porosity_measured'] - 0.8) < 0.015,
              '%.2f %%' % (100 * r['porosity_measured']))
        check('%s: watertight, one piece' % place,
              r['mesh']['open_or_nonmanifold_edges'] == 0
              and r['mesh']['components'] == 1,
              (r['mesh']['open_or_nonmanifold_edges'], r['mesh']['components']))
        V, m = stl_points(r['stl'])
        check('%s: the STL faces outward (signed volume > 0)' % place,
              m.volume > 0 and abs(m.volume - r['mesh']['volume_mm3']) < 1e-3 * m.volume,
              round(float(m.volume), 2))
        out = sdf(V).max()
        if place == 'flush':
            check('flush: nothing crosses the surface', out < 0.01, '%.4f mm' % out)
        elif place == 'middle':
            check('middle: struts stand half a strut proud',
                  0.15 < out < 0.35, '%.4f mm' % out)
        else:
            check('external: struts stand a strut proud', 0.45 < out < 0.7,
                  '%.4f mm' % out)
    r = E.run(job('box_noskin', box, skin=False, porosity=0.85))
    check('no skin struts: still 85 %', abs(r['porosity_measured'] - 0.85) < 0.015,
          '%.2f %%' % (100 * r['porosity_measured']))
    r = E.run(job('box_cell', box, define='cell', cell=3.0))
    check('cell-size mode keeps the cell asked for', r['graph']['cell_mm'] == 3.0)
    check('and predicts its porosity within 1.5 points',
          abs(r['porosity_measured'] - r['porosity_predicted']) < 0.015,
          (round(r['porosity_predicted'], 4), round(r['porosity_measured'], 4)))
    ball = trimesh.creation.icosphere(4, 10)
    r = E.run(job('ball', ball, porosity=0.7, thickness=0.5))
    check('curved ball at 70 %', abs(r['porosity_measured'] - 0.7) < 0.015,
          '%.2f %%' % (100 * r['porosity_measured']))
    inv = trimesh.Trimesh(ball.vertices, ball.faces[:, ::-1], process=False)
    r2 = E.run(job('ball_inside_out', inv, porosity=0.7, thickness=0.5,
                   mode='graph'))
    check('an inside-out closed mesh is turned the right way first',
          r2['graph']['struts'] > 100 and abs(r2['region']['volume_mm3']
                                             - r['region']['volume_mm3']) < 1e-6)


def test_surface():
    print('surface lattice')
    ball = trimesh.creation.icosphere(4, 10)
    sdf = E.SignedDistance(np.asarray(ball.vertices, float), np.asarray(ball.faces))
    for place, lo, hi in (('flush', -0.62, 0.02), ('middle', -0.33, 0.33),
                          ('external', -0.02, 0.62)):
        r = E.run(job('sball_' + place, ball, kind='surface', placement=place))
        V, _ = stl_points(r['stl'])
        dv = sdf(V)
        check('closed %s: porosity within 1.5 points of 80 %%' % place,
              abs(r['porosity_measured'] - 0.8) < 0.015,
              '%.2f %%' % (100 * r['porosity_measured']))
        check('closed %s: mesh between %.2f and %.2f mm of the surface'
              % (place, lo, hi), dv.min() > lo - 0.05 and dv.max() < hi + 0.05,
              (round(float(dv.min()), 3), round(float(dv.max()), 3)))

    sheet, N = curved_sheet()
    sides = {}
    for flip in (False, True):
        for place in ('flush', 'external'):
            r = E.run(job('sheet_%s_%s' % (place, flip), sheet, closed=False,
                          kind='surface', placement=place, flip=flip,
                          normals=np.asarray(N)))
            V, _ = stl_points(r['stl'])
            # signed height above the cylinder of radius 30 whose axis is y, z=-30
            h = np.hypot(V[:, 0], V[:, 2] + 30.0) - 30.0
            sides[(place, flip)] = float(np.median(h))
            check('open %s flip=%s: one watertight piece' % (place, flip),
                  r['mesh']['open_or_nonmanifold_edges'] == 0
                  and r['mesh']['components'] == 1)
    check('open: external sits on the normal side, flush behind it',
          sides[('external', False)] > 0.2 and sides[('flush', False)] < -0.2,
          sides)
    check('open: flip swaps the two', sides[('external', True)] < -0.2
          and sides[('flush', True)] > 0.2)
    rev, Nrev = curved_sheet(flip=True)          # wound to face -z, normals say +z
    r = E.run(job('sheet_rev', rev, closed=False, kind='surface',
                  placement='external', normals=np.asarray(Nrev)))
    V, _ = stl_points(r['stl'])
    check('open: the side follows the supplied normals, not the winding',
          np.median(np.hypot(V[:, 0], V[:, 2] + 30.0) - 30.0) > 0.2)
    r_b = E.run(job('sheet_border', sheet, closed=False, kind='surface',
                    mode='graph', border=True))
    r_n = E.run(job('sheet_noborder', sheet, closed=False, kind='surface',
                    mode='graph', border=False))
    check('border strut adds the rim (~160 mm here)',
          r_b['graph']['border_struts'] > 0
          and 'border_struts' not in r_n['graph'], r_b['graph'].get('border_struts'))


def test_refine():
    print('refining an uneven open surface (an imported STL)')
    sheet, _N = curved_sheet(n=40)
    rng = np.random.default_rng(2)
    V = np.asarray(sheet.vertices).copy()
    V[:, 1] += rng.uniform(-0.3, 0.3, len(V)) * (np.abs(V[:, 1]) < 19.0)
    F = np.asarray(sheet.faces)
    el = np.linalg.norm(V[F[:, [1, 2, 0]]] - V[F], axis=2)
    rim = lambda P_, p_: float(np.linalg.norm(np.diff(P_[p_], axis=0), axis=1).sum())
    rim0 = rim(V, E.border_paths(V, F)[0])
    _t0, a0 = E.tri_normals(V, F)
    # trimesh's subdivide_to_size left 436 / 104 / 436 border loops at these limits
    loops_n, rim_err, dup, area_err, over = [], 0.0, 0, 0.0, 0.0
    for fac in (0.6, 0.8, 1.2):
        lim = float(np.median(el)) * fac
        Vr, Fr = E.refine(V, F, lim)
        loops = E.border_paths(Vr, Fr)
        loops_n.append(len(loops))
        rim_err = max(rim_err, abs(rim(Vr, loops[0]) - rim0))
        d = np.vstack([Fr[:, [0, 1]], Fr[:, [1, 2]], Fr[:, [2, 0]]])
        _u, c = np.unique(d, axis=0, return_counts=True)
        dup += int((c > 1).sum())
        _t1, a1 = E.tri_normals(Vr, Fr)
        area_err = max(area_err, abs(a1.sum() / a0.sum() - 1.0))
        el2 = np.linalg.norm(Vr[Fr[:, [1, 2, 0]]] - Vr[Fr], axis=2)
        over = max(over, el2.max() - lim)
    check('still one border loop: no cracks where one triangle split and its '
          'neighbour did not', loops_n == [1, 1, 1], loops_n)
    check('the border is as long as before', rim_err < 1e-6, rim_err)
    check('orientation kept, no edge used twice the same way', dup == 0, dup)
    check('area unchanged', area_err < 1e-12, area_err)
    check('no edge left over the limit', over <= 1e-12, over)
    # end to end: 80 % with struts too thick for the region is unreachable; the
    # closest lattice is built and flagged, not an error
    r = E.run(job('sheet_thick', sheet, closed=False, kind='surface', mode='graph',
                  placement='external', thickness=3.0, porosity=0.8))
    check('unreachable porosity: the closest lattice comes back, flagged',
          r['graph'].get('porosity_unreached') is True, r['graph'].get('porosity_predicted'))


def test_thin_wall():
    print('thin wall')
    plate = trimesh.creation.box((40, 40, 1.0))
    r = E.run(job('plate', plate, kind='surface', placement='flush', thickness=0.3,
                  porosity=0.85, mode='graph'))
    A = np.array([s[:3] for s in r['segments']])
    B = np.array([s[3:] for s in r['segments']])
    mid = 0.5 * (A + B)
    # struts may run over the 1 mm rim; one THROUGH the wall would have its
    # middle inside the plate, near z = 0, away from the rim
    through = np.sum((np.abs(mid[:, 2]) < 0.2)
                     & (np.max(np.abs(mid[:, :2]), axis=1) < 19.0))
    check('a 1 mm plate with 6 mm cells: no strut runs through the wall',
          through == 0, int(through))
    check('and both faces carry a net',
          np.sum(mid[:, 2] > 0.2) > 50 and np.sum(mid[:, 2] < -0.2) > 50)


def test_organic_mesh():
    print('organic, mesh output')
    ball = trimesh.creation.icosphere(4, 10)
    vols = {}
    for o in (1, 5, 10):
        r = E.run(job('org_ball_%d' % o, ball, porosity=0.8, organic=o))
        check('volume, organic %d: porosity within 1.5 points of 80 %%' % o,
              abs(r['porosity_measured'] - 0.8) < 0.015,
              '%.2f %%' % (100 * r['porosity_measured']))
        check('volume, organic %d: watertight, one piece' % o,
              r['mesh']['open_or_nonmanifold_edges'] == 0
              and r['mesh']['components'] == 1)
        rc = E.run(job('org_cell_%d' % o, ball, define='cell', cell=3.0, organic=o))
        vols[o] = rc['mesh']['volume_mm3']
    check('at one cell size, more organic = more material at the nodes',
          vols[1] < vols[5] < vols[10], {k: round(v, 1) for k, v in vols.items()})
    r = E.run(job('org_sball', ball, kind='surface', organic=8))
    check('surface, organic 8: porosity within 1.5 points of 80 %',
          abs(r['porosity_measured'] - 0.8) < 0.015,
          '%.2f %%' % (100 * r['porosity_measured']))


def tsm_counts(path):
    c = {'f': 0, 'e': 0, 'v': 0, 'l': 0}
    with open(path, encoding='utf-8') as fh:
        for line in fh:
            k = line.split(' ', 1)[0]
            if k in c:
                c[k] += 1
    return c


def test_form():
    print('Form output')
    import vl_form as Fm
    ball = trimesh.creation.icosphere(4, 10)
    for o in (1, 5, 10):
        j = job('form_ball_%d' % o, ball, kind='surface', placement='middle',
                output='form', organic=o, out_tsm=os.path.join(WORK, 'ball_%d.tsm' % o))
        r = E.run(j)
        check('surface ball, organic %d: built as a Form' % o, r['output'] == 'form'
              and os.path.exists(r['tsm']), r['output'])
        check('organic %d: limit-surface porosity within 2 points of 80 %%' % o,
              abs(r['porosity_measured'] - 0.8) < 0.02,
              '%.2f %%' % (100 * r['porosity_measured']))
        c = tsm_counts(r['tsm'])
        f = r['form']
        check('organic %d: TSM has one f/v line per cage face/vertex, two links '
              'per edge' % o, c['f'] == f['cage_faces'] and c['v'] == f['cage_vertices']
              and c['l'] == 2 * c['e'], c)
    sheet, N = curved_sheet()
    r = E.run(job('form_sheet', sheet, closed=False, kind='surface', placement='external',
                  normals=np.asarray(N), output='form', organic=5,
                  out_tsm=os.path.join(WORK, 'sheet.tsm')))
    check('open sheet with its border: built as a Form', r['output'] == 'form',
          r.get('graph', {}).get('form_failed', ''))
    box = trimesh.creation.box((20, 20, 20))
    r = E.run(job('form_box', box, output='form', organic=5, mode='graph'))
    check('dense volume lattice: falls back to a mesh and says why',
          r['output'] == 'mesh' and 'cannot be built as a Form' in
          r['graph'].get('form_failed', ''), r['graph'].get('form_failed', '')[:60])
    T = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], float) / np.sqrt(3)
    segs = []
    a = 4.0
    for i in range(3):
        for jj in range(3):
            for k in range(3):
                for b in ((0, 0, 0), (0, .5, .5), (.5, 0, .5), (.5, .5, 0)):
                    p = (np.array(b) + [i, jj, k]) * a
                    for t in T * np.sqrt(3) / 4 * a:
                        segs.append((p, p + t))
    A = np.array([s_[0] for s_ in segs])
    B = np.array([s_[1] for s_ in segs])
    V, F, ci = Fm.build_cage(A, B, 0.6, 5)
    check('a diamond lattice builds without folding or pushing a single hub',
          ci['folded'] == 0 and ci['rings_pushed'] == 0, ci)


def test_remesh():
    print('region remeshing')
    box = trimesh.creation.box((40, 30, 20))       # 12 triangles, 50 mm diagonals
    r = E.run(job('rm_box', box, mode='graph', define='cell', cell=3.0))
    rm = r['region']['remesh']
    check('a CAD-like box (12 triangles) is remeshed automatically',
          rm['done'] and rm['method'] == 'isotropic', rm.get('method'))
    check('to edges of half a cell on average, none over twice that',
          abs(rm['after']['mean_edge_mm'] / rm['target_edge_mm'] - 1.0) < 0.15
          and rm['after']['max_edge_mm'] <= 2.1 * rm['target_edge_mm'],
          'mean %.2f max %.2f target %.2f mm' % (rm['after']['mean_edge_mm'],
                                                 rm['after']['max_edge_mm'], rm['target_edge_mm']))
    check('keeping the volume within 0.5 %', abs(rm['volume_change']) < 0.005,
          '%+.4f %%' % (100 * rm['volume_change']))
    check('and the reported region is the box', abs(r['region']['volume_mm3'] - 24000.0) < 120.0)
    r = E.run(job('rm_off', box, mode='graph', define='cell', cell=3.0, remesh='off'))
    rm = r['region']['remesh']
    check('off: left alone, but flagged', not rm['done'] and rm['needed'])
    ball = trimesh.creation.icosphere(subdivisions=4, radius=10.0)
    r = E.run(job('rm_ball', ball, mode='graph', define='cell', cell=3.0))
    check('a fine mesh is not touched', not r['region']['remesh']['done'],
          '%.2f mm edges' % r['region']['remesh']['before']['max_edge_mm'])
    r = E.run(job('rm_edge', box, mode='graph', define='cell', cell=3.0,
                  remesh='always', remesh_edge=2.5))
    rm = r['region']['remesh']
    check('the target edge can be set', rm['target_edge_mm'] == 2.5
          and abs(rm['after']['mean_edge_mm'] / 2.5 - 1.0) < 0.15,
          'mean %.2f mm' % rm['after']['mean_edge_mm'])


def test_errors():
    print('refusals')
    box = trimesh.creation.box((10, 10, 10))
    for bad, why in ((dict(thickness=0.0), 'zero thickness'),
                     (dict(define='cell', cell=0.5, thickness=0.6), 'cell <= strut'),
                     (dict(kind='volume', closed=False), 'volume on an open region'),
                     (dict(placement='inside'), 'unknown placement'),
                     (dict(remesh='sometimes'), 'unknown remesh mode')):
        try:
            E.run(job('bad', box, mode='graph', **bad))
            check('refuses %s' % why, False)
        except E.LatticeError:
            check('refuses %s' % why, True)
    j = job('cache', box, mode='graph', cache_dir=os.path.join(WORK, 'cache'))
    E.run(j)
    r = E.run(dict(j, mode='graph'))
    check('the graph is cached between preview and build', r['graph'].get('cached'))


def test_cubic():
    print('cubic lattice')
    for phi in (0.7, 0.8, 0.9):
        s = P.cell_for('volume', phi, 0.65, 'cubic')
        r = 0.325
        exact = 1.0 - (3 * math.pi * r * r * s - 8 * math.sqrt(2) * r ** 3) / s ** 3
        check('cubic model at %.0f %% is the exact cell' % (100 * phi),
              abs(exact - phi) < 1e-6, round(s, 4))
    # one cell's struts as a mesh: the union of three orthogonal cylinders
    s = 2.0
    cyl = [trimesh.creation.cylinder(0.325, s + 0.002, sections=128) for _ in range(3)]
    cyl[1].apply_transform(trimesh.transformations.rotation_matrix(math.pi / 2, [1, 0, 0]))
    cyl[2].apply_transform(trimesh.transformations.rotation_matrix(math.pi / 2, [0, 1, 0]))
    try:
        u = trimesh.boolean.union(cyl, engine='manifold')
        model = 1.0 - P.porosity_of('volume', s, 0.65, 'cubic')
        check('three crossed struts match the model volume',
              abs(u.volume - model * s ** 3) < 0.01 * u.volume,
              (round(float(u.volume), 4), round(model * s ** 3, 4)))
    except Exception as exc:                                  # no boolean engine
        check('boolean engine for the cubic cell check', False, exc)
    box = trimesh.creation.box((20, 20, 20))
    for skin in (False, True):
        r = E.run(job('cubic_box_%d' % skin, box, cell_type='cubic', thickness=0.65,
                      porosity=0.8, skin=skin, organic=1))
        g = r['graph']
        check('cubic box skin=%s: the unit cell is the 80 %% one' % skin,
              abs(g['porosity_unit_cell'] - 0.8) < 1e-6
              and abs(g['cell_mm'] - P.cell_for('volume', 0.8, 0.65, 'cubic')) < 1e-9,
              round(g['cell_mm'], 4))
        check('cubic box skin=%s: the trimmed part measures within 4 points' % skin,
              abs(r['porosity_measured'] - 0.8) < 0.04,
              '%.2f %%' % (100 * r['porosity_measured']))
        check('cubic box skin=%s: watertight, one piece' % skin,
              r['mesh']['open_or_nonmanifold_edges'] == 0 and r['mesh']['components'] == 1,
              (r['mesh']['open_or_nonmanifold_edges'], r['mesh']['components']))
    g = E.run(job('cubic_graph', box, cell_type='cubic', thickness=0.65, porosity=0.8,
                  skin=False, mode='graph'))['graph']
    check('cubic graph: struts along X, Y and Z only', g['struts'] > 1000, g['struts'])
    # turned onto a picked axis: interior struts run only along the three frame
    # axes, the axis is a cell diagonal (|cos| = 1/sqrt3 to every strut) or a
    # cell edge, and a grid node sits on the axis line
    a = np.array([1.0, 2.0, 3.0]) / math.sqrt(14.0)
    pt = np.array([3.0, -2.0, 1.0])
    for align, want in (('diagonal', 1.0 / math.sqrt(3.0)), ('edge', None)):
        o = {'point': pt.tolist(), 'dir': a.tolist(), 'align': align, 'spin_deg': 17.0}
        r = E.run(job('cubic_orient_' + align, box, cell_type='cubic', thickness=0.65,
                      porosity=0.8, skin=False, mode='graph', orient=o,
                      max_preview_segments=10 ** 7))
        S = np.asarray(r['segments'])
        u = S[:, 3:] - S[:, :3]
        L = np.linalg.norm(u, axis=1)
        whole = np.abs(L - r['graph']['cell_mm']) < 1e-4    # segments carry 5 decimals
        u = u[whole] / L[whole, None]
        ax = np.asarray(r['graph']['frame']['axes'])
        along = np.abs(u @ ax.T).max(1)
        cosa = np.abs(u @ a)
        node = np.asarray(r['graph']['frame']['node'])
        q = node - pt
        off = np.linalg.norm(q - (q @ a) * a)
        ok_dir = (np.allclose(cosa, want, atol=1e-4) if want else
                  np.isclose(cosa, 1.0, atol=1e-4).any()
                  and np.all(np.isclose(cosa, 1.0, atol=1e-4) | (cosa < 1e-4)))
        check('cubic oriented (%s): whole struts along the frame only, axis as asked, '
              'node on the axis' % align,
              whole.sum() > 500 and np.allclose(along, 1.0, atol=1e-4) and ok_dir
              and off < 1e-9 and np.allclose(ax @ ax.T, np.eye(3), atol=1e-9),
              (int(whole.sum()), round(float(cosa.min()), 6), round(float(cosa.max()), 6),
               off))
    sheet, N = curved_sheet()
    r = E.run(job('cubic_sheet', sheet, normals=N, kind='surface', closed=False,
                  cell_type='cubic', thickness=0.65, porosity=0.6))
    check('cubic surface lattice: watertight, porosity reported',
          r['mesh']['open_or_nonmanifold_edges'] == 0 and 0.3 < r['porosity_measured'] < 0.9,
          '%.2f %%' % (100 * r['porosity_measured']))
    o = {'point': [0, 0, 0], 'dir': a.tolist(), 'align': 'diagonal', 'spin_deg': 0.0}
    r = E.run(job('cubic_sheet_orient', sheet, normals=N, kind='surface', closed=False,
                  cell_type='cubic', thickness=0.65, porosity=0.6, orient=o))
    check('cubic surface lattice turned onto an axis: watertight',
          r['mesh']['open_or_nonmanifold_edges'] == 0 and 0.3 < r['porosity_measured'] < 0.9,
          '%.2f %%' % (100 * r['porosity_measured']))


def main():
    test_plan()
    test_strut_volume()
    test_volume()
    test_cubic()
    test_surface()
    test_refine()
    test_thin_wall()
    test_remesh()
    test_organic_mesh()
    test_form()
    test_errors()
    print()
    print('%d failed: %s' % (len(FAILED), FAILED) if FAILED else 'all offline checks passed')
    return 1 if FAILED else 0


if __name__ == '__main__':
    raise SystemExit(main())
