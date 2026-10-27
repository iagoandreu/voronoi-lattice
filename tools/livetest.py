"""Live checks inside Fusion, through the dev bridge. Each job opens a NEW
document and closes it unsaved; the active document is never touched.

    py -3.10 tools/livetest.py

What the dialog decides from the selection (closed or not) and what lands in the
design (a closed mesh body, in the right component, the right size, replaceable)
- the parts the offline tests cannot reach.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import vl as client                                          # noqa: E402

FAILED = []


def check(name, cond, detail=''):
    print('  %s %s%s' % ('ok  ' if cond else 'FAIL', name,
                         '   %s' % (detail,) if detail != '' else ''), flush=True)
    if not cond:
        FAILED.append(name)


def job(cmd, payload):
    r = client.send(cmd, payload)
    if not r.get('ok'):
        print(r.get('error'))
        raise SystemExit(1)
    return r['result']


def inside(a, b, tol):
    """bbox a within bbox b, to tol mm."""
    return all(a[0][k] >= b[0][k] - tol and a[1][k] <= b[1][k] + tol for k in range(3))


def main():
    if not client.alive():
        print('The dev bridge is not answering - is Voronoi Lattice running?')
        return 2
    print('ping', job('ping', {}))

    print('reading the selection')
    r = job('region', {'shape': 'box'})
    check('a solid box is closed, 8000 mm3, 2400 mm2',
          r['closed'] and abs(r['volume'] - 8000) < 1e-6 and abs(r['area'] - 2400) < 1e-6,
          r)
    r = job('region', {'shape': 'face'})
    check('one face of it is open, 400 mm2, and offers "whole body"',
          not r['closed'] and abs(r['area'] - 400) < 1e-6 and r['faces_of_one_solid'], r)
    r = job('region', {'shape': 'face', 'whole': True})
    check('ticking "whole body" makes it the closed solid', r['closed'], r['what'])
    r = job('region', {'shape': 'faces_all'})
    check('all six faces of a solid count as the closed solid', r['closed'], r['what'])
    r = job('region', {'shape': 'sheet'})
    check('a lone surface body is open', not r['closed'], r['what'])

    print('building')
    r = job('build', {'shape': 'box'})
    b = r['body']
    check('volume flush: a closed mesh body', b['closed'] and b['triangles'] > 1000,
          b['triangles'])
    check('inside the box (flush)', inside(b['bbox_mm'], r['source_bbox_mm'], 0.01),
          b['bbox_mm'])
    check('porosity measured near 80 %', abs(r['porosity_measured'] - 0.8) < 0.02,
          round(r['porosity_measured'], 4))
    check('Fusion volume = engine volume (units right)',
          abs(b['volume_mm3'] / r['mesh']['volume_mm3'] - 1) < 0.01,
          (round(b['volume_mm3'], 2), round(r['mesh']['volume_mm3'], 2)))
    check('the preview graph was reused by the build', r['graph_cached_on_build'])
    # MeshBody.volume is not trusted: it read -789 inside-out and 0.0 the right
    # way round. The signed volume of the body's own triangles is.
    check('the body faces outward (signed volume of its triangles > 0)',
          b['volume_mm3'] > 0, round(b['volume_mm3'], 2))
    check('the lattice got its own base feature (parametric design)',
          r['design']['parametric'] and r['design']['base_features'] == 2,
          r['design'])

    r = job('build', {'shape': 'box', 'direct': True, 'replace_check': True})
    check('direct design: built, no base feature, and replace still works',
          r['body']['closed'] and not r['design']['parametric']
          and r['replace']['mesh_bodies_after'] == 0, r['replace'])

    r = job('build', {'shape': 'box', 'params': {'placement': 'external'}})
    sb, bb = r['source_bbox_mm'], r['body']['bbox_mm']
    grow = max(max(sb[0][k] - bb[0][k], bb[1][k] - sb[1][k]) for k in range(3))
    check('external: stands about one strut outside the box', 0.4 < grow < 0.75,
          round(grow, 3))

    r = job('build', {'shape': 'ball', 'params': {'kind': 'surface',
                                                  'placement': 'middle'}})
    check('ball surface lattice: closed body',
          r['body']['closed'], r['porosity_measured'])

    r = job('build', {'shape': 'face', 'params': {'placement': 'external'}})
    bb = r['body']['bbox_mm']
    check('top face, external: lattice sits above the face (z 10 to 10.6)',
          bb[0][2] > 9.95 and bb[1][2] < 10.7, (bb[0][2], bb[1][2]))
    r = job('build', {'shape': 'face', 'params': {'placement': 'flush'}})
    bb = r['body']['bbox_mm']
    check('top face, flush: lattice sits under the face (z 9.4 to 10)',
          bb[0][2] > 9.3 and bb[1][2] < 10.05, (bb[0][2], bb[1][2]))

    r = job('build', {'shape': 'occ_box'})
    sb, bb = r['source_bbox_mm'], r['body']['bbox_mm']
    check('box in a moved occurrence: lattice goes in that component',
          r['body']['component'] != '' and r['body']['component'] != 'root'
          and 'Untitled' not in r['body']['component'], r['body']['component'])

    r = job('build', {'shape': 'box', 'replace_check': True})
    rep = r['replace']
    check('"replace" removes the lattice and nothing else',
          rep['removed'] >= 1 and rep['mesh_bodies_after'] == rep['mesh_bodies_before'] - 1,
          rep)
    print('Form output')
    for o in (1, 10):
        r = job('build', {'shape': 'ball', 'params': {
            'kind': 'surface', 'placement': 'middle', 'output': 'form', 'organic': o}})
        b = r['body']
        check('ball surface, organic %d: a Form feature and a closed solid' % o,
              r['output'] == 'form' and b['solid'] and b['form_features'] == 1,
              (r['output'], b.get('faces')))
        check('organic %d: Fusion volume = engine limit volume within 1.5 %%' % o,
              abs(b['volume_mm3'] / r['form']['volume_mm3'] - 1) < 0.015,
              (round(b['volume_mm3'], 2), round(r['form']['volume_mm3'], 2)))
        check('organic %d: porosity near 80 %%' % o,
              abs(r['porosity_measured'] - 0.8) < 0.02, round(r['porosity_measured'], 4))
    r = job('build', {'shape': 'face', 'params': {'placement': 'external',
                                                  'output': 'form', 'organic': 5}})
    bb, sb = r['body']['bbox_mm'], r['source_bbox_mm']
    check('top face, external, Form: a solid sitting on the face',
          r['output'] == 'form' and r['body']['solid']
          and bb[0][2] > 9.9 and bb[1][2] < 10.9, (bb[0][2], bb[1][2]))
    r = job('build', {'shape': 'box', 'params': {'output': 'form', 'organic': 5}})
    check('dense volume lattice asked as a Form: built as a mesh, and says so',
          r['output'] == 'mesh' and r['body']['closed']
          and 'NOT A FORM' in r['text'], r['graph'].get('form_failed', '')[:70])
    r = job('build', {'shape': 'ball', 'direct': True, 'params': {
        'kind': 'surface', 'output': 'form'}})
    check('direct design asked for a Form: a mesh, with the reason',
          r['output'] == 'mesh' and 'parametric' in (r['form_skipped'] or ''),
          r['form_skipped'])
    r = job('build', {'shape': 'ball', 'replace_check': True, 'params': {
        'kind': 'surface', 'output': 'form', 'organic': 3}})
    check('"replace" removes a Form lattice too',
          r['replace']['mesh_bodies_after'] == r['replace']['mesh_bodies_before'] - 1,
          r['replace'])
    print('organic, mesh output')
    r1 = job('build', {'shape': 'box', 'params': {'define': 'cell', 'cell': 3.0,
                                                  'organic': 1}})
    r9 = job('build', {'shape': 'box', 'params': {'define': 'cell', 'cell': 3.0,
                                                  'organic': 9}})
    check('organic 9 puts more material at the nodes than organic 1',
          r9['body']['volume_mm3'] > 1.3 * r1['body']['volume_mm3'],
          (round(r1['body']['volume_mm3'], 1), round(r9['body']['volume_mm3'], 1)))
    print()
    print('%d failed: %s' % (len(FAILED), FAILED) if FAILED else 'all live checks passed')
    return 1 if FAILED else 0


if __name__ == '__main__':
    raise SystemExit(main())
