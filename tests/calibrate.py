"""Measure vl_plan's constants: strut length and junction density of the
lattices build_graph makes, and BETA, the node overlap, from meshed volume.

    py -3.10 tests/calibrate.py

Volume: a 30 mm box. Surface: a 60 mm square sheet (open, border off, so its
own rim does not count) and a 25 mm sphere (closed, curved).
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
import numpy as np                                            # noqa: E402
import trimesh                                                # noqa: E402
import vl_engine as E                                         # noqa: E402


def region(m, closed):
    return E.prepare(np.asarray(m.vertices, float), np.asarray(m.faces), None, closed)


def main():
    box = region(trimesh.creation.box((30, 30, 30)), True)
    sheet_m = trimesh.creation.box((60, 60, 1))
    sheet_m.update_faces(sheet_m.face_normals[:, 2] > 0.9)
    sheet_m.remove_unreferenced_vertices()
    sheet = region(sheet_m, False)
    ball = region(trimesh.creation.icosphere(5, 12.5), True)
    rows = []
    for name, R, kind in (('box', box, 'volume'), ('sheet', sheet, 'surface'),
                          ('ball', ball, 'surface')):
        for s in (2.0, 3.0, 4.5):
            for d in (0.3, 0.6):
                job = dict(kind=kind, placement='middle', thickness=d, regularity=5,
                           skin=False, border=False, closed=R['closed'], seed=1)
                A, B, info = E.build_graph(R, job, s, 1)
                L, J = info['length_mm'], info['junctions']
                if kind == 'volume':
                    c, j = L * s * s / R['volume'], J * s ** 3 / R['volume']
                else:
                    c, j = L * s / R['area'], J * s * s / R['area']
                surf, _ = E.mesh_struts(A, B, d, d / 4)
                surf, _ = E.finish_mesh(surf, 10 ** 9, True)
                r = d / 2
                beta = (math.pi * r * r * L - surf.volume) / (r ** 3 * J)
                rows.append((name, s, d, c, j, beta))
                print('%-5s s %.1f d %.1f   C %.3f  J %.3f  beta %.3f'
                      % (name, s, d, c, j, beta), flush=True)
    for name in ('box', 'sheet', 'ball'):
        sel = [x for x in rows if x[0] == name]
        print('%-5s mean C %.3f  J %.3f  beta %.3f' % (
            name, np.mean([x[3] for x in sel]), np.mean([x[4] for x in sel]),
            np.mean([x[5] for x in sel])))


if __name__ == '__main__':
    main()
