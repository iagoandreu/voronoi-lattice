"""How ORGANIC changes the node term of the mesh output: BETA per level,
measured from meshed volume like tests/calibrate.py.

    py -3.10 tests/calibrate_organic.py
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
import numpy as np                                            # noqa: E402
import trimesh                                                # noqa: E402
import vl_engine as E                                         # noqa: E402
import vl_plan as P                                           # noqa: E402


def main():
    box = E.prepare(np.asarray(trimesh.creation.box((24, 24, 24)).vertices, float),
                    np.asarray(trimesh.creation.box((24, 24, 24)).faces), None, True)
    ball = trimesh.creation.icosphere(5, 12.5)
    ball = E.prepare(np.asarray(ball.vertices, float), np.asarray(ball.faces), None, True)
    for name, R, kind, s in (('box', box, 'volume', 3.0), ('ball', ball, 'surface', 3.0)):
        for d in (0.4, 0.6):
            job = dict(kind=kind, placement='middle', thickness=d, regularity=5,
                       skin=False, border=False, closed=True, seed=1)
            A, B, info = E.build_graph(R, job, s, 1)
            L, J = info['length_mm'], info['junctions']
            r = d / 2
            row = []
            for o in (1, 4, 7, 10):
                mo = P.mesh_organic(o, d)
                surf, _ = E.mesh_struts(A, B, d, d / 4, blend=mo['blend'])
                surf, _ = E.finish_mesh(surf, 10 ** 9, mo['smooth'])
                beta = (math.pi * r * r * L - E.signed_volume(E.outward(surf))) / (r ** 3 * J)
                row.append('o%-2d %6.2f' % (o, beta))
            print('%-4s d %.1f  ' % (name, d) + '   '.join(row), flush=True)


if __name__ == '__main__':
    main()
