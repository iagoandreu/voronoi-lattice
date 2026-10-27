"""Voronoi Lattice engine. Runs OUTSIDE Fusion, in a CPython with numpy, scipy,
trimesh and pyvista (py -3.10 on this machine):

    py -3.10 vl_engine.py job.json out.json

Everything is millimetres. vl_core converts to and from Fusion's centimetres.

Two stages, so the dialog can preview the cheap one:

  graph   strut centrelines. A VOLUME lattice is the 3D Voronoi diagram of
          Lloyd-relaxed seeds, clipped to the region; its skin struts (optional)
          are the same seeds' Voronoi diagram restricted to the surface, so they
          meet the clipped edges where those leave the volume. A SURFACE lattice
          is the restricted Voronoi diagram of seeds placed on the surface.
          Placement offsets every strut-centre point along the outward normal by
          -r (flush), 0 (middle) or +r (external).
  mesh    a watertight triangle mesh of the struts: a distance field to the
          centrelines sampled on a voxel grid, slab by slab, contoured by VTK's
          flying edges, Taubin-smoothed, decimated to the triangle budget.

The cell size comes from vl_plan, then is corrected against the length and
junction count of the graph actually built, until the predicted porosity is
within 1 % of the target. The porosity reported at the end is MEASURED from the
mesh volume, never the target echoed back.

The restricted Voronoi diagram: the surface is subdivided finer than the cells,
every vertex is labelled with its nearest seed, and each triangle edge whose two
ends disagree is cut where it crosses the bisector of the two seeds - exact,
because |x - si|^2 - |x - sj|^2 is linear in x. Triangles with three labels hold
a junction, solved in the triangle's plane. The label boundaries become polylines,
simplified by Douglas-Peucker. For seeds on a SURFACE only seeds facing the same
way count, so a thin wall's two sides do not steal each other's cells.
"""

import hashlib
import json
import math
import os
import sys
import time
import traceback

import numpy as np
from scipy.spatial import Voronoi, cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import vl_plan as P                                           # noqa: E402
import vl_form as FORM                                        # noqa: E402

GRAPH_VERSION = 4


class LatticeError(RuntimeError):
    pass


# --------------------------------------------------------------------------- #
# the region
# --------------------------------------------------------------------------- #

def weld(V, F, tol):
    """Merge vertices closer than tol. Fusion tessellates face by face, so a
    selection of faces arrives with every seam vertex duplicated."""
    pairs = cKDTree(V).query_pairs(tol, output_type='ndarray')
    n = len(V)
    if len(pairs):
        g = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])),
                       shape=(n, n))
        _, comp = connected_components(g, directed=False)
    else:
        comp = np.arange(n)
    _, first, inv = np.unique(comp, return_index=True, return_inverse=True)
    V2 = V[first]
    F2 = inv.reshape(-1)[F]
    keep = (F2[:, 0] != F2[:, 1]) & (F2[:, 1] != F2[:, 2]) & (F2[:, 0] != F2[:, 2])
    return V2, F2[keep]


def tri_normals(V, F):
    n = np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]])
    a = np.linalg.norm(n, axis=1)
    return n / np.maximum(a, 1e-30)[:, None], 0.5 * a


def vertex_normals(V, F):
    n = np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]])
    vn = np.zeros_like(V)
    for k in range(3):
        np.add.at(vn, F[:, k], n)
    return vn / np.maximum(np.linalg.norm(vn, axis=1), 1e-30)[:, None]


def boundary_edges(F):
    e = np.sort(np.vstack([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]]), axis=1)
    u, c = np.unique(e, axis=0, return_counts=True)
    return u[c == 1]


def prepare(V, F, N, closed):
    """Weld, orient, measure. Orientation matters: placement offsets along the
    OUTWARD normal. A closed region is made to enclose positive volume; an open
    one follows the normals Fusion supplied (the face normals), by majority."""
    import trimesh
    diag = float(np.linalg.norm(V.max(0) - V.min(0)))
    V, F = weld(V, F, max(diag * 1e-7, 1e-6))
    m = trimesh.Trimesh(V, F, process=False)
    trimesh.repair.fix_winding(m)
    V, F = np.asarray(m.vertices), np.asarray(m.faces)
    open_edges = len(boundary_edges(F))
    if closed:
        if open_edges:
            raise LatticeError('the region was reported closed but its mesh has '
                               '%d open edges' % open_edges)
        if m.volume < 0:
            F = F[:, ::-1]
    elif N is not None:
        Vn_given = N
        tn, _ = tri_normals(V, F)
        # the supplied normals belong to the UNWELDED vertices; match by position
        _, idx = cKDTree(Vn_given[:, :3]).query(V[F].mean(1))
        agree = np.einsum('ij,ij->i', tn, Vn_given[idx, 3:])
        if np.sum(agree < 0) > np.sum(agree > 0):
            F = F[:, ::-1]
    tn, ta = tri_normals(V, F)
    vol = 0.0
    if closed:
        vol = float(np.einsum('ij,ij->i', V[F[:, 0]],
                              np.cross(V[F[:, 1]], V[F[:, 2]])).sum() / 6.0)
    return {'V': V, 'F': F, 'area': float(ta.sum()), 'volume': vol,
            'open_edges': open_edges, 'closed': bool(closed), 'diag': diag}


class SignedDistance:
    """Signed distance to a closed triangle mesh, negative inside. VTK's
    vtkImplicitPolyDataDistance, evaluated over an array in one call."""

    def __init__(self, V, F):
        import vtk
        from vtk.util import numpy_support as ns
        self._ns = ns
        self._vtk = vtk
        pts = vtk.vtkPoints()
        pts.SetData(ns.numpy_to_vtk(np.ascontiguousarray(V, dtype=np.float64),
                                    deep=True))
        cells = vtk.vtkCellArray()
        flat = np.hstack([np.full((len(F), 1), 3, np.int64),
                          F.astype(np.int64)]).ravel()
        cells.SetCells(len(F), ns.numpy_to_vtkIdTypeArray(flat, deep=True))
        poly = vtk.vtkPolyData()
        poly.SetPoints(pts)
        poly.SetPolys(cells)
        self._poly = poly
        self._f = vtk.vtkImplicitPolyDataDistance()
        self._f.SetInput(poly)

    def __call__(self, X):
        X = np.ascontiguousarray(X, dtype=np.float64)
        if len(X) == 0:
            return np.zeros(0)
        out = self._vtk.vtkDoubleArray()
        out.SetNumberOfTuples(len(X))
        self._f.FunctionValue(self._ns.numpy_to_vtk(X, deep=True), out)
        return self._ns.vtk_to_numpy(out).copy()


def refine(V, F, max_edge):
    import trimesh
    if max_edge <= 0:
        return V, F
    V2, F2 = trimesh.remesh.subdivide_to_size(V, F, max_edge, max_iter=12)
    return np.asarray(V2, float), np.asarray(F2, np.int64)


# --------------------------------------------------------------------------- #
# seeds
# --------------------------------------------------------------------------- #

def surface_samples(V, F, Vn, M, rng):
    _, A = tri_normals(V, F)
    tid = rng.choice(len(F), int(M), p=A / A.sum())
    u, v = rng.random((2, int(M)))
    flip = u + v > 1.0
    u[flip], v[flip] = 1.0 - u[flip], 1.0 - v[flip]
    w = np.stack([1.0 - u - v, u, v], 1)
    T = F[tid]
    X = np.einsum('mk,mkj->mj', w, V[T])
    N = np.einsum('mk,mkj->mj', w, Vn[T])
    N /= np.maximum(np.linalg.norm(N, axis=1), 1e-30)[:, None]
    return X, N


def nearest_seed(X, seeds, Xn=None, Sn=None):
    """Index of the nearest seed; with normals, the nearest one FACING the same
    way, so a thin wall's other side is never chosen."""
    tree = cKDTree(seeds)
    if Xn is None or len(seeds) == 1:
        _, lab = tree.query(X, workers=-1)
        return lab
    k = min(8, len(seeds))
    _, idx = tree.query(X, k=k, workers=-1)
    ok = np.einsum('nkj,nj->nk', Sn[idx], Xn) > 0.0
    first = np.where(ok.any(1), ok.argmax(1), 0)
    return idx[np.arange(len(X)), first]


def lloyd(seeds, X, iters, Xn=None, Sn=None):
    """Move each seed to the centroid of the samples nearest to it, then snap it
    onto the nearest of those samples - which keeps a seed on the surface, and
    inside a non-convex volume, where the bare centroid might not be."""
    for _ in range(int(iters)):
        lab = nearest_seed(X, seeds, Xn, Sn)
        cnt = np.bincount(lab, minlength=len(seeds)).astype(float)
        live = cnt > 0
        c = np.zeros_like(seeds)
        for k in range(3):
            c[:, k] = np.bincount(lab, weights=X[:, k], minlength=len(seeds))
        c[live] /= cnt[live, None]
        c[~live] = seeds[~live]
        _, near = cKDTree(X).query(c, workers=-1)
        seeds = X[near]
        if Sn is not None:
            Sn = Xn[near]
    return seeds, Sn


def volume_samples(phi, lo, hi, M, rng):
    out = []
    got = 0
    for _ in range(60):
        X = lo + (hi - lo) * rng.random((max(int(M), 20000), 3))
        X = X[phi(X) < 0.0]
        out.append(X)
        got += len(X)
        if got >= M:
            break
    X = np.vstack(out)[:int(M)]
    if len(X) < 8:
        raise LatticeError('the region is too thin for this placement - nothing '
                           'is left inside it once the struts are kept off the '
                           'surface')
    return X


# --------------------------------------------------------------------------- #
# polylines
# --------------------------------------------------------------------------- #

def douglas_peucker(Pts, tol):
    n = len(Pts)
    if n <= 2:
        return np.arange(n)
    keep = np.zeros(n, bool)
    keep[0] = keep[-1] = True
    stack = [(0, n - 1)]
    while stack:
        i, j = stack.pop()
        if j <= i + 1:
            continue
        a, ab = Pts[i], Pts[j] - Pts[i]
        L2 = float(ab @ ab)
        seg = Pts[i + 1:j] - a
        if L2 < 1e-24:
            dist = np.linalg.norm(seg, axis=1)
        else:
            t = np.clip(seg @ ab / L2, 0.0, 1.0)
            dist = np.linalg.norm(seg - t[:, None] * ab, axis=1)
        k = int(np.argmax(dist))
        if dist[k] > tol:
            m = i + 1 + k
            keep[m] = True
            stack.append((i, m))
            stack.append((m, j))
    return np.nonzero(keep)[0]


def chains(n_nodes, E):
    """Split a graph into maximal paths between nodes whose degree is not 2;
    leftover cycles come out as closed chains (first node repeated last)."""
    adj = [[] for _ in range(n_nodes)]
    for k, (a, b) in enumerate(E):
        adj[a].append((b, k))
        adj[b].append((a, k))
    deg = np.array([len(x) for x in adj])
    used = np.zeros(len(E), bool)
    out = []

    def walk(start, nxt, k):
        path = [start, nxt]
        used[k] = True
        prev = start
        cur = nxt
        while deg[cur] == 2:
            (n1, k1), (n2, k2) = adj[cur]
            if not used[k1]:
                step, kk = n1, k1
            elif not used[k2]:
                step, kk = n2, k2
            else:
                break
            used[kk] = True
            path.append(step)
            prev, cur = cur, step
        return path

    for v in range(n_nodes):
        if deg[v] != 2:
            for nb, k in adj[v]:
                if not used[k]:
                    out.append(walk(v, nb, k))
    for k in range(len(E)):
        if not used[k]:
            a, b = E[k]
            out.append(walk(a, b, k))
    return out, deg


def polyline_segments(X, Nrm, paths, tol, shift):
    """Simplify each path, offset its points along the normal, return segments."""
    A, B = [], []
    Y = X + shift * Nrm
    for p in paths:
        p = np.asarray(p)
        closed = len(p) > 2 and p[0] == p[-1]
        if closed:
            half = len(p) // 2
            parts = [p[:half + 1], p[half:]]
        else:
            parts = [p]
        for q in parts:
            kept = q[douglas_peucker(X[q], tol)]
            A.append(Y[kept[:-1]])
            B.append(Y[kept[1:]])
    if not A:
        return np.zeros((0, 3)), np.zeros((0, 3))
    return np.vstack(A), np.vstack(B)


# --------------------------------------------------------------------------- #
# restricted Voronoi diagram on a triangle mesh
# --------------------------------------------------------------------------- #

def restricted_network(V, F, Vn, seeds, lab):
    """Nodes (positions, normals, is_junction) and edges of the label boundaries."""
    nV = len(V)
    la = lab[F]
    ea = F[:, [0, 1, 2]]
    eb = F[:, [1, 2, 0]]
    cross = la[:, [0, 1, 2]] != la[:, [1, 2, 0]]
    t_idx, e_idx = np.nonzero(cross)
    if len(t_idx) == 0:
        return (np.zeros((0, 3)), np.zeros((0, 3)), np.zeros(0, bool),
                np.zeros((0, 2), int), np.zeros(0, np.int64))
    a = ea[t_idx, e_idx]
    b = eb[t_idx, e_idx]
    key = np.minimum(a, b).astype(np.int64) * nV + np.maximum(a, b)
    ukey, inv = np.unique(key, return_inverse=True)
    ua, ub = ukey // nV, ukey % nV
    si, sj = seeds[lab[ua]], seeds[lab[ub]]
    Va, Vb = V[ua], V[ub]

    def f(x):
        return np.einsum('ij,ij->i', x - si, x - si) - np.einsum('ij,ij->i', x - sj, x - sj)

    fa, fb = f(Va), f(Vb)
    den = fa - fb
    t = np.where(np.abs(den) > 1e-30, fa / np.where(den == 0, 1, den), 0.5)
    t = np.clip(t, 0.0, 1.0)[:, None]
    X = Va + t * (Vb - Va)
    Nn = Vn[ua] + t * (Vn[ub] - Vn[ua])
    node = np.full(cross.shape, -1, np.int64)
    node[t_idx, e_idx] = inv.reshape(-1)

    ncross = cross.sum(1)
    two = np.nonzero(ncross == 2)[0]
    E2 = np.sort(node[two], axis=1)[:, 1:]          # the two non -1 entries

    three = np.nonzero(ncross == 3)[0]
    nj = len(three)
    T = F[three]
    s0, s1, s2 = seeds[la[three, 0]], seeds[la[three, 1]], seeds[la[three, 2]]
    tn, _ = tri_normals(V, T) if nj else (np.zeros((0, 3)), None)
    M = np.stack([2.0 * (s1 - s0), 2.0 * (s2 - s0), tn], 1)
    rhs = np.stack([np.einsum('ij,ij->i', s1, s1) - np.einsum('ij,ij->i', s0, s0),
                    np.einsum('ij,ij->i', s2, s2) - np.einsum('ij,ij->i', s0, s0),
                    np.einsum('ij,ij->i', tn, V[T[:, 0]])], 1)
    fallback = X[node[three]].mean(1) if nj else np.zeros((0, 3))
    J = fallback.copy()
    if nj:
        det = np.linalg.det(M)
        good = np.abs(det) > 1e-12 * np.maximum(1.0, np.abs(M).max((1, 2)) ** 3)
        if good.any():
            J[good] = np.linalg.solve(M[good], rhs[good][..., None])[..., 0]
        cen = V[T].mean(1)
        span = np.linalg.norm(V[T] - cen[:, None], axis=2).max(1)
        far = np.linalg.norm(J - cen, axis=1) > 1.5 * span
        J[far] = fallback[far]
    Jn = Vn[T].mean(1) if nj else np.zeros((0, 3))
    base = len(X)
    jid = base + np.arange(nj)
    E3 = np.vstack([np.stack([jid, node[three, k]], 1) for k in range(3)]) if nj \
        else np.zeros((0, 2), np.int64)
    Xall = np.vstack([X, J])
    Nall = np.vstack([Nn, Jn])
    Nall /= np.maximum(np.linalg.norm(Nall, axis=1), 1e-30)[:, None]
    is_j = np.zeros(len(Xall), bool)
    is_j[base:] = True
    E = np.vstack([E2, E3]).astype(np.int64)
    return Xall, Nall, is_j, E, ukey


def border_segments(V, Vn, paths, ukey, X, Nn, tol, shift):
    """The border strut, THROUGH the points where the cell net meets the
    boundary: those are crossing nodes on boundary edges, so inserting them
    makes border and net share their end points exactly - a junction, not two
    loose ends side by side (which a mesh would fuse but a Form cannot)."""
    nV = len(V)
    A, B = [], []
    for path in paths:
        pts, nrm, forced = [], [], []
        for i in range(len(path) - 1):
            u, v = path[i], path[i + 1]
            pts.append(V[u])
            nrm.append(Vn[u])
            forced.append(False)
            key = min(u, v) * nV + max(u, v)
            j = int(np.searchsorted(ukey, key))
            if j < len(ukey) and ukey[j] == key:
                pts.append(X[j])
                nrm.append(Nn[j])
                forced.append(True)
        pts.append(V[path[-1]])
        nrm.append(Vn[path[-1]])
        forced.append(False)
        P_, N_ = np.asarray(pts), np.asarray(nrm)
        Y = P_ + shift * N_
        cut = [0] + [i for i, f in enumerate(forced) if f] + [len(pts) - 1]
        cut = sorted(set(cut))
        for i0, i1 in zip(cut[:-1], cut[1:]):
            idx = np.arange(i0, i1 + 1)
            kept = idx[douglas_peucker(P_[idx], tol)]
            A.append(Y[kept[:-1]])
            B.append(Y[kept[1:]])
    if not A:
        return np.zeros((0, 3)), np.zeros((0, 3))
    return np.vstack(A), np.vstack(B)


def border_paths(V, F):
    be = boundary_edges(F)
    if len(be) == 0:
        return []
    paths, _ = chains(len(V), be)
    return paths


# --------------------------------------------------------------------------- #
# 3D Voronoi, clipped
# --------------------------------------------------------------------------- #

def ghost_box(lo, hi, s):
    """A shell of points around the region, so every cell that matters is
    bounded and its vertices are finite."""
    pad = 2.0 * s
    a, b = lo - pad, hi + pad
    axes = [np.linspace(a[k], b[k], max(2, int(math.ceil((b[k] - a[k]) / s)) + 1))
            for k in range(3)]
    g = np.stack(np.meshgrid(*axes, indexing='ij'), -1).reshape(-1, 3)
    on = np.zeros(len(g), bool)
    for k in range(3):
        on |= np.isclose(g[:, k], a[k]) | np.isclose(g[:, k], b[k])
    return g[on]


def voronoi_edges(seeds, lo, hi, s):
    n = len(seeds)
    vor = Voronoi(np.vstack([seeds, ghost_box(lo, hi, s)]))
    rp = vor.ridge_points
    pairs = []
    for k, rv in enumerate(vor.ridge_vertices):
        if rp[k, 0] >= n and rp[k, 1] >= n:
            continue
        if -1 in rv:
            continue
        m = len(rv)
        for i in range(m):
            pairs.append((rv[i], rv[(i + 1) % m]))
    E = np.unique(np.sort(np.asarray(pairs, np.int64), axis=1), axis=0)
    return vor.vertices, E


def clip_edges(Vv, E, phi):
    """Keep the parts of the edges inside phi < 0; return segments and the
    indices (into the segment list) of ends that lie on the boundary."""
    used = np.unique(E)
    ph = np.full(len(Vv), np.inf)
    ph[used] = phi(Vv[used])
    ins = ph < 0.0
    a, b = E[:, 0], E[:, 1]
    both = ins[a] & ins[b]
    one = ins[a] ^ ins[b]
    A = [Vv[a[both]]]
    B = [Vv[b[both]]]
    pi = np.where(ins[a[one]], a[one], b[one])
    pe = np.where(ins[a[one]], b[one], a[one])
    Xi, Xe = Vv[pi], Vv[pe]
    lo_t, hi_t = np.zeros(len(pi)), np.ones(len(pi))
    for _ in range(16):
        mid = 0.5 * (lo_t + hi_t)
        inside = phi(Xi + mid[:, None] * (Xe - Xi)) < 0.0
        lo_t = np.where(inside, mid, lo_t)
        hi_t = np.where(inside, hi_t, mid)
    Xc = Xi + lo_t[:, None] * (Xe - Xi)
    A.append(Xi)
    B.append(Xc)
    A, B = np.vstack(A), np.vstack(B)
    ends = np.arange(int(both.sum()), len(A))           # B of these is on the skin
    return A, B, ends


# --------------------------------------------------------------------------- #
# the graph
# --------------------------------------------------------------------------- #

def prune(A, B, r, min_len):
    """Drop strut groups that touch nothing else and are shorter than min_len -
    a clipped edge's stub on a sharp corner, a sliver of a cell. Two struts are
    joined when their axes come within r of each other: fused, in the mesh."""
    L = np.linalg.norm(B - A, axis=1)
    if len(A) == 0:
        return A, B, 0
    step = max(0.5 * r, 1e-3)
    n = np.maximum(1, np.ceil(L / step).astype(int))
    sid = np.repeat(np.arange(len(A)), n + 1)
    off = np.concatenate([[0], np.cumsum(n + 1)[:-1]])
    t = (np.arange(len(sid)) - off[sid]) / n[sid]
    S = A[sid] + t[:, None] * (B - A)[sid]
    pairs = cKDTree(S).query_pairs(r, output_type='ndarray')
    if len(pairs):
        g = coo_matrix((np.ones(len(pairs)), (sid[pairs[:, 0]], sid[pairs[:, 1]])),
                       shape=(len(A), len(A)))
        nc, comp = connected_components(g, directed=False)
    else:
        nc, comp = len(A), np.arange(len(A))
    tot = np.bincount(comp, weights=L, minlength=nc)
    keep = tot[comp] >= min_len
    return A[keep], B[keep], int(np.sum(tot < min_len))


def junction_count(A, B):
    pts = np.round(np.vstack([A, B]) * 1e6).astype(np.int64)
    _, c = np.unique(pts, axis=0, return_counts=True)
    return int(np.sum(c >= 3))


def build_graph(R, job, s, rng_seed):
    """Centrelines for cell size s. R is the prepared region."""
    rng = np.random.default_rng(int(rng_seed))
    d = float(job['thickness'])
    r = 0.5 * d
    kind = job['kind']
    shift = P.shift_for(job['placement'], d)
    if job.get('flip') and not R['closed']:
        shift = -shift
    iters = int(job.get('regularity', 5))
    V, F = R['V'], R['F']
    # the surface mesh the cells are traced on: finer than a cell and than a strut
    edge = min(s / 8.0, max(r, s / 30.0))
    edge = max(edge, math.sqrt(R['area'] / 2.5e6))
    Vf, Ff = refine(V, F, edge)
    Vn = vertex_normals(Vf, Ff)
    # chords may leave the traced curve by tol: keep it a small part of a strut
    tol = min(0.1 * r, 0.05 * s)
    info = {'cell_mm': s, 'trace_edge_mm': edge}
    parts_A, parts_B = [], []

    if kind == 'volume':
        sdf = SignedDistance(V, F)

        def phi(X):
            return sdf(X) - shift

        lo, hi = V.min(0) - abs(shift), V.max(0) + abs(shift)
        env = P.envelope('volume', job['placement'], d, R['volume'], R['area'])
        n = max(4, int(round(env / s ** 3)))
        Xs = volume_samples(phi, lo, hi, min(max(30 * n, 30000), 3_000_000), rng)
        seeds = Xs[rng.choice(len(Xs), min(n, len(Xs)), replace=False)]
        seeds, _ = lloyd(seeds, Xs, iters)
        Vv, E = voronoi_edges(seeds, lo, hi, s)
        A, B, ends = clip_edges(Vv, E, phi)
        info['seeds'] = len(seeds)
        if job.get('skin', True):
            lab = nearest_seed(Vf, seeds)
            X, Nn, is_j, Es, _ = restricted_network(Vf, Ff, Vn, seeds, lab)
            paths, _ = chains(len(X), Es)
            sA, sB = polyline_segments(X, Nn, paths, tol, shift)
            # a clipped edge leaves the volume at a triple point of the skin;
            # land it exactly on the skin junction it is aiming at
            Jp = X[is_j] + shift * Nn[is_j]
            if len(Jp) and len(ends):
                dd, k = cKDTree(Jp).query(B[ends], workers=-1)
                hit = dd < 0.4 * s
                B[ends[hit]] = Jp[k[hit]]
            parts_A.append(sA)
            parts_B.append(sB)
            info['skin_struts'] = len(sA)
        parts_A.insert(0, A)
        parts_B.insert(0, B)
    else:
        n = max(3, int(round(R['area'] / s ** 2)))
        Xs, Xn = surface_samples(Vf, Ff, Vn, min(max(40 * n, 40000), 3_000_000), rng)
        pick = rng.choice(len(Xs), min(n, len(Xs)), replace=False)
        seeds, Sn = lloyd(Xs[pick], Xs, iters, Xn, Xn[pick])
        lab = nearest_seed(Vf, seeds, Vn, Sn)
        X, Nn, is_j, Es, ukey = restricted_network(Vf, Ff, Vn, seeds, lab)
        paths, _ = chains(len(X), Es)
        A, B = polyline_segments(X, Nn, paths, tol, shift)
        parts_A.append(A)
        parts_B.append(B)
        info['seeds'] = len(seeds)
        if not R['closed'] and job.get('border', True):
            bA, bB = border_segments(Vf, Vn, border_paths(Vf, Ff), ukey, X, Nn,
                                     tol, shift)
            parts_A.append(bA)
            parts_B.append(bB)
            info['border_struts'] = len(bA)

    A = np.vstack([a for a in parts_A if len(a)]) if any(len(a) for a in parts_A) \
        else np.zeros((0, 3))
    B = np.vstack([b for b in parts_B if len(b)]) if any(len(b) for b in parts_B) \
        else np.zeros((0, 3))
    short = np.linalg.norm(B - A, axis=1) > 1e-9
    A, B = A[short], B[short]
    A, B, dropped = prune(A, B, r, max(1.5 * s, 3.0 * d))
    if len(A) == 0:
        raise LatticeError('no struts came out - the cells are larger than the '
                           'region; lower the cell size or raise the porosity')
    L = float(np.linalg.norm(B - A, axis=1).sum())
    info.update({'struts': len(A), 'length_mm': L, 'junctions': junction_count(A, B),
                 'dropped_groups': dropped, 'shift_mm': shift})
    return A, B, info


def target_solid(job, R):
    d = float(job['thickness'])
    phi = float(job.get('porosity', 0.5))
    env = P.envelope(job['kind'], job['placement'], d, R['volume'], R['area'])
    if job['kind'] == 'volume':
        return (1.0 - phi) * env, env
    return (1.0 - phi) * env * math.pi * d / 4.0, env


def porosity_from(job, R, vol):
    """Volume lattice: void fraction of the envelope. Surface lattice: open
    fraction of the surface, 1 - 4 V / (pi d A) - the footprint of struts of
    diameter d holding volume V."""
    d = float(job['thickness'])
    _, env = target_solid(job, R)
    if job['kind'] == 'volume':
        return 1.0 - vol / env
    return 1.0 - 4.0 * vol / (math.pi * d * env)


class Measure:
    """The volume of what each output will actually be, for the porosity loop:
    the mesh output's own field (organic blend and flush trim included), or the
    Form's Catmull-Clark limit surface. Neither is a model."""

    def __init__(self, job, R):
        self.job, self.R = job, R
        self.d = float(job['thickness'])
        self.h = self.d / P.RESOLUTION[job.get('resolution', 'medium')]
        self.organic = float(job.get('organic', 1))
        self.trim = self.near = None
        if job['placement'] == 'flush' and R['closed']:
            self.trim = SignedDistance(R['V'], R['F'])
            # the trim is only evaluated near the surface; "near" is judged
            # against points no further apart than a strut, NOT the region's own
            # vertices - a box has eight, and its faces went untrimmed
            self.near = refine(R['V'], R['F'], r_near(self.d))[0]

    def field(self, A, B):
        mo = P.mesh_organic(self.organic, self.d)
        return StrutField(A, B, self.d, self.h, mo['blend'], self.trim, self.near)

    def volume(self, A, B, output, levels=2):
        if output == 'form':
            V, Fc, ci = FORM.build_cage(A, B, self.d, self.organic)
            Pl, T = FORM.limit_mesh(V, Fc, levels)
            return FORM.mesh_volume(Pl, T), (V, Fc, ci, Pl)
        return self.field(A, B).solid_volume(), None


def graph_for(job, R, log=None, measure=None):
    """Build the graph and correct the cell size until the porosity of what
    will be BUILT is within 1.5 % of the target (porosity mode), or once
    (cell-size mode). A Form that cannot be built falls back to a mesh, and the
    reason travels with the result."""
    log = log or (lambda _m: None)
    d = float(job['thickness'])
    kind = job['kind']
    M = measure or Measure(job, R)
    output = job.get('output', 'mesh')
    fallback = None
    by_cell = job.get('define', 'porosity') == 'cell'
    s = float(job['cell']) if by_cell else P.cell_for(kind, float(job['porosity']), d)
    want, _ = target_solid(job, R)
    best = None
    hist = []
    for it in range(8):
        A, B, info = build_graph(R, job, s, job.get('seed', 1))
        try:
            vol, extra = M.volume(A, B, output)
        except FORM.CageError as exc:
            if output != 'form':
                raise
            fallback = str(exc)
            output = 'mesh'
            log('No Form possible - building a mesh: %s' % fallback)
            vol, extra = M.volume(A, B, output)
        info.update({'iterations': it + 1, 'output': output,
                     'volume_predicted_mm3': vol,
                     'porosity_predicted': porosity_from(job, R, vol)})
        if fallback:
            info['form_failed'] = fallback
        if extra is not None:
            info['cage'] = extra[2]
        if by_cell:
            return A, B, info
        q = vol / want
        info['predicted_over_target'] = q
        log('Cell size %.3f mm: porosity %.1f %%'
            % (s, 100.0 * info['porosity_predicted']))
        if best is None or abs(q - 1.0) < abs(best[2]['predicted_over_target'] - 1.0):
            best = (A, B, info)
        if abs(q - 1.0) < 0.015:
            break
        # V ~ s^-k: k = 2 for a volume lattice, 1 for a surface one - until the
        # nodes swell (organic webs scale with the node count, ~1/s^2 per area)
        # and a fixed k stalls: 86 % for 80 % at organic 10. So once two builds
        # exist, k is measured off them (a secant step in log-log).
        hist.append((math.log(s), math.log(vol)))
        k = 2.0 if kind == 'volume' else 1.0
        if len(hist) >= 2 and abs(hist[-1][0] - hist[-2][0]) > 1e-6:
            k = -(hist[-1][1] - hist[-2][1]) / (hist[-1][0] - hist[-2][0])
            k = min(max(k, 0.5), 4.0)
        s *= q ** (1.0 / k)
    return best


# --------------------------------------------------------------------------- #
# the mesh
# --------------------------------------------------------------------------- #

K_CAND = 6
K_BLEND = 16
# Flying edges puts its triangles on chords INSIDE a convex surface, so a strut
# came out 5 % light at 4 voxels per diameter. Contouring at +0.11 h^2/r (a
# fraction of the chord sagitta h^2/8r) brings it within +-1 % at 3, 4 and 6
# voxels per diameter, smoothing included - measured on struts of 0.4 and 1 mm
# in random directions (tests/test_engine.py keeps it honest).
ISO_BIAS = 0.11


def strut_samples(A, B, step, with_ids=False):
    L = np.linalg.norm(B - A, axis=1)
    n = np.maximum(1, np.ceil(L / step).astype(int))
    sid = np.repeat(np.arange(len(A)), n + 1)
    off = np.concatenate([[0], np.cumsum(n + 1)[:-1]])
    t = (np.arange(len(sid)) - off[sid]) / n[sid]
    S = A[sid] + t[:, None] * (B - A)[sid]
    return (S, sid) if with_ids else S


def r_near(d):
    return 0.5 * d


class StrutField:
    """Distance field to the struts: the EXACT distance to the nearest strut
    axis, minus r. A KD-tree over samples taken every h along the axes only
    proposes candidates: the K nearest samples name their segments, and the true
    point-to-segment distance is taken over those. (Using the sample distance
    itself - a union of spheres - scallops the surface by ~3 % of r between
    samples and made a strut 4 % too thin.)

    blend > 0 is ORGANIC: an exponential smooth union of the strut distances.
    trim, if given, is a signed distance (negative inside) the field is cut to,
    so a flush lattice cannot poke through the surface it is flush with.

    The mesher and the porosity loop both read this one field, so the volume the
    cell size is solved for is the volume that gets meshed."""

    def __init__(self, A, B, d, h, blend=0.0, trim=None, surf_pts=None, iso=None):
        self.A, self.B = A, B
        self.r = 0.5 * d
        self.h = h
        self.far = 2.0 * h
        self.blend = float(blend)
        self.iso = ISO_BIAS * h * h / self.r if iso is None else iso
        S, self.sid = strut_samples(A, B, h, with_ids=True)
        self.tree = cKDTree(S)
        self.AB = B - A
        self.AB2 = np.maximum(np.einsum('ij,ij->i', self.AB, self.AB), 1e-30)
        self.K = K_CAND if self.blend <= 0 else K_BLEND
        self.reach = self.r + self.far + h + 3.0 * self.blend
        self.trim = trim
        self.surf_tree = (cKDTree(surf_pts)
                          if (trim is not None and surf_pts is not None) else None)
        pad = self.r + 3.0 * h + 3.0 * self.blend
        self.lo = np.minimum(A.min(0), B.min(0)) - pad
        self.hi = np.maximum(A.max(0), B.max(0)) + pad

    def axis_distance(self, X):
        dist, j = self.tree.query(X, k=self.K, distance_upper_bound=self.reach,
                                  workers=-1)
        E_ = np.full(dist.shape, np.inf)
        S_ = np.full(dist.shape, -1, np.int64)
        for c in range(self.K):
            ok = np.isfinite(dist[:, c])
            if not ok.any():
                continue
            seg = self.sid[j[ok, c]]
            Xa = X[ok] - self.A[seg]
            t = np.clip(np.einsum('ij,ij->i', Xa, self.AB[seg]) / self.AB2[seg],
                        0.0, 1.0)
            E_[ok, c] = np.linalg.norm(Xa - t[:, None] * self.AB[seg], axis=1)
            S_[ok, c] = seg
        if self.blend <= 0:
            return E_.min(1)
        # each strut counted ONCE - two samples of one strut would fatten it
        order = np.argsort(S_, axis=1)
        S_ = np.take_along_axis(S_, order, 1)
        E_ = np.take_along_axis(E_, order, 1)
        dup = np.zeros(S_.shape, bool)
        dup[:, 1:] = (S_[:, 1:] == S_[:, :-1]) & (S_[:, 1:] >= 0)
        for c in range(self.K - 1, 0, -1):         # keep the closest of a run
            run = dup[:, c]
            E_[run, c - 1] = np.minimum(E_[run, c - 1], E_[run, c])
        E_[dup] = np.inf
        e0 = E_.min(1)
        fin = np.isfinite(e0)
        out = np.full(len(X), np.inf)
        w = np.exp(-(E_[fin] - e0[fin, None]) / self.blend)
        out[fin] = e0[fin] - self.blend * np.log(w.sum(1))
        return out

    def value(self, X):
        val = np.minimum(self.axis_distance(X) - self.r, self.far)
        if self.trim is not None:
            near = val < self.far
            if self.surf_tree is not None:
                ds, _ = self.surf_tree.query(X[near],
                                             distance_upper_bound=2 * self.r + 4 * self.h,
                                             workers=-1)
                sub = np.nonzero(near)[0][np.isfinite(ds)]
            else:
                sub = np.nonzero(near)[0]
            if len(sub):
                # the contour is taken at +iso and then smoothed, which would
                # leave a flush cut ~0.02 mm proud: cut a little inside
                val[sub] = np.maximum(val[sub],
                                      self.trim(X[sub]) + self.iso + 0.15 * self.h)
        return val

    def active_blocks(self, origin, shape, Bk):
        """Indices of the Bk-voxel blocks of a grid that can hold solid."""
        bi, bj, bk = (np.arange(0, n, Bk) for n in shape)
        BI = np.stack(np.meshgrid(bi, bj, bk, indexing='ij'), -1).reshape(-1, 3)
        cen = origin + (BI + 0.5 * (Bk - 1)) * self.h
        dd, _ = self.tree.query(cen, distance_upper_bound=self.reach
                                + math.sqrt(3) * Bk * self.h / 2, workers=-1)
        return BI[np.isfinite(dd)]

    def solid_volume(self, n=800000, seed=0):
        """Monte Carlo volume of what the mesher will produce - the contour
        level +iso, sampled only in blocks near a strut so most samples count.
        Measured against the finished mesh: +1.8 %, +0.8 %, -0.1 % at organic 1,
        5, 10 (the zero level read 2-3 % low: in the concave crotches of a
        junction the contour sits outside the design surface). Two seeds agree
        to ~0.5 %."""
        Bk = 8
        shape = (np.ceil((self.hi - self.lo) / self.h).astype(int) + 1)
        act = self.active_blocks(self.lo, shape, Bk)
        if len(act) == 0:
            return 0.0
        rng = np.random.default_rng(seed)
        pick = act[rng.integers(0, len(act), n)]
        X = self.lo + (pick + rng.random((n, 3)) * Bk) * self.h
        inside = self.value(X) < self.iso
        return float(inside.mean() * len(act) * (Bk * self.h) ** 3)


def mesh_struts(A, B, d, h, trim=None, surf_pts=None, log=None, iso=None,
                blend=0.0, field=None):
    """The field of StrutField contoured, slab by slab, by flying edges."""
    import pyvista as pv
    F = field or StrutField(A, B, d, h, blend, trim, surf_pts, iso)
    lo, far = F.lo, F.far
    nx, ny, nz = (np.ceil((F.hi - lo) / h).astype(int) + 1)
    Bk = 6
    off = np.stack(np.meshgrid(np.arange(Bk), np.arange(Bk), np.arange(Bk),
                               indexing='ij'), -1).reshape(-1, 3)
    layers = max(Bk, int(P.MAX_VOXELS_SLAB // max(nx * ny, 1)) // Bk * Bk)
    pieces = []
    k0 = 0
    solid_vox = 0
    while k0 < nz - 1:
        k1 = min(nz - 1, k0 + layers)
        snz = k1 - k0 + 1
        fld = np.full((nx, ny, snz), far, np.float32)
        origin = lo + np.array([0, 0, k0]) * h
        act = F.active_blocks(origin, (nx, ny, snz), Bk)
        for c0 in range(0, len(act), 8000):
            idx = (act[c0:c0 + 8000, None, :] + off[None]).reshape(-1, 3)
            ok = (idx[:, 0] < nx) & (idx[:, 1] < ny) & (idx[:, 2] < snz)
            idx = idx[ok]
            fld[idx[:, 0], idx[:, 1], idx[:, 2]] = F.value(origin + idx * h)
        solid_vox += (int(np.sum(fld[:, :, :-1] < 0)) if k1 < nz - 1
                      else int(np.sum(fld < 0)))
        grid = pv.ImageData(dimensions=(nx, ny, snz), spacing=(h, h, h),
                            origin=tuple(origin))
        grid.point_data['f'] = fld.ravel(order='F')
        piece = grid.contour([F.iso], scalars='f', method='flying_edges',
                             compute_scalars=False)
        if piece.n_cells:
            pieces.append(piece)
        if log:
            log('Meshing the struts: %d %%' % (100 * k1 // max(nz - 1, 1)))
        k0 = k1
    if not pieces:
        raise LatticeError('the struts produced no surface')
    surf = pieces[0].merge(pieces[1:]) if len(pieces) > 1 else pieces[0]
    surf = surf.clean(tolerance=1e-6 * h, absolute=True, **TRIS_ONLY).triangulate()
    return surf, {'grid': [int(nx), int(ny), int(nz)], 'voxel_mm': h,
                  'voxel_volume_mm3': solid_vox * h ** 3}


# clean() by default turns a collapsed triangle into a line, and one line makes
# vtkMassProperties refuse the whole mesh
TRIS_ONLY = dict(lines_to_points=False, polys_to_lines=False)


def finish_mesh(surf, max_tris, smooth=10):
    import pyvista as pv
    raw = surf.n_cells
    if smooth:
        surf = surf.smooth_taubin(n_iter=int(smooth), pass_band=0.1)
    if surf.n_cells > max_tris:
        surf = surf.decimate_pro(1.0 - max_tris / surf.n_cells, preserve_topology=True,
                                 feature_angle=60.0, splitting=False,
                                 boundary_vertex_deletion=False)
    surf = surf.clean(**TRIS_ONLY).triangulate()
    conn = surf.connectivity('all')
    reg = np.asarray(conn.cell_data['RegionId'])
    sizes = np.bincount(reg)
    small = np.nonzero(sizes < max(30, 1e-4 * surf.n_cells))[0]
    if len(small):
        conn = conn.extract_cells(~np.isin(reg, small)).extract_surface(algorithm='dataset_surface')
        surf = pv.PolyData(conn).clean(**TRIS_ONLY).triangulate()
    edges = surf.extract_feature_edges(boundary_edges=True, non_manifold_edges=True,
                                       feature_edges=False, manifold_edges=False)
    return surf, {'raw_triangles': int(raw), 'triangles': int(surf.n_cells),
                  'open_or_nonmanifold_edges': int(edges.n_cells),
                  'components': int(len(sizes) - len(small)),
                  'specks_removed': int(len(small))}


def signed_volume(surf):
    pts = np.asarray(surf.points, float)
    f = np.asarray(surf.faces).reshape(-1, 4)[:, 1:]
    return float(np.einsum('ij,ij->i', pts[f[:, 0]],
                           np.cross(pts[f[:, 1]], pts[f[:, 2]])).sum() / 6.0)


def outward(surf):
    """Flying edges winds its triangles towards the LOWER field value - the
    inside here - so the raw lattice is inside-out: Fusion read -790 mm3.
    vtkMassProperties reports |V| and hides it; the signed sum does not."""
    if signed_volume(surf) < 0:
        f = np.asarray(surf.faces).reshape(-1, 4).copy()
        f[:, 1:] = f[:, [3, 2, 1]]
        surf.faces = f.ravel()
    return surf


def write_stl(surf, path):
    pts = np.asarray(surf.points, dtype=np.float32)
    faces = np.asarray(surf.faces).reshape(-1, 4)[:, 1:]
    tri = pts[faces]
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    n /= np.maximum(np.linalg.norm(n, axis=1), 1e-30)[:, None]
    rec = np.zeros(len(faces), dtype=[('n', '<f4', 3), ('v', '<f4', (3, 3)),
                                      ('a', '<u2')])
    rec['n'] = n
    rec['v'] = tri
    tmp = path + '.part'
    with open(tmp, 'wb') as fh:
        fh.write(b'Voronoi Lattice'.ljust(80, b' '))
        fh.write(np.uint32(len(faces)).tobytes())
        fh.write(rec.tobytes())
    os.replace(tmp, path)


# --------------------------------------------------------------------------- #
# job
# --------------------------------------------------------------------------- #

GRAPH_KEYS = ('kind', 'placement', 'flip', 'skin', 'border', 'define', 'porosity',
              'cell', 'thickness', 'regularity', 'seed', 'closed', 'organic',
              'output', 'resolution')

FUSION_FORM_S_PER_FACE = 0.0042    # measured: 2195 faces 9.3 s ... 10097 faces 42 s


def load_region(job):
    m = job['mesh']
    V = np.fromfile(m['coords'], dtype='<f8').reshape(-1, 3)
    F = np.fromfile(m['tris'], dtype='<i4').reshape(-1, 3).astype(np.int64)
    N = None
    if m.get('normals') and os.path.exists(m['normals']):
        nn = np.fromfile(m['normals'], dtype='<f8').reshape(-1, 3)
        if len(nn) == len(V):
            N = np.hstack([V, nn])
    return V, F, N


def graph_key(job):
    m = job['mesh']
    h = hashlib.sha1()
    for k in ('coords', 'tris'):
        with open(m[k], 'rb') as fh:
            h.update(fh.read())
    h.update(json.dumps([GRAPH_VERSION] + [job.get(k) for k in GRAPH_KEYS],
                        default=str).encode())
    return h.hexdigest()[:20]


def validate(job):
    d = float(job.get('thickness', 0))
    if d <= 0:
        raise LatticeError('strut thickness must be positive')
    if job['kind'] not in ('volume', 'surface'):
        raise LatticeError('choose a volume or a surface lattice')
    if job['kind'] == 'volume' and not job.get('closed'):
        raise LatticeError('a volume lattice needs a closed region')
    if job['placement'] not in P.PLACEMENTS:
        raise LatticeError('placement must be flush, middle or external')
    if job.get('output', 'mesh') not in ('form', 'mesh'):
        raise LatticeError('output must be form or mesh')
    if not 1.0 <= float(job.get('organic', 1)) <= 10.0:
        raise LatticeError('organic must be between 1 and 10')
    if job.get('define', 'porosity') == 'cell':
        if float(job['cell']) <= d:
            raise LatticeError('the cell size must be larger than the strut '
                               'thickness, or the pores close')
    elif not 0.0 < float(job['porosity']) < 1.0:
        raise LatticeError('porosity must be between 0 and 100 %')


def run(job, log=None):
    t0 = time.time()
    log = log or (lambda _m: None)
    validate(job)
    V, F, N = load_region(job)
    R = prepare(V, F, N, bool(job.get('closed')))
    timing = {'prepare_s': time.time() - t0}
    M = Measure(job, R)
    cache = None
    if job.get('cache_dir'):
        os.makedirs(job['cache_dir'], exist_ok=True)
        cache = os.path.join(job['cache_dir'], 'graph_%s.npz' % graph_key(job))
    if cache and os.path.exists(cache):
        z = np.load(cache, allow_pickle=False)
        A, B = z['A'], z['B']
        info = json.loads(str(z['info']))
        info['cached'] = True
    else:
        log('Tracing the cells')
        A, B, info = graph_for(job, R, log, M)
        if cache:
            np.savez(cache, A=A, B=B, info=json.dumps(info))
    timing['graph_s'] = time.time() - t0 - timing['prepare_s']
    d = float(job['thickness'])
    _, env = target_solid(job, R)
    output = info.get('output', 'mesh')
    out = {'region': {'closed': R['closed'], 'area_mm2': R['area'],
                      'volume_mm3': R['volume'], 'envelope': env},
           'graph': info, 'timing': timing, 'output': output,
           'porosity_predicted': info['porosity_predicted']}
    if 'cage' in info:
        out['fusion_form_s'] = info['cage']['cage_faces'] * FUSION_FORM_S_PER_FACE
    if job.get('mode') == 'graph':
        cap = int(job.get('max_preview_segments', 60000))
        k = np.arange(len(A))
        if len(k) > cap:
            k = np.linspace(0, len(A) - 1, cap).astype(int)
        out['segments'] = np.hstack([A[k], B[k]]).round(5).tolist()
        out['timing']['total_s'] = time.time() - t0
        return out

    t1 = time.time()
    if output == 'form':
        log('Building the Form cage')
        Vc, Fc, ci = FORM.build_cage(A, B, d, float(job.get('organic', 1)))
        Pl, T = FORM.limit_mesh(Vc, Fc, 3)
        vol = FORM.mesh_volume(Pl, T)
        ci['volume_mm3'] = vol
        if M.trim is not None:
            # a Form cannot be trimmed: report how far its hubs bulge past the
            # surface it is meant to sit flush inside
            ci['protrudes_mm'] = max(0.0, float(M.trim(Pl).max()))
        out['porosity_measured'] = porosity_from(job, R, vol)
        out['form'] = ci
        out['fusion_form_s'] = ci['cage_faces'] * FUSION_FORM_S_PER_FACE
        tmp = job['out_tsm'] + '.part'
        with open(tmp, 'w', encoding='utf-8') as fh:
            fh.write(FORM.write_tsm(Vc, Fc))
        os.replace(tmp, job['out_tsm'])
        out['tsm'] = job['out_tsm']
        timing['form_s'] = time.time() - t1
        timing['total_s'] = time.time() - t0
        return out

    mo = P.mesh_organic(float(job.get('organic', 1)), d)
    surf, mi = mesh_struts(A, B, d, M.h, log=log, field=M.field(A, B))
    timing['field_s'] = time.time() - t1
    t1 = time.time()
    log('Smoothing and decimating %d triangles' % surf.n_cells)
    surf, fi = finish_mesh(surf, int(job.get('max_tris', 1_500_000)),
                           job.get('smooth', mo['smooth']))
    timing['finish_s'] = time.time() - t1
    mi.update(fi)
    surf = outward(surf)
    vol = signed_volume(surf)
    mi['volume_mm3'] = vol
    mi['blend_mm'] = mo['blend']
    out['porosity_measured'] = porosity_from(job, R, vol)
    out['mesh'] = mi
    log('Writing %d triangles' % surf.n_cells)
    write_stl(surf, job['out_stl'])
    out['stl'] = job['out_stl']
    timing['total_s'] = time.time() - t0
    return out


def main(argv):
    jp, op = argv[1], argv[2]
    logp = op + '.log'

    def log(m):
        with open(logp, 'a', encoding='utf-8') as fh:
            fh.write('%8.2f  %s\n' % (time.time() - T0, m))

    try:
        with open(jp, encoding='utf-8') as fh:
            job = json.load(fh)
        res = {'ok': True, 'result': run(job, log)}
    except LatticeError as exc:
        res = {'ok': False, 'error': str(exc), 'user_error': True}
    except FORM.CageError as exc:
        res = {'ok': False, 'error': str(exc), 'user_error': True}
    except P.PlanError as exc:
        res = {'ok': False, 'error': str(exc), 'user_error': True}
    except Exception:
        res = {'ok': False, 'error': traceback.format_exc()}
    tmp = op + '.part'
    with open(tmp, 'w', encoding='utf-8') as fh:
        json.dump(res, fh)
    os.replace(tmp, op)
    return 0 if res['ok'] else 1


T0 = time.time()

if __name__ == '__main__':
    raise SystemExit(main(sys.argv))
