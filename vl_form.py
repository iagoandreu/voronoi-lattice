"""Voronoi Lattice - the strut graph as a T-Spline (Form) control cage.

Runs in the engine's Python (numpy, scipy). Millimetres in, TSM (centimetres) out.

Why a cage and not a converted mesh: Fusion's Convert Mesh > Organic needs the
Product Design Extension, which this licence lacks, and a lattice mesh has
hundreds of thousands of triangles no conversion could keep. A T-Spline body is
the Catmull-Clark LIMIT surface of a control cage, and a lattice has a natural
cage: a square tube along every strut, a hub at every junction.

    tube   four quads between two square rings of the same frame, so the tube
           does not twist; the limit of a square ring of half-width a is a
           near-circle of radius 0.917a-0.943a, so the ring is sized r / 0.93
    hub    the convex hull of the rings of the struts that meet there, minus the
           rings' own faces; triangle pairs that close a flat quad are merged.
           Every ring corner must be a hull vertex - a ring buried inside the
           hull (two struts at a small angle) is pushed out along its strut
           until it is not
    end    a free end (degree 1) is capped with the ring's own quad

ORGANIC (1..10) moves the hub rings: further from the node and flared, the limit
surface webs the struts together like trabecular bone; close to it and with a
second, tight ring, the struts meet nearly as cylinders do.

The cage is checked for being a closed 2-manifold and consistently wound before
anything is written: Fusion answers a broken cage with no body AND a half-built
form feature that blocks every later one.
"""

import math

import numpy as np
from scipy.spatial import ConvexHull, cKDTree

RING_LIMIT = 0.93          # limit radius of a square ring / its half-width
T_MAX = 0.45               # a hub ring may sit at most this far along its chain

TSM_HEADER = """#TS0200

degree 3
cap-type G1CAPS
star-smoothness 0
units 1 meters
end-conditions SUBD_CREASES
star-knot-rule NURCCS"""

TSM_FOOTER = """tol 1.00000000000000008e-05
geom-tol 1.00000000000000008e-05
ver 1044
behavior-version 14.2.0
compat-version 13.3.0"""


class CageError(RuntimeError):
    pass


def organic_params(level, d):
    """level 1..10 -> hub ring distance, hub flare, and whether to add a tight
    second ring. Distances in mm."""
    o = min(max(float(level), 1.0), 10.0)
    return {'t': (0.55 + 0.2 * (o - 1.0)) * d,
            'flare': 1.0 + 0.05 * (o - 1.0),
            'tight': o <= 3.0,
            'merge': 1.5 * d}


# --------------------------------------------------------------------------- #
# the graph
# --------------------------------------------------------------------------- #

def graph_from_segments(A, B):
    P = np.vstack([A, B])
    key = np.round(P * 1e5).astype(np.int64)
    _, first, inv = np.unique(key, axis=0, return_index=True, return_inverse=True)
    X = P[first]
    n = len(A)
    E = np.stack([inv[:n], inv[n:]], 1)
    E = E[E[:, 0] != E[:, 1]]
    E = np.unique(np.sort(E, 1), axis=0)
    return X, E


def chain_graph(X, E):
    """The polylines between junctions (degree != 2). A curved strut is ONE
    chain - a tube that bends - not a string of hubs: on a surface net four
    nodes in five are only polyline points following the curvature. A loop with
    no junction on it comes out with closed=True."""
    adj = [[] for _ in range(len(X))]
    for k, (a, b) in enumerate(E):
        adj[a].append((b, k))
        adj[b].append((a, k))
    deg = np.array([len(x) for x in adj])
    used = np.zeros(len(E), bool)
    chains = []

    def walk(start, nb, k):
        used[k] = True
        path, cur = [start, nb], nb
        while deg[cur] == 2 and cur != start:
            (n1, k1), (n2, k2) = adj[cur]
            nxt, kk = (n1, k1) if not used[k1] else (n2, k2)
            if used[kk]:
                break
            used[kk] = True
            path.append(nxt)
            cur = nxt
        return path

    for v in range(len(X)):
        if deg[v] in (0, 2):
            continue
        for nb, k in adj[v]:
            if not used[k]:
                chains.append({'path': walk(v, nb, k), 'closed': False})
    for k in range(len(E)):
        if not used[k]:
            a, b = E[k]
            path = walk(a, b, k)
            closed = path[-1] == path[0]
            chains.append({'path': path[:-1] if closed else path, 'closed': closed})
    return chains


def arclength(P):
    return np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))])


def merge_junctions(X, chains, L_min):
    """Collapse junctions joined by a chain shorter than L_min into one hub at
    their mean: a Voronoi diagram is full of short edges, and a hub needs room
    along each chain for its ring. ONE pass, over connected components of the
    short chains: repeating it moved each hub, made new short chains, and
    cascaded into 33-strut clusters."""
    for _ in range(1):
        parent = {}

        def find(i):
            while parent.get(i, i) != i:
                i = parent[i]
            return i

        hit = False
        for c in chains:
            if c['closed']:
                continue
            a, b = c['path'][0], c['path'][-1]
            if a != b and arclength(X[c['path']])[-1] < L_min:
                ra, rb = find(a), find(b)
                if ra != rb:
                    parent[rb] = ra
                    hit = True
        if not hit:
            break
        groups = {}
        for c in chains:
            if c['closed']:
                continue
            for v in (c['path'][0], c['path'][-1]):
                groups.setdefault(find(v), set()).add(v)
        X = X.copy()
        for root, members in groups.items():
            if len(members) > 1:
                X[root] = X[sorted(members)].mean(0)
        out = []
        for c in chains:
            if c['closed']:
                out.append(c)
                continue
            path = list(c['path'])
            path[0], path[-1] = find(path[0]), find(path[-1])
            if path[0] == path[-1] and arclength(X[path])[-1] < 3.0 * L_min:
                continue                           # a short loop on one hub
            out.append({'path': path, 'closed': False})
        chains = dedupe(X, out)
    return X, chains


def transport(T):
    """Rotation-minimising frames along a polyline, one per segment: (p, q)
    carried from each tangent to the next by the smallest rotation, so a tube
    that follows a curve does not twist."""
    p, q = frame(T[0])
    out = [(p, q)]
    for i in range(1, len(T)):
        t0, t1 = T[i - 1], T[i]
        v = np.cross(t0, t1)
        c = float(np.dot(t0, t1))
        if c > -0.999:
            p = p * c + np.cross(v, p) + v * np.dot(v, p) / (1.0 + c)
        p = p - np.dot(p, t1) * t1
        p /= np.linalg.norm(p)
        out.append((p, np.cross(t1, p)))
    return out


class Chain:
    """One tube: polyline P, its arc length and a frame per segment."""

    def __init__(self, P, closed=False):
        if closed:
            P = np.vstack([P, P[:1]])
        self.P = P
        self.s = arclength(P)
        self.S = float(self.s[-1])
        D = np.diff(P, axis=0)
        self.T = D / np.linalg.norm(D, axis=1)[:, None]
        self.F = transport(self.T)

    def ring(self, s, half):
        """Ring corners at arc length s, and the tangent there."""
        s = min(max(s, 0.0), self.S)
        k = int(np.searchsorted(self.s, s, side='right')) - 1
        k = min(max(k, 0), len(self.T) - 1)
        u = (s - self.s[k]) / max(self.s[k + 1] - self.s[k], 1e-30)
        c = self.P[k] + u * (self.P[k + 1] - self.P[k])
        p, q = self.F[k]
        return ring(c, self.T[k], p, q, half), self.T[k]


def connect(V, faces, c0, c1):
    """Four quads from ring c0 to ring c1, corners matched by the rotation
    that pairs them closest - a hub ring and a transported frame need not
    line up corner for corner."""
    best, off = None, 0
    for o in range(4):
        dd = sum(np.linalg.norm(V[c0 + i] - V[c1 + (i + o) % 4]) for i in range(4))
        if best is None or dd < best:
            best, off = dd, o
    for i in range(4):
        j = (i + 1) % 4
        faces.append((c0 + i, c0 + j, c1 + (j + off) % 4, c1 + (i + off) % 4))


def dedupe(X, chains, max_angle_deg=15.0):
    """After a merge two chains can join the same pair of junctions and leave
    them in the same direction - parallel tubes no hub can separate (measured:
    a min angle of 0.0 at every failed hub). Keep the shorter of such twins."""
    cosmax = math.cos(math.radians(max_angle_deg))
    by_pair = {}
    for i, c in enumerate(chains):
        if c['closed']:
            continue
        a, b = c['path'][0], c['path'][-1]
        by_pair.setdefault((min(a, b), max(a, b)), []).append(i)
    drop = set()
    for (a, b), ids in by_pair.items():
        if len(ids) < 2:
            continue

        def leave(i, node):
            p = chains[i]['path']
            q = p if p[0] == node else p[::-1]
            e = X[q[1]] - X[q[0]]
            return e / (np.linalg.norm(e) + 1e-30)

        def length(i):
            return float(arclength(X[chains[i]['path']])[-1])
        ids = sorted(ids, key=length)
        for x in range(len(ids)):
            for y in range(x):
                i, j = ids[x], ids[y]
                if j in drop or i in drop:
                    continue
                if (np.dot(leave(i, a), leave(j, a)) > cosmax
                        or np.dot(leave(i, b), leave(j, b)) > cosmax):
                    drop.add(i)
    return [c for i, c in enumerate(chains) if i not in drop]


def merge_pairs(X, chains, pairs):
    """Fold each (v, w) junction pair into one hub at their mean."""
    parent = {}

    def find(i):
        while parent.get(i, i) != i:
            i = parent[i]
        return i
    for v, w in pairs:
        rv, rw = find(v), find(w)
        if rv != rw:
            parent[rw] = rv
    groups = {}
    for v in list(parent) + [p for p in parent.values()]:
        groups.setdefault(find(v), set()).add(v)
    X = X.copy()
    for root, members in groups.items():
        X[root] = X[sorted(members)].mean(0)
    out = []
    for c in chains:
        if c['closed']:
            out.append(c)
            continue
        path = list(c['path'])
        path[0], path[-1] = find(path[0]), find(path[-1])
        if path[0] == path[-1] and len(path) <= 3:
            continue
        out.append({'path': path, 'closed': False})
    return X, dedupe(X, out)


def build_cage(A, B, d, organic, log=None):
    """Points (mm) and faces of the control cage for struts A[i]-B[i].

    A junction whose rings cannot all be placed - struts too short, or meeting
    at too small an angle - is folded into its nearest neighbour across its
    shortest chain, and the cage tried again: merging only where a hub fails
    keeps junctions small, where merging every short edge up front cascaded
    into 12-strut clusters that failed themselves."""
    op = organic_params(organic, d)
    X, E = graph_from_segments(A, B)
    X, raw = merge_junctions(X, chain_graph(X, E), op['merge'])
    folded = 0
    for _round in range(13):
        chains, hubs, failed, t_end, pushed, inc, shared, blame = _hubs(X, raw, d, op)
        if not failed:
            break
        pairs = []
        for v in failed:
            best = None
            for ci, end in inc[v]:
                if ci not in blame.get(v, ()):
                    continue
                a, b, closed, ch = chains[ci]
                w = b if end == 0 else a
                if w != v and (best is None or ch.S < best[0]):
                    best = (ch.S, w)
            if best is not None:
                pairs.append((v, best[1]))
        if _round >= 12 or not pairs:
            break
        folded += len(pairs)
        X, raw = merge_pairs(X, raw, pairs)
    if failed:
        # never drop a strut to get a cage through: tried, and it cascaded until
        # a 1760 mm3 lattice came out as 535 mm3 of something else
        raise CageError('%d junctions cannot be built as a Form even after '
                        'folding %d into their neighbours: struts that short '
                        'and that acute overlap for longer than a tube can '
                        'separate them' % (len(failed), folded))
    V, faces, info = _assemble(chains, hubs, t_end, inc, d, op, pushed, folded,
                               organic, shared)
    return V, faces, info


def _hubs(X, raw, d, op):
    chains = []
    for c in raw:
        P = X[c['path']]
        keep = np.concatenate([[True], np.linalg.norm(np.diff(P, axis=0), axis=1) > 1e-6])
        P = P[keep]
        if len(P) < 2 or (c['closed'] and len(P) < 3):
            continue
        chains.append((c['path'][0], c['path'][-1], c['closed'], Chain(P, c['closed'])))
    if not chains:
        raise CageError('no struts left after merging the short ones')
    half = 0.5 * d / RING_LIMIT
    flare = half * op['flare']

    inc = {}                                  # junction -> [(chain, end)]
    for ci, (a, b, closed, ch) in enumerate(chains):
        if not closed:
            inc.setdefault(a, []).append((ci, 0))
            inc.setdefault(b, []).append((ci, 1))
    t_end = np.full((len(chains), 2), op['t'])
    # a chain too short to hold a ring at each end AND a tube between them gets
    # ONE ring at its middle, shared by both hubs - its edges then border a hub
    # face on either side and no tube is built
    gap = 0.45 * d
    shared = set()
    for ci, (a, b, closed, ch) in enumerate(chains):
        if (not closed and a != b and len(inc[a]) >= 2 and len(inc[b]) >= 2
                and ch.S < 2.0 * op['t'] + gap):
            shared.add(ci)
            t_end[ci] = 0.5 * ch.S

    def end_ring(ci, end, size):
        ch = chains[ci][3]
        if ci in shared:
            pts, tan = ch.ring(0.5 * ch.S, half)
        else:
            t = min(t_end[ci, end], T_MAX * ch.S)
            pts, tan = ch.ring(t if end == 0 else ch.S - t, size)
        return pts, (tan if end == 0 else -tan)

    pushed = 0
    hubs = {}
    failed = []
    blame = {}                # failed junction -> chains whose rings broke it
    for v, ends in inc.items():
        if len(ends) < 2:
            continue
        scale = 1.0
        for _attempt in range(40):
            pts, ring_of, ed = [], [], []
            size = max(flare * scale, 0.7 * half)
            for j, (ci, end) in enumerate(ends):
                rp, e = end_ring(ci, end, size)
                pts.append(rp)
                ed.append(e)
                ring_of += [j] * 4
            pts, ring_of, ed = np.vstack(pts), np.array(ring_of), np.array(ed)
            # a ring's square is a hull face - so its edges exist for the tube
            # to attach to - only if the whole hub lies behind its plane
            bad = []
            for j in range(len(ends)):
                c = pts[4 * j:4 * j + 4].mean(0)
                if float(((pts[ring_of != j] - c) @ ed[j]).max()) > -1e-3 * d:
                    bad.append(j)
            faces = None
            if not bad:
                try:
                    faces = hub_faces(pts, ring_of, pts.mean(0))
                except Exception:
                    faces = None
                if faces is None:
                    cosm = ed @ ed.T
                    np.fill_diagonal(cosm, -2)
                    bad = [int(np.argmax(cosm.max(1)))]
            if faces is not None:
                hubs[v] = (pts, faces)
                break
            # push out every ring that breaks the rule, not one per attempt: a
            # six-strut hub with fine angles ran out of attempts the old way
            stuck = True
            for j in bad:
                ci, end = ends[j]
                if ci not in shared and t_end[ci, end] < T_MAX * chains[ci][3].S:
                    stuck = False
                    t_end[ci, end] *= 1.25
                    pushed += 1
            if stuck:
                # nothing left to push: make this hub's rings smaller
                if size <= 0.7 * half + 1e-12:
                    failed.append(v)
                    blame[v] = [ends[j][0] for j in bad]
                    break
                scale *= 0.85
        else:
            failed.append(v)
            blame[v] = [ci for ci, _e in ends]
    return chains, hubs, failed, t_end, pushed, inc, shared, blame


def _assemble(chains, hubs, t_end, inc, d, op, pushed, folded, organic, shared):
    half = 0.5 * d / RING_LIMIT

    def t_of(ci, end):
        return min(t_end[ci, end], T_MAX * chains[ci][3].S)

    verts, faces = [], []
    ring_id = {}
    for v, (pts, hf) in hubs.items():
        local = []
        for j, (ci, end) in enumerate(inc[v]):
            other = (ci, 1 - end)
            if ci in shared and other in ring_id:
                ring_id[(ci, end)] = ring_id[other]       # the same ring
            else:
                ring_id[(ci, end)] = len(verts)
                verts.extend(pts[4 * j:4 * j + 4])
            local += [ring_id[(ci, end)] + i for i in range(4)]
        for f in hf:
            faces.append(tuple(local[i] for i in f))

    gap = 0.45 * d
    for ci, (a, b, closed, ch) in enumerate(chains):
        if ci in shared:
            continue
        if closed:
            seq = []
            for s in ch.s[:-1]:
                seq.append(len(verts))
                verts.extend(ch.ring(s, half)[0])
            seq.append(seq[0])
        else:
            for end, s_at in ((0, 0.0), (1, ch.S)):
                if (ci, end) not in ring_id:          # a free end: cap it
                    ring_id[(ci, end)] = len(verts)
                    verts.extend(ch.ring(s_at, half)[0])
                    faces.append(tuple(ring_id[(ci, end)] + i for i in range(4)))
            ta = t_of(ci, 0) if (ci, 0) in inc.get(a, []) and a in hubs else 0.0
            tb = t_of(ci, 1) if (ci, 1) in inc.get(b, []) and b in hubs else 0.0
            stations = list(ch.s[1:-1])
            if op['tight']:
                stations += [ta + 0.3 * d, ch.S - tb - 0.3 * d]
            stations = sorted(x for x in stations if ta + gap < x < ch.S - tb - gap)
            thin = []
            for x in stations:
                if not thin or x - thin[-1] > gap:
                    thin.append(x)
            seq = [ring_id[(ci, 0)]]
            for x in thin:
                seq.append(len(verts))
                verts.extend(ch.ring(x, half)[0])
            seq.append(ring_id[(ci, 1)])
        for c0, c1 in zip(seq[:-1], seq[1:]):
            connect(verts, faces, c0, c1)

    V = np.asarray(verts, float)
    faces = orient(V, faces)
    check_manifold(len(V), faces)
    info = {'junctions': len(hubs), 'chains': len(chains), 'rings_pushed': pushed,
            'cage_vertices': int(len(V)), 'cage_faces': len(faces),
            'cage_quads': sum(1 for f in faces if len(f) == 4),
            'cage_triangles': sum(1 for f in faces if len(f) == 3),
            'folded': folded, 'shared_rings': len(shared),
            'organic': float(organic)}
    return V, faces, info


# rings and hubs
# --------------------------------------------------------------------------- #

def frame(e):
    ref = np.eye(3)[int(np.argmin(np.abs(e)))]
    p = np.cross(e, ref)
    p /= np.linalg.norm(p)
    return p, np.cross(e, p)


def ring(c, e, p, q, half):
    """Four corners of a square ring centred at c, normal e, half-width half,
    counter-clockwise seen from +e."""
    return np.array([c + half * (sx * p + sy * q)
                     for sx, sy in ((1, 1), (-1, 1), (-1, -1), (1, -1))])


def hub_faces(pts, ring_of, center):
    """Hull of the hub's ring corners, minus each ring's own cap; triangles
    oriented outward. Returns faces in local indices, or None when a ring is
    buried (a corner not on the hull)."""
    hull = ConvexHull(pts)
    if len(set(hull.vertices.tolist())) < len(pts):
        return None
    tris = []
    for s in hull.simplices:
        rs = ring_of[s]
        if rs[0] == rs[1] == rs[2]:
            continue
        a, b, c = pts[s]
        n = np.cross(b - a, c - a)
        if np.dot(n, (a + b + c) / 3.0 - center) < 0:
            s = s[::-1]
        tris.append(tuple(int(i) for i in s))
    return merge_quads(pts, tris, ring_of)


def merge_quads(pts, tris, ring_of, max_dihedral_deg=8.0):
    """Pair triangles sharing an edge into a quad when they are nearly coplanar
    and the edge is not a ring edge: fewer star points, smoother limit."""
    cosmin = math.cos(math.radians(max_dihedral_deg))
    nrm = []
    for t in tris:
        a, b, c = pts[list(t)]
        n = np.cross(b - a, c - a)
        nrm.append(n / (np.linalg.norm(n) + 1e-30))
    edge = {}
    for i, t in enumerate(tris):
        for k in range(3):
            u, v = t[k], t[(k + 1) % 3]
            edge[(u, v)] = i
    used = [False] * len(tris)
    out = []
    for i, t in enumerate(tris):
        if used[i]:
            continue
        done = False
        for k in range(3):
            u, v, w = t[k], t[(k + 1) % 3], t[(k + 2) % 3]
            j = edge.get((v, u))
            if j is None or used[j] or j == i:
                continue
            if ring_of[u] == ring_of[v] and abs(u - v) in (1, 3):
                continue                       # a ring edge: keep it a boundary
            if np.dot(nrm[i], nrm[j]) < cosmin:
                continue
            x = [p for p in tris[j] if p not in (u, v)][0]
            out.append((w, u, x, v))
            used[i] = used[j] = True
            done = True
            break
        if not done:
            used[i] = True
            out.append(t)
    return out


# checks, limit surface, TSM
# --------------------------------------------------------------------------- #

def check_manifold(n_verts, faces):
    use = {}
    for f in faces:
        n = len(f)
        for k in range(n):
            u, v = f[k], f[(k + 1) % n]
            use[(u, v)] = use.get((u, v), 0) + 1
    bad = [e for e, c in use.items() if c != 1 or use.get((e[1], e[0]), 0) != 1]
    if bad:
        raise CageError('the cage is not a closed, consistently wound surface '
                        '(%d bad half-edges)' % len(bad))
    seen = np.zeros(n_verts, bool)
    for f in faces:
        seen[list(f)] = True
    if not seen.all():
        raise CageError('%d cage vertices are on no face' % int((~seen).sum()))


def orient(V, faces):
    """Consistent winding by BFS over shared edges, then outward by the sign of
    the enclosed volume of each connected piece."""
    half = {}
    for fi, f in enumerate(faces):
        n = len(f)
        for k in range(n):
            half.setdefault((f[k], f[(k + 1) % n]), []).append(fi)
    adj = [[] for _ in faces]
    for (a, b), fs in half.items():
        for fi in fs:
            for fj in half.get((b, a), ()):
                adj[fi].append((fj, False))
            for fj in fs:
                if fj != fi:
                    adj[fi].append((fj, True))
    flip = [False] * len(faces)
    comp = [-1] * len(faces)
    nc = 0
    for root in range(len(faces)):
        if comp[root] >= 0:
            continue
        comp[root] = nc
        stack = [root]
        while stack:
            fi = stack.pop()
            for fj, must in adj[fi]:
                if comp[fj] < 0:
                    comp[fj] = nc
                    flip[fj] = flip[fi] ^ must
                    stack.append(fj)
        nc += 1
    faces = [tuple(f)[::-1] if flip[i] else tuple(f) for i, f in enumerate(faces)]
    vol = np.zeros(nc)
    for i, f in enumerate(faces):
        p0 = V[f[0]]
        for k in range(1, len(f) - 1):
            vol[comp[i]] += np.dot(p0, np.cross(V[f[k]], V[f[k + 1]])) / 6.0
    return [f[::-1] if vol[comp[i]] < 0 else f for i, f in enumerate(faces)]


def subdivide(points, faces):
    """One Catmull-Clark step on a closed mesh of any face arity (vectorised
    enough for 10^5 faces). Returns new points and all-quad faces."""
    P = np.asarray(points, float)
    nv = len(P)
    nf = len(faces)
    fp = np.array([P[list(f)].mean(0) for f in faces])
    edge_id = {}
    e_faces = []
    e_verts = []
    for fi, f in enumerate(faces):
        n = len(f)
        for k in range(n):
            u, v = f[k], f[(k + 1) % n]
            key = (u, v) if u < v else (v, u)
            if key not in edge_id:
                edge_id[key] = len(e_verts)
                e_verts.append(key)
                e_faces.append([fi])
            else:
                e_faces[edge_id[key]].append(fi)
    ev = np.array(e_verts)
    ef = np.array([x[:2] for x in e_faces])
    ep = (P[ev[:, 0]] + P[ev[:, 1]] + fp[ef[:, 0]] + fp[ef[:, 1]]) / 4.0
    val = np.bincount(ev.ravel(), minlength=nv).astype(float)
    Fsum = np.zeros_like(P)
    Fcnt = np.zeros(nv)
    for fi, f in enumerate(faces):
        for i in f:
            Fsum[i] += fp[fi]
            Fcnt[i] += 1
    Rsum = np.zeros_like(P)
    mid = 0.5 * (P[ev[:, 0]] + P[ev[:, 1]])
    np.add.at(Rsum, ev[:, 0], mid)
    np.add.at(Rsum, ev[:, 1], mid)
    n = val[:, None]
    Pn = (Fsum / Fcnt[:, None] + 2.0 * Rsum / n + (n - 3.0) * P) / n
    newP = np.vstack([Pn, fp, ep])
    out = []
    for fi, f in enumerate(faces):
        m = len(f)
        for k in range(m):
            u = f[k]
            a = f[(k + 1) % m]
            b = f[k - 1]
            ka = edge_id[(u, a) if u < a else (a, u)]
            kb = edge_id[(b, u) if b < u else (u, b)]
            out.append((u, nv + nf + ka, nv + fi, nv + nf + kb))
    return newP, out


def limit_mesh(V, faces, levels=3):
    P, F = V, faces
    for _ in range(levels):
        P, F = subdivide(P, F)
    F = np.asarray(F)
    T = np.vstack([F[:, [0, 1, 2]], F[:, [0, 2, 3]]])
    return P, T


def mesh_volume(P, T):
    return float(np.einsum('ij,ij->i', P[T[:, 0]],
                           np.cross(P[T[:, 1]], P[T[:, 2]])).sum() / 6.0)


def write_tsm(points_mm, faces):
    """TSM text (centimetres). Faces counter-clockwise from outside; TSM wants
    clockwise, so they are reversed here. Format: see MeshToPart docs/TSM-FORMAT.md."""
    faces = [tuple(f)[::-1] for f in faces]
    edge_id = {}
    for q in faces:
        n = len(q)
        for k in range(n):
            key = (q[k], q[(k + 1) % n])
            key = key if key[0] < key[1] else (key[1], key[0])
            if key not in edge_id:
                edge_id[key] = len(edge_id)
    start, n_links = [], 0
    for q in faces:
        start.append(n_links)
        n_links += len(q)
    nxt = [0] * n_links
    prv = [0] * n_links
    mate = [-1] * n_links
    dest = [0] * n_links
    fidx = [0] * n_links
    eidx = [0] * n_links
    by_edge = {}
    for f, q in enumerate(faces):
        n = len(q)
        s = start[f]
        for k in range(n):
            i = s + k
            nxt[i] = s + (k + 1) % n
            prv[i] = s + (k - 1) % n
            dest[i] = q[(k + 1) % n]
            fidx[i] = f
            key = (q[k], q[(k + 1) % n])
            key = key if key[0] < key[1] else (key[1], key[0])
            eidx[i] = edge_id[key]
            by_edge.setdefault(edge_id[key], []).append(i)
    for e, ls in by_edge.items():
        if len(ls) != 2:
            raise CageError('edge %d has %d links' % (e, len(ls)))
        mate[ls[0]], mate[ls[1]] = ls[1], ls[0]
    vert_link = {}
    for i in range(n_links):
        vert_link.setdefault(dest[i], i)
    out = [TSM_HEADER]
    out += ['f %d 0' % i for i in start]
    out += ['e %d 1' % by_edge[e][0] for e in range(len(edge_id))]
    out += ['v %d NORTH' % vert_link[v] for v in range(len(points_mm))]
    out += ['l %d %d %d %d %d %d 0' % (nxt[i], prv[i], mate[i], dest[i], fidx[i],
                                       eidx[i]) for i in range(n_links)]
    out.append('')
    out += ['0g %.12g %.12g %.12g 1' % (p[0] * 0.1, p[1] * 0.1, p[2] * 0.1)
            for p in points_mm]
    out.append('')
    out.append(TSM_FOOTER)
    return '\n'.join(out) + '\n'
