# Voronoi Lattice

A Fusion 360 add-in that fills a volume with a Voronoi lattice, or covers a surface
or a set of faces with one. It follows the 3-matic lattice workflow: pick the
region, say whether the lattice goes inside or on the skin, where the struts sit
relative to the surface, then the porosity and strut thickness. The result is a
**Form (T-spline) body**, or a watertight mesh body where a Form is not
possible. Its porosity is **measured** on what was built and reported, not
assumed from the target.

Command: **Voronoi Lattice**, in *Solid > Create* and *Mesh > Modify*.

## The dialog

| Input | What it does |
|---|---|
| **Region** | A solid, surface or mesh body, or one or more faces. Clicking a solid picks a face: pick the body in the browser, or tick **Whole body**. |
| **Lattice** | Shown only when the region is **closed**, and reset to *choose* every time a closed region is picked, so the choice is always made on purpose. *Volume* fills the region with a 3D Voronoi foam. *Surface* puts a Voronoi net on its skin only. An open region is always *Surface*. |
| **Placement** | *Flush*: the struts sit inside, their outer side on the surface, and nothing pokes through (the mesh is also trimmed to the surface). *Middle*: struts centred on the surface. *External*: struts outside, their inner side on the surface. On an open surface, "inside" is the back of the face; **Flip side** swaps it. |
| **Skin struts** | Volume lattice: also run struts along the surface, closing the cut cells, so no strut ends loose on the skin. |
| **Border strut** | Open surface: run a strut along its edges. |
| **Define by** | *Porosity* (the cell size is solved for it) or *Cell size* (the porosity follows, and is shown). |
| **Strut thickness** | Strut diameter. |
| **Regularity** | Lloyd relaxation steps on the seeds: 0 is raw random Voronoi, 5 is even but irregular, 20 is close to equal cells. |
| **Random seed** | Another seed gives another lattice with the same statistics. |
| **Shape > Organic** | 1 to 10: how the struts meet. 1 is close to crisp cylinders; 10 webs them together like trabecular bone. The cell size is re-solved so the porosity still holds. |
| **Shape > Output** | *Form*: a real T-spline body, smooth and editable in the Form workspace. *Mesh*: a watertight triangle mesh, any density. A Form that cannot be built becomes a mesh, and the preview and the report both say why. |
| **Mesh > Resolution** | Voxels per strut diameter: 3, 4 or 6. |
| **Mesh > Max triangles** | The mesh is decimated (topology-preserving) down to this. |
| **Region remeshing (advanced) > Remesh the region** | Automatic (default): when the region's triangle mesh has an edge over 5 mm (or over twice the target) or more than 5 % of its area in sliver triangles - what a CAD body tessellates into - it is remeshed isotropically (MeshLab, feature edges over 30 deg kept), and the result is accepted only if it stays closed with the volume within 0.5 %; otherwise it is subdivided exactly. Always / Never. The report says what was done. |
| **Region remeshing > Target edge** | 0 = half a cell, never over 5 mm (raised only past 600 k triangles). |
| **Preview struts** | Draws the strut centrelines live, in orange. |

The panel underneath shows, live, the cell size, the pore size, the cell count, the
strut length and the expected mesh size, and warns when the pores come out
narrower than the struts.

### What "porosity" means

- **Volume lattice**: the void fraction of the lattice envelope. For *flush* that is
  the region itself. For *middle* and *external* it is the region grown by the
  amount the struts stand proud (first order: V + A x offset).
- **Surface lattice**: the open fraction of the surface, i.e. the part a light shone
  along the normal would get through. Computed as 1 - 4 V / (pi d A).

## How it works

All geometry runs in a separate Python (`vl_engine.py`) with numpy, scipy, trimesh
and pyvista. Fusion's own Python has none of them.

- **Volume**: Lloyd-relaxed seeds inside the region, offset by the placement.
  The 3D Voronoi diagram (scipy, bounded by a ghost shell) is clipped to the region
  by bisection on a VTK signed distance. Skin struts are the same seeds' Voronoi
  diagram *restricted to the surface*, so each clipped edge lands on a skin
  junction.
- **Surface**: seeds on the surface (area-weighted samples plus Lloyd). The
  surface is subdivided finer than a cell and a strut, every vertex is labelled
  with its nearest seed, and each triangle edge whose two ends disagree is cut
  exactly where it crosses the two seeds' bisector. The cut is exact because
  |x - si|^2 - |x - sj|^2 is linear in x. Triangles with three labels hold a
  junction, solved in their plane. Only seeds facing the same way compete, so a
  thin wall's two sides keep their own cells. The label boundaries are simplified
  by Douglas-Peucker (tolerance 0.1 r).
- **Placement** offsets every strut-centre point along the outward normal by
  -r, 0 or +r.
- **Mesh**: the exact distance to the strut axes, minus r, sampled on a voxel grid
  slab by slab. A KD-tree over axis samples only proposes candidate segments. The
  field is contoured by VTK flying edges, Taubin-smoothed, and decimated. Flying
  edges cuts chords inside a convex surface (a strut came out 5 % light at
  4 voxels per diameter), so it contours at +0.11 h^2/r. That puts strut volume
  within +-1 % at every resolution.
- **Porosity to cell size**: the model is V = pi r^2 L - beta r^3 J (J = junctions;
  beta = 5.1 for 3D nodes, 2.0 for surface nodes, measured by
  `tests/calibrate.py`). It gives a first cell size; then the cell size is
  corrected against the L and J of the graph actually built, until the
  prediction is within 1 %. The preview's graph is cached, so OK only meshes.

## Form output

This licence has no Product Design Extension, so Fusion's own Convert Mesh >
Organic is not available. A lattice has a natural control cage anyway, and the
body is built from it through `TSplineBodies.addByTSMDescription` (the TSM format
was decoded for Mesh to Part; see its `docs/TSM-FORMAT.md`). `vl_form.py` builds
the cage:

- **Tube**: four quads between square rings, with frames carried along a curved
  strut by parallel transport, so a strut that follows a surface does not twist.
  The limit of a square ring of half-width a is a near-circle of radius
  0.917a to 0.943a, so the ring is sized r / 0.93.
- **Hub**: at a junction, the convex hull of the rings of the struts that meet
  there, minus each ring's own square; coplanar triangle pairs become quads. A
  ring's square is a hull face only if the whole hub lies behind its plane, so
  rings are pushed out along their struts until that holds. A strut too short
  for two rings gets one shared ring at its middle.
- **Organic** moves the hub rings further out and flares them: further out means
  broader webs.

The porosity loop measures the Catmull-Clark limit surface of the actual cage
(two subdivision levels while iterating, three for the reported value). Measured
in Fusion, its volume matches the engine to 0.5 %.

**When a Form is not possible.** A Voronoi *volume* lattice of typical density
(cells under about 8 strut diameters) has short edges that merge into junctions
where two struts leave at 30-40 degrees. Two 0.6 mm struts at that angle overlap
for about 1.7 mm, longer than the struts themselves, so no tube-per-strut cage
can separate them. The builder folds such a junction into its neighbour and
tries again, but it never drops a strut to get through. A cage it cannot finish
becomes a mesh, with the reason given. Surface lattices, and regular lattices
like diamond, build without a single fold. A direct-modelling design also gets a
mesh: a Form needs the timeline.

**Size.** Fusion takes about 4.2 ms per cage face: 2 195 faces build in 9 s and
10 097 faces in 42 s. The preview shows the face count and the time to expect.

## Organic, mesh output

For the mesh output, organic is an exponential smooth union of the strut
distances (radius 0.06 (o - 1) d), with each strut counted once, plus more
Taubin smoothing. The porosity loop measures the volume of that very field by
Monte Carlo, in the blocks near a strut (within 1.8 % of the finished mesh at
organic 1, 0.1 % at 10). A fitted node constant cannot do this job: the overlap
term goes from +5 to -46 across organic levels and depends on cell and strut
size.

## Install

```
powershell -ExecutionPolicy Bypass -File install.ps1
```

Then in Fusion: *Utilities > Add-Ins*, select **VoronoiLattice**, **Run**. The
engine needs a Python with numpy, scipy, trimesh and pyvista. It is found
automatically (`py -3.10` here), or set `VORONOILATTICE_PYTHON`, or put the
interpreter's path in `python_path.txt` next to the add-in.

`VoronoiLattice.py` is a thin shell. `vl_core.py` and `vl_plan.py` hot-reload on
every command event, and `vl_engine.py` runs fresh per build, so only shell
changes need a Stop/Run.

## Tests

- `py -3.10 tests/test_engine.py` runs offline checks with known answers: exact
  capsule volumes, porosity reached within 1.5 points, flush never crosses the
  surface, each placement on its side, flip, normals winning over winding, no
  strut through a thin wall, refusals, and the graph cache.
- `py -3.10 tools/livetest.py` runs live checks through the dev bridge
  (`devbridge.on`): closed or open detection per selection kind, the mesh body in
  the right component and base feature, units, placement against the source,
  replace, Form bodies (Fusion's volume against the engine's limit volume),
  the fallbacks to mesh, and organic. Every job runs in a throwaway document.
- `py -3.10 tests/calibrate_organic.py` shows why organic is measured, not
  modelled.
- `py -3.10 tests/calibrate.py` re-measures the constants in `vl_plan.py`.

## Limits

- A Form at flush placement cannot be trimmed. The report gives how far its hubs
  bulge past the surface.

- At a sharp edge, offsetting along averaged vertex normals pulls the skin struts
  in by up to (1 - cos 45 deg) r. Flush is trimmed, so nothing pokes out, but
  middle and external sit slightly low along sharp edges.
- Labelling by straight-line distance assumes cells smaller than the surface's
  radius of curvature.
- The output is a mesh body. Converting a lattice of thousands of struts to a
  B-Rep is not practical in Fusion.
