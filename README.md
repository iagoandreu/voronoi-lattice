# Lattice

A Fusion 360 add-in that fills a volume with a lattice - an irregular Voronoi
foam or a regular cell: cubic, BCC, FCC, octet truss, diamond or Kelvin - or
covers a surface or a set of faces with a Voronoi or cubic one.
The workflow: pick the
region, say whether the lattice goes inside or on the skin, where the struts sit
relative to the surface, then the porosity and strut thickness. The result is a
**Form (T-spline) body**, or a watertight mesh body where a Form is not
possible. Its porosity is **measured** on what was built and reported, not
assumed from the target.

Command: **Lattice**, in *Solid > Create* and *Mesh > Modify*.

## How to use

### Where the button is

Design workspace, two places (the same command):

- **Solid** tab > **CREATE** panel > **Lattice**
- **Mesh** tab > **MODIFY** panel > **Lattice**

### Before you start

- You need the region to fill or cover: one solid, surface or mesh body, or one
  or more faces.
- **This add-in needs an external Python.** The geometry runs outside Fusion,
  because Fusion's own Python has no numpy. Install a Python (3.10 or newer is
  what it is built and tested on) with these packages:

  ```
  pip install numpy scipy trimesh pyvista pymeshlab
  ```

  - `numpy`, `scipy`, `trimesh`, `pyvista` are required. The add-in only accepts
    a Python in which all four import. VTK comes with pyvista.
  - `pymeshlab` is used to remesh the region (see **Remesh the region**) and,
    when present, for decimation. The search below does not check for it, but
    the region remesh fails without it, so install it too.
- **How the add-in finds that Python**, in this order. The first one in which
  `import numpy, scipy, trimesh, pyvista` succeeds is used:
  1. the environment variable **`LATTICE_PYTHON`**: the full path to
     `python.exe`;
  2. a file **`python_path.txt`** in the installed add-in folder
     (`%APPDATA%\Autodesk\Autodesk Fusion 360\API\AddIns\Lattice\`), containing
     only the full path to `python.exe`;
  3. `%LOCALAPPDATA%\Programs\Python\Python310\python.exe`, then `Python311`,
     `Python312`, `Python313`, `Python39` (the default per-user install
     locations);
  4. `py -3.10`;
  5. `python`.

  If none works, **Create** ends in a message box *Lattice failed.* with *No
  Python with numpy, scipy, trimesh and pyvista was found*, and **Preview
  struts** draws nothing (the error goes to `Lattice.log` in the add-in folder).
  Once a Python has been found it is kept until Fusion restarts, so restart
  Fusion after you change `LATTICE_PYTHON` or `python_path.txt`.

### Step by step

1. Click **Lattice**.
2. **Region**: select the body or faces.
   - Clicking a solid in the canvas picks a face. To use the whole solid, pick
     the body in the browser, or tick **Whole body** (shown when the picked
     faces all belong to one solid).
   - Select either one body or faces, not both, and one body at a time.
   - Check the first lines of the panel at the bottom: **region** names the body
     and says what it is (*solid body*, *mesh body, closed*, *surface body,
     open*, *N faces, open*), and **size** gives its volume (closed) and area.
3. **Lattice** group:

   | Input | Default (first use) | What it does |
   |---|---|---|
   | **Preset** | Cubic - cell 1 mm, strut 0.65 mm | Fills in Cell, Define by, the porosity or cell size, Strut thickness, Organic, Output and Placement. Changing any of those afterwards turns it back to *Custom*. Geometry only. See the presets table below. |
   | **Cell** | Cubic - struts along the cube edges (from the preset) | *Voronoi - irregular foam*, or a regular cell: *Cubic*, *BCC*, *FCC*, *Octet truss*, *Diamond*, *Kelvin*. See **The cells** below. On a *Surface* lattice only Voronoi and Cubic are allowed. |
   | **Grid orientation** | Design X, Y, Z | Regular cells only. *Design X, Y, Z*: struts along the design axes, a node at the centre of the region. *Along an edge or axis*: the grid is turned onto a picked axis, with a node on it. |
   | **Axis** | - | Shown for *Along an edge or axis*. A straight edge, sketch line, construction axis or cylindrical face. |
   | **Align** | Opposite cell corners on the axis | Shown for *Along an edge or axis*. *Opposite cell corners on the axis*: the cell diagonal lies on the axis. *A cell edge on the axis*: one strut direction lies on it. |
   | **Spin about the axis** | 0 deg | Shown for *Along an edge or axis*. Turns the grid round the axis. |
   | **Lattice** | - choose - | Shown only for a **closed** region, and reset to *- choose -* every time you change the region. *Volume - struts fill the inside* or *Surface - struts on the skin only*. Create stays greyed out until you choose. An open region is always a surface lattice. |
   | **Placement** | Flush - inside the surface | *Flush*: struts inside, their outer side on the surface; the mesh is trimmed so nothing pokes through. *Middle*: struts centred on the surface. *External*: struts outside, their inner side on the surface. |
   | **Flip side** | off | Open region only. Swaps which side of the surface counts as inside. |
   | **Skin struts** | on | Volume lattice only. Also runs struts along the surface, closing the cut cells. |
   | **Border strut** | on | Open region only. Runs a strut along the surface's edges. |

4. **Porosity and thickness** group:

   | Input | Default (first use) | Unit | What it does |
   |---|---|---|---|
   | **Define by** | Cell size (from the preset) | - | *Porosity*: you give the porosity, the cell size is solved. *Cell size*: you give the cell, the porosity follows and is shown. |
   | **Porosity (%)** | 80 | % | Shown for *Porosity*. 1 to 99. |
   | **Cell size** | 1 (from the preset) | mm | Shown for *Cell size*. Must be larger than the strut thickness. |
   | **Strut thickness** | 0.65 (from the preset) | mm | Strut diameter. |
   | **Regularity** | 5 | steps | Voronoi only. 0 to 20. 0 is a raw random Voronoi (uneven cells, short struts), 5 is even but irregular, 20 is close to equal cells. |
   | **Random seed** | 1 | - | Voronoi only. 1 to 99999. Another seed gives another lattice with the same statistics. |

5. **Shape** group:

   | Input | Default (first use) | What it does |
   |---|---|---|
   | **Organic** | 1 (from the preset) | Slider 1 to 10: how the struts meet. 1 is close to crisp cylinders, 10 webs them together. The cell size is re-solved so the porosity still holds. |
   | **Output** | Mesh body (from the preset) | *Form (T-spline body)*: a smooth body you can edit in the Form environment. *Mesh body*: a watertight triangle mesh, any density. |

6. **Mesh** group (collapsed; shown only for *Mesh body*):

   | Input | Default | What it does |
   |---|---|---|
   | **Resolution** | Medium (4 voxels per strut) | *Coarse (3 voxels per strut)*, *Medium (4 voxels per strut)*, *Fine (6 voxels per strut)*. Finer is smoother and slower. |
   | **Max triangles (k)** | 1500 | Thousands of triangles, 50 to 20000. The mesh is decimated down to this. |

7. **Region remeshing (advanced)** group (collapsed):

   | Input | Default | Unit | What it does |
   |---|---|---|---|
   | **Remesh the region** | Automatic - long edges or slivers | - | The lattice is traced on a triangle mesh of the region. *Automatic*: remesh it when an edge is over 5 mm (or over twice the target), or over 5 % of the area is sliver triangles. *Always*. *Never*. |
   | **Target edge (0 = auto)** | 0 | mm | Edge length of the remeshed region. 0 means half a cell, never over 5 mm. |

8. Last options:

   | Input | Default | What it does |
   |---|---|---|
   | **Replace the lattice made before** | on | Deletes **every** lattice this add-in made earlier in the design before building the new one. |
   | **Hide the region afterwards** | on | Hides the region's body after the build. |
   | **Preview struts** | off, at every opening | Draws the strut centrelines in orange. Every change re-runs the strut layout in the external Python, so Fusion waits on each change; on a large region adjust first, then preview. |

9. Read the panel (see *Checking the result*), then click **Create**. It stays
   greyed out while something is missing or invalid; the panel says what.
10. A progress window shows the engine's steps; **Cancel** stops it. For a Form,
    a second progress window shows while Fusion builds the body, with the
    expected time. At the end a message box gives the report.

The presets (all cubic, Flush, Organic 1, Mesh body):

| Preset | Sets |
|---|---|
| Custom | nothing |
| Cubic - cell 1 mm, strut 0.65 mm | cell 1 mm, strut 0.65 mm (unit cell 39 %) |
| Cubic - pore 1 mm, strut 0.65 mm | cell 1.65 mm, strut 0.65 mm (72 %) |
| Cubic - cell 1 mm, 80 % | cell 1 mm, strut 0.325 mm |
| Cubic - strut 0.65 mm, 80 % | porosity 80 %, strut 0.65 mm, cell solved |
| Cubic - strut 0.65 mm, 70 % | porosity 70 %, strut 0.65 mm, cell solved |
| Cubic - strut 0.65 mm, 90 % | porosity 90 %, strut 0.65 mm, cell solved |

The dialog remembers its last settings (saved when you click Create). If the
remembered preset is not *Custom*, it is applied again when the dialog opens.

### What you get

- One body in the **component that owns the selected region**, in that
  component's coordinates, so it moves with that component.
- It is named `<Cell> <kind> <placement> <porosity>% d<strut> o<organic>`, with
  the **measured** porosity, for example `Cubic volume flush 39% d0.65 o1`.
- **Mesh body** output: a mesh body. In a parametric design it sits inside a
  Base Feature with the same name.
- **Form** output: a Form feature with the same name, holding one T-spline body.
- A Form needs the timeline. In a direct design, or when no Form cage can be
  built for the lattice, a mesh body is made instead, and the preview and the
  report say why.
- The region is hidden if **Hide the region afterwards** is on. Nothing is
  combined with the region.

### Checking the result

Live in the panel, before you build:

- **cell size** and approximate pore size (mm);
- **porosity**: for a regular cell, the unit cell's porosity; for Voronoi with
  *Cell size*, the approximate porosity;
- **grid** (regular cells): design axes, or cell diagonal / cell edge along the axis and
  the spin;
- **cells**: approximate count and total strut length (mm);
- **mesh**: voxel size and expected triangles, or for a Form, a note that the
  cage size is given once previewed;
- warnings: *pores narrower than the struts*; *fewer than 4 cells - lower the
  cell size*; *N cells will take minutes*; for a volume lattice at Middle or
  External, how far the struts stand proud of the surface.

With **Preview struts** on, also:

- **measured**: the cell size and the porosity predicted on the actual strut
  layout;
- **Form**: number of cage faces and the expected build time in Fusion;
- *No Form possible - this will be a mesh*, with the reason;
- *N % cannot be reached*, with the reason and the closest value.

After the build, the message box reports:

- lattice type, placement, organic level, cell size, strut thickness, strut
  count and total length;
- **Porosity measured** on the built mesh or on the Form's limit surface, next to
  the value you asked for (Cubic: next to the unit cell value);
- Mesh: triangle count, voxel size, and **watertight** or the number of open or
  non-manifold edges, and the number of connected pieces;
- Form: cage faces, junctions, and how many junctions were folded into a
  neighbour; at Flush, how far the Form bulges past the surface;
- *NOT A FORM* and the reason, if a mesh was built instead;
- region remeshing done (edge lengths and volume change), or a warning that the
  region has long edges and remeshing is off;
- loose strut groups removed, and the total time.

Before you use the lattice, check: the measured porosity against your target,
**watertight** with **1 connected piece** for a mesh you will print, and the
pore size against the strut thickness.

### Where it fits

Use it for porous or lattice regions of a patient-specific implant or guide, for
example a bone-contact region filled with a cubic lattice of set porosity, or a
Voronoi net on a surface. Model the region as its own body first, then combine
the lattice with the rest of the part.

### Known limits

- **Replace the lattice made before** removes every lattice the add-in made in
  the design, not only the last one. Turn it off to keep several.
- A Form at Flush cannot be trimmed; its hubs can bulge past the surface (the
  report gives by how much). Use Mesh output when nothing may cross the surface.
- Dense Voronoi volume lattices (cells under about 8 strut diameters) often
  cannot be built as a Form and fall back to a mesh.
- Fusion takes about 4.2 ms per Form cage face (about 10 000 faces in 42 s).
- At sharp edges, Middle and External skin struts sit slightly low (Flush is
  trimmed, so it is not affected).
- Surface lattices assume cells smaller than the surface's radius of curvature.
- For a regular cell the quoted porosity is the unit cell's; the trimmed part's
  own porosity can differ and is reported after the build. With *Skin* on, the
  skin struts and the short ties to them add material, so the part is denser
  than its unit cell (a 24 mm diamond box at 80 %: about 76 %).
- BCC, FCC, octet, diamond and Kelvin fill a volume only. Their struts meet at
  up to 12 per node (octet, FCC), which a Form cage may not manage; the build then
  falls back to a mesh and says why.
- A big mesh uses several CPU cores and a few GB of memory; the engine stops
  after one hour.
- The engine writes its temporary files to a `work` folder inside the installed
  add-in folder.

## The dialog

| Input | What it does |
|---|---|
| **Region** | A solid, surface or mesh body, or one or more faces. Clicking a solid picks a face: pick the body in the browser, or tick **Whole body**. |
| **Preset** | Cubic parameter sets; choosing one fills in the cell, porosity, strut, organic, output and placement, and touching any of those turns it back to *Custom*. They set **geometry only** - cell, strut and porosity - and are not a reproduction of any published implant, its material, process or test results. A 1 mm cell with 0.65 mm struts is only a 39 % unit cell (pores 0.35 mm), so both readings of "1 mm" are listed: *cell 1 mm, strut 0.65 mm* (the default, 39 %), *pore 1 mm, strut 0.65 mm* (cell 1.65 mm, 72 %), *cell 1 mm, 80 %* (strut 0.325 mm), and *strut 0.65 mm, 70 / 80 / 90 %* with the cell solved (2.00 mm at 80 %). |
| **Cell** | *Voronoi*: an irregular foam from relaxed random seeds. The regular cells (*Cubic*, *BCC*, *FCC*, *Octet truss*, *Diamond*, *Kelvin*) tile cubes of side *cell size*; see **The cells**. |
| **Grid orientation** | Regular cells only. *Design X, Y, Z*: struts along the design axes, a node at the centre of the region. *Along an edge or axis*: pick a straight edge, sketch line, construction axis or cylindrical face; the grid is turned onto it with a node on that line. **Align** *Opposite cell corners on the axis* puts the cell diagonal along it - the axis runs through two opposite vertices of every cell it crosses; *A cell edge on the axis* puts one strut direction along it. **Spin about the axis** turns the grid round the axis. |
| **Lattice** | Shown only when the region is **closed**, and reset to *choose* every time a closed region is picked, so the choice is always made on purpose. *Volume* fills the region with the lattice. *Surface* puts a Voronoi or cubic net on its skin only. An open region is always *Surface*. |
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

### The cells

Every regular cell tiles cubes of side *cell size* (s), so the grid orientation,
the skin and the cell-size / porosity solving work the same for all of them. The
skin of a regular volume lattice follows the cube walls; the interior struts that
end on the skin are tied to it by short struts along the surface.

| Cell | Struts | Strut length per cell | Struts per node | Cell size for 80 % with 0.65 mm struts |
|---|---|---|---|---|
| Cubic | the cube edges | 3 s | 6 | 2.00 mm |
| BCC | cube centre to the 8 corners | 4 sqrt3 s = 6.93 s | 8 | 3.05 mm |
| FCC | the face diagonals | 6 sqrt2 s = 8.49 s | 4 at face centres, 12 at corners | 3.36 mm |
| Octet truss | every FCC point to its 12 nearest | 12 sqrt2 s = 16.97 s | 12 | 4.76 mm |
| Diamond | the diamond-cubic bonds, s sqrt3 / 4 long | 4 sqrt3 s = 6.93 s | 4 | 3.05 mm |
| Kelvin | the edges of truncated octahedra, s / (2 sqrt2) long | 6 sqrt2 s = 8.49 s | 4 | 3.31 mm |

BCC and diamond have the same solid fraction at every d/s - both are rods along
the cube diagonals at the same density - but diamond's nodes join only four
struts, and its struts are half as long.

### What "porosity" means

- **Cubic cell**: the porosity of the **unit cell**, solved exactly - three struts
  of length s per cell and one node where three meet:
  `V_cell = 3 pi r^2 s - 8 sqrt(2) r^3` (three Steinmetz bicylinders less the
  tricylinder; checked against a boolean of three struts to 0.1 %). That is what
  a quoted "80 %" means, and a trimmed part cannot be matched instead: in a
  20 mm box the part's own porosity jumps between 75.6 and 82.8 % with how many
  grid planes fit, nothing between. The part's measured porosity is still
  reported next to the unit cell's.
- **BCC, FCC, octet, diamond, Kelvin**: also the porosity of the unit cell, but
  their struts overlap at the nodes in ways no closed formula follows, so the
  solid fraction was **measured** once per cell for d/s from 0.02 up to solid
  (4 M points per cell against the struts; the same measurement gives the cubic
  cell's exact value to 5e-5) and is interpolated from that table. The tests
  re-measure it at d/s 0.3 with other points.

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

It copies the add-in into
`%APPDATA%\Autodesk\Autodesk Fusion 360\API\AddIns\Lattice`. Then in Fusion:
*Utilities > Add-Ins*, select **Lattice**, **Run**. The engine needs a Python
with numpy, scipy, trimesh and pyvista (and pymeshlab for region remeshing). Set
`LATTICE_PYTHON` to its `python.exe`, or put that path in `python_path.txt` in
the installed add-in folder; otherwise the default per-user Python installs,
`py -3.10` and `python` are tried, in the order given under *How to use*.

`Lattice.py` is a thin shell. `vl_core.py` and `vl_plan.py` hot-reload on
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
- The output is a Form (T-spline) body or a mesh body, never a B-Rep solid built
  strut by strut: that is not practical in Fusion for thousands of struts.
