"""Voronoi Lattice - all the Fusion-side code. RE-IMPORTED ON EVERY COMMAND EVENT.

Nothing may keep state in a module global between events; anything the dialog
decides has to be readable back off its inputs. The one exception is the region
cache below, parked in sys.modules so it survives the reload: turning a body into
triangles is the slow part of a preview, and the selection does not change while
a spinner does.

What this file does:
    read the selection     one body (solid, surface or mesh) or any set of
                           faces; whether it is CLOSED decides whether the dialog
                           asks "volume or surface"
    write it as triangles  in millimetres, component space, for vl_engine
    preview                the engine's strut centrelines as custom graphics
    build                  the engine's STL imported as a mesh body in the
                           component that owns the selection

The geometry is vl_engine's, run in a separate Python with numpy, scipy and
pyvista (Fusion's own Python has none of them); vl_plan, stdlib only, gives the
dialog its live numbers.
"""

import array
import importlib
import json
import math
import os
import subprocess
import sys
import time
import traceback
import types

import adsk.core
import adsk.fusion

import vl_plan as P

importlib.reload(P)

ADDIN_DIR = os.path.dirname(os.path.realpath(__file__))
RES_DIR = os.path.join(ADDIN_DIR, 'resources', 'VoronoiLattice')
LOG_PATH = os.path.join(ADDIN_DIR, 'VoronoiLattice.log')
SETTINGS = os.path.join(ADDIN_DIR, 'last_settings.json')
ENGINE = os.path.join(ADDIN_DIR, 'vl_engine.py')
WORK = os.path.join(ADDIN_DIR, 'work')

CMD_ID = 'VoronoiLatticeCreate'
CMD_NAME = 'Voronoi Lattice'
VERSION = '0.2.0'

CM = 10.0          # internal cm -> mm
MM = 0.1           # mm -> internal cm
ATTR = 'VoronoiLattice'
NO_WINDOW = 0x08000000 if os.name == 'nt' else 0

KINDS = [('', '- choose -'),
         ('volume', 'Volume - struts fill the inside'),
         ('surface', 'Surface - struts on the skin only')]
PLACES = [('flush', 'Flush - inside the surface'),
          ('middle', 'Middle - centred on the surface'),
          ('external', 'External - outside the surface')]
DEFINE = [('porosity', 'Porosity'), ('cell', 'Cell size')]
OUTPUTS = [('form', 'Form (T-spline body)'), ('mesh', 'Mesh body')]
RESOLUTIONS = [('coarse', 'Coarse (3 voxels per strut)'),
               ('medium', 'Medium (4 voxels per strut)'),
               ('fine', 'Fine (6 voxels per strut)')]

DEFAULTS = {
    'place': 0, 'flip': False, 'skin': True, 'border': True, 'define': 0,
    'porosity': 80.0, 'cell': 3.0, 'thickness': 0.6, 'regularity': 5, 'seed': 1,
    'resolution': 1, 'max_ktris': 1500, 'replace': True, 'hide': True,
    'preview': False, 'organic': 3, 'output': 0,
}


def log(msg):
    try:
        with open(LOG_PATH, 'a', encoding='utf-8') as fh:
            fh.write('%s  %s\n' % (time.strftime('%H:%M:%S'), msg))
    except Exception:
        pass


def short_error(tb):
    for line in reversed(tb.strip().splitlines()):
        line = line.strip()
        if line and not line.startswith(('File "', 'Traceback', '^', '~')):
            return line
    return tb.strip().splitlines()[-1]


def load_settings():
    s = dict(DEFAULTS)
    try:
        with open(SETTINGS, encoding='utf-8') as fh:
            s.update(json.load(fh))
    except Exception:
        pass
    return s


def save_settings(s):
    try:
        with open(SETTINGS, 'w', encoding='utf-8') as fh:
            json.dump(s, fh, indent=2)
    except Exception:
        pass


def _state():
    st = sys.modules.get('_vl_state')
    if st is None:
        st = types.ModuleType('_vl_state')
        st.regions = {}
        sys.modules['_vl_state'] = st
    return st


# --------------------------------------------------------------------------- #
# finding a Python with numpy, scipy and pyvista
# --------------------------------------------------------------------------- #

def find_python():
    cands = []
    if os.environ.get('VORONOILATTICE_PYTHON'):
        cands.append([os.environ['VORONOILATTICE_PYTHON']])
    txt = os.path.join(ADDIN_DIR, 'python_path.txt')
    if os.path.exists(txt):
        try:
            with open(txt, encoding='utf-8') as fh:
                line = fh.read().strip()
            if line:
                cands.append([line])
        except Exception:
            pass
    local = os.environ.get('LOCALAPPDATA', '')
    for v in ('310', '311', '312', '313', '39'):
        cands.append([os.path.join(local, 'Programs', 'Python', 'Python%s' % v,
                                   'python.exe')])
    cands.append(['py', '-3.10'])
    cands.append(['python'])
    for c in cands:
        try:
            r = subprocess.run(c + ['-c', 'import numpy, scipy, trimesh, pyvista'],
                               capture_output=True, timeout=90,
                               creationflags=NO_WINDOW)
            if r.returncode == 0:
                return c
        except Exception:
            continue
    return None


def python_cmd():
    st = _state()
    if not getattr(st, 'python', None):
        st.python = find_python() or []
    if not st.python:
        raise RuntimeError(
            'No Python with numpy, scipy, trimesh and pyvista was found. Set the '
            'VORONOILATTICE_PYTHON environment variable, or put the full path to '
            'the interpreter in python_path.txt next to the add-in.')
    return list(st.python)


# --------------------------------------------------------------------------- #
# the region
# --------------------------------------------------------------------------- #

def selected(sel_input):
    return [sel_input.selection(i).entity for i in range(sel_input.selectionCount)]


def _native(e):
    try:
        n = e.nativeObject
        return n if n else e
    except Exception:
        return e


def _context(e):
    try:
        return e.assemblyContext
    except Exception:
        return None


def read_region(entities, whole=False):
    """What the selection is. Raises ValueError with a message for the dialog.

    A body counts as CLOSED when Fusion says it is a solid (or a closed mesh);
    a set of faces when it is every face of a solid body."""
    if not entities:
        raise ValueError('Select a body or faces.')
    bodies = [e for e in entities if adsk.fusion.BRepBody.cast(e)
              or adsk.fusion.MeshBody.cast(e)]
    faces = [e for e in entities if adsk.fusion.BRepFace.cast(e)]
    if bodies and faces:
        raise ValueError('Select either one body or faces, not both.')
    if len(bodies) > 1:
        raise ValueError('Select one body at a time.')
    ctx = _context(entities[0])
    if faces:
        owners = {}
        for f in faces:
            b = _native(f).body
            owners[b.entityToken] = b
        if len(owners) == 1:
            b = list(owners.values())[0]
            if whole or (b.isSolid and len(faces) == b.faces.count):
                bodies, faces = [b], []
    if bodies:
        b = _native(bodies[0])
        mb = adsk.fusion.MeshBody.cast(b)
        if mb:
            closed = bool(mb.isClosed)
            what = 'mesh body, %s' % ('closed' if closed else 'open')
            area = mb.area * 100.0
            vol = mb.volume * 1000.0 if closed else 0.0
        else:
            closed = bool(b.isSolid)
            what = 'solid body' if closed else 'surface body, open'
            area = b.area * 100.0
            vol = b.volume * 1000.0 if closed else 0.0
        return {'entities': [b], 'closed': closed, 'what': what, 'name': b.name,
                'area': area, 'volume': vol, 'component': b.parentComponent,
                'context': ctx, 'body': b, 'faces_of_one_solid': None}
    nat = [_native(f) for f in faces]
    area = sum(f.area for f in nat) * 100.0
    owner_solid = None
    if len({f.body.entityToken for f in nat}) == 1 and nat[0].body.isSolid:
        owner_solid = nat[0].body
    return {'entities': nat, 'closed': False,
            'what': '%d face%s, open' % (len(nat), '' if len(nat) == 1 else 's'),
            'name': nat[0].body.name, 'area': area, 'volume': 0.0,
            'component': nat[0].body.parentComponent, 'context': ctx,
            'body': nat[0].body, 'faces_of_one_solid': owner_solid}


def _calc_mesh(mesh_manager, bbox):
    """Triangles of a body or face. Tolerances tied to its own size: a tight
    tolerance with no side or aspect limit has taken this Fusion down."""
    diag = bbox.minPoint.distanceTo(bbox.maxPoint) or 1.0
    calc = mesh_manager.createMeshCalculator()
    calc.surfaceTolerance = max(diag / 4000.0, 1.0e-4)
    calc.maxSideLength = diag / 15.0
    calc.maxAspectRatio = 10.0
    calc.maxNormalDeviation = math.radians(10.0)
    tm = calc.calculate()
    return (tm.nodeCoordinatesAsDouble, tm.nodeIndices, tm.normalVectorsAsDouble)


def region_triangles(reg):
    """(coords cm, indices, normals or None), in the owning component's space."""
    coords, idx, nrm = [], [], []
    have_normals = True
    for e in reg['entities']:
        mb = adsk.fusion.MeshBody.cast(e)
        if mb:
            pm = mb.mesh
            if pm and pm.triangleCount > 0 and pm.quadCount == 0 and pm.polygonCount == 0:
                c, t = pm.nodeCoordinatesAsDouble, pm.triangleNodeIndices
            else:
                dm = mb.displayMesh
                c, t = dm.nodeCoordinatesAsDouble, dm.nodeIndices
            n = None
        else:
            c, t, n = _calc_mesh(e.meshManager, e.boundingBox)
        base = len(coords) // 3
        coords.extend(c)
        idx.extend(i + base for i in t)
        if n is not None and len(n) == len(c):
            nrm.extend(n)
        else:
            have_normals = False
    return coords, idx, (nrm if have_normals else None)


def region_key(reg):
    return '|'.join(e.entityToken for e in reg['entities'])


def write_region(reg):
    """Write the region as raw little-endian arrays (mm) for the engine, once
    per selection - the preview asks again on every spinner click."""
    st = _state()
    key = region_key(reg)
    hit = st.regions.get(key)
    if hit and all(os.path.exists(p) for p in hit.values() if p):
        return hit
    coords, idx, nrm = region_triangles(reg)
    os.makedirs(WORK, exist_ok=True)
    tag = 'r%d' % (len(st.regions) % 8)
    paths = {'coords': os.path.join(WORK, tag + '_v.bin'),
             'tris': os.path.join(WORK, tag + '_f.bin'),
             'normals': os.path.join(WORK, tag + '_n.bin') if nrm else None}
    with open(paths['coords'], 'wb') as fh:
        array.array('d', [v * CM for v in coords]).tofile(fh)
    with open(paths['tris'], 'wb') as fh:
        array.array('i', idx).tofile(fh)
    if nrm:
        with open(paths['normals'], 'wb') as fh:
            array.array('d', nrm).tofile(fh)
    for k in [k for k, v in st.regions.items() if v.get('coords') == paths['coords']]:
        del st.regions[k]
    st.regions[key] = paths
    return paths


# --------------------------------------------------------------------------- #
# reading the dialog
# --------------------------------------------------------------------------- #

def _item(inputs, iid):
    return inputs.itemById(iid)


def _drop(inputs, iid, table):
    it = _item(inputs, iid).selectedItem
    return table[it.index][0] if it else table[0][0]


def region_from(inputs):
    sel = _item(inputs, 'sel')
    whole = _item(inputs, 'whole')
    return read_region(selected(sel), whole.isVisible and whole.value)


def kind_from(inputs, reg):
    if not reg['closed']:
        return 'surface'
    return _drop(inputs, 'kind', KINDS) or None


def params_from(inputs, reg):
    kind = kind_from(inputs, reg)
    return {
        'kind': kind,
        'closed': reg['closed'],
        'placement': _drop(inputs, 'place', PLACES),
        'flip': bool(_item(inputs, 'flip').value) and not reg['closed'],
        'skin': bool(_item(inputs, 'skin').value),
        'border': bool(_item(inputs, 'border').value),
        'define': _drop(inputs, 'define', DEFINE),
        'porosity': _item(inputs, 'porosity').value / 100.0,
        'cell': _item(inputs, 'cell').value * CM,
        'thickness': _item(inputs, 'thickness').value * CM,
        'regularity': int(_item(inputs, 'regularity').value),
        'seed': int(_item(inputs, 'seed').value),
        'resolution': _drop(inputs, 'resolution', RESOLUTIONS),
        'max_tris': int(_item(inputs, 'max_ktris').value) * 1000,
        'organic': int(_item(inputs, 'organic').valueOne),
        'output': _drop(inputs, 'output', OUTPUTS),
    }


def settings_from(inputs):
    def idx(iid):
        it = _item(inputs, iid).selectedItem
        return it.index if it else 0
    return {
        'place': idx('place'), 'flip': _item(inputs, 'flip').value,
        'skin': _item(inputs, 'skin').value, 'border': _item(inputs, 'border').value,
        'define': idx('define'), 'porosity': _item(inputs, 'porosity').value,
        'cell': _item(inputs, 'cell').value * CM,
        'thickness': _item(inputs, 'thickness').value * CM,
        'regularity': int(_item(inputs, 'regularity').value),
        'seed': int(_item(inputs, 'seed').value), 'resolution': idx('resolution'),
        'max_ktris': int(_item(inputs, 'max_ktris').value),
        'replace': _item(inputs, 'replace').value, 'hide': _item(inputs, 'hide').value,
        'preview': _item(inputs, 'preview').value,
        'organic': int(_item(inputs, 'organic').valueOne), 'output': idx('output'),
    }


def check_params(p):
    """The engine checks these too; checking here keeps OK greyed out."""
    if p['kind'] is None:
        raise ValueError('This region is closed: choose a Volume or a Surface '
                         'lattice.')
    if p['thickness'] <= 0:
        raise ValueError('The strut thickness must be positive.')
    if p['define'] == 'cell':
        if p['cell'] <= p['thickness']:
            raise ValueError('The cell size must be larger than the strut '
                             'thickness, or the pores close.')
        return p['cell']
    return P.cell_for(p['kind'], p['porosity'], p['thickness'])


# --------------------------------------------------------------------------- #
# the dialog
# --------------------------------------------------------------------------- #

def _esc(t):
    return str(t).replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')


def _rows(pairs):
    out = ['<table cellpadding="0" cellspacing="0">']
    for k, v in pairs:
        out.append('<tr><td width="105"><b>%s</b></td><td>%s</td></tr>'
                   % (_esc(k), v))
    out.append('</table>')
    return ''.join(out)


def on_created(app, cmd):
    s = load_settings()
    _state().regions = {}
    cmd.isRepeatable = False
    cmd.okButtonText = 'Create'
    cmd.setDialogInitialSize(380, 640)
    inputs = cmd.commandInputs
    rl = adsk.core.ValueInput.createByReal
    tl = adsk.core.DropDownStyles.TextListDropDownStyle

    sel = inputs.addSelectionInput(
        'sel', 'Region',
        'A solid, surface or mesh body, or one or more faces. Clicking a solid '
        'picks a face: pick the body in the browser, or tick "Whole body".')
    for f in ('SolidBodies', 'SurfaceBodies', 'MeshBodies', 'Faces'):
        sel.addSelectionFilter(f)
    sel.setSelectionLimits(1, 0)
    inputs.addBoolValueInput('whole', 'Whole body', True, '', False)

    g = inputs.addGroupCommandInput('gLattice', 'Lattice')
    gi = g.children
    kd = gi.addDropDownCommandInput('kind', 'Lattice', tl)
    for i, (_k, lab) in enumerate(KINDS):
        kd.listItems.add(lab, i == 0)
    kd.tooltip = ('Asked every time the region is closed: Volume fills it with a '
                  '3D Voronoi foam, Surface puts a Voronoi net on its skin only.')
    pl = gi.addDropDownCommandInput('place', 'Placement', tl)
    for i, (_k, lab) in enumerate(PLACES):
        pl.listItems.add(lab, i == int(s['place']))
    pl.tooltip = ('Where the struts sit relative to the surface. Flush: inside, '
                  'their outer side on the surface - nothing pokes through. '
                  'Middle: centred on it. External: outside, their inner side on '
                  'it. On an open surface "inside" is the back of the face.')
    gi.addBoolValueInput('flip', 'Flip side', True, '', bool(s['flip']))
    gi.addBoolValueInput('skin', 'Skin struts', True, '', bool(s['skin']))
    _item(inputs, 'skin').tooltip = ('Volume lattice: also run struts along the '
                                     'surface, closing the cut cells, so no strut '
                                     'ends loose on the skin.')
    gi.addBoolValueInput('border', 'Border strut', True, '', bool(s['border']))
    _item(inputs, 'border').tooltip = 'Open surface: run a strut along its edges.'

    q = inputs.addGroupCommandInput('gSize', 'Porosity and thickness')
    qi = q.children
    df = qi.addDropDownCommandInput('define', 'Define by', tl)
    for i, (_k, lab) in enumerate(DEFINE):
        df.listItems.add(lab, i == int(s['define']))
    qi.addFloatSpinnerCommandInput('porosity', 'Porosity (%)', '', 1.0, 99.0, 1.0,
                                   float(s['porosity']))
    _item(inputs, 'porosity').tooltip = (
        'Volume lattice: void fraction of the lattice envelope. Surface lattice: '
        'open fraction of the surface. The cell size is solved for it; the value '
        'reached is measured on the mesh and reported.')
    qi.addValueInput('cell', 'Cell size', 'mm', rl(float(s['cell']) * MM))
    qi.addValueInput('thickness', 'Strut thickness', 'mm',
                     rl(float(s['thickness']) * MM))
    _item(inputs, 'thickness').tooltip = 'Strut diameter.'
    qi.addIntegerSpinnerCommandInput('regularity', 'Regularity', 0, 20, 1,
                                     int(s['regularity']))
    _item(inputs, 'regularity').tooltip = (
        'Lloyd relaxation steps on the seeds. 0 is a raw random Voronoi (very '
        'uneven cells, short struts); 5 is even but still irregular; 20 is '
        'close to a foam of equal cells.')
    qi.addIntegerSpinnerCommandInput('seed', 'Random seed', 1, 99999, 1,
                                     int(s['seed']))

    sh = inputs.addGroupCommandInput('gShape', 'Shape')
    si = sh.children
    org = si.addIntegerSliderCommandInput('organic', 'Organic', 1, 10, False)
    org.valueOne = int(s['organic'])
    org.tooltip = ('How the struts meet. 1: nearly as cylinders do, a crisp joint. '
                   '10: webbed together like trabecular bone. The cell size is '
                   're-solved so the porosity still holds.')
    ou = si.addDropDownCommandInput('output', 'Output', tl)
    for i, (_k, lab) in enumerate(OUTPUTS):
        ou.listItems.add(lab, i == int(s['output']))
    ou.tooltip = ('Form: a real T-spline body, smooth and editable in the Form '
                  'workspace - built from the struts themselves, a tube per strut '
                  'and a hub per junction. Works for surface lattices and slender '
                  'volume ones; a lattice too dense for a cage is built as a mesh '
                  'instead, and the dialog says why. Mesh: a watertight triangle '
                  'mesh, any density.')

    m = inputs.addGroupCommandInput('gMesh', 'Mesh')
    m.isExpanded = False
    mi = m.children
    rs = mi.addDropDownCommandInput('resolution', 'Resolution', tl)
    for i, (_k, lab) in enumerate(RESOLUTIONS):
        rs.listItems.add(lab, i == int(s['resolution']))
    mi.addIntegerSpinnerCommandInput('max_ktris', 'Max triangles (k)', 50, 20000,
                                     250, int(s['max_ktris']))

    inputs.addBoolValueInput('replace', 'Replace the lattice made before', True,
                             '', bool(s['replace']))
    inputs.addBoolValueInput('hide', 'Hide the region afterwards', True, '',
                             bool(s['hide']))
    # off at every opening, whatever was used last: a preview recomputes on every change
    inputs.addBoolValueInput('preview', 'Preview struts', True, '', False)
    info = inputs.addTextBoxCommandInput('info', '', '', 12, True)
    info.isFullWidth = True
    refresh(inputs)


def refresh(inputs):
    reg = None
    try:
        reg = region_from(inputs)
    except Exception:
        pass
    closed = bool(reg and reg['closed'])
    kind = kind_from(inputs, reg) if reg else None
    _item(inputs, 'kind').isVisible = closed
    _item(inputs, 'flip').isVisible = bool(reg) and not closed
    _item(inputs, 'border').isVisible = bool(reg) and not closed
    _item(inputs, 'skin').isVisible = kind == 'volume'
    _item(inputs, 'whole').isVisible = bool(reg and reg['faces_of_one_solid'])
    cell = _drop(inputs, 'define', DEFINE) == 'cell'
    _item(inputs, 'porosity').isVisible = not cell
    _item(inputs, 'cell').isVisible = cell
    _item(inputs, 'gMesh').isVisible = _drop(inputs, 'output', OUTPUTS) == 'mesh'
    _item(inputs, 'info').formattedText = describe(inputs)


def describe(inputs):
    try:
        reg = region_from(inputs)
    except ValueError as exc:
        return _rows([('', _esc(exc))])
    except Exception:
        log(traceback.format_exc())
        return _rows([('', _esc(short_error(traceback.format_exc())))])
    pairs = [('region', '%s - %s' % (_esc(reg['name']), _esc(reg['what'])))]
    if reg['closed']:
        pairs.append(('size', 'V %.4g mm&sup3;, A %.4g mm&sup2;'
                      % (reg['volume'], reg['area'])))
    else:
        pairs.append(('size', 'A %.4g mm&sup2;' % reg['area']))
    p = params_from(inputs, reg)
    if p['kind'] is None:
        pairs.append(('', '<b>Closed region: choose Volume or Surface.</b>'))
        return _rows(pairs)
    try:
        s = check_params(p)
    except (ValueError, P.PlanError) as exc:
        pairs.append(('', '<b>%s</b>' % _esc(exc)))
        return _rows(pairs)
    d = p['thickness']
    est = P.estimate(p['kind'], s, d, reg['volume'], reg['area'], p['placement'],
                     p['resolution'])
    phi = P.porosity_of(p['kind'], s, d)
    pairs.append(('cell size', '%.3f mm, pores about %.3f mm' % (s, est['pore_mm'])))
    if p['define'] == 'cell':
        pairs.append(('porosity', 'about %.1f %%' % (phi * 100.0)))
    pairs.append(('cells', 'about %d, %.0f mm of strut' % (round(est['cells']),
                                                            est['length_mm'])))
    if p['output'] == 'mesh':
        tris = est['raw_triangles']
        pairs.append(('mesh', 'voxel %.3f mm, ~%.2g M triangles%s'
                      % (est['voxel_mm'], tris / 1e6,
                         ' (decimated to %.2g M)' % (p['max_tris'] / 1e6)
                         if tris > p['max_tris'] else '')))
    else:
        pairs.append(('form', 'cage size and Fusion time once previewed'))
    if p['organic'] > 1:
        pairs.append(('organic', '%d: the cell size grows to keep the porosity'
                      % p['organic']))
    if est['pore_mm'] < d:
        pairs.append(('', '<b>pores narrower than the struts</b>'))
    if est['cells'] < 4:
        pairs.append(('', '<b>fewer than 4 cells</b> - lower the cell size'))
    if est['cells'] > 60000:
        pairs.append(('', '<b>%d cells will take minutes</b>' % est['cells']))
    if p['placement'] != 'flush' and p['kind'] == 'volume':
        pairs.append(('', 'struts stand up to %.3g mm proud of the surface'
                      % (d if p['placement'] == 'external' else 0.5 * d)))
    return _rows(pairs)


def on_input_changed(app, inputs, changed_id):
    if changed_id in ('sel', 'whole'):
        # "always ask": a newly picked closed region starts with no choice made
        _item(inputs, 'kind').listItems.item(0).isSelected = True
    refresh(inputs)


def on_validate(inputs):
    try:
        reg = region_from(inputs)
        check_params(params_from(inputs, reg))
        return True
    except Exception:
        return False


# --------------------------------------------------------------------------- #
# running the engine
# --------------------------------------------------------------------------- #

def build_job(reg, p, mode):
    job = dict(p)
    job['mode'] = mode
    job['mesh'] = write_region(reg)
    job['cache_dir'] = os.path.join(WORK, 'cache')
    job['out_stl'] = os.path.join(WORK, 'lattice.stl')
    job['out_tsm'] = os.path.join(WORK, 'lattice.tsm')
    return job


def run_engine(job, ui=None, timeout=3600):
    """Run vl_engine on a job. With ui, show a cancellable progress dialog fed
    from the engine's log, and keep Fusion responsive while it works."""
    os.makedirs(WORK, exist_ok=True)
    jp = os.path.join(WORK, 'job_%s.json' % job['mode'])
    op = os.path.join(WORK, 'out_%s.json' % job['mode'])
    for f in (op, op + '.log'):
        if os.path.exists(f):
            os.remove(f)
    with open(jp, 'w', encoding='utf-8') as fh:
        json.dump(job, fh)
    cmd = python_cmd() + [ENGINE, jp, op]
    if ui is None:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout,
                           creationflags=NO_WINDOW)
        rc, tail = r.returncode, (r.stdout[-600:] + r.stderr[-600:])
    else:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                creationflags=NO_WINDOW)
        pd = ui.createProgressDialog()
        pd.isCancelButtonShown = True
        pd.show(CMD_NAME, 'Starting...', 0, 100, 0)
        t0 = time.time()
        try:
            while proc.poll() is None:
                adsk.doEvents()
                if pd.wasCancelled:
                    proc.kill()
                    raise RuntimeError('Cancelled.')
                if time.time() - t0 > timeout:
                    proc.kill()
                    raise RuntimeError('The lattice engine ran past %d s.' % timeout)
                msg = _last_line(op + '.log')
                if msg:
                    pd.message = '%s  (%.0f s)' % (msg, time.time() - t0)
                time.sleep(0.15)
        finally:
            pd.hide()
        out, err = proc.communicate()
        rc, tail = proc.returncode, (out[-600:] + err[-600:]).decode('utf-8', 'replace')
    if not os.path.exists(op):
        raise RuntimeError('The lattice engine produced nothing (rc=%s)\n%s' % (rc, tail))
    with open(op, encoding='utf-8') as fh:
        res = json.load(fh)
    if not res.get('ok'):
        if res.get('user_error'):
            raise ValueError(res['error'])
        raise RuntimeError('The lattice engine failed:\n%s' % res.get('error', '')[-1500:])
    return res['result']


def _last_line(path):
    try:
        with open(path, encoding='utf-8') as fh:
            lines = fh.read().strip().splitlines()
        return lines[-1].split('  ', 1)[-1].strip() if lines else ''
    except Exception:
        return ''


# --------------------------------------------------------------------------- #
# preview: strut centrelines as custom graphics
# --------------------------------------------------------------------------- #

def _affine(ctx):
    """The occurrence transform as o + x ex + y ey + z ez, read through Fusion's
    own transformBy - Matrix3D.asArray's ordering is not documented."""
    if ctx is None:
        return None
    m = ctx.transform2
    pts = []
    for v in ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)):
        p = adsk.core.Point3D.create(*v)
        p.transformBy(m)
        pts.append((p.x, p.y, p.z))
    o = pts[0]
    return o, [tuple(pts[k][j] - o[j] for j in range(3)) for k in (1, 2, 3)]


PREVIEW_ID = 'VoronoiLatticePreview'


def clear_preview(des):
    """Fusion may or may not drop graphics made in a discarded preview; never
    rely on it - remove ours by id."""
    groups = des.rootComponent.customGraphicsGroups
    for i in range(groups.count - 1, -1, -1):
        g = groups.item(i)
        if g.id == PREVIEW_ID:
            g.deleteMe()


def draw_segments(des, segments, ctx):
    clear_preview(des)
    group = des.rootComponent.customGraphicsGroups.add()
    group.id = PREVIEW_ID
    T = _affine(ctx)
    flat = []
    for s in segments:
        for p in (s[0:3], s[3:6]):
            x, y, z = p[0] * MM, p[1] * MM, p[2] * MM
            if T:
                o, (ex, ey, ez) = T
                x, y, z = (o[0] + x * ex[0] + y * ey[0] + z * ez[0],
                           o[1] + x * ex[1] + y * ey[1] + z * ez[1],
                           o[2] + x * ex[2] + y * ey[2] + z * ez[2])
            flat.extend((x, y, z))
    coords = adsk.fusion.CustomGraphicsCoordinates.create(flat)
    lines = group.addLines(coords, [], False)
    lines.color = adsk.fusion.CustomGraphicsSolidColorEffect.create(
        adsk.core.Color.create(230, 120, 30, 255))
    lines.weight = 1.5
    lines.depthPriority = 1
    return group


def on_preview(app, inputs, args):
    if not _item(inputs, 'preview').value:
        return
    reg = region_from(inputs)
    p = params_from(inputs, reg)
    check_params(p)
    res = run_engine(build_job(reg, p, 'graph'))
    # any write to an input wipes custom graphics: text first, lines last
    info = _item(inputs, 'info')
    text = describe(inputs) + _rows(preview_rows(p, res))
    if info.formattedText != text:
        info.formattedText = text
    draw_segments(adsk.fusion.Design.cast(app.activeProduct), res['segments'],
                  reg['context'])


def preview_rows(p, res):
    g = res['graph']
    rows = [('measured', 'cell %.3f mm, porosity %.1f %% (%s)'
             % (g['cell_mm'], 100.0 * res['porosity_predicted'],
                'Form limit surface' if res['output'] == 'form' else 'mesh field'))]
    if res['output'] == 'form' and 'cage' in g:
        rows.append(('Form', '%d cage faces, about %.0f s to build in Fusion'
                     % (g['cage']['cage_faces'], res.get('fusion_form_s', 0))))
    if g.get('form_failed'):
        rows.append(('', '<b>No Form possible - this will be a mesh.</b> %s'
                     % _esc(g['form_failed'])))
    return rows


# --------------------------------------------------------------------------- #
# build
# --------------------------------------------------------------------------- #

def owning_base(body):
    """baseOrFormFeature RAISES (InternalValidationError : fragment) in a direct
    design instead of returning None."""
    try:
        bf = body.baseOrFormFeature
        return bf if bf and bf.isValid else None
    except Exception:
        return None


def mesh_volume_mm3(mb):
    """From the triangles: MeshBody.volume read 0.0 on a closed lattice."""
    pm = mb.mesh
    c, t = pm.nodeCoordinatesAsDouble, pm.triangleNodeIndices
    v = 0.0
    for i in range(0, len(t), 3):
        a, b, d = 3 * t[i], 3 * t[i + 1], 3 * t[i + 2]
        ax, ay, az = c[a], c[a + 1], c[a + 2]
        bx, by, bz = c[b], c[b + 1], c[b + 2]
        dx, dy, dz = c[d], c[d + 1], c[d + 2]
        v += ax * (by * dz - bz * dy) - ay * (bx * dz - bz * dx) + az * (bx * dy - by * dx)
    return v / 6.0 * 1000.0


def remove_previous(des):
    """Delete what an earlier run made - found by the attribute it wrote. Base
    features first: deleting one takes its mesh body with it, and a mesh body
    does not report the base feature it lives in."""
    gone = 0
    for kind in (adsk.fusion.FormFeature, adsk.fusion.BaseFeature,
                 adsk.fusion.MeshBody):
        # base features first, then search again: deleting a base feature
        # LEAVES its mesh body in the design (measured), so the body goes too
        for a in list(des.findAttributes(ATTR, 'role')):
            ent = kind.cast(a.parent) if a.isValid else None
            try:
                if ent and ent.isValid:
                    ent.deleteMe()
                    gone += 1
            except Exception:
                log(traceback.format_exc())
    return gone


def import_mesh(des, comp, stl, name, summary):
    parametric = des.designType == adsk.fusion.DesignTypes.ParametricDesignType
    base = None
    if parametric:
        base = comp.features.baseFeatures.add()
        base.name = name
        base.startEdit()
    try:
        made = comp.meshBodies.add(stl, adsk.fusion.MeshUnits.MillimeterMeshUnit,
                                   base) if base else \
            comp.meshBodies.add(stl, adsk.fusion.MeshUnits.MillimeterMeshUnit)
        body = made.item(0)
        body.name = name
        body.attributes.add(ATTR, 'role', 'lattice')
        body.attributes.add(ATTR, 'params', json.dumps(summary))
    finally:
        if base:
            base.finishEdit()
    if base:
        base.attributes.add(ATTR, 'role', 'lattice')
        body = comp.meshBodies.item(comp.meshBodies.count - 1)
    return body, base


def import_form(des, comp, tsm, name, summary, ui=None, seconds=0.0):
    """A Form body from the engine's TSM. A cage Fusion refuses leaves a
    half-built form feature that makes every later formFeatures.add() fail
    (InternalValidationError : baseFeature), so it is deleted on failure."""
    with open(tsm, encoding='utf-8') as fh:
        text = fh.read()
    try:
        des.timeline.moveToEnd()
    except Exception:
        pass
    pd = None
    if ui is not None:
        pd = ui.createProgressDialog()
        pd.show(CMD_NAME, 'Fusion is building the Form body (about %.0f s)...'
                % seconds, 0, 1, 0)
        adsk.doEvents()
    ff = comp.features.formFeatures.add()
    try:
        ff.startEdit()
        try:
            tb = ff.tSplineBodies.addByTSMDescription(text)
        finally:
            ff.finishEdit()
        if tb is None:
            raise RuntimeError('Fusion returned no T-spline body')
    except Exception:
        try:
            ff.deleteMe()
        except Exception:
            log(traceback.format_exc())
        raise
    finally:
        if pd:
            pd.hide()
    ff.name = name
    ff.attributes.add(ATTR, 'role', 'lattice')
    ff.attributes.add(ATTR, 'params', json.dumps(summary))
    body = ff.bodies.item(ff.bodies.count - 1) if ff.bodies.count else None
    if body:
        body.name = name
    return body, ff


def lattice_name(p, res):
    kind = 'volume' if p['kind'] == 'volume' else 'surface'
    return 'Voronoi %s %s %.0f%% d%.2f o%d' % (kind, p['placement'],
                                               res['porosity_measured'] * 100.0,
                                               p['thickness'], p.get('organic', 1))


def build(app, ui, reg, p):
    des = adsk.fusion.Design.cast(app.activeProduct)
    check_params(p)
    p = dict(p)
    if (p.get('output') == 'form'
            and des.designType != adsk.fusion.DesignTypes.ParametricDesignType):
        p['output'] = 'mesh'
        p['form_skipped'] = ('a Form needs a parametric design (timeline on); '
                             'this one is direct')
    res = run_engine(build_job(reg, p, 'mesh'), ui)
    summary = {'version': VERSION, 'params': p,
               'porosity_measured': res['porosity_measured'],
               'cell_mm': res['graph']['cell_mm'], 'output': res['output']}
    name = lattice_name(p, res)
    if res['output'] == 'form':
        summary['cage_faces'] = res['form']['cage_faces']
        body, _ff = import_form(des, reg['component'], res['tsm'], name, summary,
                                ui, res.get('fusion_form_s', 0.0))
    else:
        summary['triangles'] = res['mesh']['triangles']
        body, _base = import_mesh(des, reg['component'], res['stl'], name, summary)
    res['form_skipped'] = p.get('form_skipped')
    return body, res


def report(p, res):
    g = res['graph']
    what = {'volume': 'Volume', 'surface': 'Surface'}[p['kind']]
    asked = ('' if p['define'] == 'cell' else
             '  (asked %.1f %%)' % (p['porosity'] * 100.0))
    lines = ['%s lattice, %s, organic %d.' % (what, p['placement'],
                                             p.get('organic', 1)),
             '',
             'Cell size %.3f mm, strut %.3f mm, %d struts (%.0f mm).'
             % (g['cell_mm'], p['thickness'], g['struts'], g['length_mm'])]
    if res['output'] == 'form':
        f = res['form']
        lines += ['Porosity measured on the Form limit surface: %.1f %%%s'
                  % (res['porosity_measured'] * 100.0, asked),
                  'Form body: %d cage faces (%d quads), %d junctions%s.'
                  % (f['cage_faces'], f['cage_quads'], f['junctions'],
                     ', %d folded into a neighbour' % f['folded'] if f['folded'] else '')]
        if f.get('protrudes_mm', 0.0) > 0.005:
            lines.append('Flush: the Form bulges up to %.3f mm past the surface at '
                         'its hubs - a Form cannot be trimmed.' % f['protrudes_mm'])
    else:
        m = res['mesh']
        if p.get('output') == 'form' or res.get('form_skipped'):
            lines += ['', 'NOT A FORM: %s' % (res.get('form_skipped')
                                             or g.get('form_failed', '')),
                      'Built as a mesh instead.', '']
        lines += ['Porosity measured on the mesh: %.1f %%%s'
                  % (res['porosity_measured'] * 100.0, asked),
                  'Mesh: %d triangles, voxel %.3f mm, %s.'
                  % (m['triangles'], m['voxel_mm'],
                     'watertight' if m['open_or_nonmanifold_edges'] == 0 else
                     '%d OPEN OR NON-MANIFOLD EDGES' % m['open_or_nonmanifold_edges']),
                  '%d connected piece%s.' % (m['components'],
                                             '' if m['components'] == 1 else 's')]
    if g.get('dropped_groups'):
        lines.append('%d loose strut groups removed.' % g['dropped_groups'])
    lines.append('Took %.1f s.' % res['timing']['total_s'])
    return '\n'.join(lines)


def on_execute(app, ui, inputs):
    des = adsk.fusion.Design.cast(app.activeProduct)
    clear_preview(des)
    reg = region_from(inputs)
    p = params_from(inputs, reg)
    save_settings(settings_from(inputs))
    if _item(inputs, 'replace').value:
        remove_previous(des)
    try:
        body, res = build(app, ui, reg, p)
    except ValueError as exc:
        ui.messageBox(str(exc), CMD_NAME)
        return None
    if _item(inputs, 'hide').value:
        try:
            reg['body'].isLightBulbOn = False
        except Exception:
            log(traceback.format_exc())
    log('%s %s' % (VERSION, json.dumps({'params': p, 'graph': res['graph'],
                                        'output': res['output'],
                                        'mesh': res.get('mesh'), 'form': res.get('form'),
                                        'porosity': res['porosity_measured']})))
    ui.messageBox(report(p, res), CMD_NAME)
    return body


def on_destroy(app):
    des = adsk.fusion.Design.cast(app.activeProduct)
    if des:
        clear_preview(des)


# --------------------------------------------------------------------------- #
# toolbar
# --------------------------------------------------------------------------- #

TOOLTIP = ('Fill a volume with a Voronoi lattice, or cover a surface or faces '
           'with one.\n\n'
           'Pick a solid, surface or mesh body, or faces. A closed region asks '
           'whether the lattice fills it (Volume) or covers its skin (Surface). '
           'Choose where the struts sit - flush inside, centred on, or outside '
           'the surface - then the porosity (or cell size) and strut thickness. '
           'The result is a watertight mesh body; its porosity is measured, not '
           'assumed.')

PANELS = ['SolidCreatePanel', 'ParaMeshModifyPanel']


def install(app, ui, created_handler_factory):
    old = ui.commandDefinitions.itemById(CMD_ID)
    if old:
        old.deleteMe()
    folder = RES_DIR if os.path.isdir(RES_DIR) and os.listdir(RES_DIR) else ''
    cmd_def = ui.commandDefinitions.addButtonDefinition(CMD_ID, CMD_NAME, TOOLTIP,
                                                        folder)
    cmd_def.commandCreated.add(created_handler_factory())
    placed = []
    for pid in PANELS:
        panel = ui.allToolbarPanels.itemById(pid)
        if panel and not panel.controls.itemById(CMD_ID):
            panel.controls.addCommand(cmd_def)
            placed.append(pid)
    return placed


def uninstall(ui):
    for pid in PANELS:
        panel = ui.allToolbarPanels.itemById(pid)
        if panel:
            c = panel.controls.itemById(CMD_ID)
            if c:
                c.deleteMe()
    d = ui.commandDefinitions.itemById(CMD_ID)
    if d:
        d.deleteMe()


# --------------------------------------------------------------------------- #
# dev bridge - everything that builds runs in a NEW document, never the user's
# --------------------------------------------------------------------------- #

def dev_ping(app, job):
    doc = app.activeDocument
    return {'version': VERSION, 'document': doc.name if doc else None,
            'python': python_cmd()}


def _in_scratch(app, job, fn):
    prev = app.activeDocument
    doc = app.documents.add(adsk.core.DocumentTypes.FusionDesignDocumentType)
    try:
        des = adsk.fusion.Design.cast(doc.products.itemByProductType('DesignProductType'))
        if job.get('direct'):
            des.designType = adsk.fusion.DesignTypes.DirectDesignType
        return fn(des)
    finally:
        if not job.get('keep'):
            doc.close(False)
            if prev:
                prev.activate()


def _add_temp(des, tmp, name, comp=None):
    comp = comp or des.rootComponent
    parametric = des.designType == adsk.fusion.DesignTypes.ParametricDesignType
    if parametric:
        base = comp.features.baseFeatures.add()
        base.startEdit()
        comp.bRepBodies.add(tmp, base)
        base.finishEdit()
    else:
        comp.bRepBodies.add(tmp)
    body = comp.bRepBodies.item(comp.bRepBodies.count - 1)
    body.name = name
    return body


def _shape(des, shape):
    """Test regions, sizes in mm: a 20 mm box, a 20 mm ball, the top face of
    the box, an open surface body (that face on its own), a box inside a moved
    occurrence."""
    tbm = adsk.fusion.TemporaryBRepManager.get()
    if shape == 'ball':
        tmp = tbm.createSphere(adsk.core.Point3D.create(0, 0, 0), 1.0)
        return [_add_temp(des, tmp, 'ball')]
    box = adsk.core.OrientedBoundingBox3D.create(
        adsk.core.Point3D.create(0, 0, 0), adsk.core.Vector3D.create(1, 0, 0),
        adsk.core.Vector3D.create(0, 1, 0), 2.0, 2.0, 2.0)
    if shape == 'occ_box':
        m = adsk.core.Matrix3D.create()
        m.translation = adsk.core.Vector3D.create(5.0, 0, 0)
        occ = des.rootComponent.occurrences.addNewComponent(m)
        body = _add_temp(des, tbm.createBox(box), 'box', occ.component)
        return [body.createForAssemblyContext(occ)]
    body = _add_temp(des, tbm.createBox(box), 'box')
    if shape == 'box':
        return [body]
    top = max(body.faces, key=lambda f: f.centroid.z)
    if shape == 'face':
        return [top]
    if shape == 'faces_all':
        return list(body.faces)
    if shape == 'sheet':
        sheet = _add_temp(des, tbm.copy(top), 'sheet')
        return [sheet]
    raise ValueError('unknown shape %r' % shape)


def dev_region(app, job):
    """What the dialog would say about a test region."""
    def run(des):
        reg = read_region(_shape(des, job['shape']), job.get('whole', False))
        out = {k: reg[k] for k in ('closed', 'what', 'area', 'volume')}
        out['faces_of_one_solid'] = bool(reg['faces_of_one_solid'])
        return out
    return _in_scratch(app, job, run)


def _bbox_mm(b):
    bb = b.boundingBox
    return [[bb.minPoint.x * CM, bb.minPoint.y * CM, bb.minPoint.z * CM],
            [bb.maxPoint.x * CM, bb.maxPoint.y * CM, bb.maxPoint.z * CM]]


def dev_build(app, job):
    """Build a lattice on a test region through the same build() the dialog
    uses, and report what landed in the design."""
    def run(des):
        ents = _shape(des, job['shape'])
        reg = read_region(ents, job.get('whole', False))
        p = dict(kind='volume' if reg['closed'] else 'surface', placement='flush',
                 flip=False, skin=True, border=True, define='porosity',
                 porosity=0.8, cell=3.0, thickness=0.6, regularity=5, seed=1,
                 resolution='coarse', max_tris=400000, closed=reg['closed'],
                 organic=1, output='mesh')
        p.update(job.get('params', {}))
        graph = run_engine(build_job(reg, p, 'graph'))
        t0 = time.time()
        body, res = build(app, None, reg, p)
        mb = adsk.fusion.MeshBody.cast(body)
        src = reg['body']
        out = {'region': {k: reg[k] for k in ('closed', 'what', 'area', 'volume')},
               'preview_segments': len(graph['segments']),
               'preview_output': graph['output'],
               'graph_cached_on_build': bool(res['graph'].get('cached')),
               'porosity_measured': res['porosity_measured'],
               'output': res['output'], 'build_s': time.time() - t0,
               'mesh': res.get('mesh'), 'form': res.get('form'),
               'graph': res['graph'], 'text': report(p, res),
               'form_skipped': res.get('form_skipped')}
        if res['output'] == 'form':
            bb = adsk.fusion.BRepBody.cast(body)
            out['body'] = {'name': bb.name, 'solid': bb.isSolid,
                           'faces': bb.faces.count,
                           'volume_mm3': bb.volume * 1000.0,
                           'bbox_mm': _bbox_mm(bb),
                           'component': bb.parentComponent.name,
                           'form_features': sum(c.features.formFeatures.count
                                                for c in des.allComponents)}
            out['source_bbox_mm'] = _bbox_mm(src)
        else:
            out['body'] = {'name': mb.name, 'closed': mb.isClosed,
                           'triangles': mb.mesh.triangleCount,
                           'volume_mm3': mesh_volume_mm3(mb),
                           'fusion_volume_mm3': mb.volume * 1000.0,
                           'bbox_mm': _bbox_mm(mb),
                           'component': mb.parentComponent.name,
                           'in_base_feature': owning_base(mb) is not None}
        parametric = des.designType == adsk.fusion.DesignTypes.ParametricDesignType
        out['design'] = {'parametric': parametric,
                         'base_features': sum(c.features.baseFeatures.count
                                              for c in des.allComponents)
                         if parametric else 0}
        out['source_bbox_mm'] = _bbox_mm(src)
        if job.get('replace_check'):
            def count():
                n = sum(c.meshBodies.count for c in des.allComponents)
                if parametric:
                    n += sum(c.features.formFeatures.count for c in des.allComponents)
                return n
            n0 = count()
            gone = remove_previous(des)
            out['replace'] = {'removed': gone, 'mesh_bodies_before': n0,
                              'mesh_bodies_after': count()}
        return out
    return _in_scratch(app, job, run)


def dev_tsm(app, job):
    """Load a TSM file as a Form body in a scratch document; time it and report
    what Fusion made of it. For measuring how large a cage Fusion will take."""
    def run(des):
        if des.designType != adsk.fusion.DesignTypes.ParametricDesignType:
            des.designType = adsk.fusion.DesignTypes.ParametricDesignType
        with open(job['path'], encoding='utf-8') as fh:
            text = fh.read()
        comp = des.rootComponent
        t0 = time.time()
        ff = comp.features.formFeatures.add()
        ff.startEdit()
        out = {}
        try:
            tb = ff.tSplineBodies.addByTSMDescription(text)
            out['tspline'] = bool(tb)
        except Exception:
            out['error'] = short_error(traceback.format_exc())
        finally:
            t1 = time.time()
            ff.finishEdit()
        t2 = time.time()
        out['add_s'] = round(t1 - t0, 2)
        out['finish_s'] = round(t2 - t1, 2)
        out['bodies'] = []
        for i in range(comp.bRepBodies.count):
            b = comp.bRepBodies.item(i)
            out['bodies'].append({'name': b.name, 'solid': b.isSolid,
                                  'faces': b.faces.count,
                                  'volume_mm3': b.volume * 1000.0 if b.isSolid else None})
        out['total_s'] = round(time.time() - t0, 2)
        return out
    return _in_scratch(app, job, run)


DEV_COMMANDS = {
    'tsm': dev_tsm,
    'ping': dev_ping,
    'region': dev_region,
    'build': dev_build,
}
