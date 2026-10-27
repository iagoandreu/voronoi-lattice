"""Dev-bridge client: send a job to the running Voronoi Lattice add-in, print the answer.

    py -3.10 tools/vl.py ping
    py -3.10 tools/vl.py region '{"shape": "faces_all"}'
    py -3.10 tools/vl.py build '{"shape": "ball", "params": {"kind": "surface"}}'

region and build run in a NEW document, closed unsaved afterwards ("keep": true
leaves it open); the active document is never touched.

Runs outside Fusion; the add-in's poller picks the job up and answers in devjobs/outbox.
"""

import json
import os
import sys
import time

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ADDIN = os.path.join(os.environ['APPDATA'],
                     r'Autodesk\Autodesk Fusion 360\API\AddIns\VoronoiLattice')
IN = os.path.join(ADDIN, 'devjobs', 'inbox')
OUT = os.path.join(ADDIN, 'devjobs', 'outbox')
BEAT = os.path.join(ADDIN, 'devjobs', 'heartbeat.json')


def alive():
    try:
        return time.time() - os.path.getmtime(BEAT) < 5.0
    except Exception:
        return False


def send(cmd, payload=None, timeout=900.0):
    os.makedirs(IN, exist_ok=True)
    os.makedirs(OUT, exist_ok=True)
    job = {'cmd': cmd}
    if payload:
        job.update(payload)
    jid = '%s_%d' % (cmd, int(time.time() * 1000))
    job['job_id'] = jid
    tmp = os.path.join(IN, jid + '.json.part')
    with open(tmp, 'w', encoding='utf-8') as fh:
        json.dump(job, fh)
    os.replace(tmp, os.path.join(IN, jid + '.json'))

    dst = os.path.join(OUT, jid + '.json')
    t0 = time.time()
    while time.time() - t0 < timeout:
        if os.path.exists(dst):
            with open(dst, encoding='utf-8') as fh:
                return json.load(fh)
        time.sleep(0.3)
    return {'ok': False, 'error': 'timeout after %.0fs (bridge alive: %s)'
            % (timeout, alive())}


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(__doc__)
        raise SystemExit(2)
    extra = json.loads(sys.argv[2]) if len(sys.argv) > 2 else None
    res = send(sys.argv[1], extra)
    text = res.get('result', {}).pop('text', None) if isinstance(res.get('result'), dict) else None
    print(json.dumps(res, indent=2, default=str))
    if text:
        print('\n' + text)
