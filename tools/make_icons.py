"""Toolbar icons: a patch of Voronoi struts, drawn with Pillow at 8x and
downsampled, in the 18 files Fusion looks for.

    py -3.10 tools/make_icons.py

A NEW icon needs a Fusion restart before the toolbar shows it.
"""

import os

import numpy as np
from PIL import Image, ImageDraw
from scipy.spatial import Voronoi

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, 'resources', 'VoronoiLattice')
SS = 8

LIGHT = (70, 70, 70, 255)
DARK = (225, 225, 225, 255)
ACCENT = (230, 120, 30, 255)
DISABLED = (160, 160, 160, 255)


def seeds():
    """A jittered 4 x 4 grid: few enough cells to read at 16 px."""
    rng = np.random.default_rng(3)
    g = np.stack(np.meshgrid(np.linspace(0.08, 0.92, 4), np.linspace(0.08, 0.92, 4)), -1)
    g = g.reshape(-1, 2) + rng.normal(0, 0.06, (16, 2))
    return g


def edges(pts):
    ghost = np.array([[x, y] for x in (-1, 0.5, 2) for y in (-1, 0.5, 2) if (x, y) != (0.5, 0.5)])
    vor = Voronoi(np.vstack([pts, ghost]))
    out = []
    for rv in vor.ridge_vertices:
        if -1 in rv:
            continue
        a, b = vor.vertices[rv[0]], vor.vertices[rv[1]]
        out.append((a, b))
    return out


def draw(size, color, accent):
    n = size * SS
    img = Image.new('RGBA', (n, n), (0, 0, 0, 0))
    d = ImageDraw.Draw(img)
    pad = 0.06 * n
    w = max(1, int(round(n * (0.085 if size <= 16 else 0.065))))
    clip = (pad, pad, n - pad, n - pad)
    d.rounded_rectangle(clip, radius=0.16 * n, outline=color, width=w)
    mask = Image.new('L', (n, n), 0)
    ImageDraw.Draw(mask).rounded_rectangle(clip, radius=0.16 * n, fill=255)
    net = Image.new('RGBA', (n, n), (0, 0, 0, 0))
    dn = ImageDraw.Draw(net)
    for a, b in edges(seeds()):
        pa = pad + a * (n - 2 * pad)
        pb = pad + b * (n - 2 * pad)
        dn.line([tuple(pa), tuple(pb)], fill=accent, width=w)
        for p in (pa, pb):
            r = w * 0.5
            dn.ellipse([p[0] - r, p[1] - r, p[0] + r, p[1] + r], fill=accent)
    img.paste(net, (0, 0), Image.composite(net, Image.new('RGBA', (n, n)), mask))
    return img.resize((size, size), Image.LANCZOS)


def main():
    os.makedirs(OUT, exist_ok=True)
    for base in (16, 32, 64):
        for scale, suffix in ((1, ''), (2, '@2x')):
            px = base * scale
            for tag, col, acc in (('', LIGHT, ACCENT), ('-dark', DARK, ACCENT),
                                  ('-disabled', DISABLED, DISABLED)):
                draw(px, col, acc).save(os.path.join(OUT, '%dx%d%s%s.png'
                                                     % (base, base, tag, suffix)))
    print('wrote', len(os.listdir(OUT)), 'icons to', OUT)


if __name__ == '__main__':
    main()
