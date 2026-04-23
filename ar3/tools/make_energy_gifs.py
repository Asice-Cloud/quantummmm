#!/usr/bin/env python3
"""Create GIFs from energy surface PNGs in results/.

It looks for files named: energy_E{band}_czz_{value:.3f}.png
and creates anim_energy_E{band}.gif cycling over czz values.
"""
from pathlib import Path
import re
import imageio.v2 as imageio

HERE = Path(__file__).resolve().parent.parent
RESULTS = HERE / 'results'

def collect():
    pattern = re.compile(r'energy_E([0-3])_czz_([-0-9.]+)\.png')
    imgs = {}
    for p in sorted(RESULTS.iterdir()):
        m = pattern.match(p.name)
        if not m:
            continue
        band = int(m.group(1))
        czz = float(m.group(2))
        imgs.setdefault(band, []).append((czz, p))
    for band in imgs:
        imgs[band].sort()
    return imgs

def make_gifs(imgs):
    for band, entries in imgs.items():
        frames = [str(p) for _, p in entries]
        out = RESULTS / f'anim_energy_E{band}.gif'
        images = [imageio.imread(fn) for fn in frames]
        # slower playback: increase duration per frame
        imageio.mimsave(out, images, duration=1.8)
        print('Wrote', out)

def main():
    imgs = collect()
    if not imgs:
        print('No energy images found in results/')
        return
    make_gifs(imgs)

if __name__ == '__main__':
    main()
