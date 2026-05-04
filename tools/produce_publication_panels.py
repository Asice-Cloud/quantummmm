#!/usr/bin/env python3
"""Create publication-quality PNG panels from existing result images.

This script searches for existing figure images in `results/`, overlays
mapping-parameter annotations (A0,B0,C0), and re-saves high-resolution
PNG files `results/pub_Fig{n}.png` for n=1..5.
"""
import os
import glob
import numpy as np
import matplotlib.pyplot as plt


def find_image_for_fig(n):
    # look for common filenames containing Fign or Fig{n}
    patterns = [f'results/**/*Fig{n}*.png', f'results/*Fig{n}*.png', f'results/**/*fig{n}*.png', f'results/*fig{n}*.png']
    for pat in patterns:
        matches = glob.glob(pat, recursive=True)
        if matches:
            return matches[0]
    # fallback: any image containing the digit n
    matches = glob.glob(f'results/**/*{n}*.png', recursive=True)
    return matches[0] if matches else None


def main():
    os.makedirs('results', exist_ok=True)
    # load mapping constants if available
    mapping_file = 'results/mapping_fit_ABC.npz'
    if os.path.exists(mapping_file):
        d = np.load(mapping_file)
        A0 = float(d.get('A0_fit', np.nan))
        B0 = float(d.get('B0_fit', np.nan))
        C0 = float(d.get('C0_fit', np.nan))
    else:
        A0 = B0 = C0 = np.nan

    saved = []
    for n in range(1, 6):
        src = find_image_for_fig(n)
        if src is None:
            print(f'No source image found for Fig{n}, skipping.')
            continue
        img = plt.imread(src)
        fig = plt.figure(figsize=(6, 4), dpi=300)
        ax = fig.add_subplot(111)
        ax.imshow(img)
        ax.axis('off')
        # annotation box
        ann = f'A0={A0:.6g}\nB0={B0:.6g}\nC0={C0:.6g}' if not np.isnan(A0) else 'mapping: N/A'
        ax.text(0.01, 0.98, ann, transform=ax.transAxes, fontsize=8,
                verticalalignment='top', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
        out = f'results/pub_Fig{n}.png'
        plt.tight_layout(pad=0)
        fig.savefig(out, dpi=300)
        plt.close(fig)
        saved.append(out)
        print('Saved', out, 'from', src)

    if not saved:
        print('No figures created.')
    else:
        print('Created publication panels:', saved)


if __name__ == '__main__':
    main()
