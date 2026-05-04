#!/usr/bin/env python3
"""Assemble existing result PNGs into paper-style figure panels for direct comparison.

This script reads the PNGs produced by `tools/reproduce_figs.py` and arranges
them into composite figures saved under `results/` as `paper_style_Fig1.png`, etc.
"""
import os
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

OUT = 'results'
os.makedirs(OUT, exist_ok=True)


def load_img(path):
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    return mpimg.imread(path)


def fig1_panel():
    src = os.path.join(OUT, 'reproduce_Fig1_ldos_panel.png')
    img = load_img(src)
    fig, ax = plt.subplots(1,1, figsize=(8,3))
    ax.imshow(img)
    ax.axis('off')
    out = os.path.join(OUT, 'paper_style_Fig1.png')
    plt.tight_layout()
    plt.savefig(out, dpi=300)
    plt.close()
    print('Saved', out)


def fig2_panel():
    Ts = [400, 450, 500]
    fig, axs = plt.subplots(2, 3, figsize=(12, 6))
    for j, T in enumerate(Ts):
        mfile = os.path.join(OUT, f'reproduce_tetron_MZM_T{T}_delta0.0.png')
        afile = os.path.join(OUT, f'reproduce_tetron_ABS_T{T}_delta0.2.png')
        axs[0, j].imshow(load_img(mfile))
        axs[0, j].axis('off')
        axs[0, j].set_title(f'MZM T={T}')
        axs[1, j].imshow(load_img(afile))
        axs[1, j].axis('off')
        axs[1, j].set_title(f'ABS T={T}')
    out = os.path.join(OUT, 'paper_style_Fig2.png')
    plt.tight_layout()
    plt.savefig(out, dpi=300)
    plt.close()
    print('Saved', out)


def fig3_panel():
    a = os.path.join(OUT, 'reproduce_Fig3_abs_eigen_vs_time.png')
    b = os.path.join(OUT, 'reproduce_Fig3_overlap_vs_T.png')
    fig, axs = plt.subplots(1, 2, figsize=(10,4))
    axs[0].imshow(load_img(a))
    axs[0].axis('off')
    axs[0].set_title('ABS eigenenergy vs time')
    axs[1].imshow(load_img(b))
    axs[1].axis('off')
    axs[1].set_title('Final overlap vs T')
    out = os.path.join(OUT, 'paper_style_Fig3.png')
    plt.tight_layout()
    plt.savefig(out, dpi=300)
    plt.close()
    print('Saved', out)


def fig4_panel():
    a = os.path.join(OUT, 'reproduce_Fig4_modulated_eigs.png')
    fig, ax = plt.subplots(1,1,figsize=(6,4))
    ax.imshow(load_img(a))
    ax.axis('off')
    ax.set_title('Fig.4 modulated eigs')
    out = os.path.join(OUT, 'paper_style_Fig4.png')
    plt.tight_layout()
    plt.savefig(out, dpi=300)
    plt.close()
    print('Saved', out)


def fig5_panel():
    a = os.path.join(OUT, 'reproduce_Fig5_modulation_amplitude.png')
    fig, ax = plt.subplots(1,1,figsize=(6,4))
    ax.imshow(load_img(a))
    ax.axis('off')
    ax.set_title('Fig.5 modulation amplitude')
    out = os.path.join(OUT, 'paper_style_Fig5.png')
    plt.tight_layout()
    plt.savefig(out, dpi=300)
    plt.close()
    print('Saved', out)


def main():
    fig1_panel()
    fig2_panel()
    fig3_panel()
    fig4_panel()
    fig5_panel()


if __name__ == '__main__':
    main()
