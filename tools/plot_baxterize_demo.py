#!/usr/bin/env python3
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def plot_baxter(npz_path):
    data = np.load(npz_path, allow_pickle=True)
    Ueff = data.get('Ueff')
    if Ueff is None:
        Ueff = data.get('Ueff', None)
    if Ueff is None:
        print('No Ueff found in', npz_path)
        return

    outdir = Path('results/plots')
    outdir.mkdir(parents=True, exist_ok=True)

    # heatmap of absolute values
    plt.figure()
    plt.imshow(np.abs(Ueff), cmap='viridis', interpolation='nearest')
    plt.colorbar(label='|U_eff|')
    plt.title('Absolute value of U_eff (2x2)')
    plt.xticks([0,1])
    plt.yticks([0,1])
    plt.tight_layout()
    p1 = outdir / (Path(npz_path).stem + '_Ueff_abs.png')
    plt.savefig(p1, dpi=150)
    plt.close()

    # real and imag bar plot
    plt.figure(figsize=(6,4))
    labels = ['00','01','10','11']
    vals_real = [Ueff[0,0].real, Ueff[0,1].real, Ueff[1,0].real, Ueff[1,1].real]
    vals_imag = [Ueff[0,0].imag, Ueff[0,1].imag, Ueff[1,0].imag, Ueff[1,1].imag]
    x = np.arange(len(labels))
    width = 0.35
    plt.bar(x - width/2, vals_real, width, label='Real')
    plt.bar(x + width/2, vals_imag, width, label='Imag')
    plt.xticks(x, labels)
    plt.ylabel('Component value')
    plt.title('U_eff components (real and imag)')
    plt.legend()
    plt.tight_layout()
    p2 = outdir / (Path(npz_path).stem + '_Ueff_components.png')
    plt.savefig(p2, dpi=150)
    plt.close()

    # phase of determinant and fidelity if available
    det = np.linalg.det(Ueff)
    try:
        from math import atan2
        phase = np.angle(det)
    except Exception:
        phase = None

    # small text summary
    summary = outdir / (Path(npz_path).stem + '_summary.txt')
    with open(summary, 'w') as f:
        f.write(f'Ueff:\n{Ueff}\n')
        f.write(f'det(Ueff) = {det}\n')
        if phase is not None:
            f.write(f'phase(det) = {phase}\n')

    print('Saved plots to', outdir)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: plot_baxterize_demo.py path/to/baxterize_demo_...npz')
        sys.exit(1)
    plot_baxter(sys.argv[1])
