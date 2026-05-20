#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.pyplot as plt

ROOT = os.path.dirname(__file__)
csv_path = os.path.join(ROOT, 'pauli_schur_compare.csv')

arr = np.genfromtxt(csv_path, delimiter=',', names=True)
e1 = arr['E1']

def save_plot(x, ys, labels, fname, ylog=False):
    plt.figure(figsize=(6,4))
    for y,lab in zip(ys, labels):
        plt.plot(x, y, '-o', label=lab)
    plt.xlabel('E1')
    plt.legend()
    if ylog:
        plt.yscale('log')
    plt.tight_layout()
    plt.savefig(os.path.join(ROOT, fname), dpi=150)
    plt.close()

# Plot K_fro and norm_fro on log scale
save_plot(e1, [arr['K_fro'], arr['norm_fro']], ['K_fro','norm_fro'], 'K_and_norm_fro.png', ylog=True)

# Plot norm_pauli and norm_fro linear
save_plot(e1, [arr['norm_pauli'], arr['norm_fro']], ['norm_pauli','norm_fro'], 'norms_compare.png')

# Plot N_pauli vs N_fro_scaled
save_plot(e1, [arr['N_pauli'], arr['N_fro_scaled']], ['N_pauli','N_fro_scaled'], 'N_pauli_vs_N_fro_scaled.png')

# Plot opdiff and scale
save_plot(e1, [arr['opdiff'], arr['scale']], ['opdiff','scale'], 'opdiff_and_scale.png')

# Save simple stats
out = os.path.join(ROOT, 'pauli_schur_diagnostics.txt')
with open(out, 'w') as f:
    # compute Pearson correlation if enough points
    try:
        corr = np.corrcoef(arr['N_pauli'], arr['N_fro_scaled'])[0,1]
    except Exception:
        corr = np.nan
    f.write('Pearson corr N_pauli vs N_fro_scaled: %g\n' % corr)
    f.write('mean scale: %g, std scale: %g\n' % (np.mean(arr['scale']), np.std(arr['scale'])))
    f.write('mean opdiff: %g, max opdiff: %g\n' % (np.mean(arr['opdiff']), np.max(arr['opdiff'])))

print('Saved diagnostics to', ROOT)

if __name__ == '__main__':
    pass
