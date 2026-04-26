#!/usr/bin/env python3
"""Plot LDOS (density) from a saved density file.
Usage: python3 tools/plot_ldos.py <density-file> [out-prefix]
"""
import sys, os
import numpy as np

def plot_density(fname, out_prefix=None):
    arr = np.loadtxt(fname)
    N = arr.size
    if out_prefix is None:
        base = os.path.splitext(os.path.basename(fname))[0]
    else:
        base = out_prefix
    try:
        import matplotlib.pyplot as plt
    except Exception:
        print('matplotlib not available')
        return
    xs = np.arange(N)
    plt.figure(figsize=(8,3))
    plt.plot(xs, arr, '-o', markersize=3)
    plt.xlabel('site')
    plt.ylabel('LDOS (density)')
    plt.title(base + ' LDOS')
    plt.tight_layout()
    out = f'results/{base}_ldos.png'
    plt.savefig(out)
    plt.close()

    # log plot
    plt.figure(figsize=(8,3))
    plt.semilogy(xs, arr + 1e-40, '-o', markersize=3)
    plt.xlabel('site')
    plt.ylabel('LDOS (log)')
    plt.title(base + ' LDOS (log)')
    plt.tight_layout()
    out2 = f'results/{base}_ldos_log.png'
    plt.savefig(out2)
    plt.close()
    print('Saved', out, out2)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: python3 tools/plot_ldos.py <density-file> [out-prefix]')
        sys.exit(1)
    fname = sys.argv[1]
    prefix = sys.argv[2] if len(sys.argv) > 2 else None
    plot_density(fname, prefix)
