#!/usr/bin/env python3
"""Compute peak index and simple 1/e decay estimate from a density file.

Usage: python3 tools/compute_density_xi.py results/xyz_demo_density0.txt
"""
import sys
import numpy as np


def compute(fname):
    arr = np.loadtxt(fname)
    N = arr.size
    peak = int(np.argmax(arr))
    peak_val = float(arr[peak])
    right = peak
    while right < N and arr[right] >= peak_val/np.e:
        right += 1
    xi_est = (right - peak) if right < N else float('nan')
    print(f'file={fname}')
    print(f'peak_idx={peak}, peak_val={peak_val:.6e}, xi_est={xi_est}')
    return peak, peak_val, xi_est


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: python3 tools/compute_density_xi.py <density-file>')
        sys.exit(1)
    compute(sys.argv[1])
