#!/usr/bin/env python3
import sys
from pathlib import Path
import numpy as np


def apply_phases(npz_in, a, b, out_suffix='_compensated'):
    p = Path(npz_in)
    data = np.load(p, allow_pickle=True)
    U = data['U'] if 'U' in data else None
    Vlow = data['Vlow0'] if 'Vlow0' in data else None
    Ueff = data['Ueff'] if 'Ueff' in data else None
    if U is None or Vlow is None:
        raise RuntimeError('Input npz must contain U and Vlow0')

    D = np.diag([np.exp(1j*a), np.exp(1j*b)])
    # lift D
    Vlift = Vlow @ D @ Vlow.conj().T
    U_new = Vlift @ U
    Ueff_new = Vlow.conj().T @ U_new @ Vlow

    outp = p.with_name(p.stem + out_suffix + p.suffix)
    np.savez(outp, U=U_new, Ueff=Ueff_new, Vlow0=Vlow)
    return outp


def main():
    if len(sys.argv) < 4:
        print('Usage: apply_phase_and_plot.py path/to/npz a b')
        sys.exit(1)
    npz = sys.argv[1]
    a = float(sys.argv[2])
    b = float(sys.argv[3])
    outp = apply_phases(npz, a, b)
    print('Saved compensated results to', outp)
    # call plotting script
    import subprocess
    subprocess.run(['python3', 'tools/plot_baxterize_demo.py', str(outp)])


if __name__ == '__main__':
    main()
