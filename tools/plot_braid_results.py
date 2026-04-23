#!/usr/bin/env python3
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def load_and_plot(npz_path):
    data = np.load(npz_path, allow_pickle=True)
    times = data.get('times')
    gaps = data.get('gaps')
    ldos_mid = data.get('ldos_mid')
    Ueff_list = data.get('Ueff_list')
    Ueff_final = data.get('Ueff_final')
    Vlow0 = data.get('Vlow0')

    outdir = Path('results/plots')
    outdir.mkdir(parents=True, exist_ok=True)

    # plot LDOS sampled at snapshot times
    # Ensure `times` is defined even if `Ueff_list` is missing.
    if times is None:
        if Ueff_list is not None:
            times = np.arange(Ueff_list.shape[0])
        else:
            # single-snapshot case: use a single time value 0
            times = np.array([0])

    # sample ldos to match times length
    if ldos_mid is not None and len(ldos_mid) > 0:
        idxs = (np.linspace(0, len(ldos_mid)-1, len(times))).round().astype(int)
        ldos_sampled = ldos_mid[idxs]
    else:
        ldos_sampled = None

    # Fidelity per snapshot
    def ideal_braid_2x2():
        sy = np.array([[0, -1j],[1j, 0]], dtype=complex)
        B = np.cos(np.pi/4.0) * np.eye(2, dtype=complex) + 1j * np.sin(np.pi/4.0) * sy
        return B

    B2 = ideal_braid_2x2()
    fid_list = []
    frob_list = []
    single_snapshot = False
    if Ueff_list is None:
        if Ueff_final is None:
            raise RuntimeError('No Ueff_list or Ueff_final found in npz; cannot compute fidelity time series')
        # single snapshot case
        single_snapshot = True
        Ueff = Ueff_final
        fid = float(np.real(np.trace(Ueff.conj().T @ B2)) / 2.0)
        frob = float(np.linalg.norm(Ueff - B2))
        fid_list = [fid]
        frob_list = [frob]
    else:
        for Ueff in Ueff_list:
            fid = np.real(np.trace(Ueff.conj().T @ B2)) / 2.0
            frob = np.linalg.norm(Ueff - B2)
            fid_list.append(fid)
            frob_list.append(frob)

    # LDOS plot
    if ldos_sampled is not None:
        plt.figure()
        plt.plot(times, ldos_sampled, '-o')
        plt.xlabel('time')
        plt.ylabel('LDOS(mid)')
        plt.title('Mid-site LDOS (sampled)')
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(outdir / 'ldos_mid.png', dpi=150)
        plt.close()

    # gap plot (per step)
    if gaps is not None:
        plt.figure()
        plt.plot(np.arange(len(gaps)), gaps)
        plt.xlabel('step')
        plt.ylabel('instantaneous gap')
        plt.title('Instantaneous gap vs step')
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(outdir / 'gaps.png', dpi=150)
        plt.close()

    # fidelity and frob over snapshots
    plt.figure()
    if single_snapshot:
        plt.plot([times[0] if hasattr(times, '__len__') and len(times)>0 else 0], fid_list, 'o', label='fid_like (single)')
        plt.plot([times[0] if hasattr(times, '__len__') and len(times)>0 else 0], frob_list, 'x', label='frob_diff (single)')
        plt.title('Projected fidelity / frob (single snapshot)')
    else:
        plt.plot(times, fid_list, label='fid_like')
        plt.plot(times, frob_list, label='frob_diff')
        plt.title('Projected fidelity / frob vs time')
    plt.xlabel('time')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(outdir / 'fidelity_vs_time.png', dpi=150)
    plt.close()

    # plot Ueff matrix components (real/imag) if we have a time series
    if Ueff_list is None:
        print('Warning: no Ueff_list found in file; skipping time-series component plot. Saved single-snapshot summaries instead.')
        # save heatmaps of abs and phase for the final Ueff
        U = Ueff_final
        plt.figure(figsize=(6,3))
        plt.subplot(1,2,1)
        plt.imshow(np.abs(U), vmin=0, vmax=np.abs(U).max())
        plt.colorbar()
        plt.title('abs(U_eff_final)')
        plt.subplot(1,2,2)
        plt.imshow(np.angle(U))
        plt.colorbar()
        plt.title('arg(U_eff_final)')
        plt.tight_layout()
        plt.savefig(outdir / 'Ueff_final_abs_phase.png', dpi=150)
        plt.close()
    else:
        Ueff_arr = np.array(Ueff_list)
        labels = ['00','01','10','11']
        plt.figure(figsize=(8,6))
        for i in range(2):
            for j in range(2):
                comp = Ueff_arr[:, i, j]
                plt.plot(times, comp.real, label=f'Re({i}{j})')
                plt.plot(times, comp.imag, '--', label=f'Im({i}{j})')
        plt.xlabel('time')
        plt.title('U_eff components (real solid / imag dashed)')
        plt.legend(ncol=2, fontsize='small')
        plt.tight_layout()
        plt.savefig(outdir / 'Ueff_components.png', dpi=150)
        plt.close()

    # final Ueff info
    np.savez(outdir / 'Ueff_final_summary.npz', Ueff_final=Ueff_final)
    print('Plots saved to', outdir)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: plot_braid_results.py path/to/braid_sim_record_...npz')
        sys.exit(1)
    load_and_plot(sys.argv[1])
