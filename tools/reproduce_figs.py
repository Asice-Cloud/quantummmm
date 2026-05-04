#!/usr/bin/env python3
"""Driver to reproduce key figures from Chen et al. (PRB 105, 054507).

This script uses existing tools in `tools/` to generate figure panels
corresponding to Fig.1..Fig.5 in the paper. It intentionally reuses the
simplified Kitaev-chain and two-level path models in this repository.
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# ensure repo root importable
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from tools import embed_kitaev
from tools import tetron_path_sim as tetron
from tools import bloch_rotation_sim as bloch
import tools.paper_params as P


def map_gates_to_links(g1, g2, g3, g4, t0, Delta0, L, mu0, VD, qd_width):
    mu = mu0 * np.ones(L)
    for i in range(qd_width):
        mu[i] = mu[i] - VD
    t_links = t0 * np.ones(L - 1)
    t_links_mod = t_links.copy()
    t_links_mod[0] = t0 * g1
    if L > 2:
        t_links_mod[1] = t0 * g3
    Delta_links = Delta0 * np.ones(L - 1)
    Delta_mod = Delta_links.copy()
    Delta_mod[0] = Delta0 * (g1 if g1 > 0 else 1e-3)
    if L > 2:
        Delta_mod[1] = Delta0 * (g3 if g3 > 0 else 1e-3)
    return mu, t_links_mod, Delta_mod


def reproduce_fig1():
    print('Reproducing Fig.1: LDOS snapshots')
    os.makedirs(P.OUTDIR, exist_ok=True)
    # call existing snapshot routine with paper params
    embed_kitaev.snapshot_ldos(L=P.L, t0=P.t0, Delta0=P.DELTA, mu0=P.mu0, VD=P.VD, qd_width=P.QD_WIDTH)
    # combine saved snapshots into a single figure for easier comparison
    imgs = [os.path.join(P.OUTDIR, f'ldos_snapshot_{name}.png') for name in ('init', 'after_step1', 'after_step2', 'after_step3')]
    fig, axs = plt.subplots(1, 4, figsize=(12,3))
    for ax, img in zip(axs, imgs):
        im = plt.imread(img)
        ax.imshow(im)
        ax.axis('off')
    out = os.path.join(P.OUTDIR, 'reproduce_Fig1_ldos_panel.png')
    plt.tight_layout()
    plt.savefig(out)
    plt.close()
    print('Saved', out)


def reproduce_fig2():
    print('Reproducing Fig.2: wavefunction evolution / tetron path')
    os.makedirs(P.OUTDIR, exist_ok=True)
    for T in P.FIG2_TS:
        # MZM case
        tetron.run_sim(T_step=T, n_per_step=400, delta=0.0, save_prefix=os.path.join(P.OUTDIR, 'reproduce_tetron_MZM'))
        # ABS case (use delta ~ 0.2 as representative)
        tetron.run_sim(T_step=T, n_per_step=400, delta=0.2, save_prefix=os.path.join(P.OUTDIR, 'reproduce_tetron_ABS'))
    print('Fig.2 simulations saved to results/')


def reproduce_fig3():
    print('Reproducing Fig.3: ABS eigenenergy during braiding and overlap vs T')
    os.makedirs(P.OUTDIR, exist_ok=True)
    T = P.FIG3_T
    tlist, step_idx, slist, dt = tetron.make_time_grid(T_step=T, n_per_step=P.FIG3_N_PER_STEP)
    Emin = np.zeros_like(tlist)
    for i in range(len(tlist)):
        step = int(step_idx[i])
        s = float(slist[i])
        g1, g2, g3, g4 = tetron.gates_at(step, s)
        mu, t_links_mod, Delta_mod = map_gates_to_links(g1, g2, g3, g4, P.t0, P.DELTA, P.L, P.mu0, P.VD, P.QD_WIDTH)
        H = embed_kitaev.build_bdg(mu, t_links_mod, Delta_mod)
        E = np.linalg.eigvalsh(H)
        # pick smallest absolute-energy eigenvalue (closest to zero)
        Emin[i] = np.min(np.abs(E))

    plt.figure()
    plt.plot(tlist / T, Emin)
    plt.xlabel('normalized time (t/T)')
    plt.ylabel('min |E| (BdG)')
    out = os.path.join(P.OUTDIR, 'reproduce_Fig3_abs_eigen_vs_time.png')
    plt.savefig(out)
    plt.close()
    print('Saved', out)

    # overlap vs T (scan a range similar to paper)
    Tscan = np.linspace(50, 1000, 20)
    overlaps = []
    for TT in Tscan:
        tetron.run_sim(T_step=TT, n_per_step=300, delta=0.2, save_prefix=os.path.join(P.OUTDIR, 'tmp_tetron_ABS'))
        data = np.load(os.path.join(P.OUTDIR, f'tmp_tetron_ABS_T{int(TT)}_delta0.2.npy'), allow_pickle=True).item()
        overlaps.append(np.abs(data['overlaps'][-1]))
    plt.figure()
    plt.plot(Tscan / P.DELTA, overlaps, '-o')
    plt.xlabel('T (in units 1/Δ)')
    plt.ylabel('final overlap |⟨ψ(0)|ψ(6T)⟩|')
    out2 = os.path.join(P.OUTDIR, 'reproduce_Fig3_overlap_vs_T.png')
    plt.savefig(out2)
    plt.close()
    print('Saved', out2)


def reproduce_fig4_and_5():
    print('Reproducing Fig.4 & Fig.5: effect of sinusoidal modulation (qualitative)')
    os.makedirs(P.OUTDIR, exist_ok=True)
    # We mimic Vx modulation by modulating the QD potential depth VD(t)
    def run_modulation(Vx0, Vx1, T=200.0):
        tlist, step_idx, slist, dt = tetron.make_time_grid(T_step=T, n_per_step=300)
        Emin = np.zeros_like(tlist)
        for i in range(len(tlist)):
            t = tlist[i]
            # sinusoidal modulation factor (mimic Zeeman effect on ABS energy)
            mod = Vx0 + Vx1 * np.cos(np.pi * t / T)
            # vary VD accordingly (rescale)
            VD_t = P.VD * (1.0 + mod)
            step = int(step_idx[i])
            s = float(slist[i])
            g1, g2, g3, g4 = tetron.gates_at(step, s)
            mu, t_links_mod, Delta_mod = map_gates_to_links(g1, g2, g3, g4, P.t0, P.DELTA, P.L, P.mu0, VD_t, P.QD_WIDTH)
            H = embed_kitaev.build_bdg(mu, t_links_mod, Delta_mod)
            E = np.linalg.eigvalsh(H)
            Emin[i] = np.min(np.abs(E))
        return tlist, Emin

    # Fig.4 case
    tlist, Emin = run_modulation(P.FIG4_VX0, P.FIG4_VX1, T=200.0)
    plt.figure()
    plt.plot(tlist / 200.0, Emin)
    plt.xlabel('normalized time (t/T)')
    plt.ylabel('min |E| (BdG)')
    out = os.path.join(P.OUTDIR, 'reproduce_Fig4_modulated_eigs.png')
    plt.savefig(out)
    plt.close()
    print('Saved', out)

    # Fig.5 amplitude comparison
    plt.figure()
    for amp in P.FIG5_VX1_OPTIONS:
        tlist, Emin = run_modulation(P.FIG4_VX0, amp, T=200.0)
        plt.plot(tlist / 200.0, Emin, label=f'Vx1={amp:.3f}')
    plt.xlabel('normalized time (t/T)')
    plt.ylabel('min |E| (BdG)')
    plt.legend()
    out2 = os.path.join(P.OUTDIR, 'reproduce_Fig5_modulation_amplitude.png')
    plt.savefig(out2)
    plt.close()
    print('Saved', out2)


def main():
    reproduce_fig1()
    reproduce_fig2()
    reproduce_fig3()
    reproduce_fig4_and_5()


if __name__ == '__main__':
    main()
