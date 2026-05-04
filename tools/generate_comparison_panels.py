#!/usr/bin/env python3
"""Generate comparison panels: BdG vs mapped Pauli prediction for Fig.2-5.

Produces PNGs in `results/compare_Fig*.png` and saves RMSE metrics.
"""
import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt

# repo root
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from tools import tetron_path_sim as tetron
from tools.reproduce_figs import map_gates_to_links
from tools import embed_kitaev
import tools.paper_params as P


def compute_bdg_minE(T_step=100.0, n_per_step=300, mod_fn=None):
    tlist, step_idx, slist, dt = tetron.make_time_grid(T_step=T_step, n_per_step=n_per_step)
    N = len(tlist)
    E_bdg = np.zeros(N)
    S = np.zeros(N)
    M = np.zeros(N)
    Theta = np.zeros(N)

    for i in range(N):
        step = int(step_idx[i])
        s = float(slist[i])
        theta = tetron.theta_from_time(step, s)
        Theta[i] = theta
        g1, g2, g3, g4 = tetron.gates_at(step, s)
        if mod_fn is None:
            VD_t = P.VD
        else:
            t = tlist[i]
            VD_t = mod_fn(t)
        mu, t_links_mod, Delta_mod = map_gates_to_links(g1, g2, g3, g4, P.t0, P.DELTA, P.L, P.mu0, VD_t, P.QD_WIDTH)
        H = embed_kitaev.build_bdg(mu, t_links_mod, Delta_mod)
        E = np.linalg.eigvalsh(H)
        E_bdg[i] = np.min(np.abs(E))
        t_left = t_links_mod[0] if len(t_links_mod) > 0 else 0.0
        t_right = t_links_mod[1] if len(t_links_mod) > 1 else 0.0
        S[i] = t_left + t_right
        M[i] = mu[0] - mu[1] if len(mu) > 1 else 0.0

    return tlist, E_bdg, S, M, Theta


def E_pred_from_params(A0, B0, C0, ts, S, M, Theta):
    return np.sqrt((A0 * S) ** 2 + (B0 * np.sin(Theta * ts)) ** 2 + (C0 * M) ** 2)


def plot_compare(tlist, E_bdg, E_pred, out, title, params, rmse):
    fig, ax = plt.subplots(figsize=(6.5, 3.5), dpi=200)
    T = max(tlist) - min(tlist)
    ax.plot(tlist / T, E_bdg, label='BdG min |E|')
    ax.plot(tlist / T, E_pred, '--', label='mapped Pauli E_pred')
    ax.set_xlabel('normalized time (t/T)')
    ax.set_ylabel('energy (Δ units)')
    ax.set_title(title)
    ax.legend()
    ann = f"A0={params[0]:.6g}\nB0={params[1]:.6g}\nC0={params[2]:.6g}\nts={params[3]:.6g}\nRMSE={rmse:.3e}"
    ax.text(0.01, 0.98, ann, transform=ax.transAxes, fontsize=8, verticalalignment='top',
            bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
    plt.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    print('Saved', out)


def main():
    os.makedirs('results', exist_ok=True)
    metrics = {}

    # Fig.2: multiple T values
    d2 = np.load('results/mapping_fit_fig2.npz')
    fig2_outs = []
    for T in P.FIG2_TS:
        key = f'T{int(T)}'
        if key not in d2.files:
            print('missing fit for', key)
            continue
        params = d2[key]
        if params.size < 4:
            params = np.concatenate([params, [1.0]])
        A0, B0, C0, ts = params[:4]
        tlist, E_bdg, S, M, Theta = compute_bdg_minE(T_step=T, n_per_step=300)
        E_pred = E_pred_from_params(A0, B0, C0, ts, S, M, Theta)
        rmse = np.sqrt(np.mean((E_bdg - E_pred) ** 2))
        out = f'results/compare_Fig2_T{int(T)}.png'
        plot_compare(tlist, E_bdg, E_pred, out, f'Fig.2 comparison T={int(T)}', (A0, B0, C0, ts), rmse)
        metrics[f'Fig2_T{int(T)}_rmse'] = rmse
        fig2_outs.append(out)

    # combine Fig2 panels into one image
    if fig2_outs:
        imgs = [plt.imread(p) for p in fig2_outs]
        h = max(im.shape[0] for im in imgs)
        total_w = sum(im.shape[1] for im in imgs)
        import numpy as _np
        canvas = _np.ones((h, total_w, 3), dtype=_np.float32)
        x = 0
        for im in imgs:
            H, W = im.shape[0], im.shape[1]
            canvas[:H, x:x+W, :] = im[:, :, :3]
            x += W
        fig_all = plt.figure(figsize=(9, 3), dpi=200)
        plt.imshow(canvas)
        plt.axis('off')
        out_all = 'results/compare_Fig2_panel.png'
        plt.tight_layout(pad=0)
        fig_all.savefig(out_all)
        plt.close(fig_all)
        print('Saved', out_all)
        metrics['Fig2_panel'] = out_all

    # Fig.3
    if os.path.exists('results/mapping_fit_fig3.npz'):
        d3 = np.load('results/mapping_fit_fig3.npz')
        A0, B0, C0, ts = float(d3['A0_fit']), float(d3['B0_fit']), float(d3['C0_fit']), float(d3['ts_fit'])
        T = P.FIG3_T
        tlist, E_bdg, S, M, Theta = compute_bdg_minE(T_step=T, n_per_step=P.FIG3_N_PER_STEP)
        E_pred = E_pred_from_params(A0, B0, C0, ts, S, M, Theta)
        rmse = np.sqrt(np.mean((E_bdg - E_pred) ** 2))
        out = 'results/compare_Fig3.png'
        plot_compare(tlist, E_bdg, E_pred, out, 'Fig.3 comparison', (A0, B0, C0, ts), rmse)
        metrics['Fig3_rmse'] = rmse

    # Fig.4
    if os.path.exists('results/mapping_fit_fig4.npz'):
        d4 = np.load('results/mapping_fit_fig4.npz')
        res = d4['res']
        A0, B0, C0, ts = [float(x) for x in res]
        T = 200.0
        Vx0 = P.FIG4_VX0
        Vx1 = P.FIG4_VX1
        def mod_fn(t):
            mod = Vx0 + Vx1 * np.cos(np.pi * t / T)
            return P.VD * (1.0 + mod)
        tlist, E_bdg, S, M, Theta = compute_bdg_minE(T_step=T, n_per_step=300, mod_fn=mod_fn)
        E_pred = E_pred_from_params(A0, B0, C0, ts, S, M, Theta)
        rmse = np.sqrt(np.mean((E_bdg - E_pred) ** 2))
        out = 'results/compare_Fig4.png'
        plot_compare(tlist, E_bdg, E_pred, out, 'Fig.4 comparison (modulated)', (A0, B0, C0, ts), rmse)
        metrics['Fig4_rmse'] = rmse

    # Fig.5 (multiple amplitudes)
    for amp in P.FIG5_VX1_OPTIONS:
        pattern = f'results/mapping_fit_fig5_amp{amp:.6g}.npz'
        if not os.path.exists(pattern):
            # try glob
            matches = glob.glob('results/mapping_fit_fig5_amp*.npz')
            matches = [m for m in matches if f'{amp:.6g}' in m]
            if matches:
                pattern = matches[0]
            else:
                continue
        d5 = np.load(pattern)
        res = d5['res']
        A0, B0, C0, ts = [float(x) for x in res]
        T = 200.0
        def mod_fn(t, Vx1=amp):
            mod = P.FIG4_VX0 + Vx1 * np.cos(np.pi * t / T)
            return P.VD * (1.0 + mod)
        tlist, E_bdg, S, M, Theta = compute_bdg_minE(T_step=T, n_per_step=300, mod_fn=lambda t,amp=amp: mod_fn(t, Vx1=amp))
        E_pred = E_pred_from_params(A0, B0, C0, ts, S, M, Theta)
        rmse = np.sqrt(np.mean((E_bdg - E_pred) ** 2))
        out = f'results/compare_Fig5_amp{amp:.6g}.png'
        plot_compare(tlist, E_bdg, E_pred, out, f'Fig.5 comparison amp={amp}', (A0, B0, C0, ts), rmse)
        metrics[f'Fig5_amp{amp:.6g}_rmse'] = rmse

    np.savez('results/compare_metrics.npz', **metrics)
    print('Saved comparison metrics to results/compare_metrics.npz')


if __name__ == '__main__':
    main()
