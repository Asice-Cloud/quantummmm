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
from tools import lindblad_twolevel as lindblad


# default LDOS etas and Lindblad cases to include in comparison overlays
DEFAULT_ETAS = [1e-2, 1e-3, 1e-4]
LINDblad_CASES = {
    'coherent': (0.0, 0.0),
    'dephasing': (0.5, 0.0),
    'relax': (0.0, 0.1),
    'both': (0.3, 0.05),
}


def compute_bdg_minE(T_step=100.0, n_per_step=300, mod_fn=None, mode='minE', eta=None, site_index=0):
    """Compute either min|E| (mode='minE') or site LDOS at E=0 (mode='ldos').

    Returns tlist, values, S, M, Theta, step_idx, slist, dt
    """
    tlist, step_idx, slist, dt = tetron.make_time_grid(T_step=T_step, n_per_step=n_per_step)
    N = len(tlist)
    vals = np.zeros(N)
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
        if mode == 'minE':
            E = np.linalg.eigvalsh(H)
            vals[i] = np.min(np.abs(E))
        elif mode == 'ldos':
            ldos, Egrid = embed_kitaev.compute_zero_ldos(H, eta=eta)
            # choose site index (default 0)
            vals[i] = ldos[site_index] if site_index < len(ldos) else np.max(ldos)
        else:
            raise ValueError('unknown mode')

        t_left = t_links_mod[0] if len(t_links_mod) > 0 else 0.0
        t_right = t_links_mod[1] if len(t_links_mod) > 1 else 0.0
        S[i] = t_left + t_right
        M[i] = mu[0] - mu[1] if len(mu) > 1 else 0.0

    return tlist, vals, S, M, Theta, step_idx, slist, dt


def E_pred_from_params(A0, B0, C0, ts, S, M, Theta):
    return np.sqrt((A0 * S) ** 2 + (B0 * np.sin(Theta * ts)) ** 2 + (C0 * M) ** 2)


def plot_compare(tlist, E_bdg, E_pred, out, title, params, rmse):
    fig, ax = plt.subplots(figsize=(6.5, 3.5), dpi=200)
    T = max(tlist) - min(tlist)
    ax.plot(tlist / T, E_bdg, label='BdG min |E|')
    ax.plot(tlist / T, E_pred, '--', label='mapped Pauli E_pred')
    # paper-style x ticks: map normalized [0,1] -> step labels 0..3
    ticks = np.linspace(0.0, 1.0, 4)
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(i) for i in range(4)])
    ax.set_xlabel(r'step $t/T$')
    ax.set_ylabel(r'Energy ($\Delta$)')
    ax.set_title(title)
    ax.legend()
    ann = f"A0={params[0]:.6g}\nB0={params[1]:.6g}\nC0={params[2]:.6g}\nts={params[3]:.6g}\nRMSE={rmse:.3e}"
    ax.text(0.01, 0.98, ann, transform=ax.transAxes, fontsize=8, verticalalignment='top',
            bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
    plt.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    print('Saved', out)


def lorentzian_ldos_from_E(E, eta):
    # unnormalized Lorentzian-like proxy for LDOS
    return eta / (E ** 2 + eta ** 2)


def build_two_level_H_list(A0, B0, C0, ts, step_idx, slist):
    H_list = []
    for step, s in zip(step_idx, slist):
        g1, g2, g3, g4 = tetron.gates_at(int(step), float(s))
        S = float(g1 + g3)
        M = float(g2 - g4)
        theta = tetron.theta_from_time(int(step), float(s))
        dx = A0 * S
        dy = B0 * np.sin(theta * ts)
        dz = C0 * M
        H = dx * lindblad.sigma_x + dy * lindblad.sigma_y + dz * lindblad.sigma_z
        H_list.append(H)
    return H_list


def main():
    os.makedirs('results', exist_ok=True)
    metrics = {}

    # Fig.2: multiple T values — include LDOS etas and Lindblad overlays
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
        tlist, E_bdg, S, M, Theta, step_idx, slist, dt = compute_bdg_minE(T_step=T, n_per_step=300, mode='minE')
        E_pred = E_pred_from_params(A0, B0, C0, ts, S, M, Theta)

        # run Lindblad on the mapped two-level H(t)
        H_list = build_two_level_H_list(A0, B0, C0, ts, step_idx, slist)
        lindblad_rnorms = {}
        for name, (g_deph, g_rel) in LINDblad_CASES.items():
            bloch = lindblad.run_lindblad_time_series(H_list, dt=dt, gamma_deph=g_deph, gamma_relax=g_rel)
            rnorm = np.linalg.norm(bloch, axis=1)
            lindblad_rnorms[name] = rnorm

        # Energy comparison with Lindblad-damped mapped predictions
        rmse = np.sqrt(np.mean((E_bdg - E_pred) ** 2))
        fig, ax = plt.subplots(figsize=(6.5, 3.5), dpi=200)
        Tnorm = max(tlist) - min(tlist)
        ax.plot(tlist / Tnorm, E_bdg, label='BdG min |E|')
        ax.plot(tlist / Tnorm, E_pred, '--', label='mapped Pauli E_pred')
        for name, rnorm in lindblad_rnorms.items():
            ax.plot(tlist / Tnorm, E_pred * rnorm, ':', label=f'Pauli damped ({name})')
        # show x-axis in paper 'step t/T' units (0..6)
        ticks = np.linspace(0.0, 1.0, 7)
        ax.set_xticks(ticks)
        ax.set_xticklabels([str(i) for i in range(7)])
        ax.set_xlabel('step t/T')
        ax.set_ylabel('energy (Δ units)')
        # energy axis in Fig.2 uses full gap units
        ax.set_ylim(-1.0, 1.0)
        ax.set_title(f'Fig.2 comparison T={int(T)} (with Lindblad overlays)')
        ax.legend(fontsize=8)
        ann = f"A0={A0:.6g}\nB0={B0:.6g}\nC0={C0:.6g}\nts={ts:.6g}\nRMSE={rmse:.3e}"
        ax.text(0.01, 0.98, ann, transform=ax.transAxes, fontsize=8, verticalalignment='top',
                bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
        plt.tight_layout()
        out = f'results/compare_Fig2_T{int(T)}.png'
        fig.savefig(out)
        plt.close(fig)
        print('Saved', out)
        metrics[f'Fig2_T{int(T)}_rmse'] = rmse
        fig2_outs.append(out)

        # LDOS comparison for a few eta values: overlay predicted Lorentzian and damped versions
        for eta in DEFAULT_ETAS:
            t2, ldos_bdg, S2, M2, Theta2, step_idx2, slist2, dt2 = compute_bdg_minE(T_step=T, n_per_step=300, mode='ldos', eta=eta, site_index=0)
            # map E_pred -> Lorentzian LDOS proxy
            ldos_pred = lorentzian_ldos_from_E(E_pred, eta)
            # scale predicted to BdG amplitude for visual comparison
            scale = (np.max(ldos_bdg) / np.max(ldos_pred)) if np.max(ldos_pred) > 0 else 1.0
            ldos_pred_scaled = ldos_pred * scale
            # damped predictions using Lindblad rnorms
            fig, ax = plt.subplots(figsize=(6.5, 3.5), dpi=200)
            ax.plot(t2 / Tnorm, ldos_bdg, label=f'BdG LDOS eta={eta:.0e}')
            ax.plot(t2 / Tnorm, ldos_pred_scaled, '--', label='Pred Lorentzian (coherent)')
            for name, rnorm in lindblad_rnorms.items():
                ax.plot(t2 / Tnorm, ldos_pred_scaled * rnorm, ':', label=f'Pred damped ({name})')
            # paper-style x ticks and LDOS label
            ticks = np.linspace(0.0, 1.0, 4)
            ax.set_xticks(ticks)
            ax.set_xticklabels([str(i) for i in range(4)])
            ax.set_xlabel(r'step $t/T$')
            ax.set_ylabel('LDOS (arb. units)')
            # set y-range roughly to match paper heatmaps when available
            ax.set_ylim(0, 405)
            ax.set_title(f'Fig.2 LDOS compare T={int(T)} eta={eta:.0e}')
            ax.legend(fontsize=8)
            plt.tight_layout()
            out_ldos = f'results/compare_Fig2_T{int(T)}_ldos_eta{int(-np.log10(eta))}.png'
            fig.savefig(out_ldos)
            plt.close(fig)
            print('Saved', out_ldos)

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
        tlist, E_bdg, S, M, Theta, step_idx, slist, dt = compute_bdg_minE(T_step=T, n_per_step=P.FIG3_N_PER_STEP, mode='minE')
        E_pred = E_pred_from_params(A0, B0, C0, ts, S, M, Theta)

        # Lindblad overlays
        H_list = build_two_level_H_list(A0, B0, C0, ts, step_idx, slist)
        lindblad_rnorms = {}
        for name, (g_deph, g_rel) in LINDblad_CASES.items():
            bloch = lindblad.run_lindblad_time_series(H_list, dt=dt, gamma_deph=g_deph, gamma_relax=g_rel)
            lindblad_rnorms[name] = np.linalg.norm(bloch, axis=1)

        rmse = np.sqrt(np.mean((E_bdg - E_pred) ** 2))
        fig, ax = plt.subplots(figsize=(6.5, 3.5), dpi=200)
        Tnorm = max(tlist) - min(tlist)
        ax.plot(tlist / Tnorm, E_bdg, label='BdG min |E|')
        ax.plot(tlist / Tnorm, E_pred, '--', label='mapped Pauli E_pred')
        for name, rnorm in lindblad_rnorms.items():
            ax.plot(tlist / Tnorm, E_pred * rnorm, ':', label=f'Pauli damped ({name})')
        # paper-style x ticks: 0..3 for the 3-step path
        ticks = np.linspace(0.0, 1.0, 4)
        ax.set_xticks(ticks)
        ax.set_xticklabels([str(i) for i in range(4)])
        ax.set_xlabel(r'step $t/T$')
        ax.set_ylabel(r'Energy ($\Delta$)')
        # Fig.3 zoomed energy panel: small energy range
        ax.set_ylim(-0.02, 0.02)
        ax.set_title('Fig.3 comparison (with Lindblad overlays)')
        ax.legend(fontsize=8)
        ann = f"A0={A0:.6g}\nB0={B0:.6g}\nC0={C0:.6g}\nts={ts:.6g}\nRMSE={rmse:.3e}"
        ax.text(0.01, 0.98, ann, transform=ax.transAxes, fontsize=8, verticalalignment='top',
                bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
        plt.tight_layout()
        out = 'results/compare_Fig3.png'
        fig.savefig(out)
        plt.close(fig)
        print('Saved', out)
        metrics['Fig3_rmse'] = rmse

        # LDOS comparisons
        for eta in DEFAULT_ETAS:
            t2, ldos_bdg, S2, M2, Theta2, step_idx2, slist2, dt2 = compute_bdg_minE(T_step=T, n_per_step=P.FIG3_N_PER_STEP, mode='ldos', eta=eta)
            ldos_pred = lorentzian_ldos_from_E(E_pred, eta)
            scale = (np.max(ldos_bdg) / np.max(ldos_pred)) if np.max(ldos_pred) > 0 else 1.0
            ldos_pred_scaled = ldos_pred * scale
            fig, ax = plt.subplots(figsize=(6.5, 3.5), dpi=200)
            ax.plot(t2 / (max(t2) - min(t2)), ldos_bdg, label=f'BdG LDOS eta={eta:.0e}')
            ax.plot(t2 / (max(t2) - min(t2)), ldos_pred_scaled, '--', label='Pred Lorentzian (coherent)')
            for name, rnorm in lindblad_rnorms.items():
                ax.plot(t2 / (max(t2) - min(t2)), ldos_pred_scaled * rnorm, ':', label=f'Pred damped ({name})')
            ticks = np.linspace(0.0, 1.0, 4)
            ax.set_xticks(ticks)
            ax.set_xticklabels([str(i) for i in range(4)])
            ax.set_xlabel(r'step $t/T$')
            ax.set_ylabel('LDOS (arb. units)')
            ax.set_ylim(0, 405)
            ax.set_title(f'Fig.3 LDOS compare eta={eta:.0e}')
            ax.legend(fontsize=8)
            plt.tight_layout()
            out_ldos = f'results/compare_Fig3_ldos_eta{int(-np.log10(eta))}.png'
            fig.savefig(out_ldos)
            plt.close(fig)
            print('Saved', out_ldos)

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
        tlist, E_bdg, S, M, Theta, step_idx, slist, dt = compute_bdg_minE(T_step=T, n_per_step=300, mod_fn=mod_fn, mode='minE')
        E_pred = E_pred_from_params(A0, B0, C0, ts, S, M, Theta)

        H_list = build_two_level_H_list(A0, B0, C0, ts, step_idx, slist)
        lindblad_rnorms = {}
        for name, (g_deph, g_rel) in LINDblad_CASES.items():
            bloch = lindblad.run_lindblad_time_series(H_list, dt=dt, gamma_deph=g_deph, gamma_relax=g_rel)
            lindblad_rnorms[name] = np.linalg.norm(bloch, axis=1)

        rmse = np.sqrt(np.mean((E_bdg - E_pred) ** 2))
        fig, ax = plt.subplots(figsize=(6.5, 3.5), dpi=200)
        Tnorm = max(tlist) - min(tlist)
        ax.plot(tlist / Tnorm, E_bdg, label='BdG min |E|')
        ax.plot(tlist / Tnorm, E_pred, '--', label='mapped Pauli E_pred')
        for name, rnorm in lindblad_rnorms.items():
            ax.plot(tlist / Tnorm, E_pred * rnorm, ':', label=f'Pauli damped ({name})')
        # paper-style x ticks: 0..3 for step t/T
        ticks = np.linspace(0.0, 1.0, 4)
        ax.set_xticks(ticks)
        ax.set_xticklabels([str(i) for i in range(4)])
        ax.set_xlabel(r'step $t/T$')
        ax.set_ylabel(r'Energy ($\Delta$)')
        ax.set_ylim(-0.02, 0.02)
        ax.set_title('Fig.4 comparison (modulated, with Lindblad overlays)')
        ax.legend(fontsize=8)
        ann = f"A0={A0:.6g}\nB0={B0:.6g}\nC0={C0:.6g}\nts={ts:.6g}\nRMSE={rmse:.3e}"
        ax.text(0.01, 0.98, ann, transform=ax.transAxes, fontsize=8, verticalalignment='top',
                bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
        plt.tight_layout()
        out = 'results/compare_Fig4.png'
        fig.savefig(out)
        plt.close(fig)
        print('Saved', out)
        metrics['Fig4_rmse'] = rmse

        for eta in DEFAULT_ETAS:
            t2, ldos_bdg, S2, M2, Theta2, step_idx2, slist2, dt2 = compute_bdg_minE(T_step=T, n_per_step=300, mod_fn=mod_fn, mode='ldos', eta=eta)
            ldos_pred = lorentzian_ldos_from_E(E_pred, eta)
            scale = (np.max(ldos_bdg) / np.max(ldos_pred)) if np.max(ldos_pred) > 0 else 1.0
            ldos_pred_scaled = ldos_pred * scale
            fig, ax = plt.subplots(figsize=(6.5, 3.5), dpi=200)
            ax.plot(t2 / (max(t2) - min(t2)), ldos_bdg, label=f'BdG LDOS eta={eta:.0e}')
            ax.plot(t2 / (max(t2) - min(t2)), ldos_pred_scaled, '--', label='Pred Lorentzian (coherent)')
            for name, rnorm in lindblad_rnorms.items():
                ax.plot(t2 / (max(t2) - min(t2)), ldos_pred_scaled * rnorm, ':', label=f'Pred damped ({name})')
            ticks = np.linspace(0.0, 1.0, 4)
            ax.set_xticks(ticks)
            ax.set_xticklabels([str(i) for i in range(4)])
            ax.set_xlabel(r'step $t/T$')
            ax.set_ylabel('LDOS (arb. units)')
            ax.set_ylim(0, 405)
            ax.set_title(f'Fig.4 LDOS compare eta={eta:.0e}')
            ax.legend(fontsize=8)
            plt.tight_layout()
            out_ldos = f'results/compare_Fig4_ldos_eta{int(-np.log10(eta))}.png'
            fig.savefig(out_ldos)
            plt.close(fig)
            print('Saved', out_ldos)

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
        tlist, E_bdg, S, M, Theta, step_idx, slist, dt = compute_bdg_minE(T_step=T, n_per_step=300, mod_fn=lambda t,amp=amp: mod_fn(t, Vx1=amp), mode='minE')
        E_pred = E_pred_from_params(A0, B0, C0, ts, S, M, Theta)

        H_list = build_two_level_H_list(A0, B0, C0, ts, step_idx, slist)
        lindblad_rnorms = {}
        for name, (g_deph, g_rel) in LINDblad_CASES.items():
            bloch = lindblad.run_lindblad_time_series(H_list, dt=dt, gamma_deph=g_deph, gamma_relax=g_rel)
            lindblad_rnorms[name] = np.linalg.norm(bloch, axis=1)

        rmse = np.sqrt(np.mean((E_bdg - E_pred) ** 2))
        out = f'results/compare_Fig5_amp{amp:.6g}.png'
        fig, ax = plt.subplots(figsize=(6.5, 3.5), dpi=200)
        Tnorm = max(tlist) - min(tlist)
        ax.plot(tlist / Tnorm, E_bdg, label='BdG min |E|')
        ax.plot(tlist / Tnorm, E_pred, '--', label='mapped Pauli E_pred')
        for name, rnorm in lindblad_rnorms.items():
            ax.plot(tlist / Tnorm, E_pred * rnorm, ':', label=f'Pauli damped ({name})')
        # paper-style x ticks: 0..3
        ticks = np.linspace(0.0, 1.0, 4)
        ax.set_xticks(ticks)
        ax.set_xticklabels([str(i) for i in range(4)])
        ax.set_xlabel(r'step $t/T$')
        ax.set_ylabel(r'Energy ($\Delta$)')
        ax.set_ylim(-0.02, 0.02)
        ax.set_title(f'Fig.5 comparison amp={amp} (with Lindblad overlays)')
        ax.legend(fontsize=8)
        ann = f"A0={A0:.6g}\nB0={B0:.6g}\nC0={C0:.6g}\nts={ts:.6g}\nRMSE={rmse:.3e}"
        ax.text(0.01, 0.98, ann, transform=ax.transAxes, fontsize=8, verticalalignment='top',
                bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
        plt.tight_layout()
        fig.savefig(out)
        plt.close(fig)
        print('Saved', out)
        metrics[f'Fig5_amp{amp:.6g}_rmse'] = rmse

        for eta in DEFAULT_ETAS:
            t2, ldos_bdg, S2, M2, Theta2, step_idx2, slist2, dt2 = compute_bdg_minE(T_step=T, n_per_step=300, mod_fn=lambda t,amp=amp: mod_fn(t, Vx1=amp), mode='ldos', eta=eta)
            ldos_pred = lorentzian_ldos_from_E(E_pred, eta)
            scale = (np.max(ldos_bdg) / np.max(ldos_pred)) if np.max(ldos_pred) > 0 else 1.0
            ldos_pred_scaled = ldos_pred * scale
            fig, ax = plt.subplots(figsize=(6.5, 3.5), dpi=200)
            ax.plot(t2 / (max(t2) - min(t2)), ldos_bdg, label=f'BdG LDOS eta={eta:.0e}')
            ax.plot(t2 / (max(t2) - min(t2)), ldos_pred_scaled, '--', label='Pred Lorentzian (coherent)')
            for name, rnorm in lindblad_rnorms.items():
                ax.plot(t2 / (max(t2) - min(t2)), ldos_pred_scaled * rnorm, ':', label=f'Pred damped ({name})')
            ticks = np.linspace(0.0, 1.0, 4)
            ax.set_xticks(ticks)
            ax.set_xticklabels([str(i) for i in range(4)])
            ax.set_xlabel(r'step $t/T$')
            ax.set_ylabel('LDOS (arb. units)')
            ax.set_ylim(0, 405)
            ax.set_title(f'Fig.5 LDOS compare amp={amp} eta={eta:.0e}')
            ax.legend(fontsize=8)
            plt.tight_layout()
            out_ldos = f'results/compare_Fig5_amp{amp:.6g}_ldos_eta{int(-np.log10(eta))}.png'
            fig.savefig(out_ldos)
            plt.close(fig)
            print('Saved', out_ldos)

    np.savez('results/compare_metrics.npz', **metrics)
    print('Saved comparison metrics to results/compare_metrics.npz')


if __name__ == '__main__':
    main()
