#!/usr/bin/env python3
"""Run ABS verification checks: PH symmetry, PH on eigenvectors, localization fit, Josephson E(phi).
Prints results and saves a small report file.
"""
import os
import argparse
import numpy as np
from numpy.linalg import eigh
# optional import for using pulse-generated chain
try:
    from tools.pulse_abs import h_from_xxz, pauli_map_to_params, build_chain_from_params
except ModuleNotFoundError:
    # when running script as `python3 tools/abs_checks.py`, package name may not be resolvable
    from pulse_abs import h_from_xxz, pauli_map_to_params, build_chain_from_params

# local builders (copied minimal from compute_abs)

def build_bdg_spatial(N, t_hop, Delta_bond, mu_site):
    A = np.zeros((N, N), dtype=complex)
    B = np.zeros((N, N), dtype=complex)
    for i in range(N):
        A[i, i] = -mu_site[i]
        if i < N - 1:
            t = t_hop[i]
            A[i, i + 1] = -t
            A[i + 1, i] = -t
            Delta = Delta_bond[i]
            B[i, i + 1] = Delta
            B[i + 1, i] = -Delta
    top = np.concatenate((A, B), axis=1)
    bottom = np.concatenate((-B.conj(), -A.T), axis=1)
    Hbdg = np.concatenate((top, bottom), axis=0)
    return Hbdg


def example_sn(HN=200, t0=1.0, Delta_L=0.3, Delta_R=0.0, mu=0.0):
    N = HN
    mid = N // 2
    t_hop = np.full(N - 1, t0, dtype=complex)
    Delta_bond = np.zeros(N - 1, dtype=complex)
    Delta_bond[:mid - 1] = Delta_L
    Delta_bond[mid - 1:] = Delta_R
    mu_site = np.full(N, mu, dtype=complex)
    H = build_bdg_spatial(N, t_hop, Delta_bond, mu_site)
    vals, vecs = eigh(H)
    return N, H, vals, vecs


def example_josephson(HN=200, t0=1.0, Delta=0.3, w=5, mu=0.0, phi=0.5):
    N = HN
    left_end = (N - w) // 2
    right_start = left_end + w
    t_hop = np.full(N - 1, t0, dtype=complex)
    Delta_bond = np.zeros(N - 1, dtype=complex)
    for i in range(0, left_end):
        Delta_bond[i] = Delta * np.exp(-1j * phi / 2)
    for i in range(right_start, N - 1):
        Delta_bond[i] = Delta * np.exp(1j * phi / 2)
    mu_site = np.full(N, mu, dtype=complex)
    H = build_bdg_spatial(N, t_hop, Delta_bond, mu_site)
    vals, vecs = eigh(H)
    return N, H, vals, vecs


def check_ph_pairs(vals):
    vals_r = np.sort(np.real_if_close(vals))
    paired = np.max(np.abs(vals_r + vals_r[::-1]))
    return paired


def check_ph_on_vecs(H, vecs, vals, ncheck=6):
    N = H.shape[0] // 2
    PH = np.block([[np.zeros((N,N)), np.eye(N)],[np.eye(N), np.zeros((N,N))]])
    order = np.argsort(np.abs(vals))
    residuals = []
    for k in range(min(ncheck, len(order))):
        i = order[k]
        E = vals[i]
        v = vecs[:, i]
        w = PH @ v.conj()
        res = np.linalg.norm(H @ w + E * w)
        residuals.append(float(res))
    return residuals


def fit_exponential_decay(rho, peak_idx, max_len=50):
    # Try one-sided fits to the right and to the left from the peak.
    N = len(rho)
    results = []
    eps = 1e-12

    def fit_side(yvals):
        mask = yvals > eps
        if mask.sum() < 3:
            return None
        x = np.arange(len(yvals))[mask]
        y = yvals[mask]
        logy = np.log(y)
        A = np.vstack([np.ones_like(x), -2 * x]).T
        sol, *_ = np.linalg.lstsq(A, logy, rcond=None)
        c, k = sol
        A0 = np.exp(c)
        xi = 1.0 / k if abs(k) > 0 else np.inf
        # compute residual norm
        pred = A @ sol
        resid = np.linalg.norm(pred - logy)
        return {'A': float(A0), 'kappa': float(k), 'xi': float(xi), 'resid': float(resid), 'n': int(mask.sum())}

    # right side (including peak)
    right_len = min(max_len, N - peak_idx)
    right_vals = rho[peak_idx: peak_idx + right_len]
    res_r = fit_side(right_vals)
    if res_r is not None:
        res_r['side'] = 'right'
        results.append(res_r)

    # left side (including peak), reverse so distance increases
    left_len = min(max_len, peak_idx + 1)
    left_vals = rho[peak_idx - left_len + 1: peak_idx + 1][::-1]
    res_l = fit_side(left_vals)
    if res_l is not None:
        res_l['side'] = 'left'
        results.append(res_l)

    if not results:
        return None
    # choose best by smallest residual, then by more points
    results.sort(key=lambda r: (r['resid'], -r['n']))
    return results[0]


def run_checks_on_H(H, N, tag='pulse'):
    os.makedirs('results', exist_ok=True)
    out_lines = []
    vals, vecs = eigh(H)
    paired_err = check_ph_pairs(vals)
    out_lines.append(f'PH pair max |E+E_rev| = {paired_err:.3e}')
    ph_res = check_ph_on_vecs(H, vecs, vals, ncheck=8)
    out_lines.append('PH on vectors residuals: ' + ', '.join([f'{r:.3e}' for r in ph_res]))
    mags = np.abs(vals)
    order = np.argsort(mags)
    modes_info = []
    for k in range(4):
        i = order[k]
        E = vals[i]
        modes_info.append((i, float(E)))
    out_lines.append('Lowest 4 modes (index, E): ' + ', '.join([f'{a}:{b:.3e}' for a,b in modes_info]))
    for kk in range(1,3):
        i = order[kk]
        E = vals[i]
        psi = vecs[:, i]
        rho = np.abs(psi[:N])**2 + np.abs(psi[N:])**2
        peak = np.argmax(rho)
        fit = fit_exponential_decay(rho, peak, max_len=80)
        if fit:
            out_lines.append(f'mode{kk} E={E:.3e} peak={peak} xi={fit["xi"]:.3e} kappa={fit["kappa"]:.3e} side={fit.get("side")} resid={fit.get("resid"):.3e}')
        else:
            out_lines.append(f'mode{kk} E={E:.3e} peak={peak} fit_failed')
    # Josephson scan using example_josephson for comparison
    phi_vals = np.linspace(0, np.pi, 13)
    lowest = []
    for phi in phi_vals:
        N2, H2, vals2, vecs2 = example_josephson(phi=phi)
        idx = np.argmin(np.abs(vals2))
        lowest.append(float(np.abs(vals2[idx])))
    out_lines.append('Josephson lowest |E| per phi: ' + ', '.join([f'{v:.3e}' for v in lowest]))
    report = '\n'.join(out_lines)
    print(report)
    report_path = f'results/abs_checks_{tag}_report.txt'
    with open(report_path,'w') as f:
        f.write(report + '\n')
    return report


def run_checks():
    # default behavior: original example-based checks
    os.makedirs('results', exist_ok=True)
    out_lines = []
    # S-N
    N, Hsn, vals_sn, vecs_sn = example_sn()
    paired_err = check_ph_pairs(vals_sn)
    out_lines.append(f'PH pair max |E+E_rev| = {paired_err:.3e}')
    ph_res = check_ph_on_vecs(Hsn, vecs_sn, vals_sn, ncheck=8)
    out_lines.append('PH on vectors residuals: ' + ', '.join([f'{r:.3e}' for r in ph_res]))
    # pick lowest nonzero mode
    mags = np.abs(vals_sn)
    order = np.argsort(mags)
    modes_info = []
    for k in range(4):
        i = order[k]
        E = vals_sn[i]
        modes_info.append((i, float(E)))
    out_lines.append('Lowest 4 modes (index, E): ' + ', '.join([f'{a}:{b:.3e}' for a,b in modes_info]))
    # localization fit for mode 1 and 2
    for kk in range(1,3):
        i = order[kk]
        E = vals_sn[i]
        psi = vecs_sn[:, i]
        rho = np.abs(psi[:N])**2 + np.abs(psi[N:])**2
        peak = np.argmax(rho)
        fit = fit_exponential_decay(rho, peak, max_len=80)
        if fit:
            out_lines.append(f'mode{kk} E={E:.3e} peak={peak} xi={fit["xi"]:.3e} kappa={fit["kappa"]:.3e}')
        else:
            out_lines.append(f'mode{kk} E={E:.3e} peak={peak} fit_failed')

    # Josephson scan
    phi_vals = np.linspace(0, np.pi, 13)
    lowest = []
    for phi in phi_vals:
        N2, H2, vals2, vecs2 = example_josephson(phi=phi)
        idx = np.argmin(np.abs(vals2))
        lowest.append(float(np.abs(vals2[idx])))
    out_lines.append('Josephson lowest |E| per phi: ' + ', '.join([f'{v:.3e}' for v in lowest]))

    report = '\n'.join(out_lines)
    print(report)
    with open('results/abs_checks_report.txt','w') as f:
        f.write(report + '\n')

    return report

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--use-pulse', action='store_true', help='Run checks on the pulse_abs constructed chain (eta=0.6, Delta_left=0.3)')
    args = parser.parse_args()
    if args.use_pulse:
        # reproduce pulse_abs parameters and build chain
        eta = 0.6
        mapping = h_from_xxz(eta)
        t, Delta, mu, U = pauli_map_to_params(mapping)
        # force left pairing as in pulse_abs demo
        Delta_left = 0.3
        N = 200
        H, t_hop, Delta_bond, mu_site = build_chain_from_params(N, t, Delta_left, mu, t, 0.0, mu, interface_width=4)
        run_checks_on_H(H, N, tag='pulse')
    else:
        run_checks()
