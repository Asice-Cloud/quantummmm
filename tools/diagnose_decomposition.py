#!/usr/bin/env python3
"""Compute and save U, U_geom, U_dyn for a given delta.
"""
import os
import math
import argparse
import numpy as np

# import local helpers from the same directory
import verify_from_R as vr


def save_decomposition(delta=0.015, N=600, alphas=None):
    os.makedirs('results', exist_ok=True)
    us, du, H4s = vr.compute_Hs_from_R(delta, N=N)
    H2_list, ds = vr.build_H2_list_from_H4s(H4s)
    dt = du
    Evals = np.array([np.linalg.norm(d) for d in ds])
    I0 = float(np.sum(Evals) * dt)
    alpha_target = math.pi / I0 if I0 > 0 else float('nan')

    if alphas is None:
        alphas = [1.0]
    if not math.isnan(alpha_target) and alpha_target not in alphas:
        alphas.append(alpha_target)

    out_data = {'delta': delta, 'dt': dt, 'I0': I0, 'alpha_target': alpha_target, 'us': us}

    for alpha in alphas:
        Hs = [alpha * H for H in H2_list]
        U_final, Ulist = vr.compute_U_from_Hlist(Hs, dt)
        U_geom, U_dyn, theta, gamma = vr.adiabatic_decomposition(Hs, dt)
        prod = U_geom @ U_dyn
        res_norm = float(np.linalg.norm(U_final - prod))
        axis, angle = vr.rot_axis_angle_from_U(U_final)
        out_data[f'alpha_{alpha:.6g}'] = {
            'U_final': U_final,
            'U_geom': U_geom,
            'U_dyn': U_dyn,
            'res_norm': res_norm,
            'axis': axis,
            'angle': angle,
            'theta': theta,
            'gamma': gamma,
        }
        # save a comparison figure
        outfig = f'results/compare_decomp_delta{delta:.3g}_alpha{alpha:.6g}.png'
        try:
            vr.plot_U_comparison(U_final, U_geom, U_dyn, outfig)
        except Exception as e:
            print('Failed to plot comparison:', e)

    np.savez(f'results/decompose_delta{delta:.3g}.npz', **out_data)
    print('Saved decomposition results to', f'results/decompose_delta{delta:.3g}.npz')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--delta', type=float, default=0.015)
    p.add_argument('--N', type=int, default=600)
    p.add_argument('--alpha', type=float, nargs='*', default=None)
    args = p.parse_args()
    save_decomposition(delta=args.delta, N=args.N, alphas=args.alpha)
