#!/usr/bin/env python3
"""Analyze tetron path simulation results and produce diagnostic plots.
"""
import glob
import os
import re
import numpy as np
import matplotlib.pyplot as plt


def load_result(path):
    data = np.load(path, allow_pickle=True)
    if isinstance(data, np.ndarray) and data.shape == ():
        data = data.item()
    return data


def collect_metrics(pattern):
    files = sorted(glob.glob(pattern))
    Ts = []
    final_overlaps = []
    final_p0 = []
    final_p1 = []
    for f in files:
        m = re.search(r"_T(\d+)_", f)
        if not m:
            continue
        T = int(m.group(1))
        d = load_result(f)
        overlaps = d['overlaps']
        probs = d['probs']
        Ts.append(T)
        final_overlaps.append(abs(overlaps[-1]))
        final_p0.append(probs[-1, 0])
        final_p1.append(probs[-1, 1])
    return np.array(Ts), np.array(final_overlaps), np.array(final_p0), np.array(final_p1)


def plot_overlap_vs_T():
    os.makedirs('results', exist_ok=True)
    T_m, ov_m, p0_m, p1_m = collect_metrics('results/tetron_MZM_T*_delta0.0.npy')
    T_a, ov_a, p0_a, p1_a = collect_metrics('results/tetron_ABS_T*_delta0.2.npy')

    plt.figure()
    if T_m.size:
        idx = np.argsort(T_m)
        plt.plot(T_m[idx], ov_m[idx], 'o-', label='MZM final overlap')
    if T_a.size:
        idx = np.argsort(T_a)
        plt.plot(T_a[idx], ov_a[idx], 's--', label='ABS final overlap')
    plt.xlabel('T_step')
    plt.ylabel('final overlap |<psi0|psi(T)>|')
    plt.legend()
    plt.grid(True)
    out = 'results/analysis_final_overlap_vs_T.png'
    plt.savefig(out)
    print('Saved', out)


def plot_timeseries_T(T=400):
    os.makedirs('results', exist_ok=True)
    f_m = f'results/tetron_MZM_T{T}_delta0.0.npy'
    f_a = f'results/tetron_ABS_T{T}_delta0.2.npy'
    if os.path.exists(f_m):
        d = load_result(f_m)
        t = d['t']
        probs = d['probs']
        plt.figure(figsize=(8,4))
        plt.plot(t, probs[:,0], label='MZM p(site1)')
        plt.plot(t, probs[:,1], label='MZM p(site2)')
        plt.xlabel('time')
        plt.ylabel('probability')
        plt.title(f'MZM timeseries T={T}')
        plt.legend()
        plt.tight_layout()
        out = f'results/analysis_timeseries_MZM_T{T}.png'
        plt.savefig(out)
        print('Saved', out)
    if os.path.exists(f_a):
        d = load_result(f_a)
        t = d['t']
        probs = d['probs']
        plt.figure(figsize=(8,4))
        plt.plot(t, probs[:,0], label='ABS p(site1)')
        plt.plot(t, probs[:,1], label='ABS p(site2)')
        plt.xlabel('time')
        plt.ylabel('probability')
        plt.title(f'ABS timeseries T={T}')
        plt.legend()
        plt.tight_layout()
        out = f'results/analysis_timeseries_ABS_T{T}.png'
        plt.savefig(out)
        print('Saved', out)


def main():
    plot_overlap_vs_T()
    plot_timeseries_T(T=400)


if __name__ == '__main__':
    main()
