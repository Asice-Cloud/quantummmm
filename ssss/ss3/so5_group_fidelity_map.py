"""SO(5) group-space fidelity map for the PRB111 Fig. 1(d) protocol.

The plotted quantity is F_group = F_surv * F_rot_cond, as defined in n2.md.
Units: x=tau/(100/meV), so the dimensionless integration time is q=t/(100/meV).
The 3-step protocol is repeated twice (total time 6 tau), matching the paper's
two successive swaps.
"""
from pathlib import Path
import os
import numpy as np
import matplotlib.pyplot as plt

TC = 0.3       # meV
E0 = 0.3       # meV
E1 = float(os.environ.get('SO5_E1', '0.01'))  # meV; text says 0.001, caption says 0.01
SO5_NORMALIZATION = 2.0  # dot(gamma)=2 h gamma; the factor 2 from model.md
N_X, N_Y = 121, 101
XMAX = 12.0
YMIN, YMAX = -1.0, 1.0
DQ = 0.01      # q units; midpoint/RK4 step (phase step <= 0.072 for tc)

B23 = np.array([[1., 0., 0.], [0., 0., -1.], [0., 1., 0.]])
B = B23 @ B23                 # two exchanges at t=6 tau

def env(u, rising):
    return (1.0 - np.cos(np.pi*u))/2 if rising else (1.0 + np.cos(np.pi*u))/2

def omega_batch(t2, t3, ed, t1):
    """Return batch of 5x5 antisymmetric generators in meV."""
    n = len(t1)
    o = np.zeros((n, 5, 5), dtype=float)
    s = SO5_NORMALIZATION
    o[:, 0, 1] =  s*E1; o[:, 1, 0] = -s*E1
    o[:, 3, 1] =  s*t2; o[:, 1, 3] = -s*t2
    o[:, 4, 0] = -s*t1; o[:, 0, 4] =  s*t1
    o[:, 3, 2] = -s*t3; o[:, 2, 3] =  s*t3
    o[:, 3, 4] =  s*ed; o[:, 4, 3] = -s*ed
    return o

def integrate_one_tau(x, ratios):
    """Integrate all t1/E1 ratios for one tau value in parallel."""
    ny = len(ratios)
    t1amp = E1 * 10.0**ratios
    R = np.broadcast_to(np.eye(5), (ny, 5, 5)).copy()
    # q is time measured in 100/meV; each segment has length x.
    for cyc in range(2):
        for seg in range(3):
            if x == 0:
                continue
            nstep = max(1, int(np.ceil(x/DQ)))
            h = x/nstep
            for k in range(nstep):
                u0 = k/nstep
                u1 = (k+0.5)/nstep
                u2 = (k+1.0)/nstep
                def vals(u):
                    if seg == 0:
                        return TC*env(u, True), 0., E0*env(u, False), env(u, True)
                    if seg == 1:
                        return TC*env(u, False), TC*env(u, True), 0., env(u, False)
                    return 0., TC*env(u, False), E0*env(u, True), 0.
                a = vals(u0); b = vals(u1); c = vals(u2)
                O1 = 100.0*omega_batch(a[0], a[1], a[2], t1amp*a[3])
                O2 = 100.0*omega_batch(b[0], b[1], b[2], t1amp*b[3])
                O3 = O2
                O4 = 100.0*omega_batch(c[0], c[1], c[2], t1amp*c[3])
                K1 = np.einsum('nij,njk->nik', O1, R)
                K2 = np.einsum('nij,njk->nik', O2, R + 0.5*h*K1)
                K3 = np.einsum('nij,njk->nik', O3, R + 0.5*h*K2)
                K4 = np.einsum('nij,njk->nik', O4, R + h*K3)
                R += h*(K1 + 2*K2 + 2*K3 + K4)/6.0
    # Remove accumulated roundoff drift with the polar orthogonal factor.
    u, _, vh = np.linalg.svd(R)
    R = np.einsum('nij,njk->nik', u, vh)
    # Keep SO(5), correcting the rare numerical det=-1 branch.
    det = np.linalg.det(R)
    bad = det < 0
    if np.any(bad):
        u[bad, :, -1] *= -1
        R[bad] = np.einsum('nij,njk->nik', u[bad], vh[bad])
    M = R[:, :3, :3]
    N = R[:, 3:, :3]
    Fsurv = np.sum(M*M, axis=(1,2))/3.0
    # Proper orthogonal polar factor of each 3x3 projected block.
    u, _, vh = np.linalg.svd(M)
    Q = np.einsum('nij,njk->nik', u, vh)
    detq = np.linalg.det(Q)
    bad = detq < 0
    if np.any(bad):
        u[bad, :, -1] *= -1
        Q[bad] = np.einsum('nij,njk->nik', u[bad], vh[bad])
    Frot = (1.0 + np.einsum('ij,nji->n', B, Q))/4.0
    return Fsurv, Frot, Fsurv*Frot

def main():
    xs = np.linspace(0.0, XMAX, N_X)
    ys = np.linspace(YMIN, YMAX, N_Y)
    fs = np.empty((N_Y, N_X)); surv = np.empty_like(fs); rot = np.empty_like(fs)
    for j, x in enumerate(xs):
        surv[:, j], rot[:, j], fs[:, j] = integrate_one_tau(x, ys)
        if j % 10 == 0:
            print(f'{j+1}/{N_X}: tau/(100/meV)={x:.2f}', flush=True)
    np.savez('so5_group_fidelity_map.npz', x=xs, y=ys, F_group=fs,
             F_surv=surv, F_rot_cond=rot, E1=E1, tc=TC, E0=E0)
    fig, ax = plt.subplots(figsize=(7.0, 4.8), constrained_layout=True)
    im = ax.pcolormesh(xs, ys, fs, shading='auto', cmap='viridis', vmin=0, vmax=1)
    cs = ax.contour(xs, ys, fs, levels=[0.2, 0.4, 0.6, 0.8], colors='w', linewidths=.35, alpha=.7)
    ax.clabel(cs, inline=True, fontsize=7, fmt='%.1f')
    ax.set_xlabel(r'$\tau\ /(100\,\mathrm{meV}^{-1})$')
    ax.set_ylabel(r'$\log_{10}(t_1/E_1)$')
    ax.set_title(r'$SO(5)$ group-space braiding fidelity, $E_1=0.01$ meV')
    fig.colorbar(im, ax=ax, label=r'$F_{\rm group}=F_{\rm surv}F_{\rm rot}^{\rm cond}$')
    fig.savefig('so5_group_fidelity_map.png', dpi=220)
    print('wrote so5_group_fidelity_map.png and so5_group_fidelity_map.npz')

if __name__ == '__main__':
    main()
