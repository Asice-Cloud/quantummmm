#!/usr/bin/env python3
"""Fig 1(d) reproduction — high-res smoothed with contour lines."""
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.ndimage import gaussian_filter, zoom
import time; pi = np.pi

# ── Parameters ──
tc = E0 = 0.3;  E1_FIXED = 0.001
# τ in units of 100/meV (as in paper Fig 1d)
# tau_c = tau_p * 100 is the actual τ in meV⁻¹

# ── JW operators ──
sp = np.array([[0.,0.],[1.,0.]], complex); sm = np.array([[0.,1.],[0.,0.]], complex)
sz = np.array([[1.,0.],[0.,-1.]], complex); I2 = np.eye(2, dtype=complex)
c1  = np.kron(np.kron(sm,I2),I2); c1d = np.kron(np.kron(sp,I2),I2)
c2  = np.kron(np.kron(sz,sm),I2); c2d = np.kron(np.kron(sz,sp),I2)
c3  = np.kron(np.kron(sz,sz),sm); c3d = np.kron(np.kron(sz,sz),sp)
g1=c1+c1d; g2=1j*(c1-c1d); g3=c2+c2d; ga=c3+c3d; gb=1j*(c3-c3d)
def hm(M): return (M+M.conj().T)*0.5
H_Ed=hm(1j*(ga@gb)); H_E1=hm(1j*(g1@g2)); H_t2=hm(1j*(ga@g2))
H_t1=hm(-1j*(gb@g1)); H_t3=hm(-1j*(ga@g3))
print("Pre-computed ✓")

# ── Protocol ──
def pr(s,tau): return 0.5*(1-np.cos(pi*s/tau))
def pf(s,tau): return 0.5*(1+np.cos(pi*s/tau))
def cp(t,tau,t1m):
    T=6*tau; tm=t%T; k=int(tm//tau); s=tm-k*tau
    if k in (0,3): return (E0*pf(s,tau), t1m*pr(s,tau), tc*pr(s,tau), 0.)
    elif k in (1,4): return (0., t1m*pf(s,tau), tc*pf(s,tau), tc*pr(s,tau))
    else: return (E0*pr(s,tau), 0., 0., tc*pf(s,tau))

def evolve(psi,tau,t1m,ns):
    dt=tau/ns; psi=psi.copy()
    for i in range(6*ns):
        t=i*dt
        e0,t10,t20,t30=cp(t,tau,t1m)
        em,t1m2,t2m2,t3m2=cp(t+.5*dt,tau,t1m)
        e1,t11,t21,t31=cp(t+dt,tau,t1m)
        H0=e0*H_Ed+E1_FIXED*H_E1+t10*H_t1+t20*H_t2+t30*H_t3
        Hm=em*H_Ed+E1_FIXED*H_E1+t1m2*H_t1+t2m2*H_t2+t3m2*H_t3
        H1=e1*H_Ed+E1_FIXED*H_E1+t11*H_t1+t21*H_t2+t31*H_t3
        k1=dt*(-1j*(H0@psi))
        k2=dt*(-1j*(Hm@(psi+.5*k1)))
        k3=dt*(-1j*(Hm@(psi+.5*k2)))
        k4=dt*(-1j*(H1@(psi+k3)))
        psi=psi+(k1+2*k2+2*k3+k4)/6.
        psi/=np.sqrt(np.real(psi.conj()@psi))
    return psi

# ── Init & fidelity ──
psi0=np.zeros(8,complex); psi0[4]=1.  # |100>
P0=np.zeros((8,8),complex)
for i in range(4): P0[i,i]=1.
def fid(psi): return np.real(psi.conj()@P0@psi)

# ── Sanity ──
E1s=E1_FIXED; E1_FIXED=0.
print(f"Sanity: {fid(evolve(psi0,500.,0.,200)):.8f}")
E1_FIXED=E1s

# ── Sweep ──
N_TAU, N_T1 = 80, 30
tv = np.linspace(0.2, 12.0, N_TAU); td = tv
t1v = E1_FIXED * 10**np.linspace(-1.0, 1.0, N_T1)  # log scale: t₁/E₁ ∈ [0.1, 10]
lg_vals = np.linspace(-1.0, 1.0, N_T1)  # for axis display
grid = np.zeros((N_T1, N_TAU))
print(f"\n{N_T1}×{N_TAU} sweep..."); t0 = time.time()
for j, tt in enumerate(t1v):
    for i, tau_p in enumerate(tv):
        tau = tau_p * 100.0
        ns = max(80, min(int(tau / 1.0), 1500))
        grid[j, i] = fid(evolve(psi0, tau, tt, ns))
    if (j + 1) % 8 == 0:
        print(f"  {j+1}/{N_T1} {time.time()-t0:.1f}s")
print(f"Done: {time.time()-t0:.1f}s")

# ── Smooth + upsample ──
grid_s = gaussian_filter(grid, sigma=0.3, mode='nearest')
ZOOM = 4
grid_z = zoom(grid_s, ZOOM, order=3)
tau_z = np.linspace(0.2, 12.0, N_TAU * ZOOM)
lg_z  = np.linspace(-1.0, 1.0, N_T1 * ZOOM)

# ── Plot ──
colors = ['#0d0887', '#46039f', '#7201a8', '#9711a1', '#b4298c',
          '#c94d71', '#d76e56', '#de8d3e', '#e8ab31', '#f0c92b', '#fae724']
cmap = LinearSegmentedColormap.from_list('paper', colors, N=256)

fig, ax = plt.subplots(figsize=(8, 6))
Tz, Lz = np.meshgrid(tau_z, lg_z)

# Filled contours
levels = np.linspace(0, 1, 25)
cf = ax.contourf(Tz, Lz, grid_z, levels=levels, cmap=cmap, extend='both')

# White contour lines at 0.1, 0.3, 0.5, 0.7, 0.9
cs = ax.contour(Tz, Lz, grid_z, levels=[0.1, 0.3, 0.5, 0.7, 0.9],
                colors='white', linewidths=0.5, alpha=0.6)
ax.clabel(cs, inline=True, fontsize=7, fmt='%.1f')

cb = plt.colorbar(cf, ax=ax, label=r'$|\langle\psi_{1-}(6\tau)|\psi_{1+}(0)\rangle|$',
                   pad=0.02)
cb.set_ticks([0, 0.25, 0.5, 0.75, 1])

ax.set_xlabel(r'$\tau$ $(100/\mathrm{meV})$', fontsize=13)
ax.set_ylabel(r'$\lg(t_1/E_1)$', fontsize=13)
ax.set_title(f'Fig 1(d): $E_1={E1_FIXED}$ meV, $t_c=E_0={tc}$ meV', fontsize=12)
ax.set_xlim(0.2, 12.0); ax.set_ylim(-1.0, 1.0)
plt.tight_layout()
plt.savefig('fig1d_final.png', dpi=200)
plt.savefig('fig1d_final.pdf')
print("\nSaved.")
