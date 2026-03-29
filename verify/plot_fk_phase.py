import numpy as np
import matplotlib.pyplot as plt
from math import sqrt

# Honeycomb delta vectors (conventional choices)
delta_x = np.array([0.0, -1.0])
delta_y = np.array([sqrt(3)/2, 0.5])
delta_z = np.array([-sqrt(3)/2, 0.5])

def f_of_k(kx, ky, Jx, Jy, Jz):
    ex = np.exp(1j*(kx*delta_x[0] + ky*delta_x[1]))
    ey = np.exp(1j*(kx*delta_y[0] + ky*delta_y[1]))
    ez = np.exp(1j*(kx*delta_z[0] + ky*delta_z[1]))
    return Jx*ex + Jy*ey + Jz*ez

# Compute min |f(k)| on a k-grid
def min_abs_f(Jx, Jy, Jz, nk=201):
    kxs = np.linspace(-np.pi, np.pi, nk)
    kys = np.linspace(-np.pi, np.pi, nk)
    KX, KY = np.meshgrid(kxs, kys, indexing='xy')
    F = f_of_k(KX, KY, Jx, Jy, Jz)
    return np.min(np.abs(F)), np.abs(F)

# Barycentric plotting helpers
Vx = np.array([0.0, 0.0])
Vy = np.array([1.0, 0.0])
Vz = np.array([0.5, sqrt(3)/2])

def barycentric_point(Jx, Jy, Jz):
    return Jx*Vx + Jy*Vy + Jz*Vz

# Phase diagram over simplex
def compute_phase_diagram(nr=101, nk=121, tol=1e-3):
    rs = np.linspace(0,1,nr)
    ss = np.linspace(0,1,nr)
    pts = []
    gapless = []
    for r in rs:
        for s in ss:
            if r + s <= 1.0:
                Jx = r
                Jy = s
                Jz = 1.0 - r - s
                m, _ = min_abs_f(Jx, Jy, Jz, nk=nk)
                pts.append(barycentric_point(Jx, Jy, Jz))
                gapless.append(m < tol)
    P = np.array(pts)
    G = np.array(gapless)
    return P, G

if __name__ == '__main__':
    # Example |f(k)| for isotropic point
    Jx = Jy = Jz = 1.0/3.0
    m, F = min_abs_f(Jx, Jy, Jz, nk=301)
    print('Isotropic min|f(k)| =', m)
    plt.figure(figsize=(6,5))
    kxs = np.linspace(-np.pi, np.pi, F.shape[0])
    kys = np.linspace(-np.pi, np.pi, F.shape[1])
    plt.contourf(kxs, kys, F, levels=100, cmap='magma')
    plt.colorbar(label='|f(k)|')
    plt.title(f'|f(k)| at J=(1/3,1/3,1/3), min={m:.3e}')
    plt.xlabel('kx')
    plt.ylabel('ky')
    plt.tight_layout()
    plt.savefig('verify/f_isotropic.png', dpi=200)
    print('Saved verify/f_isotropic.png')

    # Phase diagram
    print('Computing phase diagram (this may take a few seconds)...')
    P, G = compute_phase_diagram(nr=121, nk=121, tol=1e-3)
    plt.figure(figsize=(6,5))
    plt.scatter(P[~G,0], P[~G,1], c='tab:blue', s=6, label='gapped')
    plt.scatter(P[G,0], P[G,1], c='tab:orange', s=8, label='gapless')
    # plot triangle vertices and labels
    plt.plot([Vx[0], Vy[0], Vz[0], Vx[0]], [Vx[1], Vy[1], Vz[1], Vx[1]], 'k-')
    plt.text(Vx[0]-0.03, Vx[1]-0.03, 'Jx=1')
    plt.text(Vy[0]+0.02, Vy[1]-0.03, 'Jy=1')
    plt.text(Vz[0], Vz[1]+0.02, 'Jz=1', ha='center')

    # analytic boundary lines: J_alpha = 1/2
    # For each line fix J_alpha=0.5 and let the other two vary summing to 0.5
    nline = 201
    # Jx = 1/2 line
    line_pts = [barycentric_point(0.5, y, 0.5 - y) for y in np.linspace(0,0.5,nline)]
    L = np.array(line_pts)
    plt.plot(L[:,0], L[:,1], 'k--', lw=1.2, label='analytic boundary')
    # Jy = 1/2 line
    line_pts = [barycentric_point(x, 0.5, 0.5 - x) for x in np.linspace(0,0.5,nline)]
    L = np.array(line_pts)
    plt.plot(L[:,0], L[:,1], 'k--', lw=1.2)
    # Jz = 1/2 line
    line_pts = [barycentric_point(x, 0.5 - x, 0.5) for x in np.linspace(0,0.5,nline)]
    L = np.array(line_pts)
    plt.plot(L[:,0], L[:,1], 'k--', lw=1.2)

    plt.title('Phase diagram (blue=gapped, orange=gapless)')
    plt.axis('equal')
    plt.axis('off')
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig('verify/phase_diagram.png', dpi=200)
    print('Saved verify/phase_diagram.png')
