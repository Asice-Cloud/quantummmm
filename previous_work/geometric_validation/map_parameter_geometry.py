#!/usr/bin/env python3
"""Map finite E1/eta geometry using forward states and backward targets.

Run from the project root:
  .venv/bin/python -B geometric_validation/map_parameter_geometry.py

Full ten-component state equations are integrated in parallel over parameter
points. The spherical interpretation is evaluated afterwards, not used in
propagation. Refine contours and sampled ridge maxima with the full equations.
"""

import json
import os
from pathlib import Path
import sys

# Keep large flattened ODE vectors from triggering excessive BLAS threading.
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
os.environ.setdefault('MPLCONFIGDIR', '/tmp/quantumsss-matplotlib')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import linear_sum_assignment

from verify_bottleneck import B_INITIAL, INDICES, PSI_INITIAL, PSI_TARGET, controls


ROOT = Path(__file__).resolve().parent
TAU = 100.
RTOL, ATOL = 2e-11, 2e-13
E_AXIS = np.linspace(0, .03, 61)
ETA_AXIS = np.linspace(0, .05, 101)


def state_rhs(stage, v, y, E, eta):
    d, x, a, b, u, c, w, h, j, k = np.moveaxis(y, -1, 0)
    p, r, ed, gate = controls(stage, v)
    q = eta*gate
    return 2*TAU*np.stack((
        -p*a-q*w,
        E*u+r*a-q*j,
        E*c-r*x+ed*b+p*d-q*k,
        E*w-ed*a,
        -E*x+r*c+p*h,
        -E*a-r*u+ed*w,
        -E*b-ed*c+q*d-p*k,
        ed*j-p*u,
        r*k-ed*h+q*x,
        -r*j+q*a+p*w,
    ), axis=-1)


def fields_matrix(y):
    # State order: d, ReZ3, ReZa, ReZb, ImZ3, ImZa, ImZb, h3a, h3b, hab.
    indices = ((0, 1), (0, 2), (0, 3), (0, 4), (1, 2),
               (1, 3), (1, 4), (2, 3), (2, 4), (3, 4))
    B = np.zeros(y.shape[:-1]+(5, 5))
    for pos, (i, j) in enumerate(indices):
        B[..., i, j], B[..., j, i] = y[..., pos], -y[..., pos]
    return B


def evolve(initial, E, eta, stages, backward=False, sensitivity=False):
    y = initial.copy()
    shape = y.shape
    def rhs(stage, v, flat):
        rows = flat.reshape(shape)
        if not sensitivity:
            return state_rhs(stage, v, rows, E, eta).ravel()
        result = state_rhs(stage, v, rows, E[:, None], eta[:, None])
        d, x, a, b, u, c, w, h, j, k = rows[:, 0].T
        z = np.zeros_like(d)
        result[:, 1] += 2*TAU*np.stack((z, u, c, w, -x, -a, -b, z, z, z), axis=-1)
        gate = controls(stage, v)[3]
        result[:, 2] += 2*TAU*gate*np.stack((-w, -j, -k, z, z, z, d, z, x, a), axis=-1)
        return result.ravel()
    for stage in stages:
        sol = solve_ivp(lambda v, flat: rhs(stage, v, flat),
                        (1., 0.) if backward else (0., 1.), y.ravel(),
                        method='DOP853', rtol=RTOL, atol=ATOL, max_step=.025)
        assert sol.success, sol.message
        y = sol.y[:, -1].reshape(shape)
    return y


def evaluate(E, eta, sensitivity=False):
    E, eta = np.broadcast_arrays(np.atleast_1d(E), np.atleast_1d(eta))
    E, eta = E.ravel(), eta.ravel()
    initial = np.zeros((len(E), 3, 10) if sensitivity else (len(E), 10))
    if sensitivity:
        initial[:, 0, 0] = 1
    else:
        initial[:, 0] = 1
    b3 = evolve(initial, E, eta, (0, 1, 2), sensitivity=sensitivity)
    q3 = evolve(-initial, E, eta, (4, 3), backward=True, sensitivity=sensitivity)
    if sensitivity:
        b, q = b3[:, 0], q3[:, 0]
    else:
        b, q = b3, q3
    B, Q = fields_matrix(b), fields_matrix(q)
    P, Pq = -B@B, -Q@Q
    product = np.einsum('nij,nij->n', P, Pq)
    inner = np.sum(b*q, axis=-1)
    F, d = (product+2*inner)/4, -inner
    bi, qi = b[:, [0, 1, 4]], q[:, [0, 1, 4]]
    beta, chi = np.linalg.norm(bi, axis=1), np.linalg.norm(qi, axis=1)
    cos_theta = np.clip(np.sum(bi*qi, axis=1)/(beta*chi), -1, 1)
    sphere_F = (1+cos_theta)**2/4
    azimuth = np.arctan2(bi[:, 2], bi[:, 1])-np.arctan2(qi[:, 2], qi[:, 1])
    azimuth = np.arctan2(np.sin(azimuth), np.cos(azimuth))
    k, kq = np.clip(P[:, 2, 2], 0, 1), np.clip(Pq[:, 2, 2], 0, 1)
    cap = (1+np.sqrt((1-k)*(1-kq))+np.sqrt(k*kq))**2/4
    result = {'F': F, 'd': d, 'B3': B, 'Q3': Q, 'beta': beta, 'chi': chi,
              'cos_theta': cos_theta, 'theta_over_pi': np.arccos(cos_theta)/np.pi,
              'sphere_F': sphere_F, 'sphere_error': abs(F-sphere_F),
              'ancilla_bivector_norm': np.sqrt(np.maximum(0, 1-np.minimum(beta, chi)**2)),
              'k': k, 'k_required': kq, 'weight_cap': cap,
              'azimuth_over_pi': azimuth/np.pi,
              'actual_latitude': np.arcsin(np.clip(bi[:, 0]/beta, -1, 1)),
              'required_latitude': np.arcsin(np.clip(qi[:, 0]/chi, -1, 1)),
              'position_loss': np.sum((P-Pq)**2, axis=(1, 2))/8,
              'orientation_loss': np.sum((B-Q)**2, axis=(1, 2))/8}
    result['latitude_factor'] = np.cos(result['actual_latitude'])**4
    result['phase_factor'] = np.cos(azimuth/2)**4
    result['factorized_F'] = result['latitude_factor']*result['phase_factor']
    if sensitivity:
        Sb, Sq = fields_matrix(b3[:, 1:]), fields_matrix(q3[:, 1:])
        dP = -(Sb@B[:, None]+B[:, None]@Sb)
        dPq = -(Sq@Q[:, None]+Q[:, None]@Sq)
        d_inner = np.sum(b3[:, 1:]*q[:, None]+b[:, None]*q3[:, 1:], axis=-1)
        d_prod = (np.sum(dP*Pq[:, None], axis=(-1, -2))
                  +np.sum(P[:, None]*dPq, axis=(-1, -2)))
        result['F_gradient'] = (d_prod+2*d_inner)/4
        result['d_gradient'] = -d_inner
        db, dq = b3[:, 1:], q3[:, 1:]
        result['phase_gradient'] = (
            (b[:, None, 1]*db[:, :, 4]-b[:, None, 4]*db[:, :, 1])
            /(b[:, None, 1]**2+b[:, None, 4]**2)
            -(q[:, None, 1]*dq[:, :, 4]-q[:, None, 4]*dq[:, :, 1])
            /(q[:, None, 1]**2+q[:, None, 4]**2))
    return result


def contour_brackets(grid, peaks):
    brackets = []
    for key, level in (('d', 0.), ('F', .9), ('F', .99)):
        values = grid[key]
        for col, E in enumerate(E_AXIS):
            for row in range(len(ETA_AXIS)-1):
                low, high = values[row, col]-level, values[row+1, col]-level
                if low*high < 0 or low == 0:
                    brackets.append((E, ETA_AXIS[row], ETA_AXIS[row+1], key, level))
    # A narrow high-fidelity island can peak above a threshold between samples.
    # Use refined local peaks to supply both brackets if the grid missed it.
    for peak in peaks:
        col = np.argmin(abs(E_AXIS-peak['E']))
        cell = np.clip(np.searchsorted(ETA_AXIS, peak['eta'])-1, 0, len(ETA_AXIS)-2)
        for level in (.9, .99):
            if peak['F'] > level and max(grid['F'][cell:cell+2, col]) < level:
                brackets.extend(((peak['E'], ETA_AXIS[cell], peak['eta'], 'F', level),
                                 (peak['E'], peak['eta'], ETA_AXIS[cell+1], 'F', level)))
    phase = np.sin(np.pi*grid['azimuth_over_pi'])
    for col, E in enumerate(E_AXIS):
        for row in range(len(ETA_AXIS)-1):
            if phase[row, col]*phase[row+1, col] < 0:
                brackets.append((E, ETA_AXIS[row], ETA_AXIS[row+1], 'phase_alignment', 0.))
    return brackets


def refine_contours(grid, peaks):
    brackets = contour_brackets(grid, peaks)
    E = np.array([b[0] for b in brackets])
    lo, hi = np.array([b[1] for b in brackets]), np.array([b[2] for b in brackets])
    index = np.array([{'d': 0, 'F': 1, 'phase_alignment': 2}[b[3]] for b in brackets])
    level = np.array([b[4] for b in brackets])
    initial = evaluate(E, lo)
    def residual(result):
        return np.where(index == 0, result['d'],
            np.where(index == 1, result['F'],
                     np.sin(np.pi*result['azimuth_over_pi'])))-level
    left = residual(initial)
    for _ in range(19):
        mid = (lo+hi)/2
        value = residual(evaluate(E, mid))
        replace = left*value > 0
        lo, left = np.where(replace, mid, lo), np.where(replace, value, left)
        hi = np.where(replace, hi, mid)
    eta = (lo+hi)/2
    result = evaluate(E, eta, sensitivity=True)
    entries = []
    for i, bracket in enumerate(brackets):
        if index[i] == 2 and np.cos(np.pi*result['azimuth_over_pi'][i]) < 0:
            continue
        key = ('d_gradient', 'F_gradient', 'phase_gradient')[index[i]]
        grad = result[key][i]
        entries.append({'E': float(E[i]), 'eta': float(eta[i]),
            'quantity': bracket[3], 'level': bracket[4],
            'residual': float(residual(result)[i]),
            'gradient_E': float(grad[0]), 'gradient_eta': float(grad[1]),
            'deta_dE': float(-grad[0]/grad[1]) if abs(grad[1]) > 1e-8 else None,
            'eta_bracket_width': float(hi[i]-lo[i])})
    return entries


def refine_ridges(grid):
    candidates = []
    values = grid['F']
    for col, E in enumerate(E_AXIS):
        for row in range(1, len(ETA_AXIS)-1):
            if values[row, col] >= max(values[row-1, col], values[row+1, col]):
                candidates.append((E, ETA_AXIS[row-1], ETA_AXIS[row+1]))
    E = np.array([c[0] for c in candidates])
    lo, hi = np.array([c[1] for c in candidates]), np.array([c[2] for c in candidates])
    ratio = (np.sqrt(5)-1)/2
    a, b = hi-ratio*(hi-lo), lo+ratio*(hi-lo)
    fa, fb = evaluate(E, a)['F'], evaluate(E, b)['F']
    for _ in range(25):
        choose = fa < fb
        lo, hi = np.where(choose, a, lo), np.where(choose, hi, b)
        a, b = hi-ratio*(hi-lo), lo+ratio*(hi-lo)
        fa, fb = evaluate(E, a)['F'], evaluate(E, b)['F']
    eta = (lo+hi)/2
    exact = evaluate(E, eta)
    peaks = [{'E': float(E[i]), 'eta': float(eta[i]), 'F': float(exact['F'][i]),
              'd': float(exact['d'][i]), 'theta_over_pi': float(exact['theta_over_pi'][i])}
             for i in range(len(E))]
    best = []
    for col, e in enumerate(E_AXIS):
        endpoints = [{'E': float(e), 'eta': float(ETA_AXIS[row]),
                      'F': float(values[row, col]), 'd': float(grid['d'][row, col]),
                      'theta_over_pi': float(grid['theta_over_pi'][row, col])}
                     for row in (0, len(ETA_AXIS)-1)]
        best.append(max([p for p in peaks if p['E'] == e]+endpoints, key=lambda p: p['F']))
    return peaks, best


def wave_reference(E, eta):
    current = np.broadcast_to(PSI_INITIAL, (len(E), 5)).copy()
    for stage in range(5):
        def rhs(v, flat):
            psi = flat.reshape(-1, 5)
            p, r, ed, gate = controls(stage, v)
            q = eta*gate
            return (2*TAU*np.stack((E*psi[:, 1]+q*psi[:, 4],
                -E*psi[:, 0]-p*psi[:, 3], r*psi[:, 3],
                p*psi[:, 1]-r*psi[:, 2]+ed*psi[:, 4],
                -q*psi[:, 0]-ed*psi[:, 3]), axis=-1)).ravel()
        solution = solve_ivp(rhs, (0, 1), current.ravel(), method='DOP853',
                             rtol=5e-13, atol=5e-15, max_step=.01)
        assert solution.success
        current = solution.y[:, -1].reshape(-1, 5)
    return abs(current@PSI_TARGET.conj())**2, 2*np.imag(current[:, 0]*current[:, 1].conj())


def verify(grid, roots, best):
    rng = np.random.default_rng(741)
    points = [(float(E_AXIS[c]), float(ETA_AXIS[r]))
              for r, c in zip(rng.integers(0, len(ETA_AXIS), 16),
                               rng.integers(0, len(E_AXIS), 16))]
    points += [(p['E'], p['eta']) for p in best[::10]]
    points += [(p['E'], p['eta']) for p in roots[::max(1, len(roots)//12)]]
    points += [(.01, .01), (.01, .023), (.02, .023)]
    E, eta = np.array(points).T
    geometry = evaluate(E, eta)
    Fw, dw = wave_reference(E, eta)
    checks = {'wave_F': float(np.max(abs(geometry['F']-Fw))),
              'wave_d': float(np.max(abs(geometry['d']-dw))),
              'root_residual': max(abs(p['residual']) for p in roots),
              'loss_identity': float(np.max(abs(
                  grid['position_loss']+grid['orientation_loss']-(1-grid['F']))))}
    for name in ('B3', 'Q3'):
        B = grid[name]
        checks[name+'_plane_constraint'] = float(np.max(abs(B@B@B+B)))
        checks[name+'_norm'] = float(np.max(abs(np.sum(B**2, axis=(-1, -2))-2)))
    for key, value in checks.items():
        assert value < (5e-7 if key == 'root_residual' else 2e-8), (key, value)
    # Exact parameter derivatives on regular contours and on the sampled ridge.
    rows = [best[i] for i in (0, 10, 20, 30, 40, 50, 60)]
    grad = evaluate([p['E'] for p in rows], [p['eta'] for p in rows], sensitivity=True)
    for i, p in enumerate(rows):
        p['F_gradient_E'] = float(grad['F_gradient'][i, 0])
        p['F_gradient_eta'] = float(grad['F_gradient'][i, 1])
    # Independent direct-wave sensitivities from the previous chain tool.
    from verify_second_cycle_chain import direct_wave_sensitivity
    gradient_error = 0.
    for i in (2, 4):
        case = {'E1': rows[i]['E'], 'eta': rows[i]['eta'], 'tau': TAU}
        _, Fg, dg = direct_wave_sensitivity(case)
        gradient_error = max(gradient_error,
            float(np.max(abs(Fg-grad['F_gradient'][i]))),
            float(np.max(abs(dg-grad['d_gradient'][i]))))
    checks['independent_wave_gradient'] = gradient_error
    assert gradient_error < 2e-7
    return {'checked_parameter_points': len(points), 'errors': checks}


def high_fidelity_curves(ax, roots):
    """Join refined slice roots; do not interpolate coarse eta-grid contours."""
    for level, color in ((.9, '#ffb000'), (.99, '#ff4d4d')):
        curves, active = [], []
        for E in E_AXIS:
            points = sorted([p['eta'] for p in roots
                             if p['quantity'] == 'F' and p['level'] == level
                             and abs(p['E']-E) < 1e-12])
            new_active, assigned = [], set()
            if active and points:
                distance = abs(np.array([curves[i][-1][1] for i in active])[:, None]
                               -np.array(points)[None, :])
                left, right = linear_sum_assignment(distance)
                for i, j in zip(left, right):
                    if distance[i, j] < .004:
                        curve_id = active[i]
                        curves[curve_id].append((E, points[j]))
                        new_active.append(curve_id)
                        assigned.add(j)
            for j, eta in enumerate(points):
                if j not in assigned:
                    curves.append([(E, eta)])
                    new_active.append(len(curves)-1)
            active = new_active
        for curve in curves:
            values = np.array(curve)
            ax.plot(values[:, 0], values[:, 1], color=color, lw=1.2,
                    marker='.' if len(curve) == 1 else None, ms=2)


def figures(grid, best, roots):
    x, y = E_AXIS, ETA_AXIS
    extent = (x[0], x[-1], y[0], y[-1])
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), layout='constrained')
    panels = [
        ('F', '(a) Full-state fidelity', 'viridis', 0, 1),
        ('d', '(b) Final orientation', 'RdBu_r', -1, 1),
        ('theta_over_pi', r'(c) Boundary angle $\theta/\pi$', 'viridis_r', 0, 1),
        ('ancilla_bivector_norm', '(d) Bivector outside modes (1,2,3)', 'magma', None, None),
        ('sphere_error', r'(e) Error $|F-F_{S^2}|$', 'magma', None, None),
        ('azimuth_over_pi', r'(f) Relative azimuth $\Delta\varphi/\pi$', 'twilight', -1, 1),
    ]
    for ax, (key, title, cmap, lo, hi) in zip(axes.ravel(), panels):
        opts = {'vmin': lo, 'vmax': hi} if lo is not None else {
            'norm': LogNorm(vmin=max(1e-12, float(np.min(grid[key][grid[key]>0]))),
                            vmax=float(np.max(grid[key])))}
        im = ax.imshow(np.maximum(grid[key], 1e-14) if lo is None else grid[key],
                       origin='lower', extent=extent, aspect='auto', cmap=cmap, **opts)
        ax.set(xlabel=r'$E_1$ (meV)', ylabel=r'$\eta=|t_1|_{\rm amplitude}$ (meV)', title=title)
        fig.colorbar(im, ax=ax, shrink=.85)
    for ax in axes.ravel()[:3]:
        ax.contour(x, y, grid['d'], levels=[0], colors='white', linewidths=.8)
    high_fidelity_curves(axes[0, 0], roots)
    jumps = np.flatnonzero(abs(np.diff([p['eta'] for p in best])) > .003)+1
    for i, piece in enumerate(np.split(np.arange(len(best)), jumps)):
        axes[0, 0].plot([best[j]['E'] for j in piece], [best[j]['eta'] for j in piece],
                       color='black', ls='--', lw=1,
                       label='Best detected ridge' if i == 0 else None)
    phase_roots = [p for p in roots if p['quantity'] == 'phase_alignment']
    axes[1, 2].scatter([p['E'] for p in phase_roots], [p['eta'] for p in phase_roots],
                       s=3, c='black', label=r'$\Delta\varphi=0$ mod $2\pi$')
    axes[1, 2].legend(fontsize=7, loc='upper right')
    for E, eta, label in ((.01, .01, 'L'), (.01, .023, 'H'), (.02, .023, 'E')):
        axes[0, 0].plot(E, eta, 'o', mec='white', mfc='black', ms=5)
        axes[0, 0].annotate(label, (E, eta), xytext=(4, 4), textcoords='offset points', color='white')
    axes[0, 0].legend(loc='upper right', fontsize=8)
    fig.suptitle(r'Finite parameter geometry: $\tau=100$/meV, $t_c=E_0=0.3$ meV; '
                 r'white: $d=0$, orange: $F=0.9$, red: $F=0.99$', fontsize=12)
    fig.savefig(ROOT/'parameter_geometry.png', dpi=180)
    fig.savefig(ROOT/'parameter_geometry.pdf')
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), layout='constrained')
    for ax, key, title, cmap, vmin, vmax in (
            (axes[0], 'F', 'Final squared fidelity', 'viridis', 0, 1),
            (axes[1], 'd', 'Final orientation', 'RdBu_r', -1, 1)):
        im = ax.imshow(grid[key], origin='lower', extent=extent, aspect='auto',
                       cmap=cmap, vmin=vmin, vmax=vmax)
        ax.contour(x, y, grid['d'], levels=[0], colors='white', linewidths=.9)
        ax.set(xlabel=r'$E_1$ (meV)', ylabel=r'$\eta$ (meV)', title=title)
        fig.colorbar(im, ax=ax, shrink=.9)
    high_fidelity_curves(axes[0], roots)
    for E, eta, label in ((.01, .01, 'L'), (.01, .023, 'H'), (.02, .023, 'E')):
        axes[0].plot(E, eta, 'o', mec='white', mfc='black', ms=5)
        axes[0].annotate(label, (E, eta), xytext=(4, 4), textcoords='offset points', color='white')
    fig.suptitle(r'$\tau=100$/meV; white: $d=0$, orange: $F=0.9$, red: $F=0.99$')
    fig.savefig(ROOT/'parameter_fidelity_map.png', dpi=180)
    fig.savefig(ROOT/'parameter_fidelity_map.pdf')
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), layout='constrained')
    for E in (.005, .01, .015, .02, .025, .03):
        col = int(round(E/(E_AXIS[1]-E_AXIS[0])))
        axes[0].plot(ETA_AXIS, grid['F'][:, col], label=f'E={E:.3f}')
    axes[0].axhline(.9, color='.5', ls=':', lw=1)
    axes[0].set(xlabel=r'$\eta$ (meV)', ylabel='F', title='Fidelity sections', ylim=(-.02, 1.02))
    axes[0].legend(fontsize=8)
    axes[1].plot([p['E'] for p in best], [p['F'] for p in best], color='black')
    axes[1].axhline(.9, color='.5', ls=':', lw=1)
    axes[1].set(xlabel=r'$E_1$ (meV)', ylabel='Best F in the scanned eta interval',
                title='Upper envelope of detected ridge maxima', ylim=(-.02, 1.02))
    for ax in axes:
        ax.grid(alpha=.2)
    fig.savefig(ROOT/'parameter_sections.png', dpi=180)
    fig.savefig(ROOT/'parameter_sections.pdf')
    plt.close(fig)


def main():
    Em, etam = np.meshgrid(E_AXIS, ETA_AXIS)
    if '--reuse-grid' in sys.argv:
        with np.load(ROOT/'parameter_geometry.npz') as cache:
            assert np.array_equal(cache['E_axis'], E_AXIS)
            assert np.array_equal(cache['eta_axis'], ETA_AXIS)
            grid = {key: cache[key] for key in cache.files if key not in ('E_axis', 'eta_axis')}
        grid['latitude_factor'] = np.cos(grid['actual_latitude'])**4
        grid['phase_factor'] = np.cos(np.pi*grid['azimuth_over_pi']/2)**4
        grid['factorized_F'] = grid['latitude_factor']*grid['phase_factor']
    else:
        flat = evaluate(Em.ravel(), etam.ravel())
        grid = {key: value.reshape(Em.shape+value.shape[1:]) for key, value in flat.items()}
    print('State/target grid completed:', Em.shape, flush=True)
    np.savez_compressed(ROOT/'parameter_geometry.npz', E_axis=E_AXIS, eta_axis=ETA_AXIS, **grid)
    peaks, best = refine_ridges(grid)
    print('Refined ridge candidates:', len(peaks), flush=True)
    roots = refine_contours(grid, peaks)
    print('Refined contour crossings:', len(roots), flush=True)
    validation = verify(grid, roots, best)
    result = {'tau': TAU, 'tc': .3, 'E0': .3, 'hbar': 1,
        'E_range': [float(E_AXIS[0]), float(E_AXIS[-1])],
        'eta_range': [float(ETA_AXIS[0]), float(ETA_AXIS[-1])],
        'grid_shape_eta_E': list(Em.shape), 'point_count': int(Em.size),
        'sphere_error_max': float(np.max(grid['sphere_error'])),
        'sphere_error_quantiles': np.quantile(grid['sphere_error'], [.5, .9, .99]).tolist(),
        'factorized_F_error_max': float(np.max(abs(grid['F']-grid['factorized_F']))),
        'outside_three_mode_norm_max': float(np.max(grid['ancilla_bivector_norm'])),
        'minimum_weight_cap': float(np.min(grid['weight_cap'])),
        'latitude_opposition_error_max': float(np.max(abs(
            grid['actual_latitude']+grid['required_latitude']))),
        'refined_contours': roots, 'local_maxima': peaks,
        'best_ridge_per_E_slice': best, 'validation': validation}
    (ROOT/'parameter_geometry.json').write_text(
        json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    figures(grid, best, roots)
    print('Sphere error max:', result['sphere_error_max'], flush=True)
    print('Validation:', validation, flush=True)
    print('Selected ridge points:', best[::10], flush=True)


if __name__ == '__main__':
    main()
