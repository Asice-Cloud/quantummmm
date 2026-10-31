"""Exact parameter jets and SO(5) channels of the unified cosine function.

The base evolution is never expanded in E or eta.  A finite jet computes
exact derivatives at a parameter point by a closed triangular lift of the
same differential equation.  Taylor continuation has a numerical order,
independent of the number of parameter derivatives requested.
"""

import math
import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

if __package__:
    from .full_model import LE, LH, segment_coefficients
else:
    from full_model import LE, LH, segment_coefficients


VARIABLES = ("tau", "E", "eta")
INDICES = ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1),
           (2, 0, 0), (1, 1, 0), (1, 0, 1),
           (0, 2, 0), (0, 1, 1), (0, 0, 2))
INDEX = {a: n for n, a in enumerate(INDICES)}
LABELS = ("value", "tau", "E", "eta", "tau_tau", "tau_E", "tau_eta",
          "E_E", "E_eta", "eta_eta")
ZETA = np.array([1., -1j, 0., 0., 0.])/np.sqrt(2.)

_left, _right, _weights, _starts = [], [], [], []
for alpha in INDICES:
    _starts.append(len(_left))
    for left_index, beta in enumerate(INDICES):
        gamma = tuple(a-b for a, b in zip(alpha, beta))
        if min(gamma) < 0 or gamma not in INDEX:
            continue
        _left.append(left_index)
        _right.append(INDEX[gamma])
        _weights.append(math.prod(math.comb(a, b) for a, b in zip(alpha, beta)))
_LEFT = np.array(_left)
_RIGHT = np.array(_right)
_WEIGHTS = np.array(_weights, dtype=float)[:, None, None]
_STARTS = np.array(_starts)


def jet_product(left, right):
    """Leibniz product of raw derivatives, not factorial-normalized jets."""
    products = (left[_LEFT]@right[_RIGHT])*_WEIGHTS
    return np.add.reduceat(products, _STARTS, axis=0)


def conditioning_scales(tau):
    """Coordinate units for intermediate derivative arithmetic only."""
    return np.array([1., 1./max(1., abs(tau)), 1./max(1., abs(tau))])


def _powers(scales):
    return np.prod(scales[None, :]**np.array(INDICES), axis=1)


def coefficient_jets(kind, tau, eta, tc=.3, E0=.3, E=.01, scales=None):
    """Derivatives of P=tau*X and Q=tau*Y in scaled local coordinates."""
    scales = conditioning_scales(tau) if scales is None else np.asarray(scales)
    x, y = segment_coefficients(kind, eta, tc, E0, E)
    x_eta = LH if kind in (1, 2) else np.zeros((5, 5))
    y_eta = (-LH if kind == 1 else LH) if kind in (1, 2) else np.zeros((5, 5))
    p, q = np.zeros((10, 5, 5), complex), np.zeros((10, 5, 5), complex)
    p[0], q[0] = tau*x, tau*y
    p[1], q[1] = x, y
    p[2], p[3], q[3] = 2*tau*LE, tau*x_eta, tau*y_eta
    p[5], p[6], q[6] = 2*LE, x_eta, y_eta
    factors = _powers(scales)[:, None, None]
    return p*factors, q*factors


def _lift(jet):
    """Triangular matrix implementing multiplication by a coefficient jet."""
    block = np.zeros((50, 50), dtype=complex)
    for output, alpha in enumerate(INDICES):
        for column, beta in enumerate(INDICES):
            gamma = tuple(a-b for a, b in zip(alpha, beta))
            if min(gamma) < 0 or gamma not in INDEX:
                continue
            weight = math.prod(math.comb(a, b) for a, b in zip(alpha, beta))
            block[5*output:5*output+5, 5*column:5*column+5] = weight*jet[INDEX[gamma]]
    return block


def segment_response(kind, tau, eta, tc=.3, E0=.3, E=.01,
                     order=28, local_phase=2., scales=None):
    """One segment and all parameter derivatives of total order <=2.

    Uses the rational-z recurrence of full_model, with each
    coefficient in the finite jet algebra.  This is the first block column
    of a 50-dimensional cosine connection, without dense 50D products.
    Returned derivatives use ``scales`` (unscaled by protocol_response).
    """
    if tau < 0 or order < 4 or local_phase <= 0:
        raise ValueError("nonnegative tau, order>=4, positive local_phase required")
    p, q = coefficient_jets(kind, tau, eta, tc, E0, E, scales)
    scale = np.linalg.norm(_lift(p), 2)+np.linalg.norm(_lift(q), 2)
    panels = max(18, int(np.ceil(scale/local_phase)))
    b2, b1 = q/(2j*np.pi), p/(1j*np.pi)
    identity = np.zeros_like(p)
    identity[0] = np.eye(5)
    out = identity.copy()
    points = np.exp(1j*np.linspace(0., np.pi, panels+1))
    for z0, z1 in zip(points[:-1], points[1:]):
        step = z1-z0
        d0 = b2*(z0*z0+1.)+b1*z0
        d1 = 2*b2*z0+b1
        t, previous, previous2 = identity.copy(), np.zeros_like(p), np.zeros_like(p)
        local = t.copy()
        for n in range(order):
            following = step*(jet_product(d0, t)-2*n*z0*t)
            if n:
                following += step**2*(jet_product(d1, previous)-(n-1)*previous)
            if n > 1:
                following += step**3*jet_product(b2, previous2)
            following /= (n+1)*z0*z0
            local += following
            previous2, previous, t = previous, t, following
        out = jet_product(local, out)
    return out, panels


def protocol_response(tau, eta, tc=.3, E0=.3, E=.01, **options):
    """Return jets of W, R=W^2, B and F in physical parameter units."""
    scales = conditioning_scales(tau)
    segments, panels = [], []
    for kind in (1, 2, 3):
        value, n = segment_response(kind, tau, eta, tc, E0, E, scales=scales, **options)
        segments.append(value)
        panels.append(n)
    cycle = jet_product(segments[2], jet_product(segments[1], segments[0]))
    rotation = jet_product(cycle, cycle)
    bivector = np.zeros_like(rotation)
    bivector[0] = LE
    bivector = jet_product(rotation, jet_product(bivector, rotation.swapaxes(-1, -2)))
    amplitude = np.einsum("i,nij,j->n", ZETA, rotation, ZETA)
    fidelity = jet_product(amplitude[:, None, None], amplitude.conj()[:, None, None])[:, 0, 0]
    factors = _powers(scales)
    return {"W": cycle/factors[:, None, None],
            "R": rotation/factors[:, None, None],
            "B": bivector/factors[:, None, None],
            "amplitude": amplitude/factors, "F": fidelity.real/factors,
            "panels": panels, "imaginary_scaled_R": float(abs(rotation.imag).max())}


def so5_channels(cycle, tolerance=1e-7):
    """Two rotation-plane channels, only on a separated spectral chart.

    Angles are principal values, not unwrapped dynamical phases.  Labels
    may exchange at degeneracy; the full matrix amplitude is always valid.
    """
    w = np.asarray(cycle).real
    h = (w+w.T)/2
    total = (np.trace(w)-1)/2
    square_sum = (np.trace(w@w)+3)/4
    gap2 = 2*square_sum-total**2
    if gap2 < -1e-9:
        raise ValueError("input fails SO(5) spectral invariants")
    gap = np.sqrt(max(0., gap2))
    c = np.array([(total+gap)/2, (total-gap)/2])
    if gap < tolerance or np.min(abs(1-c)) < tolerance:
        raise ValueError("degenerate spectral chart; use full-matrix readout")
    identity = np.eye(5)
    projectors = [((h-c[0]*identity)@(h-c[1]*identity))/((1-c[0])*(1-c[1]))]
    for j in range(2):
        other = 1-j
        projectors.append(((h-identity)@(h-c[other]*identity))
                          /((c[j]-1)*(c[j]-c[other])))
    projectors = np.array(projectors)
    weights = np.einsum("i,nij,j->n", ZETA, projectors, ZETA)
    channels = -2*weights[1:]*(1-c*c)
    amplitude = np.einsum("i,ij,j", ZETA, w@w, ZETA)
    diagonal = abs(channels)**2
    interference = 2*np.real(channels[0]*channels[1].conj())
    return {"cosines": c, "angles": np.arccos(np.clip(c, -1., 1.)),
            "weights": weights, "projectors": projectors,
            "channels": channels, "amplitude": amplitude,
            "F_diagonal": diagonal, "F_interference": float(interference),
            "F_reconstructed": float(diagonal.sum()+interference),
            "amplitude_error": float(abs(amplitude-channels.sum()))}


def angle_response(cycle, cycle_derivative):
    """Principal-angle response on the regular chart used by so5_channels."""
    w, dw = np.asarray(cycle).real, np.asarray(cycle_derivative).real
    channels = so5_channels(w)
    c = channels["cosines"]
    total = c.sum()
    dtotal = np.trace(dw)/2
    dgap2 = np.trace(w@dw)-2*total*dtotal
    dc = np.array([dtotal/2+dgap2/(4*(c[0]-c[1])),
                   dtotal/2-dgap2/(4*(c[0]-c[1]))])
    sine = np.sqrt(np.maximum(0., 1-c*c))
    if np.min(sine) < 1e-7:
        raise ValueError("principal angle has a singular coordinate at 0 or pi")
    return -dc/sine


def peak_response(response):
    """Implicit-function slopes of a nondegenerate stationary F(tau)."""
    f = response["F"]
    if abs(f[4]) < 1e-14:
        raise ValueError("zero curvature; a local peak branch is not determined")
    return {"F": float(f[0]), "F_tau": float(f[1]), "F_tau_tau": float(f[4]),
            "tau_E": float(-f[5]/f[4]), "tau_eta": float(-f[6]/f[4]),
            "height_E": float(f[2]), "height_eta": float(f[3])}
