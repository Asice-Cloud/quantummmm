"""A normalized matrix function for U'=tau*(X+Y*cos(pi*u))*U.

It is evaluated using a finite-width Taylor recurrence in z=exp(i*pi*u),
continued along the upper unit semicircle. No fit to an image, Magnus
truncation, or expansion in the unknown parameters is used.
"""

import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
from scipy.integrate import solve_ivp


def L(i, j):
    out = np.zeros((5, 5))
    out[i, j], out[j, i] = 1., -1.
    return out


LE, LH, LP, LR, LD = L(0, 1), L(0, 4), L(1, 3), L(2, 3), L(3, 4)


def segment_coefficients(kind, eta, tc=.3, E0=.3, E=.01):
    """Physical-time K(u)=X+Y*cos(pi*u), in the original 5D representation."""
    if kind == 1:
        x = 2*E*LE+eta*LH-tc*LP+E0*LD
        y = -eta*LH+tc*LP+E0*LD
    elif kind == 2:
        x = 2*E*LE+eta*LH-tc*LP+tc*LR
        y = eta*LH-tc*LP-tc*LR
    elif kind == 3:
        x = 2*E*LE+tc*LR+E0*LD
        y = tc*LR-E0*LD
    else:
        raise ValueError("kind must be 1, 2, or 3")
    return x, y


def cosine_connection(x, y, taus, order=24, local_phase=1.):
    """Return the segment function at u=1 for a vector of durations.

    Coefficients are scaled by the local continuation step to avoid overflow.
    Analytic convergence is local to |z-z0|<|z0|. Finite precision accuracy
    must also be checked with independent propagation and order refinement.
    """
    taus = np.atleast_1d(np.asarray(taus, dtype=float))
    if np.any(taus < 0) or order < 4 or local_phase <= 0:
        raise ValueError("nonnegative durations, order>=4, positive phase required")
    dim = x.shape[0]
    identity = np.eye(dim, dtype=complex)
    scale = float(taus.max(initial=0))*(np.linalg.norm(x, 2)+np.linalg.norm(y, 2))
    panels = max(18, int(np.ceil(scale/local_phase)))
    b2 = taus[:, None, None]*y/(2j*np.pi)
    b1 = taus[:, None, None]*x/(1j*np.pi)
    b0 = b2
    out = np.broadcast_to(identity, (len(taus), dim, dim)).copy()
    zvalues = np.exp(1j*np.linspace(0., np.pi, panels+1))
    for z0, z1 in zip(zvalues[:-1], zvalues[1:]):
        step = z1-z0
        d0 = b2*z0*z0+b1*z0+b0
        d1 = 2*b2*z0+b1
        t = np.broadcast_to(identity, out.shape).copy()
        previous, previous2 = np.zeros_like(t), np.zeros_like(t)
        local = t.copy()
        for k in range(order):
            following = step*(d0@t-2*k*z0*t)
            if k:
                following += step**2*(d1@previous-(k-1)*previous)
            if k > 1:
                following += step**3*(b2@previous2)
            following /= (k+1)*z0*z0
            local += following
            previous2, previous, t = previous, t, following
        out = local@out
    return out, panels


def protocol_function(taus, eta, tc=.3, E0=.3, E=.01, **series_options):
    stages, panels = [], []
    for kind in (1, 2, 3):
        u, n = cosine_connection(
            *segment_coefficients(kind, eta, tc, E0, E), taus, **series_options)
        stages.append(u)
        panels.append(n)
    cycle = stages[2]@stages[1]@stages[0]
    return cycle@cycle, panels


def direct_protocol(tau, eta, tc=.3, E0=.3, E=.01):
    """Independent adaptive real-time propagation, for validation only."""
    out = np.eye(5)
    for kind in (1, 2, 3):
        x, y = segment_coefficients(kind, eta, tc, E0, E)
        def rhs(u, flat):
            return (tau*(x+y*np.cos(np.pi*u))@flat.reshape(5, 5)).ravel()
        solution = solve_ivp(rhs, (0., 1.), out.ravel(), method="DOP853",
                             rtol=3e-12, atol=3e-14)
        if not solution.success:
            raise RuntimeError(solution.message)
        out = solution.y[:, -1].reshape(5, 5)
    return out@out


def readout(rotation):
    """Original target-wavefunction overlap; keep the two readouts distinct."""
    r = np.asarray(rotation).real
    overlap = ((r[..., 0, 0]-r[..., 1, 1])
               -1j*(r[..., 0, 1]+r[..., 1, 0]))/2
    fidelity = np.abs(overlap)**2
    d = r[..., 0, 0]*r[..., 1, 1]-r[..., 0, 1]*r[..., 1, 0]
    legacy_linear = (1-d)/2
    return np.sqrt(fidelity), fidelity, legacy_linear


def ensure_full_model_cuts(data_dir):
    """Reuse validated reference arrays, or regenerate them from this recurrence.

    The reference does not use the scalar low-energy formula or image fitting.
    Fresh generation includes independent real-time checks at three durations.
    """
    import json
    from pathlib import Path
    data_dir = Path(data_dir)
    data_dir.mkdir(parents=True, exist_ok=True)
    path = data_dir / "full_model_cuts.npz"
    report_path = data_dir / "full_model_cuts.json"
    parameters = {"E1": .01, "tc": .3, "E0": .3, "K": "2A", "cycles": 2}
    taus = np.linspace(0., 1200., 601)
    ratios = (.1, 1., 2.3446, 10.)
    if path.exists() and report_path.exists():
        report = json.loads(report_path.read_text())
        assert report["parameters"] == parameters
        with np.load(path) as saved:
            assert np.array_equal(saved["tau"], taus)
            for ratio in ratios:
                assert saved[f"f_ratio_{ratio}"].shape == taus.shape
        return path
    arrays = {"tau": taus}
    report = {"parameters": parameters, "fit_parameters": 0,
              "function": "Matrix cosine function, continued by local Taylor recurrence",
              "curves": [], "independent_checks": []}
    for ratio in ratios:
        rotation, panels = protocol_function(taus, ratio*.01, order=28, local_phase=2.)
        f, squared, legacy = readout(rotation)
        arrays[f"f_ratio_{ratio}"] = f
        arrays[f"F_ratio_{ratio}"] = squared
        orthogonality = float(np.max(abs(np.swapaxes(rotation, -1, -2)@rotation-np.eye(5))))
        imaginary = float(np.max(abs(rotation.imag)))
        assert orthogonality < 1e-8 and imaginary < 1e-8
        report["curves"].append({"ratio": ratio, "series_order": 28, "local_phase": 2.,
            "continuation_panels": panels, "orthogonality_max_abs": orthogonality,
            "imaginary_max_abs": imaginary})
        for tau in (150., 300., 600.):
            reference = direct_protocol(tau, ratio*.01)
            index = int(round(tau/2))
            error = float(np.max(abs(rotation[index]-reference)))
            assert error < 2e-8
            report["independent_checks"].append({"ratio": ratio, "tau": tau,
                                                "max_matrix_abs": error})
        print(f"Full-model reference complete: r={ratio:g}", flush=True)
    np.savez_compressed(path, **arrays)
    report_path.write_text(json.dumps(report, indent=2)+"\n")
    return path


if __name__ == "__main__":
    from pathlib import Path
    print(ensure_full_model_cuts(Path(__file__).resolve().parents[1]/"data"))
