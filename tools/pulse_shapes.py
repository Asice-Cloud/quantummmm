#!/usr/bin/env python3
"""Pulse shapes and modulation helpers (Chen et al. style).

Provide callable pulse/modulation generators for use in Bloch/BdG simulations.
"""
import numpy as np


def d_pulse_middle(T_step, amplitude=0.3):
    """Return a pulse function active during the middle (Step2) interval.

    Usage: d = d_pulse_middle(T_step, amplitude); value = d(t)
    where t is in same units as T_step used by the tetron/bloch simulators.
    """
    def f(t):
        if (t >= T_step) and (t < 2.0 * T_step):
            return float(amplitude)
        return 0.0

    return f


def d_const(amplitude=0.3):
    """Constant hybridization function d(t)=amplitude."""
    return lambda t: float(amplitude)


def d_gaussian(t0, sigma, amplitude):
    """Gaussian pulse centered at t0 with width sigma."""
    def f(t):
        return float(amplitude) * np.exp(-0.5 * ((t - t0) / float(sigma)) ** 2)

    return f


def VD_modulation(T, Vx0, Vx1, base_VD):
    """Return a function giving time-dependent VD(t) used in reproduce_figs.

    VD_t = base_VD * (1 + Vx0 + Vx1 * cos(pi * t / T))
    """
    def f(t):
        mod = float(Vx0) + float(Vx1) * np.cos(np.pi * t / float(T))
        return float(base_VD) * (1.0 + mod)

    return f


def smooth_step(s, width=0.08):
    """Smooth step in [0,1] (useful to smooth gate functions). s in [0,1]."""
    x = (s - 0.5) / float(width)
    return 0.5 * (1.0 + np.tanh(x))


if __name__ == '__main__':
    # quick textual demo
    dmid = d_pulse_middle(200.0, amplitude=0.25)
    times = [0.0, 100.0, 200.0, 300.0, 400.0]
    print('d_pulse_middle samples:', [dmid(t) for t in times])
