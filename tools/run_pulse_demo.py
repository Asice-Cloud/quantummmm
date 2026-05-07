#!/usr/bin/env python3
"""Run a short demo using pulse_shapes with the Bloch rotation simulator.

Generates a few example plots saved under `results/`.
"""
import os
import sys

# ensure repo root importable
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from tools import bloch_rotation_sim as bloch
import tools.pulse_shapes as pulses


def main():
    os.makedirs('results', exist_ok=True)
    # choose a modest T list for a short demo
    T_list = [50, 100, 200]

    # Case A: MZM (no splitting)
    d_none = lambda t: 0.0
    bloch.scan_T_and_plot(T_list, d_none, label='Demo MZM (d=0)', outprefix='demo_bloch_MZM')

    # Case B: ABS constant hybridization
    d_const = pulses.d_const(0.3)
    bloch.scan_T_and_plot(T_list, d_const, label='Demo ABS const', outprefix='demo_bloch_ABS_const')

    # Case C: ABS pulsed in the middle step (assume T_step=200 for the pulse placement)
    T_step_demo = 200.0
    d_mid = pulses.d_pulse_middle(T_step_demo, amplitude=0.3)
    bloch.scan_T_and_plot(T_list, d_mid, label='Demo ABS middle pulse', outprefix='demo_bloch_ABS_pulse')

    print('Pulse demo finished. See results/demo_bloch_*.png')


if __name__ == '__main__':
    main()
