#!/usr/bin/env python3
"""Construct the ideal braid unitary for a single fermionic mode
built from two Majoranas: B = exp((pi/4) gamma1 gamma2).

In the fermionic occupation basis {|0>,|1>} this equals
U = diag(e^{-i pi/4}, e^{+i pi/4}).
"""
import numpy as np


def ideal_braid_unitary():
    phi0 = np.exp(-1j * np.pi / 4.0)
    phi1 = np.exp(+1j * np.pi / 4.0)
    U = np.diag([phi0, phi1])
    return U


if __name__ == '__main__':
    U = ideal_braid_unitary()
    print('Ideal braid unitary (2x2) in occupation basis:')
    print(U)
