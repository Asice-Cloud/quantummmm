#!/usr/bin/env python3
import numpy as np

# grid
N = 400
re = np.linspace(-2, 2, N)
im = np.linspace(-2, 2, N)
Re, Im = np.meshgrid(re, im)
Xgrid = Re + 1j * Im

b_y_list = [0.0, 0.3, 1.0]

for b_y in b_y_list:
    Y = np.exp(2j * b_y)
    X = Xgrid
    denom = X + Y

    # safe denom
    eps = 1e-12
    denom_safe = denom.copy()
    mask_small = np.abs(denom_safe) < eps
    denom_safe[mask_small] = eps

    Z_minus = (X * Y - 1.0) / denom_safe
    Z_plus = (1.0 - X * Y) / denom_safe

    diff_minus = np.abs(np.abs(Z_minus) - 1.0)
    diff_plus = np.abs(np.abs(Z_plus) - 1.0)

    min_minus = diff_minus.min()
    idx_minus = np.unravel_index(np.argmin(diff_minus), diff_minus.shape)
    x_min_minus = Xgrid[idx_minus]
    z_at_min_minus = Z_minus[idx_minus]

    min_plus = diff_plus.min()
    idx_plus = np.unravel_index(np.argmin(diff_plus), diff_plus.shape)
    x_min_plus = Xgrid[idx_plus]
    z_at_min_plus = Z_plus[idx_plus]

    print(f"b_y={b_y}:\n  F_- min | |Z|-1 | = {min_minus:.3e} at X={x_min_minus} (grid index {idx_minus}), Z={z_at_min_minus}")
    print(f"  F_+ min | |Z|-1 | = {min_plus:.3e} at X={x_min_plus} (grid index {idx_plus}), Z={z_at_min_plus}\n")

    # report if denom small region exists
    if mask_small.any():
        coords = np.column_stack(np.where(mask_small))
        print(f"  Note: denom ~0 at {coords.shape[0]} grid points (sample index {coords[0]})")

print('Done')
