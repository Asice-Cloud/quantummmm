import numpy as np
import matplotlib.pyplot as plt
from wkb2 import V, find_energy, construct_wavefunction, get_turning_points, m, hbar, a, V0
from verify_wkb import solve_schrodinger_numerically
import matplotlib

# Use non-interactive backend
matplotlib.use('Agg')

def main():
    print("=== Comparing WKB and Finite Difference Solutions at Turning Points ===\n")

    # 1. Setup Grid
    # Use high resolution to see details at turning points
    x_max = a * 2.5
    N_grid = 5000 
    x_grid = np.linspace(-x_max, x_max, N_grid)
    dx = x_grid[1] - x_grid[0]

    # 2. Exact Solution (Finite Difference Method)
    # We calculate the ground state (n=0)
    # Note: In double well, n=0 and n=1 are nearly degenerate.
    vals, vecs = solve_schrodinger_numerically(x_grid, 2)
    E_exact = vals[0]
    psi_exact = vecs[:, 0] / np.sqrt(dx) # Normalize

    # 3. WKB Solution
    # WKB ground state corresponds to n=0 in single well approximation
    E_wkb = find_energy(0)
    psi_wkb = construct_wavefunction(x_grid, E_wkb, parity='even')

    # Align phases (ensure both are positive at the peak)
    # Find index of max amplitude in the left well
    left_well_mask = x_grid < 0
    idx_max = np.argmax(np.abs(psi_exact * left_well_mask))
    
    if np.sign(psi_exact[idx_max]) != np.sign(psi_wkb[idx_max]):
        psi_exact = -psi_exact

    # 4. Identify Turning Points
    tp = get_turning_points(E_wkb)
    # tp = (-x_out, -x_in, x_in, x_out)
    x_out_left = tp[0]
    x_in_left = tp[1]

    print(f"Ground State Comparison:")
    print(f"Exact Energy (FDM): {E_exact:.6f}")
    print(f"WKB Energy:         {E_wkb:.6f}")
    print(f"Turning Points (Left Well): Outer={x_out_left:.4f}, Inner={x_in_left:.4f}")

    # 5. Plotting
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))

    # Plot 1: Full Wavefunction
    ax = axes[0]
    # Scale wavefunction for visibility against potential
    scale = 2.0
    ax.plot(x_grid, V(x_grid), 'k-', alpha=0.2, label='Potential V(x)')
    ax.plot(x_grid, psi_exact * scale + E_exact, 'b-', label='Exact (FDM)', linewidth=2, alpha=0.7)
    ax.plot(x_grid, psi_wkb * scale + E_wkb, 'r--', label='WKB (Airy)', linewidth=1.5)
    ax.set_title("Full Wavefunction (Ground State)")
    ax.set_ylim(0, V0 * a**4 * 0.6)
    ax.set_xlabel("x")
    ax.set_ylabel("Energy")
    ax.legend()

    # Plot 2: Zoom at Outer Turning Point (Left side of left well)
    ax = axes[1]
    zoom_width = 0.8
    mask_out = (x_grid > x_out_left - zoom_width) & (x_grid < x_out_left + zoom_width)
    
    ax.plot(x_grid[mask_out], psi_exact[mask_out], 'b-', label='Exact', linewidth=2.5, alpha=0.6)
    ax.plot(x_grid[mask_out], psi_wkb[mask_out], 'r--', label='WKB', linewidth=2)
    ax.axvline(x_out_left, color='k', linestyle=':', label='Turning Point')
    ax.set_title(f"Zoom: Outer Turning Point\n(x ≈ {x_out_left:.3f})")
    ax.set_xlabel("x")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 3: Zoom at Inner Turning Point (Right side of left well)
    ax = axes[2]
    mask_in = (x_grid > x_in_left - zoom_width) & (x_grid < x_in_left + zoom_width)
    
    ax.plot(x_grid[mask_in], psi_exact[mask_in], 'b-', label='Exact', linewidth=2.5, alpha=0.6)
    ax.plot(x_grid[mask_in], psi_wkb[mask_in], 'r--', label='WKB', linewidth=2)
    ax.axvline(x_in_left, color='k', linestyle=':', label='Turning Point')
    ax.set_title(f"Zoom: Inner Turning Point\n(x ≈ {x_in_left:.3f})")
    ax.set_xlabel("x")
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    output_file = 'wkb_turning_point_comparison.png'
    plt.savefig(output_file)
    print(f"\nDetailed comparison plot saved to '{output_file}'")

if __name__ == "__main__":
    main()
