import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh
import matplotlib

# Use non-interactive backend
matplotlib.use('Agg')

# Import your WKB implementation
from wkb2 import V, find_energy, construct_wavefunction, m, hbar, a, V0

def solve_schrodinger_numerically(x_grid, n_eigen):
    """
    Solve the 1D Schrodinger equation numerically using Finite Difference Method.
    H = -h^2/2m * d^2/dx^2 + V(x)
    """
    dx = x_grid[1] - x_grid[0]
    N = len(x_grid)
    
    # Kinetic energy matrix (Finite Difference 2nd derivative)
    # Stencil: [1, -2, 1] / dx^2
    factor = -hbar**2 / (2 * m * dx**2)
    diagonals = [np.ones(N-1), -2*np.ones(N), np.ones(N-1)]
    # offsets: -1 (lower), 0 (main), 1 (upper)
    T = diags(diagonals, [-1, 0, 1], shape=(N, N)) * factor
    
    # Potential energy matrix (Diagonal)
    U = diags([V(x_grid)], [0], shape=(N, N))
    
    # Hamiltonian
    H = T + U
    
    # Solve eigenvalue problem for the smallest eigenvalues (algebraic)
    # 'SA' = Smallest Algebraic
    vals, vecs = eigsh(H, k=n_eigen, which='SA')
    
    return vals, vecs

def main():
    print("=== Verifying WKB Results against Numerical Exact Solution ===\n")
    
    # 1. Setup Grid for Numerical Solver
    # Needs to be dense enough for accuracy
    x_max = a * 2.5
    N_grid = 2000
    x_grid = np.linspace(-x_max, x_max, N_grid)
    dx = x_grid[1] - x_grid[0]
    
    # 2. Get Exact Numerical Solution
    # We calculate first 4 levels (n=0,1,2,3)
    # In double well, these appear as pairs: (0,1) and (2,3)
    n_exact_levels = 4
    exact_energies, exact_vecs = solve_schrodinger_numerically(x_grid, n_exact_levels)
    
    # Normalize numerical wavefunctions: sum |psi|^2 dx = 1
    # eigsh returns normalized vectors such that sum |v|^2 = 1
    # So psi(x) = v / sqrt(dx)
    exact_wavefunctions = exact_vecs / np.sqrt(dx)
    
    # 3. Compare with WKB
    # WKB code calculates "single well" levels.
    # WKB n=0 corresponds to Exact n=0 & n=1 (Ground state doublet)
    # WKB n=1 corresponds to Exact n=2 & n=3 (First excited doublet)
    
    print(f"{'Level (Exact)':<15} | {'Exact E':<10} | {'WKB E (approx)':<15} | {'Diff':<10}")
    print("-" * 60)
    
    plt.figure(figsize=(12, 8))
    # Plot Potential
    plt.plot(x_grid, V(x_grid), 'k-', alpha=0.2, label='Potential V(x)', linewidth=2)
    
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
    
    for n_exact in range(n_exact_levels):
        # Map exact level to WKB single-well level
        n_wkb = n_exact // 2
        
        # Get WKB Energy
        E_wkb = find_energy(n_wkb)
        
        E_exact = exact_energies[n_exact]
        
        print(f"n={n_exact:<12} | {E_exact:.4f}     | {E_wkb:.4f}          | {abs(E_exact - E_wkb):.4f}")
        
        # Plotting
        psi_exact = exact_wavefunctions[:, n_exact]
        
        # Generate WKB wavefunction for comparison
        # Even parity for n_exact 0, 2; Odd for 1, 3
        parity = 'even' if n_exact % 2 == 0 else 'odd'
        psi_wkb = construct_wavefunction(x_grid, E_wkb, parity=parity)
        
        # Align phases (sign is arbitrary in eigensolvers)
        # Find index of max amplitude to check sign
        idx_max = np.argmax(np.abs(psi_exact))
        if np.sign(psi_exact[idx_max]) != np.sign(psi_wkb[idx_max]):
            psi_exact = -psi_exact
            
        # Offset for plotting
        scale = 3.0
        offset = E_exact
        
        label_exact = f'Exact n={n_exact}'
        label_wkb = f'WKB (n_well={n_wkb})' if n_exact % 2 == 0 else None # Avoid duplicate labels
        
        plt.plot(x_grid, psi_exact * scale + offset, '-', color=colors[n_exact], alpha=0.5, linewidth=2, label=label_exact)
        plt.plot(x_grid, psi_wkb * scale + offset, '--', color=colors[n_exact], linewidth=1.5, label=label_wkb)

    plt.title("Verification: Numerical Exact vs WKB Approximation")
    plt.xlabel("Position x")
    plt.ylabel("Energy / Wavefunction Amplitude")
    plt.ylim(0, V0 * a**4 * 0.8) # Focus on the wells
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=4)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    output_file = 'wkb_verification.png'
    plt.savefig(output_file)
    print(f"\nComparison plot saved to '{output_file}'")
    print("Note: WKB energies are for a single well and do not capture the tunneling splitting,")
    print("so they should lie roughly between the exact doublet energies.")

if __name__ == "__main__":
    main()
