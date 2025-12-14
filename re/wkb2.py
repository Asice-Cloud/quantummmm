import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import airy

# Constants
m = 1.0
hbar = 1.0
V0 = 1.0
a = 2.0

def V(x):
    """Double well potential"""
    return V0 * (x**2 - a**2)**2

def get_turning_points(E):
    """Return the four turning points for E < barrier height"""
    term = np.sqrt(E/V0)
    x_out = np.sqrt(a**2 + term)
    x_in = np.sqrt(a**2 - term)
    return -x_out, -x_in, x_in, x_out

def momentum(x, E):
    """Classical momentum k(x) or decay constant kappa(x)"""
    val = 2 * m * (E - V(x))
    return np.sqrt(np.abs(val))

def action_integral(x1, x2, E):
    """Compute action integral between x1 and x2"""
    return quad(lambda x: momentum(x, E), x1, x2, limit=200)[0]

def quantization_condition(E, n):
    """Bohr-Sommerfeld quantization condition for a single well"""
    tp = get_turning_points(E)
    # Left well: tp[0] to tp[1]
    integral = action_integral(tp[0], tp[1], E)
    return integral - (n + 0.5) * np.pi * hbar

def find_energy(n):
    """Find the nth energy level using single well approximation"""
    barrier_height = V0 * a**4
    # Search for root
    try:
        return brentq(quantization_condition, 1e-6, barrier_height * 0.99, args=(n,))
    except ValueError:
        print(f"Could not find energy for n={n}")
        return None

def airy_uniform_approximation(x, xt, E, direction):
    """
    Uniform Airy approximation near a turning point xt.
    direction: +1 if allowed region is to the right (x > xt)
               -1 if allowed region is to the left (x < xt)
    """
    k = momentum(x, E)
    if k < 1e-12: k = 1e-12
    
    # Calculate action S from turning point
    S = abs(action_integral(xt, x, E))
    
    # Zeta variable
    # zeta = (1.5 * S / hbar)^(2/3)
    zeta_mag = (1.5 * S / hbar)**(2/3)
    
    # Determine sign of argument for Airy function
    # If in allowed region, arg is negative.
    # If in forbidden region, arg is positive.
    is_allowed = (E > V(x))
    
    if is_allowed:
        arg = -zeta_mag
    else:
        arg = zeta_mag
        
    # Amplitude factor: (zeta / k^2)^(1/4) -> (S / k^3)^(1/6) * const
    # Using the form: sqrt(pi) * (zeta/k^2)^(1/4) * Ai(arg)
    # But we need to match the WKB amplitude 1/sqrt(k).
    # WKB: 1/sqrt(k) * cos(...)
    # Airy asymptotic: 1/sqrt(pi) * (-arg)^(-1/4) * sin(...)
    # We want to match these.
    
    # Let's use the robust formula:
    # amp = sqrt( 2 * pi * zeta_mag / k ) ? No.
    
    # Let's stick to the one that works in wkb.py which is derived from standard texts
    # amp = sqrt( zeta_mag / k ) * sqrt(pi) ?
    # wkb.py used: amp = np.sqrt(abs(zeta) / k_abs)
    # Let's use that.
    
    amp = np.sqrt(zeta_mag / k)
    
    ai_val, _, _, _ = airy(arg)
    
    # We need to handle the sign/phase correctly.
    # For a simple turning point, the solution is proportional to Ai.
    # The normalization will be handled globally later.
    
    return amp * ai_val

def construct_wavefunction(x_grid, E, parity='even'):
    """
    Construct the wavefunction for the double well.
    We construct the solution for the left well and then symmetrize.
    """
    tp = get_turning_points(E)
    x1, x2, x3, x4 = tp
    
    psi = np.zeros_like(x_grid)
    
    # Midpoint of left well
    x_mid_left = (x1 + x2) / 2
    
    for i, x in enumerate(x_grid):
        # We only compute for x < 0 (left side) and mirror
        if x > 0:
            continue
            
        # Determine which turning point to use for the approximation
        if x < x_mid_left:
            # Closer to x1 (outer left)
            # Allowed region is to the right of x1
            val = airy_uniform_approximation(x, x1, E, direction=1)
            
            # We need to ensure the sign is correct to match the other side
            # But Airy function handles the oscillation.
            # We might need a phase factor if we were matching WKB, but uniform approx handles it.
            psi[i] = val
        else:
            # Closer to x2 (inner left)
            # Allowed region is to the left of x2
            val = airy_uniform_approximation(x, x2, E, direction=-1)
            
            # There might be a sign mismatch between the two patches at x_mid_left.
            # We need to check continuity at x_mid_left.
            # However, for this simple demo, we can just compute.
            # A more robust way is to calculate the sign factor to match at x_mid_left.
            
            psi[i] = val

    # Fix continuity at x_mid_left
    # Calculate values just left and right of x_mid_left to find the ratio
    # This is a bit hacky but ensures continuity
    idx_mid = np.searchsorted(x_grid, x_mid_left)
    if 0 < idx_mid < len(x_grid):
        # We need to evaluate the functions exactly at x_mid_left to find the scaling factor
        val1 = airy_uniform_approximation(x_mid_left, x1, E, direction=1)
        val2 = airy_uniform_approximation(x_mid_left, x2, E, direction=-1)
        
        if abs(val2) > 1e-10:
            scale_factor = val1 / val2
            # Apply scaling to the second part (x > x_mid_left)
            # But wait, we are iterating.
            # Let's do it in a second pass or just apply it now.
            # Since we fill the array, we can just scale the part we just filled?
            # No, we need to scale the part connected to x2.
            pass

    # Let's redo the loop with scaling
    # Compute scale factor first
    val1_mid = airy_uniform_approximation(x_mid_left, x1, E, direction=1)
    val2_mid = airy_uniform_approximation(x_mid_left, x2, E, direction=-1)
    scale = val1_mid / val2_mid if abs(val2_mid) > 1e-10 else 1.0
    
    for i, x in enumerate(x_grid):
        if x > 0: break
        
        if x < x_mid_left:
            psi[i] = airy_uniform_approximation(x, x1, E, direction=1)
        else:
            psi[i] = airy_uniform_approximation(x, x2, E, direction=-1) * scale

    # Mirror to right side
    for i, x in enumerate(x_grid):
        if x > 0:
            # Find corresponding index for -x
            # Assuming symmetric grid
            idx_mirror = len(x_grid) - 1 - i
            if 0 <= idx_mirror < len(x_grid):
                if parity == 'even':
                    psi[i] = psi[idx_mirror]
                else:
                    psi[i] = -psi[idx_mirror]
                    
    # Normalize
    norm = np.trapezoid(psi**2, x_grid)
    psi = psi / np.sqrt(norm)
    
    return psi

def main():
    # Find first few energy levels
    n_levels = 2
    energies = []
    for n in range(n_levels):
        E = find_energy(n)
        if E is not None:
            energies.append(E)
            print(f"Level {n}: E = {E:.4f}")
            
    # Plot
    x_grid = np.linspace(-a*1.5, a*1.5, 1000)
    plt.figure(figsize=(10, 6))
    
    # Plot potential
    plt.plot(x_grid, V(x_grid), 'k-', label='V(x)', alpha=0.5)
    plt.ylim(-1, V0*a**4 * 1.2)
    
    # Plot wavefunctions
    scale = 5.0 # Scale for visibility
    
    for n, E in enumerate(energies):
        # Plot even parity
        psi_even = construct_wavefunction(x_grid, E, parity='even')
        plt.plot(x_grid, scale * psi_even + E, label=f'n={n} (Even)')
        
        # Plot odd parity (degenerate in this approximation)
        # psi_odd = construct_wavefunction(x_grid, E, parity='odd')
        # plt.plot(x_grid, scale * psi_odd + E, '--', label=f'n={n} (Odd)')
        
        # Draw energy level
        plt.axhline(E, color='gray', linestyle=':', alpha=0.5)

    plt.title("WKB Approximation for Double Well Potential with Airy Function Patching")
    plt.xlabel("x")
    plt.ylabel("Energy / Wavefunction")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig('wkb_double_well.png')
    print("Plot saved to wkb_double_well.png")

if __name__ == "__main__":
    main()
