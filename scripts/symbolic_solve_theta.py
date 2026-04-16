#!/usr/bin/env python3
"""Symbolically analyze the representative symmetric polynomial

P = e^{i(θx+θy)} - e^{i(θx+θz)} + e^{i(θy+θz)} - 1 = 0
We factor to obtain a simpler relation and enumerate analytic solution branches.
"""
from sympy import symbols, exp, I, sin, cos, simplify, factor, Eq
from sympy import pi

# define real theta symbols
θx, θy, θz = symbols('θx θy θz', real=True)

# define alpha,beta,gamma for clarity
α = θx + θy
β = θx + θz
γ = θy + θz

# original exponential polynomial
U = exp(I*α) - exp(I*β) + exp(I*γ) - 1

# factor U analytically: derive bracket form
# We compute the inner bracket B such that U = 2*I*exp(I*γ/2)*B
B = simplify(U / (2*I*exp(I*γ/2)))

with open('scripts/symbolic_solve_theta_out.txt', 'w') as f:
    f.write('Original U = e^{i(θx+θy)} - e^{i(θx+θz)} + e^{i(θy+θz)} - 1\n')
    f.write('Factorization: U = 2 i e^{i γ/2} * ( e^{i θx} sin((θy-θz)/2) + sin((θy+θz)/2) )\n\n')
    f.write('Derived B (should equal e^{i θx} sin((θy-θz)/2) + sin((θy+θz)/2)):\n')
    f.write(str(simplify(B)) + '\n\n')

    f.write('Set B = 0 gives the central relation:\n')
    f.write('  e^{i θx} * sin((θy-θz)/2) + sin((θy+θz)/2) = 0\n\n')

    f.write('Branch analysis:\n')

    # Branch 1: sin((θy-θz)/2) = 0
    f.write('Branch 1: sin((θy-θz)/2) = 0 => θy - θz = m π.\n')
    f.write('  Then equation reduces to sin((θy+θz)/2)=0 => θy + θz = n π.\n')
    f.write('  Solve: θy = θz and 2 θy = n π => θy = n π/2, θz = θy.\n')
    f.write('  In J-variables (θ = 2 J): Jy,Jz ∈ {n π/4} with Jy ≡ Jz and 2 Jy ∈ {0, π/2 (mod π)}.\n\n')

    # Branch 2: sin((θy-θz)/2) != 0, then
    f.write('Branch 2: sin((θy-θz)/2) != 0 => e^{i θx} = - sin((θy+θz)/2) / sin((θy-θz)/2).\n')
    f.write('  For RHS to have unit modulus (|e^{i θx}|=1) we need |sin(a)| = |sin(b)| where\n')
    f.write('    a=(θy+θz)/2, b=(θy-θz)/2.\n')
    f.write('  Thus sin(a)=± sin(b). Two subcases:\n')
    f.write('   (i) sin(a)=sin(b) -> a ≡ b (mod 2π) or a ≡ π - b (mod 2π).\n')
    f.write('       a ≡ b -> θz ≡ 0 (mod 2π). Leads to simpler reductions (θz=0).\n')
    f.write('       a ≡ π - b -> gives θy or θz equal to π (mod 2π) etc.\n')
    f.write('   (ii) sin(a) = - sin(b) -> a ≡ -b (mod 2π) or a ≡ π + b (mod 2π).\n')
    f.write('       These again reduce to linear congruences on θy,θz (e.g. θy or θz multiples of π).\n\n')

    f.write('Concrete consequences (collected):\n')
    f.write('- Many solutions occur at θx,θy,θz ∈ multiples of π or π/2 -> equivalently J ∈ multiples of π/2 or π/4.\n')
    f.write('- Special subfamilies: (θy=θz=0) gives arbitrary θx (since B reduces to sin 0 + 0 = 0),\n')
    f.write('  which corresponds to Jy=Jz=0 and any Jx — but recall original factorization included global prefactors so check numerically.\n')
    f.write('\nVerification note:\n')
    f.write('  The symbolic reduction shows the problem reduces to the single complex equation B=0.\n')
    f.write('  Practical solving proceeds by enumerating the elementary trigonometric branches above and\n')
    f.write('  validating each branch against the original exponential polynomial to discard spurious roots.\n')

print('Wrote scripts/symbolic_solve_theta_out.txt')
