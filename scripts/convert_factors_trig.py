from sympy import symbols, I, Matrix, simplify, expand, factor, exp, cos, sin
from sympy import kronecker_product as kron

# Pauli
sigma0 = Matrix([[1,0],[0,1]])
sigma1 = Matrix([[0,1],[1,0]])
sigma2 = Matrix([[0,-I],[I,0]])
sigma3 = Matrix([[1,0],[0,-1]])

Jx, Jy, Jz, c00 = symbols('Jx Jy Jz c00')

# build K and R
K = c00 * kron(sigma0, sigma0) + Jx * kron(sigma1, sigma1) + Jy * kron(sigma2, sigma2) + Jz * kron(sigma3, sigma3)
R = (I * K).exp()
I2 = sigma0
r12 = kron(R, I2)
r23 = kron(I2, R)
# permutation P23
P23 = Matrix([[0]*8 for _ in range(8)])
for a in range(2):
    for b in range(2):
        for c in range(2):
            src = a*4 + b*2 + c
            dst = a*4 + c*2 + b
            P23[dst, src] = 1
r13 = P23 * r12 * P23
expr = simplify(r12*r13*r23 - r23*r13*r12)

factors = []
for i in range(8):
    for j in range(8):
        e = expr[i, j]
        if e != 0:
            fe = factor(simplify(e))
            factors.append(fe)

# deduplicate by structural equality (small list)
uniq_exprs = []
for fe in factors:
    if not any(fe.equals(u) for u in uniq_exprs):
        uniq_exprs.append(fe)

print(f"Found {len(uniq_exprs)} distinct factor expressions (structural).\n")

# mapping for exp(i*theta) -> cos+I*sin using xreplace to avoid heavy expansion
map_exp = {exp(I*Jx): cos(Jx) + I*sin(Jx), exp(I*Jy): cos(Jy) + I*sin(Jy), exp(I*Jz): cos(Jz) + I*sin(Jz)}

out_lines = []
for idx, fe in enumerate(uniq_exprs):
    # only perform targeted replacement; avoid full expand/simplify unless small
    fe_trig = fe.xreplace(map_exp)
    # try lightweight simplify of real/imag
    try:
        realp, imagp = fe_trig.as_real_imag()
    except Exception:
        realp = fe_trig
        imagp = 0

    print(f"--- factor {idx} ---")
    print("orig:", fe)
    print("trig (approx):", fe_trig)
    print("real part:", realp)
    print("imag part:", imagp)
    print()

    out_lines.append(f"--- factor {idx} ---\n")
    out_lines.append("orig: " + str(fe) + "\n")
    out_lines.append("trig: " + str(fe_trig) + "\n")
    out_lines.append("real: " + str(realp) + "\n")
    out_lines.append("imag: " + str(imagp) + "\n\n")

with open('scripts/converted_factors.txt', 'w') as f:
    f.writelines(out_lines)

print('Wrote scripts/converted_factors.txt')
