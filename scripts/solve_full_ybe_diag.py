from sympy import symbols, I, Matrix, simplify, factor, expand
from sympy import kronecker_product as kron

# Pauli matrices
sigma0 = Matrix([[1,0],[0,1]])
sigma1 = Matrix([[0,1],[1,0]])
sigma2 = Matrix([[0,-I],[I,0]])
sigma3 = Matrix([[1,0],[0,-1]])
paulis = [sigma0, sigma1, sigma2, sigma3]

# symbols
Jx, Jy, Jz, c00 = symbols('Jx Jy Jz c00')

# build K = c00*I\otimes I + Jx X\otimes X + Jy Y\otimes Y + Jz Z\otimes Z
K = c00 * kron(sigma0, sigma0) + Jx * kron(sigma1, sigma1) + Jy * kron(sigma2, sigma2) + Jz * kron(sigma3, sigma3)

# matrix exponential R = exp(i K)
R = (I * K).exp()

# build r12,r23 and r13 via permutation P23
I2 = sigma0
r12 = kron(R, I2)
r23 = kron(I2, R)

# permutation P23 on three qubits swaps factors 2 and 3: acts on 8-dim basis
P23 = Matrix([[0] * 8 for _ in range(8)])
for a in range(2):
    for b in range(2):
        for c in range(2):
            src = a * 4 + b * 2 + c
            dst = a * 4 + c * 2 + b
            P23[dst, src] = 1

r13 = P23 * r12 * P23

expr = simplify(r12 * r13 * r23 - r23 * r13 * r12)

# collect unique nonzero simplified equations
eqs = []
for i in range(8):
    for j in range(8):
        e = simplify(expand(expr[i, j]))
        if e != 0:
            eqs.append(e)

uniq = []
for e in eqs:
    if e not in uniq:
        uniq.append(e)

print(f"Nonzero entries: {len(uniq)}")

# factor and print all unique equations (they will be trigonometric expressions)
with open('full_ybe_diag_out.txt', 'w') as f:
    f.write(f"Nonzero entries: {len(uniq)}\n\n")
    for idx, e in enumerate(uniq):
        fe = factor(e)
        print(f"entry{idx}: {fe}\n")
        f.write(f"entry{idx}: {str(fe)}\n\n")

print('Wrote full_ybe_diag_out.txt')
