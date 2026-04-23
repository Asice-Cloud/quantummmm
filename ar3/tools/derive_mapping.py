from sympy import symbols, I, Rational, simplify, Matrix

# Define symbolic coefficients for Pauli two-site terms
c_xx, c_yy, c_xy, c_yx = symbols('c_xx c_yy c_xy c_yx')

# Majorana labels: gamma1..gamma4
# Use integers 1..4 to denote them in monomials

# Helper: represent linear combination of monomials as dict: monomial(tuple)->coeff
# monomial is tuple of gamma indices in order

def add_term(d, mono, coeff):
    if mono in d:
        d[mono] = simplify(d[mono] + coeff)
    else:
        d[mono] = simplify(coeff)

# Multiply two sums represented as dicts
def mul(sumA, sumB):
    out = {}
    for a_mono, a_coeff in sumA.items():
        for b_mono, b_coeff in sumB.items():
            mono = a_mono + b_mono
            coeff = simplify(a_coeff * b_coeff)
            # canonicalize mono using anticommutation: gamma_i gamma_j = - gamma_j gamma_i, gamma_i^2 = 1
            # We'll reduce monomial by reordering with sign and eliminating squares
            mono_list = list(mono)
            sign = 1
            # bubble sort to reorder indices ascending while tracking sign flips
            # For fermionic operators, swapping two different gammas introduces a -1 sign
            for i in range(len(mono_list)):
                for j in range(i+1, len(mono_list)):
                    if mono_list[i] > mono_list[j]:
                        # swap
                        mono_list[i], mono_list[j] = mono_list[j], mono_list[i]
                        sign *= -1
            # eliminate squares: gamma_i gamma_i = 1 (and complete removal from tuple)
            i = 0
            while i < len(mono_list)-1:
                if mono_list[i] == mono_list[i+1]:
                    # remove both and multiply coeff by 1 (no change)
                    del mono_list[i:i+2]
                    # after deletion, restart scanning
                    i = 0
                    continue
                i += 1
            mono_t = tuple(mono_list)
            final_coeff = simplify(sign * coeff)
            add_term(out, mono_t, final_coeff)
    return out

# Define gamma expressions for c1, c1d, c2, c2d in terms of gamma operators
# c1 = (g1 - i g2)/2, c1d = (g1 + i g2)/2
# c2 = (g3 - i g4)/2, c2d = (g3 + i g4)/2

def gamma_linear_comb(coeff_pairs):
    # coeff_pairs: list of (coeff, index) representing coeff * gamma_index
    d = {}
    for coeff, idx in coeff_pairs:
        add_term(d, (idx,), simplify(coeff))
    return d

half = Rational(1,2)

c1 = gamma_linear_comb([(half,1), ( -I*half, 2)])
c1d = gamma_linear_comb([(half,1), ( I*half, 2)])

c2 = gamma_linear_comb([(half,3), ( -I*half, 4)])
c2d = gamma_linear_comb([(half,3), ( I*half, 4)])

# Identity monomial
ID = {(): 1}

# Map Pauli products (nearest neighbor) to fermionic bilinears per JW results
# σ^x_i σ^x_{i+1} -> c1d*c2d + c1d*c2 + c1*c2d + c1*c2
term_xx = {}
for mono, coeff in mul(mul(c1d, c2d), ID).items(): add_term(term_xx, mono, coeff)
for mono, coeff in mul(mul(c1d, c2), ID).items(): add_term(term_xx, mono, coeff)
for mono, coeff in mul(mul(c1, c2d), ID).items(): add_term(term_xx, mono, coeff)
for mono, coeff in mul(mul(c1, c2), ID).items(): add_term(term_xx, mono, coeff)

# σ^y_i σ^y_{i+1} -> -c1d*c2d + c1d*c2 + c1*c2d - c1*c2
term_yy = {}
for mono, coeff in mul(mul(c1d, c2d), ID).items(): add_term(term_yy, mono, -coeff)
for mono, coeff in mul(mul(c1d, c2), ID).items(): add_term(term_yy, mono, coeff)
for mono, coeff in mul(mul(c1, c2d), ID).items(): add_term(term_yy, mono, coeff)
for mono, coeff in mul(mul(c1, c2), ID).items(): add_term(term_yy, mono, -coeff)

# σ^x_i σ^y_{i+1} -> -i c1d*c2d + i c1d*c2 - i c1*c2d + i c1*c2
term_xy = {}
for mono, coeff in mul(mul(c1d, c2d), ID).items(): add_term(term_xy, mono, -I*coeff)
for mono, coeff in mul(mul(c1d, c2), ID).items(): add_term(term_xy, mono, I*coeff)
for mono, coeff in mul(mul(c1, c2d), ID).items(): add_term(term_xy, mono, -I*coeff)
for mono, coeff in mul(mul(c1, c2), ID).items(): add_term(term_xy, mono, I*coeff)

# σ^y_i σ^x_{i+1} -> -i c1d*c2d - i c1d*c2 + i c1*c2d + i c1*c2
term_yx = {}
for mono, coeff in mul(mul(c1d, c2d), ID).items(): add_term(term_yx, mono, -I*coeff)
for mono, coeff in mul(mul(c1d, c2), ID).items(): add_term(term_yx, mono, -I*coeff)
for mono, coeff in mul(mul(c1, c2d), ID).items(): add_term(term_yx, mono, I*coeff)
for mono, coeff in mul(mul(c1, c2), ID).items(): add_term(term_yx, mono, I*coeff)

# Now build total H_P (only the quadratic-contributing Pauli combinations)
HP = {}
for mono, coeff in term_xx.items(): add_term(HP, mono, c_xx*coeff)
for mono, coeff in term_yy.items(): add_term(HP, mono, c_yy*coeff)
for mono, coeff in term_xy.items(): add_term(HP, mono, c_xy*coeff)
for mono, coeff in term_yx.items(): add_term(HP, mono, c_yx*coeff)

# We want to express HP in the basis of bilinears i*gamma_a*gamma_b (a<b)
# Extract coefficients for monomials of length 2
bilin_keys = [(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)]
coeffs = {k: simplify(HP.get(k,0)) for k in bilin_keys}

# For convenience, factor out i and present A_ab where HP = i * sum_{a<b} A_ab * gamma_a gamma_b
A = {}
for k,v in coeffs.items():
    # coefficients currently correspond to gamma_a gamma_b; HP = sum v * gamma_a gamma_b
    # We want A_ab such that HP = i * sum A_ab * gamma_a gamma_b -> A_ab = -I * v
    A[k] = simplify(-I * v)

# Print results
print("Mapping HP -> i * sum_{a<b} A_ab gamma_a gamma_b  (A_ab in terms of c_xx,c_yy,c_xy,c_yx):\n")
for k in bilin_keys:
    print(f"A_{k} = {A[k]}")

# Also print HP expanded
print('\nExpanded HP (monomial -> coeff):')
for mono, coeff in sorted(HP.items(), key=lambda x: (len(x[0]), x[0])):
    print(f"{mono}: {simplify(coeff)}")

# Now express HP in the fermionic operator basis
# candidate fermion monomials: n1=c1d*c1, n2=c2d*c2, hop1=c1d*c2, hop2=c1*c2d,
# pair1=c1d*c2d, pair2=c1*c2
ferm_monoms = [ ('n1', (c1d, c1)), ('n2', (c2d, c2)), ('hop1', (c1d, c2)), ('hop2', (c1, c2d)), ('pair1', (c1d, c2d)), ('pair2', (c1, c2)) ]

# Build basis (bilinear gamma monomials) from bilin_keys
basis = bilin_keys

# Build matrix columns: each fermion monomial expands into basis coefficients
M_cols = []
for name, (L, R) in ferm_monoms:
    exp = mul(L, R)  # dict gamma_mono -> coeff
    col = [ simplify(exp.get(k, 0)) for k in basis ]
    M_cols.append(col)

# HP vector in same basis
HP_vec = [ simplify(HP.get(k, 0)) for k in basis ]

Mat = Matrix(M_cols).T
vec = Matrix(HP_vec)

try:
    sol = Mat.LUsolve(vec)
except Exception:
    sol = Mat.gauss_jordan_solve(vec) if hasattr(Mat, 'gauss_jordan_solve') else None

print('\nRepresentation of HP in fermionic bilinear basis (coefficients for [n1,n2,hop1,hop2,pair1,pair2]):')
if sol is None:
    print('Could not solve linear system for fermionic coefficients')
else:
    for (name, _), val in zip(ferm_monoms, sol):
        print(f"{name}: {simplify(val)}")

# Interpret common BdG parameters (conventions may vary):
# hopping term ~ coef_hop + coef_hop2 (hermitian), pairing ~ coef_pair + coef_pair2, onsite mu ~ coef_n
if sol is not None:
    coeffs = { name: simplify(val) for (name,_), val in zip(ferm_monoms, sol) }
    print('\nInterpreting BdG-like parameters (raw coefficients):')
    print('mu1 (n1):', coeffs['n1'])
    print('mu2 (n2):', coeffs['n2'])
    print('hop c1d*c2 (hop1):', coeffs['hop1'])
    print('hop c1*c2d (hop2):', coeffs['hop2'])
    print('pair c1d*c2d (pair1):', coeffs['pair1'])
    print('pair c1*c2 (pair2):', coeffs['pair2'])

    print('\nNote: mapping to standard parameters t, \u03b4 may require factoring conventions:')
    print('For example, a kinetic term -t(c1d*c2 + c2d*c1) has coefficient -t for c1d*c2 here.')
