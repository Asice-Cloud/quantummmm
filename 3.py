from sympy import symbols, Matrix, eye, simplify, expand
from sympy import I as symI
from sympy import KroneckerProduct, solve

a,b,c,d = symbols('a b c d', real=True)
sx = Matrix([[0,1],[1,0]])
sy = Matrix([[0,-symI],[symI,0]])
sz = Matrix([[1,0],[0,-1]])
I2 = eye(2)
kron = lambda A,B: KroneckerProduct(A,B).as_mutable()
R = a*kron(I2,I2) + b*kron(sx,sx) + c*kron(sy,sy) + d*kron(sz,sz)
kron3 = lambda A,B,C: KroneckerProduct(KroneckerProduct(A,B),C).as_mutable()
R12 = KroneckerProduct(R, I2).as_mutable()
R23 = KroneckerProduct(I2, R).as_mutable()
R13 = a*kron3(I2,I2,I2) + b*kron3(sx,I2,sx) + c*kron3(sy,I2,sy) + d*kron3(sz,I2,sz)
D = simplify(R12*R13*R23 - R23*R13*R12)
eqs = []
for i in range(8):
    for j in range(8):
        entry = simplify(D[i,j])
        if entry!=0:
            eqs.append(expand(entry))
# remove duplicates
uniq = []
for e in eqs:
    if e not in uniq:
        uniq.append(e)
eqs = uniq
print('num_equations=', len(eqs))
# simplify eqs
from sympy import factor
for i,e in enumerate(eqs,1):
    print(i, factor(e))
# Solve by computing common factorizations manually: derive simplified constraints
# We'll solve using solve for cases by scanning logical cases
sol = []
# Try solve generally (may be heavy)
try:
    sol = solve(eqs, [a,b,c,d], dict=True)
    print('solve general ->', sol)
except Exception as ex:
    print('general solve failed:', ex)
# Try solving by cases from algebraic deductions
print('\nAnalytic deductions:')
# From factoring representative equations
from sympy import cancel
E1 = factor(eqs[0])
print('E1=',E1)