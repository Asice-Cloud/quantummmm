from sympy import symbols, Matrix, eye, simplify, expand
from sympy import I as symI
from sympy import KroneckerProduct

a,b,c,d = symbols('a b c d')
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
            entry_r = expand(entry.as_real_imag()[0])
            entry_i = expand(entry.as_real_imag()[1])
            if entry_r!=0 and entry_r not in eqs:
                eqs.append(entry_r)
            if entry_i!=0 and entry_i not in eqs:
                eqs.append(entry_i)
with open('ybe_eqs.txt','w') as f:
    f.write('num_equations=%d\n'%len(eqs))
    for k,eq in enumerate(eqs,1):
        f.write('Eq%d: %s\n'%(k,str(eq)))
print('wrote ybe_eqs.txt')