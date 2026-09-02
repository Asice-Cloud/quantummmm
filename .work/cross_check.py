import numpy as np
# quaternion -> complex 2x2 : q = q0 + q1 i + q2 j + q3 k
def q2c(q0,q1,q2,q3):
    return np.array([[q0+1j*q1, q2+1j*q3],[-q2+1j*q3, q0-1j*q1]],complex)
I=np.eye(2,dtype=complex); sx=np.array([[0,1],[1,0]],complex)
sy=np.array([[0,-1j],[1j,0]],complex); sz=np.diag([1.,-1.]).astype(complex)
# model.md Gamma's (2x2 quaternion), convert to 4x4 complex (block: top row of 2x2 quat = two complex rows)
def H_from_quat(M):  # M: 2x2 array of quaternions (each (q0,q1,q2,q3)), complexified as 2x2 blocks
    out=np.zeros((4,4),complex)
    # quaternion entry (row r,col c) -> 2x2 complex block at rows 2r:2r+2, cols 2c:2c+2
    for r in range(2):
        for c in range(2):
            q0,q1,q2,q3=M[r][c]
            out[2*r:2*r+2,2*c:2*c+2]=q2c(q0,q1,q2,q3)
    return out
G1=H_from_quat([[ (0,0,0,0),(1,0,0,0) ],[ (1,0,0,0),(0,0,0,0)]])
G2=H_from_quat([[ (0,0,0,0),(0,0,-1,0) ],[ (0,0,1,0),(0,0,0,0)]])  # 0,-j / j,0
G3=H_from_quat([[ (0,0,0,0),(0,0,0,-1) ],[ (0,0,0,1),(0,0,0,0)]])  # 0,-k / k,0
G4=H_from_quat([[ (1,0,0,0),(0,0,0,0) ],[ (0,0,0,0),(-1,0,0,0)]])  # diag(1,-1)
# Gamma5 = Gamma1 Gamma2 Gamma3 Gamma4 (Clifford chirality) -- model has 5 Gamma's but lists 4 + we derive 5th
# Actually model lists Gamma_1..Gamma_5 as: [0,1],[1,0]; [0,-i],[i,0]; [0,-j],[j,0]; [0,-k],[k,0]; [1,0],[0,-1]
G5=H_from_quat([[ (1,0,0,0),(0,0,0,0) ],[ (0,0,0,0),(-1,0,0,0)]])
Gammas=[G1,G2,G3,G4,G5]
# check Clifford
ok=True
for i in range(5):
    for j in range(5):
        err=np.max(np.abs(Gammas[i]@Gammas[j]+Gammas[j]@Gammas[i]-2*(i==j)*np.eye(4)))
        if err>1e-9: ok=False
print("model Gamma Clifford ok:",ok)
# Now build my JW gammas (gamma1..gamma5 <-> g0..g4). Verify JW vs model relate by unitary
g1=np.kron(sx,I); g2=np.kron(-sy,I); ga=np.kron(sz,sx); gb=np.kron(sz,-sy); g3=ga@gb@g1@g2
JW=[g1,g2,g3,ga,gb]   # gamma1,gamma2,gamma3,gammaA,gammaB
MOD=[G1,G2,G3,G4,G5]  # gamma1,gamma2,gamma3,gammaA(=G4),gammaB(=G5)
# compare traces of products Xij=i g_i g_j - two inequivalent 4-dim complex Cl5 reps differ.
def trace_profile(gs):
    prof=[]
    for i in range(5):
        for j in range(i+1,5):
            prof.append(np.trace(1j*gs[i]@gs[j]))
    return np.array(prof)
print("JW Xij trace :",np.round(trace_profile(JW),4))
print("MOD Xij trace:",np.round(trace_profile(MOD),4))
# find unitary V: solve M V = V G by diagonalizing generators? Try: choose V from aligning a generic Hermitian combo.
# Use two commuting generators: X12 and X45 (or g1g2 and ga gb) diagonal structures.
H12_JW=1j*JW[0]@JW[1]; H45_JW=1j*JW[3]@JW[4]
H12_MOD=1j*MOD[0]@MOD[1]; H45_MOD=1j*MOD[3]@MOD[4]
print("JW diag X12:",np.diag(H12_JW).real)
print("MOD diag X12:",np.diag(H12_MOD).real)
print("JW diag X45:",np.diag(H45_JW).real)
print("MOD diag X45:",np.diag(H45_MOD).real)
