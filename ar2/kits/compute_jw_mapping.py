#!/usr/bin/env python3
import numpy as np
from itertools import product

# Pauli matrices
sx = np.array([[0,1],[1,0]],dtype=complex)
sy = np.array([[0,-1j],[1j,0]],dtype=complex)
sz = np.array([[1,0],[0,-1]],dtype=complex)
si = np.eye(2,dtype=complex)

paulis = {'0':si,'x':sx,'y':sy,'z':sz}

def kron_multi(mats):
    res = np.array([1],dtype=complex)
    for m in mats:
        res = np.kron(res,m)
    return res

def op_on_site(op, site, L):
    mats = [si]*L
    mats[site]=op
    return kron_multi(mats)

def two_site_op(op1, i, op2, j, L):
    mats = [si]*L
    mats[i]=op1
    mats[j]=op2
    return kron_multi(mats)

def build_full_H(cmu, L):
    H = np.zeros((2**L,2**L),dtype=complex)
    # two-site terms
    for i in range(L-1):
        for mu,nu in product(['0','x','y','z'],repeat=2):
            c = cmu.get((mu,nu),0.0)
            if abs(c)==0: continue
            op1 = paulis[mu]
            op2 = paulis[nu]
            H += c * two_site_op(op1,i,op2,i+1,L)
    # optionally add single-site terms c_{mu,0} and c_{0,nu} already included above when mu or nu is '0'
    return H.real if np.allclose(H.imag,0,atol=1e-12) else H

def build_H0(cmu, L):
    # H0: include only two-site terms with mu,nu in {x,y} and single-site z terms (c_{z0}, c_{0z})
    H0 = np.zeros((2**L,2**L),dtype=complex)
    for i in range(L-1):
        for mu,nu in product(['x','y'],repeat=2):
            c = cmu.get((mu,nu),0.0)
            if abs(c)==0: continue
            H0 += c * two_site_op(paulis[mu],i,paulis[nu],i+1,L)
        # single-site z contributions from c_{z0} and c_{0z}
        c_z0 = cmu.get(('z','0'),0.0)
        c_0z = cmu.get(('0','z'),0.0)
        if abs(c_z0)!=0:
            H0 += c_z0 * op_on_site(sz,i,L)
        if abs(c_0z)!=0:
            H0 += c_0z * op_on_site(sz,i+1,L)
    return H0.real if np.allclose(H0.imag,0,atol=1e-12) else H0

def spectral_gap(H):
    w, _ = np.linalg.eigh(H)
    # sort
    w = np.sort(np.real(w))
    gap = w[1]-w[0] if len(w)>1 else 0.0
    return gap, w

def op_norm(H):
    # operator norm (spectral norm) for Hermitian is max abs eigenvalue
    w = np.linalg.eigvalsh(H)
    return max(abs(w))

def analyze_sample(cmu, L=6):
    H = build_full_H(cmu,L)
    H0 = build_H0(cmu,L)
    V = H - H0
    gap, eigs = spectral_gap(H0)
    normV = op_norm(V)
    eta = normV / gap if gap>1e-12 else np.inf
    # crude fidelity estimate between ground states via perturbation bound: 1 - (normV/gap)^2
    fid_est = max(0.0, 1.0 - (eta**2)) if np.isfinite(eta) else 0.0
    return {'gap':gap,'normV':normV,'eta':eta,'fid_est':fid_est,'H0_eigs':eigs[:6]}

if __name__=='__main__':
    L=6
    samples = {}
    # Sample A: diagonal
    cmuA = {('x','x'):1.0,('y','y'):0.8,('z','z'):0.2}
    samples['diag'] = cmuA
    # Sample B: small string terms
    cmuB = {('x','x'):1.0,('y','y'):0.8,('z','z'):0.2,('x','z'):0.1,('z','x'):0.1}
    samples['small_string'] = cmuB
    # Sample C: larger string
    cmuC = {('x','x'):1.0,('y','y'):0.8,('z','z'):0.2,('x','z'):0.5,('z','x'):0.5}
    samples['large_string'] = cmuC

    for name, cmu in samples.items():
        res = analyze_sample(cmu,L=L)
        print(f"Sample {name} -> gap(H0)={res['gap']:.6f}, ||V||={res['normV']:.6f}, eta={res['eta']:.6f}, fid_est≈{res['fid_est']:.6f}")
        print('  lowest H0 eigenvalues:', np.round(res['H0_eigs'],6))
        print()
