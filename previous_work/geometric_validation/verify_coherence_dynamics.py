#!/usr/bin/env python3
"""Check the exact coherence channels against full-wavefunction evolution.

Run from the project root:
  .venv/bin/python -B geometric_validation/verify_coherence_dynamics.py

Uses the full-wavefunction paths saved by verify_bottleneck.py, verifies the
scalar memory equations on independent dense wave solutions, and tests the
first-segment invariant at eta=tc. Units: hbar=1, energies in meV.
"""

import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp
import sympy as sp

from verify_bottleneck import CASES, PSI_INITIAL, controls, wave_rhs


ROOT = Path(__file__).resolve().parent


def symbolic_checks():
    E, q, p, r, ed = sp.symbols("E q p r ed", real=True)
    d, x, y, z, u, v, w, h, j, k = sp.symbols(
        "d x y z u v w h j k", real=True
    )
    A = sp.Matrix([[0,E,0,0,q], [-E,0,0,-p,0], [0,0,0,r,0],
                   [0,p,-r,0,ed], [-q,0,0,-ed,0]])
    B = sp.Matrix([[0,d,x,y,z], [-d,0,u,v,w], [-x,-u,0,h,j],
                   [-y,-v,-h,0,k], [-z,-w,-j,-k,0]])
    DB = A*B-B*A
    Z3, Za, Zb = x+sp.I*u, y+sp.I*v, z+sp.I*w
    expressions = [
        (-q*w-p*y, DB[0,1]),
        (-sp.I*E*Z3+r*Za+sp.I*p*h-q*j, DB[0,2]+sp.I*DB[1,2]),
        (-sp.I*E*Za-r*Z3+ed*Zb+p*d-q*k, DB[0,3]+sp.I*DB[1,3]),
        (-sp.I*E*Zb-ed*Za+sp.I*q*d-sp.I*p*k, DB[0,4]+sp.I*DB[1,4]),
        (ed*j-p*u, DB[2,3]),
        (r*k-ed*h+q*x, DB[2,4]),
        (-r*j+q*y+p*w, DB[3,4]),
    ]
    assert all(sp.expand(a-b) == 0 for a, b in expressions)
    for sign in (1, -1):
        X, C = d+sign*k, Za+sign*sp.I*Zb
        V, Delta = p-sign*q, E+sign*ed
        dX = DB[0,1]+sign*DB[3,4]
        dC = (DB[0,3]+sp.I*DB[1,3]
              +sign*sp.I*(DB[0,4]+sp.I*DB[1,4]))
        assert sp.expand(dX+V*sp.re(C)+sign*r*j) == 0
        assert sp.expand(dC+sp.I*Delta*C-V*X+r*Z3) == 0
        assert sp.expand((2*X*dX+2*sp.re(sp.conjugate(C)*dC)).subs(r,0)) == 0
    return {"ten_component_equations": True, "channel_equations": True,
            "sphere_norms_when_t3_zero": True}


def matrix_from_fields(values):
    d = values[..., 0]
    z = values[..., 1:4]+1j*values[..., 4:7]
    h = values[..., 7:10]
    B = np.zeros(values.shape[:-1]+(5,5))
    B[...,0,1] = d
    B[...,0,2:] = z.real
    B[...,1,2:] = z.imag
    B[...,2,3], B[...,2,4], B[...,3,4] = np.moveaxis(h,-1,0)
    return B-np.swapaxes(B,-1,-2)


def coherence_rhs(stage, v, values, case):
    d = values[0]
    z3, za, zb = values[1:4]+1j*values[4:7]
    h3a, h3b, hab = values[7:10]
    p, r, ed, gate = controls(stage, v)
    E, q = case["E1"], case["eta"]*gate
    dz = np.array([
        -1j*E*z3+r*za+1j*p*h3a-q*h3b,
        -1j*E*za-r*z3+ed*zb+p*d-q*hab,
        -1j*E*zb-ed*za+1j*q*d-1j*p*hab,
    ])
    dh = np.array([ed*h3b-p*z3.imag,
                   r*hab-ed*h3a+q*z3.real,
                   -r*h3b+q*za.real+p*zb.imag])
    dd = -q*zb.imag-p*za.real
    return 2*case["tau"]*np.concatenate(([dd], dz.real, dz.imag, dh))


def third_segment_rotation(tau):
    def rhs(v, values):
        _, r, ed, _ = controls(2, v)
        D = np.array([[0,r,0],[-r,0,ed],[0,-ed,0]])
        return (2*tau*D@values.reshape(3,3)).ravel()
    solution = solve_ivp(rhs,(0,1),np.eye(3).ravel(),method="DOP853",
                         rtol=5e-13,atol=5e-15)
    assert solution.success
    return solution.y[:,-1].reshape(3,3)


def full_path_check(case, wave_data, Q3):
    current = np.zeros(10)
    current[0] = 1
    path = []
    grid = np.linspace(0,1,801)
    for stage in range(6):
        solution = solve_ivp(
            lambda v,y: coherence_rhs(stage,v,y,case), (0,1), current,
            method="DOP853",rtol=2e-12,atol=2e-14,t_eval=grid,
        )
        assert solution.success
        path.append(solution.y.T if stage == 0 else solution.y.T[1:])
        current = solution.y[:,-1]
    fields = np.concatenate(path)
    B = matrix_from_fields(fields)
    name = case["name"]
    psi = wave_data[f"{name}_psi"]
    B_wave = 2*np.imag(psi[..., :,None]*psi.conj()[...,None,:])
    d = fields[:,0]
    z = fields[:,1:4]+1j*fields[:,4:7]
    F = (d*d-d)/2+np.sum(abs(z)**2,axis=1)/4
    Cplus, Cminus = z[:,1]+1j*z[:,2], z[:,1]-1j*z[:,2]
    Xplus, Xminus = d+fields[:,9], d-fields[:,9]
    errors = {
        "full_wave_coherence_matrix": float(np.max(abs(B-B_wave))),
        "full_wave_fidelity": float(np.max(abs(F-wave_data[f"{name}_F"]))),
    }
    for segment in (0,3):
        sl = slice(segment*800,(segment+1)*800+1)
        for label, X, C in (("plus",Xplus,Cplus),("minus",Xminus,Cminus)):
            norms = X[sl]**2+abs(C[sl])**2
            errors[f"sphere_norm_segment_{segment+1}_{label}"] = float(
                np.max(abs(norms-norms[0]))
            )
    for segment in (2,5):
        start,end = segment*800,(segment+1)*800
        predicted = np.exp(-2j*case["E1"]*case["tau"])*z[start]@Q3.T
        errors[f"third_segment_phase_{segment+1}"] = float(
            np.max(abs(predicted-z[end]))
        )
        errors[f"third_segment_d_constant_{segment+1}"] = float(
            np.max(abs(d[start:end+1]-d[start]))
        )
        errors[f"third_segment_F_constant_{segment+1}"] = float(
            np.max(abs(F[start:end+1]-F[start]))
        )
    assert max(errors.values()) < 2e-9, errors
    return {**case,"F_end":float(F[-1]),"d_end":float(d[-1]),
            "d_at_segment_boundaries":[float(d[i*800]) for i in range(7)],
            "errors":errors}


def first_wave_solution(case):
    solution = solve_ivp(lambda v,y:wave_rhs(0,v,y,case),(0,1),PSI_INITIAL,
                         method="DOP853",rtol=8e-13,atol=8e-15,
                         dense_output=True,max_step=.01)
    assert solution.success
    return solution


def wave_bivector(state):
    return 2*np.imag(np.outer(state,state.conj()))


def memory_check(case):
    solution = first_wave_solution(case)
    scale = 2*case["tau"]
    def X(v,sign):
        B = wave_bivector(solution.sol(v))
        return B[0,1]+sign*B[3,4]
    checks = []
    for v in (.25,.6,.95):
        state = solution.sol(v)
        derivative = wave_rhs(0,v,state,case)
        dB = 2*np.imag(np.outer(derivative,state.conj())
                       +np.outer(state,derivative.conj()))
        for sign in (1,-1):
            def integrand(s):
                phase = scale*(case["E1"]*(v-s)
                        +sign*.3*((v-s)/2+(np.sin(np.pi*v)-np.sin(np.pi*s))/(2*np.pi)))
                return (1-np.cos(np.pi*s))/2*np.cos(phase)*X(s,sign)
            integral,error = quad(integrand,0,v,epsabs=1e-13,
                                  epsrel=5e-11,limit=500)
            prefactor = scale**2*(.3-sign*case["eta"])**2*(1-np.cos(np.pi*v))/2
            predicted = -prefactor*integral
            actual = dB[0,1]+sign*dB[3,4]
            checks.append({"v":v,"sign":sign,
                           "derivative_error":float(abs(predicted-actual)),
                           "estimated_quadrature_error":float(prefactor*error)})
    worst = max(check["derivative_error"] for check in checks)
    assert worst < 2e-8, checks
    return {**case,"max_derivative_error":worst,"checks":checks}


def invariant_checks():
    results = []
    for E in (0.,.01,.08):
        case = {"E1":E,"eta":.3,"tau":100.}
        solution = first_wave_solution(case)
        states = solution.sol(np.linspace(0,1,1001)).T
        B = 2*np.imag(states[..., :,None]*states.conj()[...,None,:])
        error = float(np.max(abs(B[:,0,1]+B[:,3,4]-1)))
        minimum_d = float(np.min(B[:,0,1]))
        assert error < 1e-10 and minimum_d > -1e-10
        results.append({**case,"Xplus_equals_one_error":error,
                        "minimum_sampled_d":minimum_d})
    return results


def main():
    result = {"symbolic":symbolic_checks()}
    rotations = {tau:third_segment_rotation(tau) for tau in {c["tau"] for c in CASES}}
    with np.load(ROOT/"bottleneck.npz") as wave_data:
        result["full_paths"] = [full_path_check(c,wave_data,rotations[c["tau"]])
                                for c in CASES]
    result["memory_equations"] = [memory_check(c) for c in CASES
                                  if c["name"] in ("high_fidelity","long_protocol")]
    result["eta_equals_tc_first_segment_invariant"] = invariant_checks()
    (ROOT/"coherence_dynamics.json").write_text(
        json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8"
    )
    for case in result["full_paths"]:
        print(case["name"],"d at boundaries:",case["d_at_segment_boundaries"],
              "max error:",max(case["errors"].values()))
    print("Scalar memory equation errors:",
          [c["max_derivative_error"] for c in result["memory_equations"]])
    print("eta=tc invariant errors:",
          [c["Xplus_equals_one_error"] for c in result["eta_equals_tc_first_segment_invariant"]])


if __name__ == "__main__":
    main()
