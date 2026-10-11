#!/usr/bin/env python3
"""Verify the finite-parameter second-stage block and memory equations.

Use direct complex-wavefunction evolution as an independent reference. Check
both second stages of the double cycle, the explicit Volterra kernel, and the
endpoint memory budget. Controlled input/operator swaps diagnose the existing
finite-parameter examples; they are not additional full-protocol scans.

Run from the project root:
  .venv/bin/python -B geometric_validation/verify_second_stage.py
"""

import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad_vec, solve_ivp
import sympy as sp

from verify_bottleneck import CASES, PSI_TARGET, controls, wave_rhs


ROOT = Path(__file__).resolve().parent
TC = .3


def symbolic_checks():
    E, p, q, r = sp.symbols("E p q r", real=True)
    d, x, y, z, u, v, w, h, j, k = sp.symbols(
        "d x y z u v w h j k", real=True
    )
    A = sp.Matrix([[0,E,0,0,q], [-E,0,0,-p,0], [0,0,0,r,0],
                   [0,p,-r,0,0], [-q,0,0,0,0]])
    B = sp.Matrix([[0,d,x,y,z], [-d,0,u,v,w], [-x,-u,0,h,j],
                   [-y,-v,-h,0,k], [-z,-w,-j,-k,0]])
    DB = A*B-B*A
    X = sp.Matrix([y,u,w,j])
    Y = sp.Matrix([d,x,z,-v,-h,k])
    L = sp.Matrix([[-p,0,-q,0], [r,E,0,-q], [0,0,E,0],
                   [E,r,0,0], [0,p,0,0], [q,0,p,-r]])
    DX = sp.Matrix([DB[0,3],DB[1,2],DB[1,4],DB[2,4]])
    DY = sp.Matrix([DB[0,1],DB[0,2],DB[0,4],
                    -DB[1,3],-DB[2,3],DB[3,4]])
    assert sp.expand(DX+L.T*Y) == sp.zeros(4,1)
    assert sp.expand(DY-L*X) == sp.zeros(6,1)
    assert sp.expand((X.T*DX+Y.T*DY)[0]) == 0
    ps, qs, rs = sp.symbols("p_s q_s r_s", real=True)
    row = (L*L.subs({p:ps,q:qs,r:rs}).T)[0,:]
    expected = sp.Matrix([[p*ps+q*qs,-p*rs,-E*q,-E*p,0,-p*qs-q*ps]])
    assert sp.expand(row-expected) == sp.zeros(1,6)
    ordering = [0,3,1,2,4]  # (1,a) | (2,3,b)
    reordered = A.extract(ordering,ordering)
    M = sp.Matrix([[E,0,q],[p,-r,0]])
    assert reordered == sp.BlockMatrix(
        [[sp.zeros(2),M],[-M.T,sp.zeros(3)]]
    ).as_explicit()
    dark_axis = sp.Matrix([0,q*r,p*q,0,-E*r])
    assert sp.expand(A*dark_axis) == sp.zeros(5,1)
    gram = M*M.T
    assert gram == sp.Matrix([[E**2+q**2,E*p],[E*p,p**2+r**2]])
    assert sp.expand(gram.det()-(dark_axis.T*dark_axis)[0]) == 0
    lam = sp.symbols("lambda")
    spectrum = lam*(lam**4+(E**2+q**2+p**2+r**2)*lam**2
                   +E**2*r**2+q**2*(p**2+r**2))
    assert sp.expand(A.charpoly(lam).as_expr()-spectrum) == 0
    return {"bipartite_generator": True, "block_equations": True,
            "norm_conservation": True, "first_kernel_row": True,
            "instantaneous_dark_axis": True,
            "rotation_frequencies": True}


def wave_bivector(state):
    return 2*np.imag(state[..., :,None]*state.conj()[...,None,:])


def fields_from_matrix(B):
    X = np.stack([B[...,0,3],B[...,1,2],B[...,1,4],B[...,2,4]],axis=-1)
    Y = np.stack([B[...,0,1],B[...,0,2],B[...,0,4],
                  -B[...,1,3],-B[...,2,3],B[...,3,4]],axis=-1)
    return X, Y


def matrix_from_fields(values):
    X, Y = values[..., :4], values[..., 4:]
    B = np.zeros(values.shape[:-1]+(5,5))
    for i,j,component in [(0,3,X[...,0]),(1,2,X[...,1]),
                          (1,4,X[...,2]),(2,4,X[...,3]),
                          (0,1,Y[...,0]),(0,2,Y[...,1]),
                          (0,4,Y[...,2]),(1,3,-Y[...,3]),
                          (2,3,-Y[...,4]),(3,4,Y[...,5])]:
        B[...,i,j] = component
    return B-np.swapaxes(B,-1,-2)


def block(v, case):
    p, r, _, gate = controls(1,v)
    E, q = case["E1"], case["eta"]*gate
    return np.array([[-p,0,-q,0],[r,E,0,-q],[0,0,E,0],
                     [E,r,0,0],[0,p,0,0],[q,0,p,-r]])


def block_rhs(v, values, case):
    X, Y = values[:4], values[4:]
    L = block(v,case)
    return 2*case["tau"]*np.concatenate((-L.T@Y,L@X))


def full_segment_check(case, stage, wave_data):
    name = case["name"]
    start = stage*800
    psi = wave_data[f"{name}_psi"][start:start+801]
    X0,Y0 = fields_from_matrix(wave_bivector(psi[0]))
    solution = solve_ivp(
        lambda v,y:block_rhs(v,y,case),(0,1),np.concatenate((X0,Y0)),
        method="DOP853",rtol=2e-12,atol=2e-14,
        t_eval=np.linspace(0,1,801),max_step=.01,
    )
    assert solution.success
    values = solution.y.T
    B = matrix_from_fields(values)
    X,Y = values[:,:4],values[:,4:]
    d = Y[:,0]
    F = (d*d-d)/2+(np.sum(X[:,:3]**2,axis=1)
                       +np.sum(Y[:,1:4]**2,axis=1))/4
    norms = np.sum(values**2,axis=1)
    errors = {
        "full_wave_bivector":float(np.max(abs(B-wave_bivector(psi)))),
        "full_wave_fidelity":float(np.max(abs(F-abs(psi@PSI_TARGET.conj())**2))),
        "conserved_norm":float(np.max(abs(norms-norms[0]))),
        "simple_plane_constraint":float(np.max(abs(B@B@B+B))),
    }
    assert max(errors.values()) < 2e-9, errors
    return {**case,"segment":stage+1,"d_in":float(d[0]),
            "d_out":float(d[-1]),"F_out":float(F[-1]),"errors":errors}


def dense_wave_solution(case, state):
    solution = solve_ivp(
        lambda v,y:wave_rhs(4,v,y,case),(0,1),state,
        method="DOP853",rtol=8e-13,atol=8e-15,
        dense_output=True,max_step=.005,
    )
    assert solution.success
    return solution


def memory_check(case, solution):
    X0,Y0 = fields_from_matrix(wave_bivector(solution.sol(0)))
    scale = 2*case["tau"]
    checks = []
    for v in (.2,.6,.95):
        state = solution.sol(v)
        derivative = wave_rhs(4,v,state,case)
        B = wave_bivector(state)
        DB = 2*np.imag(np.outer(derivative,state.conj())
                       +np.outer(state,derivative.conj()))
        X,Y = fields_from_matrix(B)
        _,DY = fields_from_matrix(DB)

        def integrand(s):
            _,Ys = fields_from_matrix(wave_bivector(solution.sol(s)))
            return block(s,case).T@Ys

        integral,quad_error = quad_vec(integrand,0,v,epsabs=1e-13,
                                       epsrel=5e-12,limit=600)
        predicted_X = X0-scale*integral
        predicted_DY = scale*block(v,case)@X0-scale**2*block(v,case)@integral

        def d_integrand(s):
            Bs = wave_bivector(solution.sol(s))
            _,rs,_,gs = controls(1,s)
            return (gs*((TC**2+case["eta"]**2)*Bs[0,1]
                        -2*TC*case["eta"]*Bs[3,4])
                    -TC*rs*Bs[0,2]
                    +case["E1"]*(TC*Bs[1,3]-case["eta"]*Bs[0,4]))

        d_integral,d_quad_error = quad_vec(
            d_integrand,0,v,epsabs=1e-13,epsrel=5e-12,limit=600,
        )
        _,_,_,gv = controls(1,v)
        predicted_dd = (-scale*gv*(TC*X0[0]+case["eta"]*X0[2])
                        -scale**2*gv*d_integral)
        checks.append({
            "v":v,"X_reconstruction_error":float(np.max(abs(predicted_X-X))),
            "Y_derivative_error":float(np.max(abs(predicted_DY-DY))),
            "d_memory_row_error":float(abs(predicted_dd-DY[0])),
            "estimated_integral_error":float(quad_error),
            "estimated_d_integral_error":float(d_quad_error),
        })
    worst = max(max(c[key] for key in ("X_reconstruction_error",
                "Y_derivative_error","d_memory_row_error")) for c in checks)
    assert worst < 2e-8, checks
    return {**case,"max_error":worst,"checks":checks}


def endpoint_budget(case, solution):
    Bin = wave_bivector(solution.sol(0))
    Bout = wave_bivector(solution.sol(1))
    scale,eta,E = 2*case["tau"],case["eta"],case["E1"]

    def contributions(v):
        B = wave_bivector(solution.sol(v))
        _,r,_,g = controls(1,v)
        # G/τ is the remaining integral of g over normalized local time.
        tail = (1-v)/2-np.sin(np.pi*v)/(2*np.pi)
        return np.array([
            -scale**2*(TC**2+eta**2)*tail*g*B[0,1],
            2*scale**2*TC*eta*tail*g*B[3,4],
            scale**2*TC*tail*r*B[0,2],
            -scale**2*E*tail*(TC*B[1,3]-eta*B[0,4]),
            -scale*g*(TC*B[0,3]+eta*B[1,4]),
        ])

    terms,quad_error = quad_vec(contributions,0,1,epsabs=1e-11,
                                epsrel=2e-12,limit=800)
    input_term = -scale*.5*(TC*Bin[0,3]+eta*Bin[1,4])
    budget = {"input":float(input_term),"d_memory":float(terms[0]),
              "ab_memory":float(terms[1]),"mode_3_memory":float(terms[2]),
              "E_quadrature_memory":float(terms[3])}
    change = float(Bout[0,1]-Bin[0,1])
    residual = abs(sum(budget.values())-change)
    driver_residual = abs(float(terms[4])-change)
    assert max(residual,driver_residual) < 2e-8, (case,budget,change,residual)
    return {**case,"d_in":float(Bin[0,1]),"d_out":float(Bout[0,1]),
            "d_change":change,"initial_d_velocity":float(
                -2*(TC*Bin[0,3]+eta*Bin[1,4])),
            "terms":budget,"sum_residual":float(residual),
            "integrated_driver_residual":float(driver_residual),
            "estimated_quadrature_error":float(quad_error)}


def swap_diagnostics(cases, wave_data, actual_solutions):
    pairs = [("eta_low","high_fidelity"),("high_fidelity","E_changed")]
    results = []
    for pair in pairs:
        entries = []
        for incoming in pair:
            for operator in pair:
                c = cases[operator]
                if incoming == operator:
                    solution = actual_solutions[incoming]
                else:
                    solution = dense_wave_solution(
                        c,wave_data[f"{incoming}_psi"][3200],
                    )
                state = solution.sol(1)
                entries.append({
                    "incoming_state":incoming,"second_stage_parameters":operator,
                    "d_out":float(wave_bivector(state)[0,1]),
                    "F_out":float(abs(np.vdot(PSI_TARGET,state))**2),
                })
        results.append({"cases":list(pair),"entries":entries})
    return results


def main():
    result = {"units":"hbar=1; energy in meV; tau in meV^-1",
              "symbolic":symbolic_checks()}
    cases = {c["name"]:c for c in CASES}
    with np.load(ROOT/"bottleneck.npz") as wave_data:
        result["full_segments"] = [
            full_segment_check(c,stage,wave_data)
            for c in CASES for stage in (1,4)
        ]
        solutions = {
            c["name"]:dense_wave_solution(c,wave_data[f"{c['name']}_psi"][3200])
            for c in CASES
        }
        result["memory_equations"] = [
            memory_check(c,solutions[c["name"]])
            for c in CASES if c["name"] in ("high_fidelity","long_protocol")
        ]
        result["endpoint_budgets"] = [
            endpoint_budget(c,solutions[c["name"]]) for c in CASES
        ]
        result["input_operator_swaps"] = swap_diagnostics(cases,wave_data,solutions)
    (ROOT/"second_stage.json").write_text(
        json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8",
    )
    print("Max full-segment error:",max(max(c["errors"].values())
          for c in result["full_segments"]))
    print("Memory-equation errors:",[c["max_error"]
          for c in result["memory_equations"]])
    print("Max endpoint-budget residual:",max(c["sum_residual"]
          for c in result["endpoint_budgets"]))
    for row in result["endpoint_budgets"]:
        print(row["name"],"d_in/out:",row["d_in"],row["d_out"],
              "initial d velocity:",row["initial_d_velocity"],
              "budget:",row["terms"])
    for group in result["input_operator_swaps"]:
        print("Input/operator swaps:",group)


if __name__ == "__main__":
    main()
