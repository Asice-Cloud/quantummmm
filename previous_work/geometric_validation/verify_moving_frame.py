#!/usr/bin/env python3
"""Check an exact moving frame and a proposed second-stage reduction.

The reduced model retains both accumulated rotation phases and the full
three-dimensional slow motion. It omits only the connection between the fast
rotation plane and its three-dimensional complement. These are local segment
tests with exact incoming wavefunctions, not full-protocol parameter scans.

Run from the project root:
  .venv/bin/python -B geometric_validation/verify_moving_frame.py
"""

import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad_vec, solve_ivp

from verify_bottleneck import CASES, PSI_TARGET


ROOT = Path(__file__).resolve().parent
TC = .3


def instantaneous_frame(v, E, eta):
    """Return Q, dQ/dv, mu_plus, mu_minus in original mode order.

    Column order: fast U,V; slow U,V; null direction. All derivatives are
    analytic. This chart requires nonzero singular values and distinct
    singular values, and is used here for the four positive E,eta examples.
    """
    gate = (1+np.cos(np.pi*v))/2
    dgate = -np.pi*np.sin(np.pi*v)/2
    p,q,r = TC*gate,eta*gate,TC*(1-gate)
    dp,dq,dr = TC*dgate,eta*dgate,-TC*dgate
    M = np.array([[E,0,q],[p,-r,0]])
    DM = np.array([[0,0,dq],[dp,-dr,0]])
    gram = M@M.T
    dgram = DM@M.T+M@DM.T
    a,b,d = gram[0,0],gram[0,1],gram[1,1]
    da,db,dd = dgram[0,0],dgram[0,1],dgram[1,1]
    splitting = np.hypot(a-d,2*b)
    theta = .5*np.arctan2(2*b,a-d)
    dtheta = ((a-d)*db-b*(da-dd))/splitting**2
    uplus = np.array([np.cos(theta),np.sin(theta)])
    uminus = np.array([-np.sin(theta),np.cos(theta)])
    duplus,duminus = dtheta*uminus,-dtheta*uplus
    lambda_plus = (a+d+splitting)/2
    # Use the determinant for the smaller root to avoid cancellation.
    lambda_minus = (E**2*r**2+q**2*(p**2+r**2))/lambda_plus
    muplus,muminus = np.sqrt(lambda_plus),np.sqrt(lambda_minus)
    dmuplus = (uplus@dgram@uplus)/(2*muplus)
    dmuminus = (uminus@dgram@uminus)/(2*muminus)
    vplus,vminus = M.T@uplus/muplus,M.T@uminus/muminus
    dvplus = (DM.T@uplus+M.T@duplus-vplus*dmuplus)/muplus
    dvminus = (DM.T@uminus+M.T@duminus-vminus*dmuminus)/muminus
    null = np.array([q*r,p*q,-E*r])
    dnull = np.array([dq*r+q*dr,dp*q+p*dq,-E*dr])
    n = null/np.linalg.norm(null)
    dn = (dnull-n*(n@dnull))/np.linalg.norm(null)

    Q,DQ = np.zeros((5,5)),np.zeros((5,5))
    Q[:2,0],Q[:2,2] = uplus,uminus
    Q[2:,1],Q[2:,3],Q[2:,4] = vplus,vminus,n
    DQ[:2,0],DQ[:2,2] = duplus,duminus
    DQ[2:,1],DQ[2:,3],DQ[2:,4] = dvplus,dvminus,dn
    # Fix the null-axis sign so the full frame is in SO(5).
    if np.linalg.det(Q) < 0:
        Q[:,4] *= -1
        DQ[:,4] *= -1
    inverse_order = np.argsort([0,3,1,2,4])  # Undo (1,a)|(2,3,b).
    return Q[inverse_order],DQ[inverse_order],muplus,muminus


def instantaneous_generator(v, case):
    gate = (1+np.cos(np.pi*v))/2
    E,p,q,r = case["E1"],TC*gate,case["eta"]*gate,TC*(1-gate)
    return np.array([[0,E,0,0,q],[-E,0,0,-p,0],[0,0,0,r,0],
                     [0,p,-r,0,0],[-q,0,0,0,0]])


def omega(muplus, muminus, scale):
    matrix = np.zeros((5,5))
    matrix[0,1],matrix[1,0] = scale*muplus,-scale*muplus
    matrix[2,3],matrix[3,2] = scale*muminus,-scale*muminus
    return matrix


def frame_diagnostics(case):
    scale = 2*case["tau"]
    E,eta = case["E1"],case["eta"]
    fast,slow = [],[]
    errors = {"orthogonal_frame":0.,"special_orthogonal_frame":0.,
              "skew_connection":0.,"instantaneous_rotation_blocks":0.,
              "analytic_frame_derivative":0.,"slow_connection_structure":0.}
    minimum_fast_gap = float("inf")
    for v in np.linspace(0,1,1001):
        Q,DQ,muplus,muminus = instantaneous_frame(v,E,eta)
        connection = Q.T@DQ
        A = instantaneous_generator(v,case)
        minimum_fast_gap = min(minimum_fast_gap,muplus-muminus)
        fast.append(np.linalg.norm(connection[:2,2:],2)
                    /(scale*(muplus-muminus)))
        slow.append(np.linalg.norm(connection[2:4,4],2)/(scale*muminus))
        errors["orthogonal_frame"] = max(errors["orthogonal_frame"],
                                         float(np.max(abs(Q.T@Q-np.eye(5)))))
        errors["special_orthogonal_frame"] = max(
            errors["special_orthogonal_frame"],float(abs(np.linalg.det(Q)-1)),
        )
        errors["skew_connection"] = max(errors["skew_connection"],
                                        float(np.max(abs(connection+connection.T))))
        errors["instantaneous_rotation_blocks"] = max(
            errors["instantaneous_rotation_blocks"],
            float(np.max(abs(Q.T@A@Q-omega(muplus,muminus,1)))),
        )
        slow_block = connection[2:,2:]
        expected_slow = np.zeros((3,3))
        expected_slow[1,2] = connection[3,4]
        expected_slow[2,1] = -connection[3,4]
        errors["slow_connection_structure"] = max(
            errors["slow_connection_structure"],
            float(np.max(abs(slow_block-expected_slow))),
        )
    # Independently check the analytic derivative without using it in the
    # finite-difference reference; do not differentiate at chart endpoints.
    step = 1e-5
    for v in (.1,.35,.6,.85,.95):
        finite_difference = (
            -instantaneous_frame(v+2*step,E,eta)[0]
            +8*instantaneous_frame(v+step,E,eta)[0]
            -8*instantaneous_frame(v-step,E,eta)[0]
            +instantaneous_frame(v-2*step,E,eta)[0]
        )/(12*step)
        errors["analytic_frame_derivative"] = max(
            errors["analytic_frame_derivative"],
            float(np.max(abs(finite_difference-instantaneous_frame(v,E,eta)[1]))),
        )
    assert max(errors.values()) < 2e-8, (case,errors)
    phases,_ = quad_vec(
        lambda v: scale*np.array(instantaneous_frame(v,E,eta)[2:]),
        0,1,epsabs=1e-11,epsrel=2e-12,
    )
    return {**case,"max_fast_mixing_indicator":float(max(fast)),
            "max_slow_null_mixing_indicator":float(max(slow)),
            "minimum_mu_plus_minus_mu_minus":minimum_fast_gap,
            "fast_phase":float(phases[0]),"slow_phase":float(phases[1]),
            "errors":errors}


def second_stage_solution(case, incoming_state, grid, model):
    E,eta,scale = case["E1"],case["eta"],2*case["tau"]
    initial = instantaneous_frame(0,E,eta)[0].T@incoming_state
    if model == "exact_frame":
        def rhs(v,state):
            Q,DQ,muplus,muminus = instantaneous_frame(v,E,eta)
            return (omega(muplus,muminus,scale)-Q.T@DQ)@state
        solution = solve_ivp(rhs,(0,1),initial,method="DOP853",
                             rtol=8e-13,atol=8e-15,t_eval=grid,max_step=.005)
        coordinates = solution.y.T
    else:
        # Evolve one fast phase and a complex three-component slow state.
        # The fast block's internal connection is zero in this frame.
        reduced_initial = np.concatenate(([0j],initial[2:]))
        def rhs(v,state):
            Q,DQ,muplus,muminus = instantaneous_frame(v,E,eta)
            gamma = (Q.T@DQ)[3,4] if model == "slow_geometry" else 0.
            slow_generator = np.array([[0,scale*muminus,0],
                                       [-scale*muminus,0,-gamma],[0,gamma,0]])
            return np.concatenate(([scale*muplus],slow_generator@state[1:]))
        solution = solve_ivp(rhs,(0,1),reduced_initial,method="DOP853",
                             rtol=8e-13,atol=8e-15,t_eval=grid,max_step=.005)
        phases = solution.y[0].real
        fast = np.stack([np.cos(phases)*initial[0]+np.sin(phases)*initial[1],
                         -np.sin(phases)*initial[0]+np.cos(phases)*initial[1]],axis=1)
        coordinates = np.concatenate((fast,solution.y[1:].T),axis=1)
    assert solution.success
    return np.array([instantaneous_frame(v,E,eta)[0]@state
                     for v,state in zip(grid,coordinates)])


def slow_sphere_check(case, incoming_state, reduced_result):
    E,eta,scale = case["E1"],case["eta"],2*case["tau"]
    assert 0 < E < TC
    initial = instantaneous_frame(0,E,eta)[0].T@incoming_state
    slow = initial[2:]
    slow_norm = float(np.sum(abs(slow)**2))
    fast_norm = float(np.sum(abs(initial[:2])**2))
    null_amplitude = float(abs(np.sum(slow*slow)))
    B = 2*np.imag(slow[:,None]*slow.conj()[None,:])
    b0 = np.array([B[1,2],B[2,0],B[0,1]])
    radius = float(np.linalg.norm(b0))

    def rhs(v,state):
        Q,DQ,_,muminus = instantaneous_frame(v,E,eta)
        gamma = (Q.T@DQ)[3,4]
        generator = np.array([[0,scale*muminus,0],
                              [-scale*muminus,0,-gamma],[0,gamma,0]])
        return np.concatenate((generator@state[:3],[scale*muminus]))

    solution = solve_ivp(rhs,(0,1),np.concatenate((b0,[0.])),
                         method="DOP853",rtol=8e-13,atol=8e-15,
                         dense_output=True,max_step=.005)
    assert solution.success
    b = solution.y[:3,-1]
    d_from_sphere = float(b[2])
    F_from_sphere = float(slow_norm*(1-b[2]/radius)**2/4)
    F_error = abs(F_from_sphere-reduced_result["F_out"])
    F_bound = null_amplitude/4+max(0.,slow_norm-radius)/2
    d_error = abs(d_from_sphere-reduced_result["d_out"])
    assert d_error < 2e-9
    assert F_error <= F_bound+2e-10, (case,F_error,F_bound)
    assert abs(radius**2-(slow_norm**2-null_amplitude**2)) < 2e-10
    assert null_amplitude <= fast_norm+2e-10

    # Verify the scalar memory equation for b_z, keeping its full slow phase.
    memory_checks = []
    for v in (.2,.6,.95):
        state = solution.sol(v)
        phase = state[3]
        Q,DQ,_,_ = instantaneous_frame(v,E,eta)
        gamma = (Q.T@DQ)[3,4]

        def integrand(s):
            Qs,DQs,_,_ = instantaneous_frame(s,E,eta)
            gammas = (Qs.T@DQs)[3,4]
            sphere = solution.sol(s)
            return gammas*np.cos(phase-sphere[3])*sphere[2]

        integral,_ = quad_vec(integrand,0,v,epsabs=1e-12,epsrel=2e-11)
        predicted = gamma*(b0[1]*np.cos(phase)-b0[0]*np.sin(phase)-integral)
        error = float(abs(predicted-gamma*state[1]))
        assert error < 2e-9, (case,v,error)
        memory_checks.append({"v":v,"derivative_error":error})
    return {"fast_input_weight":fast_norm,"slow_input_weight":slow_norm,
            "slow_null_amplitude":null_amplitude,"sphere_radius":radius,
            "initial_b":b0.tolist(),"d_out":d_from_sphere,
            "F_out":F_from_sphere,"d_vs_slow_model_error":float(d_error),
            "F_vs_slow_model_error":float(F_error),
            "F_input_null_error_bound":float(F_bound),
            "memory_checks":memory_checks}


def segment_checks(case, stage, wave_data):
    name = case["name"]
    reference = wave_data[f"{name}_psi"][stage*800:(stage+1)*800+1]
    grid = np.linspace(0,1,801)
    d_reference = 2*np.imag(reference[:,0]*reference[:,1].conj())
    F_reference = abs(reference@PSI_TARGET.conj())**2
    models = {}
    for model in ("exact_frame","slow_geometry","phases_only"):
        states = second_stage_solution(case,reference[0],grid,model)
        d = 2*np.imag(states[:,0]*states[:,1].conj())
        F = abs(states@PSI_TARGET.conj())**2
        models[model] = {
            "d_out":float(d[-1]),"F_out":float(F[-1]),
            "d_endpoint_error":float(abs(d[-1]-d_reference[-1])),
            "F_endpoint_error":float(abs(F[-1]-F_reference[-1])),
            "d_max_path_error":float(np.max(abs(d-d_reference))),
            "F_max_path_error":float(np.max(abs(F-F_reference))),
            "wave_max_error":float(np.max(abs(states-reference))),
            "wave_norm_error":float(np.max(abs(np.sum(abs(states)**2,axis=1)-1))),
            "wave_null_error":float(np.max(abs(np.sum(states*states,axis=1)))),
        }
        assert max(models[model][key] for key in ("wave_norm_error","wave_null_error")) < 2e-9
    exact = models["exact_frame"]
    assert max(exact[key] for key in ("wave_max_error","d_max_path_error","F_max_path_error")) < 2e-9
    # Approximation errors are reported, not used as universal guarantees.
    sphere = slow_sphere_check(case,reference[0],models["slow_geometry"])
    return {**case,"segment":stage+1,"d_reference":float(d_reference[-1]),
            "F_reference":float(F_reference[-1]),"models":models,
            "slow_sphere":sphere}


def main():
    cases = CASES[1:]  # E=eta=0 needs a different frame at degeneracies.
    result = {"units":"hbar=1; energies in meV; tau in meV^-1",
              "scope":"Local second-stage tests with exact incoming states",
              "diagnostics":[frame_diagnostics(c) for c in cases]}
    with np.load(ROOT/"bottleneck.npz") as wave_data:
        result["segments"] = [segment_checks(c,stage,wave_data)
                              for c in cases for stage in (1,4)]
    (ROOT/"moving_frame.json").write_text(
        json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8",
    )
    for c in result["diagnostics"]:
        print(c["name"],"fast/slow mixing indicators:",
              c["max_fast_mixing_indicator"],c["max_slow_null_mixing_indicator"],
              "phases:",c["fast_phase"],c["slow_phase"])
    for segment in result["segments"]:
        print(segment["name"],"segment",segment["segment"],
              "reference d/F:",segment["d_reference"],segment["F_reference"])
        for model,values in segment["models"].items():
            print(model,"d/F:",values["d_out"],values["F_out"],
                  "endpoint errors:",values["d_endpoint_error"],values["F_endpoint_error"])
        sphere = segment["slow_sphere"]
        print("sphere d/F:",sphere["d_out"],sphere["F_out"],
              "fast input weight:",sphere["fast_input_weight"],
              "max scalar-memory error:",max(c["derivative_error"]
              for c in sphere["memory_checks"]))


if __name__ == "__main__":
    main()
