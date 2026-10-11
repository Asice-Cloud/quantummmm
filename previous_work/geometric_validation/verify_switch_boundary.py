#!/usr/bin/env python3
"""Predict and independently verify a fixed-input orientation boundary.

The incoming state is eta_low at 4*tau from bottleneck.npz. Keep that same
physical state while changing only the second-stage coupling or duration.
No fit to the complete wavefunction is used for the predictions. The sphere
keeps the full slow phase; a Dyson response in the omitted fast/slow connection
provides a correction while retaining all finite energy-dependent phases.

Run from the project root:
  .venv/bin/python -B geometric_validation/verify_switch_boundary.py
"""

import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/quantumsss-matplotlib")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from verify_bottleneck import PSI_INITIAL, PSI_TARGET, wave_rhs
from verify_moving_frame import instantaneous_frame, omega


ROOT = Path(__file__).resolve().parent
E = .01
ETA_LOW, ETA_HIGH = .01, .023
STANDARD = {"rtol":2e-11,"atol":2e-13,"max_step":.01}
STRICT = {"rtol":6e-13,"atol":6e-15,"max_step":.003}


def state_observables(state):
    return {"d":float(2*np.imag(state[0]*state[1].conj())),
            "F":float(abs(np.vdot(PSI_TARGET,state))**2),
            "ell":float(np.sum(abs(state[2:])**2)),
            "norm_error":float(abs(np.sum(abs(state)**2)-1)),
            "null_error":float(abs(np.sum(state*state)))}


def verify_incoming_state(incoming):
    source_case = {"E1":E,"eta":ETA_LOW,"tau":100.}
    state = PSI_INITIAL.copy()
    for stage in range(4):
        solution = solve_ivp(lambda v,z:wave_rhs(stage,v,z,source_case),
                             (0,1),state,method="DOP853",**STRICT)
        assert solution.success
        state = solution.y[:,-1]
    error = float(np.max(abs(state-incoming)))
    assert error < 2e-10
    return {"source_case":source_case,"source_time_over_tau":4.,
            "fixed_physical_input":True,"recomputed_input_error":error,
            "observables":state_observables(incoming),
            "real":incoming.real.tolist(),"imag":incoming.imag.tolist()}


class BoundaryStudy:
    def __init__(self, incoming, tau, tolerances=None):
        self.incoming = incoming.copy()
        self.tau = tau
        self.tolerances = dict(STANDARD if tolerances is None else tolerances)
        self.cache = {}

    def endpoint(self, eta, model):
        key = (float(eta),model)
        if key in self.cache:
            return self.cache[key]
        case = {"E1":E,"eta":float(eta),"tau":self.tau}
        scale = 2*self.tau
        if model == "full":
            solution = solve_ivp(lambda v,z:wave_rhs(4,v,z,case),
                                 (0,1),self.incoming,method="DOP853",
                                 **self.tolerances)
            assert solution.success
            result = state_observables(solution.y[:,-1])
            self.cache[key] = result
            return result

        Q0 = instantaneous_frame(0,E,eta)[0]
        initial = Q0.T@self.incoming
        if model == "sphere":
            slow = initial[2:]
            B = 2*np.imag(slow[:,None]*slow.conj()[None,:])
            b0 = np.array([B[1,2],B[2,0],B[0,1]])
            radius = float(np.linalg.norm(b0))
            slow_weight = float(np.sum(abs(slow)**2))

            def rhs(v,b):
                Q,DQ,_,mu = instantaneous_frame(v,E,eta)
                gamma = (Q.T@DQ)[3,4]
                return np.array([scale*mu*b[1],
                                 -scale*mu*b[0]-gamma*b[2],gamma*b[1]])

            solution = solve_ivp(rhs,(0,1),b0,method="DOP853",**self.tolerances)
            assert solution.success
            b = solution.y[:,-1]
            result = {"d":float(b[2]),
                      "F":float(slow_weight*(1-b[2]/radius)**2/4),
                      "sphere_norm_error":float(abs(np.linalg.norm(b)-radius)),
                      "fast_input_weight":float(np.sum(abs(initial[:2])**2)),
                      "slow_input_null_amplitude":float(abs(np.sum(slow*slow))),
                      "slow_input_weight":slow_weight,"sphere_radius":radius}
            self.cache[key] = result
            return result

        assert model in ("first_order","second_order")
        # Exact moving-frame generator = G0 + V. G0 preserves the two
        # separate blocks; V is only the omitted cross-block connection.
        def rhs(v,values):
            Q,DQ,muplus,muminus = instantaneous_frame(v,E,eta)
            connection = Q.T@DQ
            V = np.zeros((5,5))
            V[:2,2:] = -connection[:2,2:]
            V[2:,:2] = -connection[2:,:2]
            G0 = omega(muplus,muminus,scale)-connection-V
            z0,z1,z2 = values.reshape(3,5)
            return np.stack([G0@z0,G0@z1+V@z0,G0@z2+V@z1]).ravel()

        values0 = np.concatenate((initial,np.zeros(10,dtype=complex)))
        solution = solve_ivp(rhs,(0,1),values0,method="DOP853",**self.tolerances)
        assert solution.success
        responses = solution.y[:,-1].reshape(3,5)
        Qend = instantaneous_frame(1,E,eta)[0]
        for order,label in ((1,"first_order"),(2,"second_order")):
            state = Qend@np.sum(responses[:order+1],axis=0)
            result = state_observables(state)
            result["response_norms"] = np.linalg.norm(responses,axis=1).tolist()
            # No renormalization: report the truncated response's error.
            self.cache[(float(eta),label)] = result
        return self.cache[key]

    def root(self, model, bracket):
        eta = brentq(lambda x:self.endpoint(x,model)["d"],*bracket,
                     xtol=5e-14,rtol=2e-13)
        return {"eta":float(eta),"endpoint":self.endpoint(eta,model)}


def mixing_diagnostics(eta, tau):
    values, gaps = [],[]
    for v in np.linspace(0,1,501):
        Q,DQ,muplus,muminus = instantaneous_frame(v,E,eta)
        gap = 2*tau*(muplus-muminus)
        values.append(np.linalg.norm((Q.T@DQ)[:2,2:],2)/gap)
        gaps.append(gap)
    return {"max_fast_mixing_indicator":float(max(values)),
            "minimum_dimensionless_fast_gap":float(min(gaps))}


def root_comparison(study, bracket):
    # Compute predictions before the full-reference root. Nothing is fitted.
    sphere = study.root("sphere",bracket)
    corrected = study.root("second_order",bracket)
    first_order = study.root("first_order",bracket)
    full = study.root("full",bracket)
    eta_full = full["eta"]
    return {"tau":study.tau,"bracket":list(bracket),
            "sphere":sphere,"first_order":first_order,
            "second_order":corrected,"full":full,
            "sphere_absolute_eta_error":float(abs(sphere["eta"]-eta_full)),
            "sphere_relative_eta_error":float(abs(sphere["eta"]-eta_full)/eta_full),
            "first_order_absolute_eta_error":float(abs(first_order["eta"]-eta_full)),
            "second_order_absolute_eta_error":float(abs(corrected["eta"]-eta_full)),
            "full_at_sphere_prediction":study.endpoint(sphere["eta"],"full"),
            "full_at_corrected_prediction":study.endpoint(corrected["eta"],"full"),
            "mixing":mixing_diagnostics(eta_full,study.tau)}


def convergence_check(incoming, standard):
    study = BoundaryStudy(incoming,100.,STRICT)
    result = {model:study.root(model,(ETA_LOW,ETA_HIGH))
              for model in ("sphere","second_order","full")}
    shifts = {model:float(abs(result[model]["eta"]-standard[model]["eta"]))
              for model in result}
    assert max(shifts.values()) < 2e-10, shifts
    return {"standard_tolerances":STANDARD,"strict_tolerances":STRICT,
            "strict_roots":result,"eta_shifts":shifts}


def finite_slope(study, model, eta):
    step = 1e-6
    slope = (study.endpoint(eta+step,model)["d"]
             -study.endpoint(eta-step,model)["d"])/(2*step)
    half = step/2
    refined = (study.endpoint(eta+half,model)["d"]
               -study.endpoint(eta-half,model)["d"])/(2*half)
    return {"d_deta":float(refined),"step_convergence_difference":float(abs(refined-slope))}


def plot_comparison(grid, full_d, sphere_d, full_F, sphere_F, zoom, roots):
    fig,axes = plt.subplots(1,3,figsize=(14.4,4.1))
    full_color,sphere_color,corrected_color = "#222222","#0072B2","#D55E00"
    axes[0].plot(grid,full_d,color=full_color,lw=1.7,label="Full wavefunction")
    axes[0].plot(grid,sphere_d,color=sphere_color,lw=1.3,ls="--",label="Slow sphere")
    axes[0].axhline(0,color=".6",lw=.8)
    axes[0].set_ylabel(r"$d_{\rm out}$")
    axes[0].set_title("Fixed-input orientation")
    axes[1].plot(grid,full_F,color=full_color,lw=1.7)
    axes[1].plot(grid,sphere_F,color=sphere_color,lw=1.3,ls="--")
    axes[1].axhline(.25,color=".6",lw=.8)
    axes[1].set_ylabel(r"$F_{\rm out}$")
    axes[1].set_title("Squared fidelity")
    center = roots["full"]["eta"]
    zoom_x = (np.array(zoom["eta"])-center)*1e6
    axes[2].plot(zoom_x,zoom["full_d"],color=full_color,lw=1.7)
    axes[2].plot(zoom_x,zoom["sphere_d"],color=sphere_color,ls="--",lw=1.3)
    axes[2].plot(zoom_x,zoom["second_order_d"],color=corrected_color,
                 ls=":",lw=2,label="Fast/slow correction, order 2")
    axes[2].axhline(0,color=".6",lw=.8)
    axes[2].axvspan(0,(roots["sphere"]["eta"]-center)*1e6,
                   color=sphere_color,alpha=.1)
    axes[2].set_xlabel(r"$(\eta-\eta_c^{\rm full})\times10^6$ [meV]")
    axes[2].set_ylabel(r"$d_{\rm out}$")
    axes[2].set_title("Boundary shift")
    for ax in axes[:2]:
        ax.set_xlabel(r"$\eta=|t_1|$ amplitude [meV]")
    for ax in axes:
        ax.grid(alpha=.18)
    handles,labels = axes[0].get_legend_handles_labels()
    h2,l2 = axes[2].get_legend_handles_labels()
    fig.legend(handles+h2,labels+l2,loc="lower center",ncol=3,
               bbox_to_anchor=(.5,-.015),frameon=False)
    fig.suptitle(r"$E_1=0.01$ meV, $\tau=100$ meV$^{-1}$; same physical incoming state")
    fig.tight_layout(rect=(0,.075,1,.95))
    fig.savefig(ROOT/"switch_boundary.png",dpi=180,bbox_inches="tight")
    fig.savefig(ROOT/"switch_boundary.pdf",bbox_inches="tight")
    plt.close(fig)


def main():
    with np.load(ROOT/"bottleneck.npz") as wave_data:
        incoming = wave_data["eta_low_psi"][3200].copy()
    result = {"scope":"Same physical input; change only stage-2 eta or duration",
              "units":"hbar=1; energy in meV; tau in meV^-1",
              "input":verify_incoming_state(incoming)}
    study = BoundaryStudy(incoming,100.)
    main_roots = root_comparison(study,(ETA_LOW,ETA_HIGH))
    result["main_boundary"] = main_roots
    print("Main boundary roots:",
          {m:main_roots[m]["eta"] for m in ("sphere","first_order","second_order","full")},
          flush=True)
    result["convergence"] = convergence_check(incoming,main_roots)
    result["boundary_slopes"] = {
        model:finite_slope(study,model,main_roots[model]["eta"])
        for model in ("sphere","full")
    }

    grid = np.linspace(ETA_LOW,ETA_HIGH,41)
    full = [study.endpoint(eta,"full") for eta in grid]
    sphere = [study.endpoint(eta,"sphere") for eta in grid]
    full_d,sphere_d = np.array([r["d"] for r in full]),np.array([r["d"] for r in sphere])
    full_F,sphere_F = np.array([r["F"] for r in full]),np.array([r["F"] for r in sphere])
    result["main_scan"] = {
        "sample_count":len(grid),"eta":grid.tolist(),
        "full_d":full_d.tolist(),"sphere_d":sphere_d.tolist(),
        "full_F":full_F.tolist(),"sphere_F":sphere_F.tolist(),
        "max_d_error":float(np.max(abs(full_d-sphere_d))),
        "max_F_error":float(np.max(abs(full_F-sphere_F))),
        "full_sampled_decreasing":bool(np.all(np.diff(full_d)<0)),
        "sphere_sampled_decreasing":bool(np.all(np.diff(sphere_d)<0)),
    }
    # A coarse scan is not a proof that other roots cannot exist.
    width = max(8e-6,2*main_roots["sphere_absolute_eta_error"])
    zoom_eta = np.linspace(main_roots["full"]["eta"]-width,
                           main_roots["sphere"]["eta"]+width,19)
    result["boundary_zoom"] = {"eta":zoom_eta.tolist()}
    for model in ("full","sphere","second_order"):
        result["boundary_zoom"][f"{model}_d"] = [
            study.endpoint(eta,model)["d"] for eta in zoom_eta
        ]
    result["boundary_zoom"]["max_second_order_d_error"] = float(np.max(abs(
        np.array(result["boundary_zoom"]["second_order_d"])
        -np.array(result["boundary_zoom"]["full_d"])
    )))

    durations = []
    # Broaden the shorter segment's bracket because its crossing moves beyond
    # eta=.023; the physical input stays fixed even when the duration changes.
    for tau,bracket in ((75.,(.01,.035)),(150.,(.01,.023))):
        local = root_comparison(BoundaryStudy(incoming,tau),bracket)
        durations.append(local)
        print("Duration",tau,"roots:",{m:local[m]["eta"]
              for m in ("sphere","second_order","full")},flush=True)
    result["duration_comparisons"] = durations
    # The exact reference and corrected response are not fitted to each other.
    full_errors = max(max(r["norm_error"],r["null_error"]) for r in full)
    assert full_errors < 2e-9
    result["full_scan_max_constraint_error"] = full_errors
    (ROOT/"switch_boundary.json").write_text(
        json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8",
    )
    np.savez_compressed(ROOT/"switch_boundary.npz",eta=grid,
                        full_d=full_d,sphere_d=sphere_d,full_F=full_F,sphere_F=sphere_F)
    plot_comparison(grid,full_d,sphere_d,full_F,sphere_F,result["boundary_zoom"],main_roots)
    print("Root convergence shifts:",result["convergence"]["eta_shifts"])
    print("Main-scan max d/F errors:",result["main_scan"]["max_d_error"],
          result["main_scan"]["max_F_error"])
    print("Corrected boundary zoom max d error:",
          result["boundary_zoom"]["max_second_order_d_error"])


if __name__ == "__main__":
    main()
