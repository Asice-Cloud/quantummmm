#!/usr/bin/env python3
"""Plot the verified complete bivector PDE using the paper's axes and palette.

E1=.01 in the parameter map; tau, tc and E0 remain controlled design choices.
All evolution uses K=2A, the original six cosine segments, and existing solver
settings. Sampling density controls the figures, not the integrator settings.
"""

import argparse
import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path

from analyze_peak_structure import (
    ATOL, MAX_STEP, PAIRS, RTOL, controls, independent_wave, propagate,
    readout, state_action,
)
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import numpy as np
from scipy.integrate import solve_ivp


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT/'state_evolution'
E1, TC, E0 = .01, .3, .3
COMPONENTS = [r'$d=B_{12}$', r'$x_1=B_{13}$', r'$x_2=B_{1a}$',
              r'$x_3=B_{1b}$', r'$y_1=B_{23}$', r'$y_2=B_{2a}$',
              r'$y_3=B_{2b}$', r'$h_1=B_{3a}$', r'$h_2=B_{3b}$', r'$h_3=B_{ab}$']
TARGET_LABEL = r'$|\langle\psi_1^+(0)|\psi_1^-(t)\rangle|$'
INITIAL_LABEL = r'$|\langle\psi_1^-(0)|\psi_1^-(t)\rangle|$'
PARAMETERS = r'$E_1=0.01\ \mathrm{meV},\quad t_c=E_0=0.3\ \mathrm{meV},\quad K=2A$'


def palette():
    root = ET.parse(ROOT/'paper_alignment/paper_contours.svg').getroot()
    rgb = []
    for node in root.iter():
        fill = node.get('fill', '')
        if fill.startswith('rgb('):
            rgb.append([float(x)/100 for x in re.findall(r'[\d.]+', fill)])
    assert len(rgb) == 10
    return ListedColormap(rgb, name='paper_fig1d')


def save(fig, stem):
    fig.savefig(OUTPUT/f'{stem}.png', dpi=190)
    fig.savefig(OUTPUT/f'{stem}.pdf')
    plt.close(fig)


def invariants(b):
    B = np.zeros(b.shape[:-1]+(5,5))
    B[..., PAIRS[0], PAIRS[1]] = b
    B -= np.swapaxes(B, -1, -2)
    return {
        'max_abs_B_cubed_plus_B': float(np.max(abs(B@B@B+B))),
        'max_abs_component_norm_squared_minus_one': float(np.max(abs(np.sum(b*b, axis=-1)-1))),
        'F_min': float(np.min(readout(b)['F'])),
        'F_max': float(np.max(readout(b)['F'])),
    }


def time_history(tau=100.):
    # The ideal reference matches Fig. 1(b); the other columns span Fig. 1(d).
    energies = np.array([0., E1, E1, E1])
    etas = np.array([0., .001, .01, .1])
    b = np.zeros((4,10))
    b[:,0] = 1.
    histories, times = [], []
    local = np.linspace(0.,1.,501)
    for stage in range(6):
        def rhs(v, flat):
            p, r, ed, gate = controls(stage, v, TC, E0)
            return (2*tau*state_action(flat.reshape(4,10), energies, etas,
                                       p, r, ed, gate)).ravel()
        solution = solve_ivp(rhs, (0.,1.), b.ravel(), method='DOP853',
                             rtol=RTOL, atol=ATOL, max_step=MAX_STEP, t_eval=local)
        if not solution.success:
            raise RuntimeError(solution.message)
        samples = solution.y.T.reshape(-1,4,10)
        start = int(stage > 0)
        histories.append(samples[start:])
        times.append(stage+local[start:])
        b = samples[-1]
    values = np.concatenate(histories).transpose(1,0,2)
    F = readout(values)['F']
    survival = F+values[...,0]
    return {'zeta': np.concatenate(times), 'tau': np.array(tau), 'E1': energies,
            'eta': etas, 'b': values, 'F': F, 'survival_F': survival}


def plot_history(data, cmap):
    signed_cmap = LinearSegmentedColormap.from_list('paper_signed',cmap.colors,N=256)
    titles = [r'Ideal reference: $E_1=|t_1|=0$',
              r'$E_1=0.01,\quad |t_1|/E_1=0.1$',
              r'$E_1=0.01,\quad |t_1|/E_1=1$',
              r'$E_1=0.01,\quad |t_1|/E_1=10$']
    fig, axes = plt.subplots(3,4,figsize=(16.2,9.5),layout='constrained',
                             gridspec_kw={'height_ratios':[1.6,1,1]},sharex=True)
    for col in range(4):
        top, squared, amplitude = axes[:,col]
        mesh = top.pcolormesh(data['zeta'], np.arange(10), data['b'][col].T,
                              cmap=signed_cmap, vmin=-1, vmax=1, shading='nearest', rasterized=True)
        top.set(yticks=np.arange(10), yticklabels=COMPONENTS if col == 0 else [],
                ylim=(9.5,-.5), title=titles[col])
        for ax, exponent in [(squared,1.), (amplitude,.5)]:
            target = np.maximum(data['F'][col],0)**exponent
            survival = np.maximum(data['survival_F'][col],0)**exponent
            power = r'$^2$' if exponent == 1. else ''
            ax.plot(data['zeta'], target, color='blue', lw=1.6, label=TARGET_LABEL+power)
            ax.plot(data['zeta'], survival, color='magenta', lw=1.3, label=INITIAL_LABEL+power)
            ax.set(ylim=(-.035,1.055),yticks=[0,.5,1.])
            ax.grid(axis='y',alpha=.14)
        for ax in axes[:,col]:
            ax.set(xlim=(0,6),xticks=np.arange(7))
            for stage in range(1,6):
                ax.axvline(stage,color='0.35',lw=.65,ls='--',alpha=.65)
            if col > 0 and ax is not top:
                ax.tick_params(labelleft=False)
        amplitude.set_xlabel(r'step ($t/\tau$)')
    axes[0,0].set_ylabel(r'Full state $\mathcal{B}$')
    axes[1,0].set_ylabel(r'Squared fidelity $F$')
    axes[2,0].set_ylabel(r'Overlap modulus $f=\sqrt{F}$')
    axes[1,0].legend(loc='upper left',fontsize=7.2,framealpha=.82)
    axes[2,0].legend(loc='upper left',fontsize=7.2,framealpha=.82)
    fig.colorbar(mesh,ax=axes[0,:],shrink=.86,pad=.012,ticks=[-1,-.5,0,.5,1],label=r'$B_{ij}$')
    fig.suptitle(r'Complete state and fidelity during two exchanges: $\tau=100/\mathrm{meV}$, $t_c=E_0=0.3$ meV, $K=2A$',fontsize=13)
    save(fig,'time_evolution')

    fig, axes = plt.subplots(2,2,figsize=(11.,8.),layout='constrained',sharex=True)
    for row, ax in enumerate(axes.flat):
        mesh=ax.pcolormesh(data['zeta'],np.arange(10),data['b'][row].T,cmap=signed_cmap,
                           vmin=-1,vmax=1,shading='nearest',rasterized=True)
        for stage in range(1,6):
            ax.axvline(stage,color='0.3',lw=.65,ls='--')
        ax.set(xlim=(0,6),ylim=(9.5,-.5),xticks=np.arange(7),yticks=np.arange(10),
               yticklabels=COMPONENTS,title=titles[row],xlabel=r'step ($t/\tau$)')
    fig.colorbar(mesh,ax=axes.ravel().tolist(),pad=.01,shrink=.85,
                  ticks=[-1,-.5,0,.5,1],label=r'$B_{ij}$')
    fig.suptitle(r'Full state $\mathcal{B}(t)$: $\tau=100/\mathrm{meV}$, $t_c=E_0=0.3$ meV, $K=2A$')
    save(fig,'B_time_evolution')

    # A separate, readable F-only figure with the exact probabilities labelled.
    fig, axes = plt.subplots(2,2,figsize=(10.5,6.7),layout='constrained',sharex=True,sharey=True)
    for row, ax in enumerate(axes.flat):
        ax.plot(data['zeta'],data['F'][row],color='blue',lw=1.7,label=TARGET_LABEL+r'$^2$')
        ax.plot(data['zeta'],data['survival_F'][row],color='magenta',lw=1.4,label=INITIAL_LABEL+r'$^2$')
        for stage in range(1,6):
            ax.axvline(stage,color='0.5',lw=.6,ls='--')
        ax.set(xlim=(0,6),ylim=(-.035,1.055),xticks=np.arange(7),yticks=[0,.5,1],
               title=titles[row],xlabel=r'step ($t/\tau$)',ylabel=r'$F$')
        ax.grid(axis='y',alpha=.15)
    axes[0,0].legend(fontsize=9)
    fig.suptitle(r'Squared fidelity: $\tau=100/\mathrm{meV}$, $t_c=E_0=0.3$ meV, $K=2A$')
    save(fig,'fidelity_evolution')


def endpoint_map():
    tau_axis, log_axis = np.linspace(0.,1200.,241), np.linspace(-1.,1.,101)
    b = np.empty((len(log_axis),len(tau_axis),10))
    # Split by duration so short protocols do not inherit a long-protocol step.
    for start in range(0,len(tau_axis),20):
        stop = min(start+20,len(tau_axis))
        tau, logarithm = np.meshgrid(tau_axis[start:stop],log_axis)
        b[:,start:stop] = propagate(tau,E1,E1*10**logarithm,response=False)['b']
        print(f'Endpoint map: {stop}/{len(tau_axis)} duration columns',flush=True)
    return {'tau':tau_axis,'log_ratio':log_axis,'b':b,**readout(b)}


def map_axes(ax):
    ax.set(xlim=(0,12),ylim=(-1,1),xticks=np.arange(0,13,2),yticks=[-1,-.5,0,.5,1],
           xlabel=r'$\tau\ (100/\mathrm{meV})$',ylabel=r'$\lg(|t_1|/E_1)$')


def map_panel(ax, fig, data, values, cmap, limits, label, title):
    # A band centered on zero prevents tiny signed numerical residuals from
    # appearing as strong alternating colors in B components.
    levels = np.linspace(*limits,12 if limits[0] < 0 else 11)
    mesh = ax.contourf(data['tau']/100,data['log_ratio'],np.clip(values,*limits),
                       levels=levels,cmap=cmap)
    map_axes(ax)
    ax.set_title(title)
    fig.colorbar(mesh,ax=ax,pad=.018,ticks=np.linspace(*limits,6),label=label)


def plot_maps(data,cmap):
    signed_cmap = LinearSegmentedColormap.from_list('paper_signed',cmap.colors,N=256)
    fig, axes = plt.subplots(1,3,figsize=(15.6,4.9),layout='constrained')
    map_panel(axes[0],fig,data,data['d'],signed_cmap,(-1,1),r'$d=B_{12}(6\tau)$',
              r'(a) State orientation $d$')
    map_panel(axes[1],fig,data,data['F'],cmap,(0,1),r'$F(6\tau)$',
              r'(b) Squared fidelity $F$')
    map_panel(axes[2],fig,data,data['f'],cmap,(0,1),
              r'$|\langle\psi_1^+(0)|\psi_1^-(6\tau)\rangle|$',
              r'(c) Paper readout $\sqrt{F}$')
    fig.suptitle(PARAMETERS,fontsize=14)
    save(fig,'endpoint_maps')

    fig, axes = plt.subplots(2,5,figsize=(17.,7.),layout='constrained')
    levels=np.linspace(-1,1,12)
    for index, ax in enumerate(axes.flat):
        mesh=ax.contourf(data['tau']/100,data['log_ratio'],np.clip(data['b'][...,index],-1,1),
                         levels=levels,cmap=signed_cmap)
        map_axes(ax)
        ax.set_title(COMPONENTS[index])
        if index%5:
            ax.set_ylabel('')
            ax.tick_params(labelleft=False)
    fig.colorbar(mesh,ax=axes.ravel().tolist(),pad=.008,shrink=.83,ticks=[-1,-.5,0,.5,1],label=r'$B_{ij}(6\tau)$')
    fig.suptitle(r'All ten components of $\mathcal{B}(6\tau)$; '+PARAMETERS,fontsize=13)
    save(fig,'B_endpoint_components')


def check_results(history, endpoints):
    report={'history_invariants':invariants(history['b']),
            'endpoint_invariants':invariants(endpoints['b'])}
    wave_errors=[]
    for row in range(4):
        wave=independent_wave(float(history['tau']),history['E1'][row],history['eta'][row])
        wave_errors.append({'E1':float(history['E1'][row]),'eta':float(history['eta'][row]),
                            'B_max_error':float(np.max(abs(wave['b']-history['b'][row,-1]))),
                            'F_error':float(abs(wave['F']-history['F'][row,-1]))})
    report['trajectory_endpoint_wave_checks']=wave_errors
    residual=history['F']-(1-readout(history['b'])['ell']-history['b'][...,0])/2
    report['geometric_F_identity_residual']=float(np.max(abs(residual)))
    errors=[]
    for stage in [2,5]:
        mask=(history['zeta']>=stage)&(history['zeta']<=stage+1)
        errors.append(float(np.max(np.ptp(history['F'][:,mask],axis=1))))
    report['third_segment_F_variations']=errors
    # Compare selected map nodes to the existing independent wave solver.
    report['map_wave_checks']=[]
    for row,col in [(0,20),(50,60),(100,200)]:
        wave=independent_wave(endpoints['tau'][col],E1,E1*10**endpoints['log_ratio'][row])
        report['map_wave_checks'].append({
            'tau':float(endpoints['tau'][col]),'log_ratio':float(endpoints['log_ratio'][row]),
            'B_max_error':float(np.max(abs(wave['b']-endpoints['b'][row,col]))),
            'F_error':float(abs(wave['F']-endpoints['F'][row,col]))})
    for key in ['history_invariants','endpoint_invariants']:
        assert report[key]['max_abs_B_cubed_plus_B']<2e-8,report[key]
        assert report[key]['max_abs_component_norm_squared_minus_one']<2e-8,report[key]
        assert -2e-8<=report[key]['F_min']<=report[key]['F_max']<=1+2e-8,report[key]
    for row in wave_errors+report['map_wave_checks']:
        assert row['B_max_error']<2e-8 and row['F_error']<2e-8,row
    assert max(errors)<2e-8
    return report


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--replot',action='store_true',help='Render the saved arrays without reintegrating.')
    args=parser.parse_args()
    OUTPUT.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size':10,'axes.titlesize':11,'mathtext.fontset':'stix',
                         'axes.linewidth':.8,'xtick.direction':'in','ytick.direction':'in',
                         'pdf.fonttype':42})
    cmap=palette()
    if args.replot:
        with np.load(OUTPUT/'time_history.npz') as stored:
            history={key:stored[key] for key in stored.files}
        with np.load(OUTPUT/'endpoint_map.npz') as stored:
            endpoints={key:stored[key] for key in stored.files}
    else:
        history=time_history()
        np.savez_compressed(OUTPUT/'time_history.npz',**history)
        plot_history(history,cmap)
        print('Saved full time trajectories and fidelity figures.',flush=True)
        endpoints=endpoint_map()
        np.savez_compressed(OUTPUT/'endpoint_map.npz',**endpoints)
        checks=check_results(history,endpoints)
        report={'E1_map_meV':E1,'tc_meV':TC,'E0_meV':E0,'generator':'K=2A',
                'physical_end_time':'6*tau','map_grid':[101,241],
                'map_tau_range':[0,1200],'map_log10_ratio_range':[-1,1],
                'trajectory_tau':100.,'trajectory_points':len(history['zeta']),
                'solver':{'method':'DOP853','rtol':RTOL,'atol':ATOL,'max_local_step':MAX_STEP},
                'palette':'Ten fill colors extracted from the original Fig. 1(d) vector paths.',
                'color_readout':'F is squared fidelity; sqrt(F) has the same definition as the paper colorbar. Signed B components use [-1,1].',
                'scope':'Finite sampling of the complete state PDE; tiny rapid structures may require a finer plotting grid.',
                'checks':checks}
        (OUTPUT/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(checks,indent=2),flush=True)
    plot_history(history,cmap)
    plot_maps(endpoints,cmap)
    print('Saved all evolution and parameter-map figures.',flush=True)


if __name__=='__main__':
    main()
