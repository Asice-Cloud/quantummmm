#!/usr/bin/env python3
import numpy as np
import json
import matplotlib.pyplot as plt
from build_micro_BdG import build_micro_BdG, extract_effective_params_from_micro
from majorana_braid_check import linear_fit

# parameter grids
bs = np.linspace(-0.5, 0.5, 21)
cs = np.linspace(-0.5, 0.5, 21)

b_plus_c = []
b_minus_c = []
t_meas = []
delta_meas = []

N = 8
for b in bs:
    for c in cs:
        H = build_micro_BdG(N, a=0.0, b=b, c=c, d=0.0, mu0=0.0)
        eff = extract_effective_params_from_micro(H)
        b_plus_c.append(b + c)
        b_minus_c.append(b - c)
        t_meas.append(eff['t_eff'])
        delta_meas.append(eff['Delta_eff'])

# fits
m_t, c_t, r2_t = linear_fit(b_plus_c, t_meas)
m_d, c_d, r2_d = linear_fit(b_minus_c, delta_meas)

# plot
xs = np.linspace(min(b_plus_c), max(b_plus_c), 200)
plt.figure(figsize=(6,4))
plt.scatter(b_plus_c, t_meas, s=8)
plt.plot(xs, m_t*xs + c_t, 'r--')
plt.xlabel('b + c')
plt.ylabel('t_eff')
plt.title('t_eff vs (b+c)')
plt.tight_layout()
plt.savefig('micro_fit_t_vs_bpc.png', dpi=200)
plt.close()

xs2 = np.linspace(min(b_minus_c), max(b_minus_c), 200)
plt.figure(figsize=(6,4))
plt.scatter(b_minus_c, delta_meas, s=8)
plt.plot(xs2, m_d*xs2 + c_d, 'r--')
plt.xlabel('b - c')
plt.ylabel('Delta_eff')
plt.title('Delta_eff vs (b-c)')
plt.tight_layout()
plt.savefig('micro_fit_Delta_vs_bmc.png', dpi=200)
plt.close()

# now sweep a,d to fit mu
as_ = np.linspace(-0.2,0.2,21)
ds = np.linspace(-0.2,0.2,21)
ads = []
mu_meas = []
for a in as_:
    for d in ds:
        H = build_micro_BdG(N, a=a, b=0.1, c=0.05, d=d, mu0=0.0)
        eff = extract_effective_params_from_micro(H)
        ads.append((a,d))
        mu_meas.append(eff['mu_eff'])

d_vals = [ad[1] for ad in ads]
m_mu, c_mu, r2_mu = linear_fit(d_vals, mu_meas)

results = {
    't_fit': {'slope':float(m_t), 'intercept':float(c_t), 'r2':float(r2_t)},
    'Delta_fit': {'slope':float(m_d), 'intercept':float(c_d), 'r2':float(r2_d)},
    'mu_fit': {'slope':float(m_mu), 'intercept':float(c_mu), 'r2':float(r2_mu)}
}

with open('mapping_from_micro.json','w') as f:
    json.dump(results, f, indent=2)

print('Saved mapping_from_micro.json and PNGs')
