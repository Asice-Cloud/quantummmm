#!/usr/bin/env python3
import numpy as np
import json
import matplotlib.pyplot as plt
from majorana_braid_check import build_bdg_two_site, extract_gap_and_zero_modes, linear_fit

# sweep (b,c) keeping a,d,mu0 fixed
a = 0.0
d = 0.0
mu0 = 0.0

bs = np.linspace(-0.5,0.5,21)
cs = np.linspace(-0.5,0.5,21)

results = {'bc_scan': [], 'ad_scan': []}

# sweep b,c diagonally (pairs) to test t vs b+c and Delta vs b-c
pairs = []
for b in bs:
    for c in cs:
        pairs.append((b,c))

t_meas = []
delta_meas = []
bc_plus = []
bc_minus = []
for b,c in pairs:
    H, params = build_bdg_two_site(a=a,b=b,c=c,d=d,mu0=mu0)
    # measured values: use the encoded mapping (since microscopic -> Kitaev mapping)
    t_meas.append(params['t'])
    delta_meas.append(params['Delta'])
    bc_plus.append(b + c)
    bc_minus.append(b - c)

# fits
m_t, c_t, r2_t = linear_fit(bc_plus, t_meas)
m_d, c_d, r2_d = linear_fit(bc_minus, delta_meas)
results['bc_scan'] = {'t_fit': {'slope':float(m_t),'intercept':float(c_t),'r2':float(r2_t)},
                      'Delta_fit': {'slope':float(m_d),'intercept':float(c_d),'r2':float(r2_d)}}

# plot results
plt.figure(figsize=(6,5))
plt.scatter(bc_plus, t_meas, s=8)
xs = np.linspace(min(bc_plus), max(bc_plus), 200)
plt.plot(xs, m_t*xs + c_t, 'r--')
plt.xlabel('b + c')
plt.ylabel('measured t')
plt.title('t vs (b+c)')
plt.tight_layout()
plt.savefig('fit_t_vs_bpc.png', dpi=200)
plt.close()

plt.figure(figsize=(6,5))
plt.scatter(bc_minus, delta_meas, s=8)
xs = np.linspace(min(bc_minus), max(bc_minus), 200)
plt.plot(xs, m_d*xs + c_d, 'r--')
plt.xlabel('b - c')
plt.ylabel('measured Delta')
plt.title('Delta vs (b-c)')
plt.tight_layout()
plt.savefig('fit_Delta_vs_bmc.png', dpi=200)
plt.close()

# sweep (a,d) to test mu vs 4*d + mu0 linearity
as_ = np.linspace(-0.2,0.2,21)
ds = np.linspace(-0.2,0.2,21)
ads = []
mu_meas = []
for a_ in as_:
    for d_ in ds:
        H, params = build_bdg_two_site(a=a_,b=0.1,c=0.05,d=d_,mu0=mu0)
        mu_meas.append(params['mu'])
        ads.append((a_, d_))
# fit mu vs d (we expect mu = 4*d + mu0)
d_vals = [ad[1] for ad in ads]
m_mu, c_mu, r2_mu = linear_fit(d_vals, mu_meas)
results['ad_scan'] = {'mu_fit': {'slope':float(m_mu),'intercept':float(c_mu),'r2':float(r2_mu)}}

# save results
with open('mapping_fit_results.json','w') as f:
    json.dump(results, f, indent=2)

print('Saved mapping_fit_results.json and PNGs')
