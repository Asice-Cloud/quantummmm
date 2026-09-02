import numpy as np, sys
sys.path.insert(0,'.work')
from dyn import step_prop
tc=0.3; E0=0.3
I2=np.eye(2); sx=np.array([[0,1],[1,0]],complex); sy=np.array([[0,-1j],[1j,0]],complex); sz=np.diag([1.,-1.]).astype(complex)

# ---- FULL JW engine (physical observable) ----
def full_prop(tau,E1,t1c,n=1500):
    U=np.eye(4,dtype=complex)
    for ph in (0,1,2,0,1,2): U=step_prop(tau,E1,t1c,tc,E0,ph,n=n)@U
    return U
def full_W(tau,E1,t1c,n=1500):
    U=full_prop(tau,E1,t1c,n)
    return 0.5*sum(abs(U[1*2+a,0*2+b])**2 for a in range(2) for b in range(2))

# ---- Reduced code-space model built from model.md protocol ----
# Paper's reduced form: U(6tau)=T exp{-i ∫ [E1,eff(t)(g1g2)+t1,eff(t)(g3g1)+ vtheta(t)(g2g3)]}
# In code space (gamma1,gamma2,gamma3 Pauli embedding chosen consistent w/ paper: 
#  sigx=-i g2g3, sigy=-i g3g1, sigz=-i g1g2  =>  g1g2 = i sigz, g3g1 = i sigy, g2g3 = i sigx
#  -i∫[...] with g1g2=i sigz -> E1 eff term gives ... but paper wrote exp(-i dt ϑ g2g3)=exp(π/2 σx)
#  since g2g3 = i sigx : -i ϑ (i sigx)= ϑ sigx. So reduced H_R = E1,eff sigz + t1,eff sigy + (ϑdot) sigx? no...
# Let's just follow paper: reduced Hamiltonian (Hermitian) driving U as exp(-i ∫ H dt):
#    H_red(t) = -E1,eff(t)*sigz  ... careful
# Instead derive from H_EM's code-sector generators mapped to code Pauli:
#   E1 term  E1*(i g1g2):   i g1g2 = -i*(g1g2)?? no. compute in JW: (i g1g2) restricted... 
# We'll fit instead: assume H_red(t)= a(t)*E1*sigz + b(t)*t1c*sigy with a(t),b(t) from protocol windows,
# geometric NOT applied as explicit σx rotation at end of each double-braid? 
# SIMPLER: The geometric braid = σx NOT (ϑ=π/2) is τ-independent. Reduced dynamics over windows:
#   step1&2: E1 on (full), t1 ramps; step3: E1 on, t1 off
#   Over full 6tau: phase from E1 ~ E1*6tau*(1/2 window avg) sigz-ish etc.
# CRUDE but parameterized: U_red = NOT * exp(-i (θE sigz + θt sigy)) with θE=kE*E1*tau_phys, θt=kt*t1*tau_phys
# Compare envelopes vs full engine for ratios.
def U_red(tau,E1,t1c,beta):
    # beta=[kE,kt] free scale (window fractions), NOT geometric
    thE=beta[0]*E1*tau; tht=beta[1]*t1c*tau
    Urot=np.linalg.matrix_power(np.eye(2,dtype=complex),0)
    # apply exp(-i θE sz) then exp(-i θt sy) (noncommuting, sequential approx of T-ordering)
    A=(1j)*0  # placeholder
    from scipy.linalg import expm
    U=expm(-1j*thE*sz)@expm(-1j*tht*sy)
    # geometric NOT on code: sigx with ϑ=π/2 -> exp(-i (π/2) sigx)? paper: NOT=exp(π/2 σx)=iσx
    NOT=1j*sx  # exp(pi/2 sigx)=i sigx (up to our sign convention); its |0>->|1> effect same
    # W measured between |psi1-> (code 1) and |psi1+ (code 0): U_total maps 0->... 
    # In MZM limit U->NOT then W=1. So final code amp from |0>: |out>=U_rot NOT|0>? or NOT U_rot? order matters.
    # paper: braiding and dynamic simultaneous: U(6τ)=Texp{-i∫[E1 g1g2 + t1 g3g1 + ϑ g2g3]} all three terms
    #   => NOT and E1/t1 phases ARE the whole integral, not separate. So:
    v0=np.array([1,0],complex)
    # full reduced propagator with piecewise: approximate T-order as U=Urot then geometric?
    Ufull=U  # ignore NOT first
    return Ufull

# empirical calibration: find kE such that pure-E1 (t1=0) envelope period matches full engine
print("CALIBRATE pure E1, t1=0: full engine W(tau) oscillations")
taus=np.linspace(200,1200,51)
E1=0.01
for tag,E1 in [("E1=.001",0.001),("E1=.01",0.01)]:
    Ws=[full_W(t,E1,0.0) for t in taus]
    # count peaks
    peaks=sum(1 for i in range(1,len(Ws)-1) if Ws[i]>Ws[i-1] and Ws[i]>Ws[i+1])
    print(tag,"E1*dtau_range=",E1*(taus[-1]-taus[0])," peaks_in_window:",peaks, " W range",round(min(Ws),3),round(max(Ws),3))
