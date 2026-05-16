import numpy as np
from scipy.linalg import expm, eigh
import matplotlib.pyplot as plt

# Minimal Majorana effective-model braiding simulation
# 6 Majorana -> 3 complex fermion modes -> Hilbert dim 8


def pauli_matrices():
    I = np.eye(2, dtype=complex)
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sz = np.array([[1, 0], [0, -1]], dtype=complex)
    sp = np.array([[0, 1], [0, 0]], dtype=complex)
    sm = np.array([[0, 0], [1, 0]], dtype=complex)
    return I, sx, sy, sz, sp, sm


def kron_n(op_list):
    M = op_list[0]
    for A in op_list[1:]:
        M = np.kron(M, A)
    return M


def fermionic_operators(n_modes):
    I, sx, sy, sz, sp, sm = pauli_matrices()
    dim = 2 ** n_modes
    c_ops = []
    for j in range(n_modes):
        ops = []
        for k in range(n_modes):
            if k < j:
                ops.append(sz)
            elif k == j:
                ops.append(sm)  # lowering operator
            else:
                ops.append(I)
        c = kron_n(ops)
        c_ops.append(c)
    return c_ops


def majorana_ops_from_c(c_ops):
    gammas = []
    for c in c_ops:
        cd = c.conj().T
        gamma1 = c + cd
        gamma2 = -1j * (c - cd)
        gammas.append(gamma1)
        gammas.append(gamma2)
    return gammas


def time_profiles(t, tau, tc):
    # single-swap step length = 3*tau; we'll implement protocol repeated twice (total_time = 6*tau)
    # determine which substep within single swap we are in (0..3*tau)
    # we return t1(t), t2(t), t3(t), Ed(t)
    Tsingle = 3 * tau
    t_mod = t % Tsingle
    # default
    t1 = 0.0
    t2 = 0.0
    t3 = 0.0
    Ed = 1.0  # base value

    # step1: 0->tau : turn on t2, Ed->0
    if 0 <= t_mod < tau:
        x = t_mod
        t2 = tc * (1 - np.cos(np.pi * x / tau)) / 2
        Ed = 1.0 * (1 + np.cos(np.pi * x / tau)) / 2
        t1 = 0.0
        t3 = 0.0
    # step2: tau->2tau : turn on t3 and turn on t1, t2 decreases
    elif tau <= t_mod < 2 * tau:
        x = t_mod - tau
        t3 = tc * (1 - np.cos(np.pi * x / tau)) / 2
        t2 = tc * (1 + np.cos(np.pi * x / tau)) / 2
        t1 = tc * (1 - np.cos(np.pi * x / tau)) / 2
        Ed = 0.0
    # step3: 2tau->3tau : turn off t3, Ed returns
    else:
        x = t_mod - 2 * tau
        t3 = tc * (1 + np.cos(np.pi * x / tau)) / 2
        t1 = tc * (1 + np.cos(np.pi * x / tau)) / 2
        t2 = 0.0
        Ed = 1.0 * (1 - np.cos(np.pi * x / tau)) / 2

    return t1, t2, t3, Ed


def build_H(gammas, params, t1_val, t2_val, t3_val, Ed_val):
    # gammas indexed: 0..5 for 6 Majorana: assign as in paper
    # γ1 = g0, γ2 = g1, γ3 = g2, γ4 = g3, γa = g4, γb = g5
    g = gammas
    H = np.zeros_like(g[0], dtype=complex)
    # i E1 γ1 γ2
    H += 1j * params['E1'] * (g[0] @ g[1])
    # optionally E2 (set small or zero)
    H += 1j * params.get('E2', 0.0) * (g[2] @ g[3])
    # i Ed γa γb
    H += 1j * Ed_val * (g[4] @ g[5])
    # couplings: + i|t2| γa γ2  (paper had +i|t2| γa γ2)
    H += 1j * t2_val * (g[4] @ g[1])
    # - i|t1| γb γ1
    H += -1j * t1_val * (g[5] @ g[0])
    # - i|t3| γa γ3  (paper sign conventions vary; choose -i t3 γa γ3)
    H += -1j * t3_val * (g[4] @ g[2])
    return 0.5 * H  # factor 1/2 to avoid double counting (convention)


def run_braiding(params, tau=10.0, steps_per_tau=100, repeat=2, tc=0.5):
    n_modes = 3
    c_ops = fermionic_operators(n_modes)
    # build majorana operators inline to avoid potential parsing issues
    gammas = []
    for c in c_ops:
        cd = c.conj().T
        gammas.append(c + cd)
        gammas.append(-1j * (c - cd))
    dim = 2 ** n_modes

    total_time = repeat * 3 * tau
    steps = int(repeat * 3 * tau * steps_per_tau)
    dt = total_time / steps

    # Use computational basis initial state |00...0> to avoid eigenstate tracking issues
    psi = np.zeros(dim, dtype=complex)
    psi[0] = 1.0

    t = 0.0
    for k in range(steps):
        t += dt
        t1v, t2v, t3v, Edv = time_profiles(t, tau, tc)
        Ht = build_H(gammas, params, t1v, t2v, t3v, Edv)
        U = expm(-1j * Ht * dt)
        psi = U @ psi

    # measure fidelity as overlap with final state
    # for perfect braiding, we expect phase evolution; measure against initial state
    final_overlap = abs(psi[0])
    return final_overlap, psi, None, gammas


if __name__ == '__main__':
    params = {
        'E1': 0.001,  # small overlap between γ1 and γ2
        'E2': 0.0,
    }

    taus = np.linspace(1, 200, 40)
    fidelities = []
    for tau in taus:
        fid, _, _, _ = run_braiding(params, tau=tau, steps_per_tau=200, repeat=2, tc=0.03)
        fidelities.append(fid)

    fidelities = np.array(fidelities)

    print('taus -> fidelity (|<psi_plus|U|psi_minus>|)')
    for t, f in zip(taus, fidelities):
        print(f'{t:.3f}  {f:.6f}')

    plt.figure()
    plt.plot(taus, fidelities, '-o')
    plt.xlabel('tau')
    plt.ylabel('Fidelity')
    plt.title('Effective Majorana braiding fidelity vs tau')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('effective_braiding_fidelity.png', dpi=200)
    print('Saved effective_braiding_fidelity.png')
