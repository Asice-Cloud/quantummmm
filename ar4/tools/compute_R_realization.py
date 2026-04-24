import numpy as np
from scipy.linalg import logm, expm

# Pauli matrices
I = np.array([[1,0],[0,1]], dtype=complex)
X = np.array([[0,1],[1,0]], dtype=complex)
Y = np.array([[0,-1j],[1j,0]], dtype=complex)
Z = np.array([[1,0],[0,-1]], dtype=complex)
PAULIS = {'0':I, 'x':X, 'y':Y, 'z':Z}

def kron(a,b):
    return np.kron(a,b)

def decompose_to_pauli(H):
    coeffs = {}
    for mu, A in PAULIS.items():
        for nu, B in PAULIS.items():
            P = kron(A,B)
            # coefficient in basis: c = Tr(P^ H)/4 (since Tr(P^ P)=4)
            c = np.trace(P.conj().T.dot(H)) / 4.0
            coeffs[f'c_{mu}{nu}'] = c
    return coeffs

# mapping to Kitaev params (from notes)
def map_to_kitaev(c):
    c_xx = c['c_xx']
    c_yy = c['c_yy']
    c_xy = c['c_xy']
    c_yx = c['c_yx']
    c_zz = c['c_zz']
    c_z0 = c['c_z0']
    c_0z = c['c_0z']
    c_00 = c['c_00']
    t = (c_xx + c_yy)
    Delta = (c_xx - c_yy)
    U = 4.0 * c_zz
    mu_site = 4.0*c_zz - 2.0*(c_z0 + c_0z)
    Ebond = c_zz - (c_z0 + c_0z) + c_00
    return {'t':t, 'Delta':Delta, 'U':U, 'mu':mu_site, 'Ebond':Ebond}

def suggest_pulse(coeffs):
    # Look for dominant XX+YY combination
    c_xx = coeffs['c_xx']
    c_yy = coeffs['c_yy']
    score = abs(c_xx) + abs(c_yy)
    suggestion = {}
    if score.real > 1e-8:
        # treat H = J*(XX+YY) -> in fermion picture this gives Majorana bilinear
        J = (c_xx + c_yy)/2.0
        # choose tau to realize angle pi/4: J * tau = pi/4
        if abs(J) < 1e-12:
            suggestion['note'] = 'XX+YY coefficient small; no simple XX+YY pulse.'
        else:
            tau = (np.pi/4.0) / (J.real if abs(J.real)>1e-12 else J)
            suggestion['pulse'] = {
                'term': 'XX+YY',
                'J (suggested)': float(np.real(J)),
                'tau (such that J*tau=pi/4)': complex(tau)
            }
    else:
        suggestion['note'] = 'No dominant XX/YY; consider decomposing H into local generators and Trotterizing.'
    return suggestion

if __name__ == '__main__':
    # Example: construct a test R from a simple local H
    H_test = 0.5 * (kron(X,X) + kron(Y,Y))
    # build a unitary R corresponding to exp(i * angle * H_test)
    angle = np.pi/4.0
    R = expm(1j * angle * H_test)

    print('Test R constructed as exp(i * angle * H_test) with angle=', angle)

    # compute Hermitian generator H_eff via log
    L = logm(R)
    H_eff = (-1j) * L
    # enforce Hermitian
    H_eff = 0.5 * (H_eff + H_eff.conj().T)

    coeffs = decompose_to_pauli(H_eff)

    print('\nDecomposed coefficients (Pauli basis):')
    for k,v in coeffs.items():
        # only print if magnitude significant
        if abs(v) > 1e-8:
            print(f'{k}: {v}')

    params = map_to_kitaev(coeffs)
    print('\nMapped Kitaev-like parameters:')
    for k,v in params.items():
        print(f'{k}: {v}')

    suggestion = suggest_pulse(coeffs)
    print('\nPulse suggestion:')
    print(suggestion)

    # Save numeric output to file for inspection
    import json
    out = {'coeffs': {k: (float(np.real(v)), float(np.imag(v))) for k,v in coeffs.items()}, 'params': {k: (float(np.real(v)), float(np.imag(v))) for k,v in params.items()}, 'suggestion': suggestion}
    with open('tools/R_decomposition_output.json','w') as f:
def json_serial(obj):
        if isinstance(obj, (complex, np.complex64, np.complex128)):
            return [float(obj.real), float(obj.imag)]
        raise TypeError("Type %s not serializable" % type(obj))
json.dump(out, f, indent=2, default=json_serial)
print('\nWrote tools/R_decomposition_output.json')
