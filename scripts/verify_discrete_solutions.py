#!/usr/bin/env python3
import numpy as np
from math import isclose

def U(theta_x, theta_y, theta_z):
    return np.exp(1j*(theta_x+theta_y)) - np.exp(1j*(theta_x+theta_z)) + np.exp(1j*(theta_y+theta_z)) - 1

def main():
    S = [0.0, np.pi/2, np.pi, 3*np.pi/2]
    sols = []
    for tx in S:
        for ty in S:
            for tz in S:
                val = U(tx,ty,tz)
                if isclose(abs(val), 0.0, rel_tol=0, abs_tol=1e-12):
                    sols.append((tx,ty,tz))
    fn = 'scripts/discrete_theta_solutions.txt'
    with open(fn, 'w') as f:
        f.write('Solutions (theta_x,theta_y,theta_z) in {0,pi/2,pi,3pi/2} that make U=0:\n')
        for s in sols:
            f.write(str(s) + '\n')
    print('Wrote', fn)
    print('Found', len(sols), 'solutions')

if __name__ == "__main__":
    main()
