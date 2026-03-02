#!/usr/bin/env python3
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def rodrigues(v, k, theta):
    k = k / np.linalg.norm(k)
    return v * np.cos(theta) + np.cross(k, v) * np.sin(theta) + k * np.dot(k, v) * (1 - np.cos(theta))


def plot_bloch_path(save_path='bloch_exchange_path.png'):
    # Sphere mesh
    u = np.linspace(0, 2 * np.pi, 200)
    v = np.linspace(0, np.pi, 100)
    x = np.outer(np.cos(u), np.sin(v))
    y = np.outer(np.sin(u), np.sin(v))
    z = np.outer(np.ones_like(u), np.cos(v))

    fig = plt.figure(figsize=(7,7))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(x, y, z, rstride=6, cstride=6, color='lightgray', alpha=0.14, linewidth=0)

    # Initial Bloch vector (north pole)
    v0 = np.array([0.0, 0.0, 1.0])

    # Define a sequence of effective doublet operators mapping to Pauli axes
    # e.g. iγ1γ2 -> Z, iγ2γ3 -> X, iγ3γ4 -> Y (axes are unit vectors)
    seq = [
        (np.array([0.0, 0.0, 1.0]), np.pi/2, 'iγ1γ2 → Z'),
        (np.array([1.0, 0.0, 0.0]), np.pi/2, 'iγ2γ3 → X'),
        (np.array([0.0, 1.0, 0.0]), np.pi/2, 'iγ3γ4 → Y'),
    ]

    points = [v0]
    labels = ['start']
    for axis, angle, name in seq:
        v_prev = points[-1]
        v_next = rodrigues(v_prev, axis, angle)
        points.append(v_next)
        labels.append(name)

    points = np.vstack(points)

    # Draw axes for reference
    def draw_axis(k, label):
        ax.quiver(0,0,0, k[0], k[1], k[2], length=1.0, color='gray', linewidth=1)
        ax.text(k[0]*1.05, k[1]*1.05, k[2]*1.05, label, color='gray')

    draw_axis(np.array([1.0,0.0,0.0]), 'X')
    draw_axis(np.array([0.0,1.0,0.0]), 'Y')
    draw_axis(np.array([0.0,0.0,1.0]), 'Z')

    # Plot points and connecting arcs
    colors = ['C0', 'C1', 'C2', 'C3']
    for i, p in enumerate(points):
        ax.quiver(0,0,0, p[0], p[1], p[2], length=1.0, color=colors[i%len(colors)], linewidth=2)
        ax.text(p[0]*1.08, p[1]*1.08, p[2]*1.08, labels[i], color=colors[i%len(colors)])

    # Draw arcs for each rotation
    for i, (axis, angle, name) in enumerate(seq):
        v_start = points[i]
        thetas = np.linspace(0, angle, 80)
        arc = np.array([rodrigues(v_start, axis, t) for t in thetas])
        ax.plot(arc[:,0], arc[:,1], arc[:,2], color='C3', linewidth=2, label=name if i==0 else '')

    # Draw full path connecting final back to start (optional visual)
    ax.plot(points[:,0], points[:,1], points[:,2], color='k', linestyle='--', linewidth=1)

    # Styling
    ax.set_box_aspect([1,1,1])
    ax.set_xlim([-1,1]); ax.set_ylim([-1,1]); ax.set_zlim([-1,1])
    ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
    ax.view_init(elev=20, azim=30)
    ax.set_title('Bloch sphere: sequence of Majorana-mediated rotations')

    # Legend and save
    ax.legend(loc='upper left')
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    print(f"Saved Bloch path figure to {save_path}")


if __name__ == '__main__':
    plot_bloch_path('bloch_exchange_path.png')
