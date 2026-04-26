"""
First Order Control Barrier Function.
The system has the form
\dot{x} = u.
The barrier function is defined as the norm between the system position
and the obstacle position minus a raduis of 0.5.
\norm(p_system - p_obstacle) - radius.


This is based on the results shown in:
Comparative Analysis of Control Barrier Functions and 
Artificial Potential Fields for Obstacle Avoidance
 http://ames.caltech.edu/singletary2020comparative.pdf

"""

import numpy as np
import matplotlib.pyplot as plt

# System setup
I = np.eye(2)
A = np.zeros((2, 2))
B = I
K = 1

dt = 0.001
t_final = 20
s = np.arange(0, t_final + dt, dt)

# States
xp = np.zeros((2, len(s)))
xp[:, 0] = np.array([0.0, 0.0])

e = np.zeros((2, len(s)))
r = np.zeros((2, len(s)))

# Simulation loop
for i in range(len(s) - 1):
    # Reference
    r[:, i] = np.array([3.0, 5.0])
    e[:, i] = xp[:, i] - r[:, i]

    # Barrier parameters
    alpha1 = 5
    gamma = np.eye(2)
    psi = np.zeros(2)

    pos = xp[:, i]

    # Obstacles
    obs1 = np.array([1.0, 2.0])
    obs2 = np.array([2.5, 3.0])

    def compute_barrier(pos, obs):
        diff = pos - obs
        norm_diff = np.linalg.norm(diff)
        h = norm_diff - 0.5                 # Barrier function
        b = diff.T @ gamma / norm_diff
        phi = diff.T @ psi / norm_diff + alpha1 * h
        return h, b, phi

    hx1, b1, phi1 = compute_barrier(pos, obs1)
    hx2, b2, phi2 = compute_barrier(pos, obs2)

    # Select closest (most restrictive) barrier
    if hx1 < hx2:
        h, b, phi = hx1, b1, phi1
    else:
        h, b, phi = hx2, b2, phi2

    # Nominal controller
    u_nominal = -K * e[:, i]

    # Barrier condition
    cond = phi + b @ u_nominal

    if cond < 0:
        u_safe = -b.T * cond / (np.linalg.norm(b) ** 2)
    else:
        u_safe = np.zeros(2)

    # Final control
    u = u_nominal + u_safe

    # Euler integration
    xp[:, i + 1] = xp[:, i] + dt * (A @ xp[:, i] + B @ u)

# --- Plot Errors ---
plt.figure()
plt.plot(s, e[0, :], label='e1')
plt.plot(s, e[1, :], label='e2')
plt.xlabel("Time")
plt.ylabel("Errors")
plt.legend()
plt.grid()

# --- Trajectory Plot ---
plt.figure()

# Draw obstacles
theta = np.linspace(0, 2*np.pi, 100)
plt.fill(1 + 0.5*np.sin(theta), 2 + 0.5*np.cos(theta), 'k')
plt.fill(2.5 + 0.5*np.sin(theta), 3 + 0.5*np.cos(theta), 'k')

# Goal
plt.plot(3, 5, 'ro')

# Trajectory
plt.plot(xp[0, :], xp[1, :], 'b')

plt.xlim([-1, 5])
plt.ylim([-1, 6])
plt.gca().set_aspect('equal', adjustable='box')

plt.xlabel("x")
plt.ylabel("y")
plt.grid()

plt.show()