# Imports
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from sympy.printing.latex import latex
from IPython.display import display, Math

# Symbolic variables
R, l, g = sp.symbols('R l g', positive=True)
psi, phi, theta, psi_dot, phi_dot, theta_dot = sp.symbols('psi phi theta psi_dot phi_dot theta_dot', real=True)

# Positions and velocitys
h = -l * sp.sin(psi)
r_s = sp.Matrix([l, 0, 0])
Omega = sp.Matrix([- sp.sin(psi) * phi_dot, psi_dot, sp.cos(psi) * phi_dot])
omega = sp.Matrix([theta_dot, 0, 0]) + Omega
Traegheitsmoment_s = 0.25 * R**2 * sp.Matrix([[2, 0, 0], [0, 1, 0], [0, 0, 1]])
v_s = Omega.cross(r_s)

# Energies
E_kin = 0.5 * v_s.dot(v_s) + 0.5 * omega.dot(Traegheitsmoment_s * omega)
E_pot = g * h
Lagrangian = E_kin - E_pot
Lagrangian = sp.simplify(Lagrangian)

# Partial Derivatives
L_psi = sp.diff(Lagrangian, psi_dot)
L_phi = sp.diff(Lagrangian, phi_dot)
L_theta = sp.diff(Lagrangian, theta_dot)

# Express generalized velocities and Lagrangian in terms of generalized momenta
Pi_theta, Pi_phi, Pi_psi = sp.symbols('Pi_theta Pi_phi Pi_psi', real=True)

dot_sub_eqns = [sp.Eq(Pi_psi, L_psi), sp.Eq(Pi_phi, L_phi), sp.Eq(Pi_theta, L_theta)]
sol_dot_sub = sp.solve(dot_sub_eqns, [psi_dot, phi_dot, theta_dot], dict=True)[0]
[psi_dot_sub, phi_dot_sub, theta_dot_sub] = [sol_dot_sub[psi_dot], sol_dot_sub[phi_dot], sol_dot_sub[theta_dot]]
Lagrangian_sub = Lagrangian.subs({psi_dot: psi_dot_sub, phi_dot: phi_dot_sub, theta_dot: theta_dot_sub})

# Hamiltonian
H = psi_dot_sub*Pi_psi + phi_dot_sub*Pi_phi + theta_dot_sub*Pi_theta - Lagrangian_sub
H = sp.simplify(H)
H_latex = latex(H)

# Define a function to compute Poisson brackets for scalar functions A and B
def poisson_bracket(A, B, q, p):
    return sum(sp.diff(A, q_i)*sp.diff(B, p_i) - sp.diff(A, p_i)*sp.diff(B, q_i) for q_i, p_i in zip(q, p))

# Define a function to compute Poisson brackets for a vector function A
def poisson_bracket_vector(A, B, q, p):
    return sp.Matrix([poisson_bracket(A_i, B, q, p) for A_i in A])

# Compute equations of motion using Poisson brackets
q = [psi, phi, theta]
p = [Pi_psi, Pi_phi, Pi_theta]
qp = [psi, phi, theta, Pi_psi, Pi_phi, Pi_theta]

H_eqn = poisson_bracket_vector(qp , H, q, p)
H_eqn = sp.simplify(H_eqn)

# Parameters
R0 = 0.125 # m
l0 = R0 # m
g0 = 9.81 # m/s^2
psi_0 = 0 # rad
phi_0 = 0 # rad
theta_0 = 0 # rad
psi_dot_0 = 0 # rad/s
phi_dot_0 = 0 # rad/s
theta_dot_0 = 10 * 2 * np.pi # rad/s
T = 1 # total time in seconds
t_span = (0, T)

param_subs = {R: R0, l: l0, g: g0}

# Initial phase space coordinates
[L_psi_sub, L_phi_sub, L_theta_sub] = [L_psi.subs(param_subs), L_phi.subs(param_subs), L_theta.subs(param_subs)]
[L_psi_num, L_phi_num, L_theta_num] = [sp.lambdify((psi, psi_dot, phi_dot, theta_dot), L_psi_sub, 'numpy'),
                                          sp.lambdify((psi, psi_dot, phi_dot, theta_dot), L_phi_sub, 'numpy'),
                                          sp.lambdify((psi, psi_dot, phi_dot, theta_dot), L_theta_sub, 'numpy')]
[Pi_psi_0, Pi_phi_0, Pi_theta_0] = [L_psi_num(psi_0, psi_dot_0, phi_dot_0, theta_dot_0),
                                    L_phi_num(psi_0, psi_dot_0, phi_dot_0, theta_dot_0),
                                    L_theta_num(psi_0, psi_dot_0, phi_dot_0, theta_dot_0)]
qp_0 = np.array([psi_0, phi_0, theta_0, Pi_psi_0, Pi_phi_0, Pi_theta_0], dtype = float)

# Numerical solution of Hamiltonian equations
H_eqn_sub = H_eqn.subs(param_subs)
H_num = sp.lambdify(qp, H_eqn_sub, 'numpy')
def H_ode(t, y):
    return np.squeeze(H_num(y[0], y[1], y[2], y[3], y[4], y[5]))
sol_H = solve_ivp(H_ode, t_span, qp_0, t_eval = np.linspace(0, T, 500))

# Extract solutions
t = sol_H.t
psi_sol = sol_H.y[0]
phi_sol= sol_H.y[1]
theta_sol = sol_H.y[2]
Pi_psi_sol = sol_H.y[3]
Pi_phi_sol = sol_H.y[4]
Pi_theta_sol = sol_H.y[5]

# Plot results

const = (np.pi/2) * np.ones_like(t)

# Psi plot
plt.figure(figsize=(7,4))
plt.plot(t, psi_sol, linewidth=2)
plt.plot(t, Pi_psi_sol, linewidth=2)
plt.plot(t, const, 'r')
plt.xlabel('t [s]')
plt.ylabel(r'$\psi(t)\ [\mathrm{rad}],\  \Pi_{\psi}(t)\ [\mathrm{rad/s}]$')
plt.title(r'Nutation')
plt.legend([r'$\psi(t)$', r'$\Pi_{\psi}(t)$', r'$\psi = \frac{\pi}{2}$'], loc='best')
plt.grid(True)

# Phi plot
plt.figure(figsize=(7,4))
plt.plot(t, phi_sol, linewidth=2)
plt.plot(t, Pi_phi_sol, linewidth=2)
plt.plot(t, 2*const, 'r')
plt.xlabel('t [s]')
plt.ylabel(r'$\phi(t)\ [\mathrm{rad}],\  \Pi_{\phi}(t)\ [\mathrm{rad/s}]$')
plt.title(r'Precession')
plt.legend([r'$\phi(t)$', r'$\Pi_{\phi}(t)$', r'$\phi = \pi$'], loc='best')
plt.grid(True)

# Theta plot
plt.figure(figsize=(7,4))
plt.plot(t, theta_sol, linewidth=2)
plt.plot(t, Pi_theta_sol, linewidth=2)
plt.plot(t, 4*const, 'r')
plt.xlabel('t [s]')
plt.ylabel(r'$\theta(t)\ [\mathrm{rad}],\  \Pi_{\theta}(t)\ [\mathrm{rad/s}]$')
plt.title(r'Spin')
plt.legend([r'$\theta(t)$', r'$\Pi_{\theta}(t)$', r'$\theta = 2\pi$'], loc='best')
plt.grid(True)

plt.show()

# Phase space portraits

# Psi plot
plt.figure(figsize=(7,4))
plt.plot(psi_sol, Pi_psi_sol, linewidth=2)
plt.xlabel(r'$\psi,\ [\mathrm{rad}]$')
plt.ylabel(r'$\dot{\psi},\ [\mathrm{rad/s}]$')
plt.title(r'Phase of $\psi$')
plt.grid(True)

# Phi plot
plt.figure(figsize=(7,4))
plt.plot(phi_sol, Pi_phi_sol, linewidth=2)
plt.xlabel(r'$\phi,\ [\mathrm{rad}]$')
plt.ylabel(r'$\dot{\phi},\ [\mathrm{[rad/s}]$')
plt.title(r'Phase of $\phi$')
plt.grid(True)

# Theta plot
plt.figure(figsize=(7,4))
plt.plot(theta_sol, Pi_theta_sol, linewidth=2)
plt.xlabel(r'$\theta,\ [\mathrm{rad}]$')
plt.ylabel(r'$\dot{\theta},\ [\mathrm{rad/s}]$')
plt.title(r'Phase of $\theta$')
plt.grid(True)

plt.show()