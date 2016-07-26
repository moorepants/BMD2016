#!/usr/bin/env python

from sympy import Symbol, symbols, Matrix, zeros, eye, init_printing

init_printing()

# build the bicycle state space
delta, phi = symbols('delta, phi', real=True)
delta_dot, phi_dot = symbols('\dot{\delta}, \dot{\phi}', real=True)
t_phi, t_delta, f_b = symbols('T_\phi, T_\delta, F_B', real=True)
t_delta_dot = symbols('\dot{T}_\delta', real=True)

b = {}
b['x'] = Matrix([phi, delta, phi_dot, delta_dot])
b['u'] = Matrix([t_phi, t_delta, f_b])
b['y'] = Matrix([delta, phi, phi_dot])

num_states = len(b['x'])
num_inputs = len(b['u'])
num_outputs = len(b['y'])

# Make the A and B matrices generic as they are just essentially submatrices in
# the closed loop system. In this case the C only returns a subset of the
# states.

b['A'] = zeros(num_states / 2).row_join(
    eye(num_states / 2)).col_join(
        Matrix(num_states / 2, num_states, lambda i, j:
               Symbol('a_{}{}'.format(i + 2, j))))
b['B'] = Matrix(num_states, num_inputs,
                lambda i, j: Symbol('b_{}{}'.format(i, j)) if i > 1 else 0)
b['C'] = eye(num_outputs).row_join(zeros(num_outputs, 1))
b['D'] = Matrix(num_outputs, num_inputs, lambda i, j: 0)


# Write the controller output as a function of the gains and the commanded
# lateral deviation.
k_delta, k_phi_dot = symbols('k_\delta, k_{\dot{\phi}}', real=True)
phi_dot_c = symbols('\dot{\phi}_c', real=True)

delta_c = k_phi_dot * (phi_dot_c - phi_dot)
Unm = k_delta * (delta_c - delta)

Unm = Unm.expand()
controller = Matrix([Unm.coeff(phi),
                     Unm.coeff(delta),
                     Unm.coeff(phi_dot),
                     Unm.coeff(phi_dot_c)]).T

# Define the neuromuscular state space.
omega = Symbol('omega_nm', real=True, positive=True)
zeta = Symbol('zeta_nm', real=True, positive=True)
n = {}
n['A'] = Matrix([[0, 1],
                 [-omega**2, -2 * zeta * omega]])
n['B'] = Matrix([[0],
                 [omega**2]])
n['C'] = Matrix([1, 0]).T
n['D'] = Matrix([0])

n['x'] = Matrix([[t_delta],
                 [t_delta_dot]])
n['u'] = Matrix([Unm])
n['y'] = Matrix([t_delta])

# Compute the differential equations for the closed loop system.
x_dot = b['A'] * b['x'] + b['B'] * b['u']
temp_1 = n['A'] * n['x']
temp_2 = (b['C'] * b['x']).col_join(Matrix([phi_dot_c]))
temp_3 = n['B'] * controller * temp_2
x_dot = x_dot.col_join(temp_1 + temp_3)

# Build the closed loop state space matrices.
c = {}
c['x'] = Matrix([phi, delta, phi_dot, delta_dot, t_delta, t_delta_dot])
c['u'] = Matrix([t_phi, f_b, phi_dot_c])
c['y'] = Matrix([phi, delta, phi_dot, t_delta])


def mat_coeff(equations, variables, i, j):
    c = equations[i].expand().coeff(variables[j])
    if c is None:
        c = 0
    return c

c['A'] = Matrix(len(c['x']), len(c['x']),
                lambda i, j: mat_coeff(x_dot, c['x'], i, j))
c['B'] = Matrix(len(c['x']), len(c['u']),
                lambda i, j: mat_coeff(x_dot, c['u'], i, j))
c['C'] = Matrix(len(c['y']), len(c['x']),
                lambda i, j: mat_coeff(c['y'], c['x'], i, j))
