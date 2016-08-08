#!/usr/bin/env python

# TODO : An easy way to get the A and B matrices for all closed loops would be
# to compute the A and B for the five loop closures and then set the outer loop
# gains to zero to get the various inner loops.

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
b['y'] = Matrix([delta, phi_dot])  # in order of the loop closures

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
b['C'] = Matrix([[1 if s == o else 0 for s in b['x']] for o in b['y']])
b['D'] = zeros(num_outputs, num_inputs)

# Define the neuromuscular state space.
omega = Symbol('\omega_{nm}', real=True, positive=True)
zeta = Symbol('\zeta_{nm}', real=True, positive=True)
n = {}
n['A'] = Matrix([[0, 1],
                 [-omega**2, -2 * zeta * omega]])
n['B'] = Matrix([0, omega**2])
n['C'] = Matrix([1, 0]).T
n['D'] = Matrix([0])

n['x'] = Matrix([[t_delta],
                 [t_delta_dot]])
n['y'] = Matrix([t_delta])

# First loop controller.
k_delta, delta_c = symbols('k_\delta, delta_c', real=True)
Unm = k_delta * (delta_c - delta)
Unm = Unm.expand()
controller = Matrix([Unm.coeff(var) for var in b['y'].col_join(Matrix([delta_c]))]).T

# Compute the differential equations for the closed loop system.
x_dot_bicycle = b['A'] * b['x'] + b['B'] * b['u']
x_dot_neuro = n['A'] * n['x'] + n['B'] * controller * b['y'].col_join(Matrix([delta_c]))
x_dot = x_dot_bicycle.col_join(x_dot_neuro)

# Build the closed loop state space matrices.
c_delta = {}
c_delta['x'] = Matrix([phi, delta, phi_dot, delta_dot, t_delta, t_delta_dot])
c_delta['u'] = Matrix([t_phi, f_b, delta_c])
c_delta['y'] = Matrix([phi, delta, phi_dot, t_delta])


def mat_coeff(equations, variables, i, j):
    c = equations[i].expand().coeff(variables[j])
    if c is None:
        c = 0
    return c

c_delta['A'] = Matrix(len(c_delta['x']), len(c_delta['x']),
                      lambda i, j: mat_coeff(x_dot, c_delta['x'], i, j))
c_delta['B'] = Matrix(len(c_delta['x']), len(c_delta['u']),
                      lambda i, j: mat_coeff(x_dot, c_delta['u'], i, j))
c_delta['C'] = Matrix(len(c_delta['y']), len(c_delta['x']),
                      lambda i, j: mat_coeff(c_delta['y'], c_delta['x'], i, j))

# Write the controller output as a function of the gains and the commanded
# lateral deviation.
k_delta, k_phi_dot = symbols('k_\delta, k_{\dot{\phi}}', real=True)
phi_dot_c = symbols('\dot{\phi}_c', real=True)

delta_c = k_phi_dot * (phi_dot_c - phi_dot)
Unm = k_delta * (delta_c - delta)

Unm = Unm.expand()
controller = Matrix([Unm.coeff(var)
                     for var in b['y'].col_join(Matrix([phi_dot_c]))]).T

# Compute the differential equations for the closed loop system.
x_dot_bicycle = b['A'] * b['x'] + b['B'] * b['u']
x_dot_neuro = n['A'] * n['x'] + n['B'] * controller * b['y'].col_join(Matrix([phi_dot_c]))
x_dot = x_dot_bicycle.col_join(x_dot_neuro)

# Build the closed loop state space matrices.
c_phidot = {}
c_phidot['x'] = Matrix([phi, delta, phi_dot, delta_dot, t_delta, t_delta_dot])
c_phidot['u'] = Matrix([t_phi, f_b, phi_dot_c])
c_phidot['y'] = Matrix([phi, delta, phi_dot, t_delta])
c_phidot['A'] = Matrix(len(c_phidot['x']), len(c_phidot['x']),
                       lambda i, j: mat_coeff(x_dot, c_phidot['x'], i, j))
c_phidot['B'] = Matrix(len(c_phidot['x']), len(c_phidot['u']),
                       lambda i, j: mat_coeff(x_dot, c_phidot['u'], i, j))
c_phidot['C'] = Matrix(len(c_phidot['y']), len(c_phidot['x']),
                       lambda i, j: mat_coeff(c_phidot['y'], c_phidot['x'], i, j))
