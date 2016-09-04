#!/usr/bin/env python

# TODO : An easy way to get the A and B matrices for all closed loops would be
# to compute the A and B for the five loop closures and then set the outer loop
# gains to zero to get the various inner loops.

from sympy import Symbol, symbols, Matrix, zeros, eye, init_printing

init_printing()


def mat_coeff(equations, variables, i, j):
    c = equations[i].expand().coeff(variables[j])
    if c is None:
        c = 0
    return c

#####################
# Bicycle state space
#####################
delta, phi = symbols('delta, phi', real=True)
delta_dot, phi_dot = symbols('\dot{\delta}, \dot{\phi}', real=True)
t_phi, t_delta, f_b = symbols('T_\phi, T_\delta, F_B', real=True)
t_delta_dot = symbols('\dot{T}_\delta', real=True)

b = {}
b['x'] = Matrix([phi, delta, phi_dot, delta_dot])
b['u'] = Matrix([t_phi, t_delta, f_b])
b['y'] = Matrix([delta, phi_dot, phi])  # in order of the loop closures

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

###########################
# Neuromuscular state space
###########################
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

#######################
# Steer loop controller
#######################
k_delta, delta_c = symbols('k_\delta, delta_c', real=True)
Unm = k_delta * (delta_c - delta)
Unm = Unm.expand()
new_y = b['y'].col_join(Matrix([delta_c]))
controller = Matrix([Unm.coeff(var) for var in new_y]).T

# Compute the differential equations for the closed loop system.
x_dot_bicycle = b['A'] * b['x'] + b['B'] * b['u']
x_dot_neuro = n['A'] * n['x'] + n['B'] * controller * new_y
x_dot = x_dot_bicycle.col_join(x_dot_neuro)

# Build the closed loop state space matrices.
c_delta = {}
c_delta['x'] = Matrix([phi, delta, phi_dot, delta_dot, t_delta, t_delta_dot])
c_delta['u'] = Matrix([t_phi, f_b, delta_c])
c_delta['y'] = Matrix([phi, delta, phi_dot, t_delta])
c_delta['A'] = Matrix(len(c_delta['x']), len(c_delta['x']),
                      lambda i, j: mat_coeff(x_dot, c_delta['x'], i, j))
c_delta['B'] = Matrix(len(c_delta['x']), len(c_delta['u']),
                      lambda i, j: mat_coeff(x_dot, c_delta['u'], i, j))
c_delta['C'] = Matrix(len(c_delta['y']), len(c_delta['x']),
                      lambda i, j: mat_coeff(c_delta['y'], c_delta['x'], i, j))

###########################
# Roll rate loop controller
###########################
k_phi_dot, phi_dot_c = symbols('k_{\dot{\phi}}, \dot{\phi}_c', real=True)

delta_c = k_phi_dot * (phi_dot_c - phi_dot)
Unm = k_delta * (delta_c - delta)
new_y = b['y'].col_join(Matrix([phi_dot_c]))
Unm = Unm.expand()
controller = Matrix([Unm.coeff(var) for var in new_y]).T

# Compute the differential equations for the closed loop system.
x_dot_bicycle = b['A'] * b['x'] + b['B'] * b['u']
x_dot_neuro = n['A'] * n['x'] + n['B'] * controller * new_y
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

######################
# Roll loop controller
######################
k_phi, phi_c = symbols('k_{\phi}, \phi_c', real=True)

phi_dot_c = k_phi * (phi_c - phi)
delta_c = k_phi_dot * (phi_dot_c - phi_dot)
Unm = k_delta * (delta_c - delta)
new_y = b['y'].col_join(Matrix([phi_c]))
Unm = Unm.expand()
Unm_phi = Unm
controller = Matrix([Unm.coeff(var) for var in new_y]).T

# Compute the differential equations for the closed loop system.
x_dot_bicycle = b['A'] * b['x'] + b['B'] * b['u']
x_dot_neuro = n['A'] * n['x'] + n['B'] * controller * new_y
x_dot = x_dot_bicycle.col_join(x_dot_neuro)

# Build the closed loop state space matrices.
c_phi = {}
c_phi['x'] = Matrix([phi, delta, phi_dot, delta_dot, t_delta, t_delta_dot])
c_phi['u'] = Matrix([t_phi, f_b, phi_c])
c_phi['y'] = Matrix([phi, delta, phi_dot, t_delta])
c_phi['A'] = Matrix(len(c_phi['x']), len(c_phi['x']),
                    lambda i, j: mat_coeff(x_dot, c_phi['x'], i, j))
c_phi['B'] = Matrix(len(c_phi['x']), len(c_phi['u']),
                    lambda i, j: mat_coeff(x_dot, c_phi['u'], i, j))
c_phi['C'] = Matrix(len(c_phi['y']), len(c_phi['x']),
                    lambda i, j: mat_coeff(c_phi['y'], c_phi['x'], i, j))

#########################
# Handling Quality Metric
#########################

U_M_expr = -(delta / k_phi_dot + phi_dot) / k_phi
U_M = symbols('U_M')

# Build the closed loop state space matrices.
c_hqm = {}
c_hqm['x'] = c_phi['x']
c_hqm['u'] = Matrix([phi_c])
c_hqm['y'] = Matrix([U_M])
c_hqm['A'] = c_phi['A']
c_hqm['B'] = Matrix(len(c_hqm['x']), len(c_hqm['u']),
                    lambda i, j: mat_coeff(x_dot, c_hqm['u'], i, j))
c_hqm['C'] = Matrix(len(c_hqm['y']), len(c_hqm['x']),
                    lambda i, j: mat_coeff(Matrix([U_M_expr]), c_hqm['x'], i, j))

s = symbols('s')
hqm = c_hqm['C'] * (s * eye(c_hqm['A'].shape[0]) - c_hqm['A']).inv() * c_hqm['B']

# TODO : Need to define A an B matrix for bicycle with the extra states: psi,
# yq.
#####################
# Yaw loop controller
#####################
###k_psi, psi_c = symbols('k_\psi, \psi_c', real=True)
###
###phi_c = k_psi * (psi_c - psi)
###phidot_c = k_phi * (phi_c - phi)
###delta_c = k_phi_dot * (phi_dot_c - phi_dot)
###Unm = k_delta * (delta_c - delta)
###
###new_y = b['y'].col_join(Matrix([psi_c]))
###
###Unm = Unm.expand()
###controller = Matrix([Unm.coeff(var) for var in new_y]).T
###
#### Compute the differential equations for the closed loop system.
###x_dot_bicycle = b['A'] * b['x'] + b['B'] * b['u']
###x_dot_neuro = n['A'] * n['x'] + n['B'] * controller * new_y
###x_dot = x_dot_bicycle.col_join(x_dot_neuro)
###
#### Build the closed loop state space matrices.
###c_psi = {}
###c_psi['x'] = Matrix([phi, delta, phi_dot, delta_dot, t_delta, t_delta_dot])
###c_psi['u'] = Matrix([t_phi, f_b, phi_c])
###c_psi['y'] = Matrix([phi, delta, phi_dot, t_delta])
###c_psi['A'] = Matrix(len(c_psi['x']), len(c_psi['x']),
                    ###lambda i, j: mat_coeff(x_dot, c_psi['x'], i, j))
###c_psi['B'] = Matrix(len(c_psi['x']), len(c_psi['u']),
                    ###lambda i, j: mat_coeff(x_dot, c_psi['u'], i, j))
###c_psi['C'] = Matrix(len(c_psi['y']), len(c_psi['x']),
                    ###lambda i, j: mat_coeff(c_psi['y'], c_psi['x'], i, j))
