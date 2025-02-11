import math
import pandas as pd
import numpy as np

############### Initial Conditions ######################
g = 9.81  # gravity constant
h = 0.01  # smoothing length
Umax = 1
c = 10*Umax  # Speed of sound in water
# c = 1500  # Speed of sound in water
gama=7
rho_0 = 1
P_0 = rho_0 * c*c/gama
C_cfl = 1.15
dt = C_cfl * h / c  # delta time from CFL condition
height = 2
width = 2
k=75
no_of_layers = 2
kai = X =500000
n = 3600
# n = math.ceil(height * width / math.pi / h**2)  # Number of particles [2D]
rho_init = 1  # constant rho
mu = 0.0010016  # at temp 20
nu = mu / rho_init  # kinematic viscosity

V_in_lam = 0
# Umax = V_in_lam = 2300 * mu / (height * rho)  # laminar flow velocity
v0 = np.zeros((k,3))
v0[:,0][0] = V_in_lam # v0 =  np.array([V_in_lam, 0, 0]) # initial velocity
fg = -g * np.array([0, 1, 0])  # gravity force vector
total_mass_of_the_system = rho_init * height * width
mj = total_mass_of_the_system / n  # mass of each particle
t=0
current_step =0

eta = 0.01 * h
length_L = 2
velocity_of_lid_U = 1
rho_0 = 1 
volume = math.pow(h,2)
no_of_particle = math.pow(length_L,2) / math.pow(h,2)
nu = 0.1 
mu = nu * rho_0
v_at_t_0 = p_at_t_0 = 0
va = va_wall = np.zeros((k,3))
va_wall[:,0][0] =  velocity_of_lid_U 
delta_t = (C_cfl * h) / c
number_of_steps = 1 / delta_t
v_prime = []  # results of step 1
vij = v0.copy()  # initial velocity
vi = v0.copy()  # initial velocity
vij = np.zeros((n, 3))  # For storing velocity differences
# Initial coordinates at time = 0
delta_y = (height - 4 * h) / (n * 2) # the y distance between the particles 
dx = delta_y
coordinates = np.zeros((n, 3))  # initialize particles' positions in 3D space
for i in range(n // 2):  # create the initial coordinates
    if i == 0:
        coordinates[i] = np.array([i*dx, delta_y * i, 0])
    else:
        coordinates[i] = np.array([i*dx, delta_y * i, 0])
for i in range(n // 2, n):  # create the initial coordinates
    coordinates[i] = np.array([0, -delta_y * (i - n // 2), 0])

print(coordinates)
# print(coordinates)
# ############### End of Initial Conditions ################
