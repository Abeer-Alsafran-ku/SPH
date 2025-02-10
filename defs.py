from ic import *
import math
import pandas as pd
import numpy as np
from sklearn.neighbors import NearestNeighbors
# ################# DEFs ###################
def distance_ij(point, coordinates, k):
    # Calculate pairwise distances using NearestNeighbors
    nbrs = NearestNeighbors(n_neighbors=k, algorithm='auto').fit(coordinates)
    distances, indices = nbrs.kneighbors(point)
    return distances, indices

def calc_rj(r, k):
    # Calculates the neighbors' positions for each particle
    rj = np.zeros((k,3))
    for i in range(k):
        rj,_ = distance_ij(np.array([[0.0, r[:,1][i], 0.0]]),r, k)
    return rj

def calc_rij(r,rj ,k):
    rij = np.zeros(k) 
    for i in range(k):
        rij[i] = r[:,1][i] - rj[:,1][i]
    return rij

def get_neighbours_indx(point ,coordinates,k):
    # Initialize the NearestNeighbors model
    nbrs = NearestNeighbors(n_neighbors=k, algorithm='auto').fit(coordinates) # it takes 2D array as an argument
    _, indices = nbrs.kneighbors(point) #indices 
    return indices

def calc_vj(v,indx): ## here need to specify the formula for calculating vj
    vj = np.zeros((k,3))
    # call indicies (distance) -> []
    # the vj is the velocity of the neighbours particles 
    # step 1: get the indecies of the neighbours particles 
    # current_point = v # to be seted 
    indinces = get_neighbours_indx(np.array([coordinates[indx]]),coordinates, 75) # sending the current point (particle) 
    # step 2: 
    for i in range(len(indinces[0])):
        vj[:,0][i] = indinces[0][i]#v[indinces[i]] #vj[] = v[indicies]
    return vj

def calc_vij(vi, vj):
    vij = np.zeros((k,3))
    for i in range(k):
        vij[:,0][i] = vi[:,0][i] - vj[:,0][i]
    return vij

def calc_del_W(rij, h):
    # Calculate the gradient of the smoothing function
    norm_rij = np.linalg.norm(rij)
    s = norm_rij / h
    if s > 0:
        del_W = 15 / (7 * math.pi * h**4) * (-9 / 4 + 19 / 8 * s - 5 / 8 * s**2) * rij
    else:
        del_W = np.zeros_like(rij)
    return del_W

def calc_viscous_i(rij, vij, h):
    del_Wij = calc_del_W(rij, h)  # Make sure rij is a vector
    term_r_delW = np.dot(rij, del_Wij)  # dot product of rij and del_Wij
    rij_n = np.linalg.norm(rij)  # distance between particles
    term_viscous_i = mj * mu * vij * term_r_delW / (rho**2) / (rij_n**2 + eta**2)  # viscous force term
    return term_viscous_i

def calc_viscous_f(rij, vij, h):
    nu_del2vi = np.zeros_like(rij)  # Shape (n, 3), assuming rij is (n, 3)
    for i in range(len(rij[0])):
        nu_del2vi[i] = calc_viscous_i(rij[i], vij[i], h)  # Each rij[i] and vij[i] should be vectors
    return nu_del2vi

def step_1_i(v_i_n, rij, vij, h, fg):
    fv = calc_viscous_f(rij, vij, h)  # Make sure rij and vij are arrays of vectors (shape (n, 3))
    F = fg + fv  # F is the total force
    v_i_prime = v_i_n + F * delta_t  # update velocity using force
    return v_i_prime

def calc_W(r,h): # this function will calculate W: smoothing particle for a specific r 
    w = np.zeros((k,3))
    for i in range(k):
        r = np.append(r,i)
        s = np.linalg.norm(r[i]) / h # = q
        if s >=0 and s <=1:
            w[:,1][i] = 15 / 7 / math.pi / math.pow(h,2) * (2/3 - 9/8 * math.pow(s,2) + 19/ 24 * math.pow(s,3) - 5/32 * math.pow(s,4))
        else: # s is not in the range
            w[:,1][i] = 0
    return w

def calc_rho_star(mj,w):
    rho_star = np.zeros((k,3))
    for i in range(len(w)):
        rho_star[:,1][i] = mj * w[:,1][i]
    return rho_star 

def calc_rho_ij(rho_star,indx): #to be calc # rho_star == rho_i(mj * w(ri)) summation 
    rho_i = rho_star.copy() #  
    rho_j = np.zeros((k,3))
    indc = get_neighbours_indx(np.array([[0.0, rho_i[:,1][indx], 0.0]]),rho_i,k)
    for i in range (len(indc[0])):
        rho_j[i] = rho_i[indc[0][i]]
    return rho_i,rho_j  

def calc_init_p(rho_0,c,gama):
    init_p = rho_0*c / math.pow(gama,2)
    return init_p

def step_2_i(r, v_prime,delta_t):
    ri_0 = np.copy(r) # to be calculated
    ri = [] # final result of step 2
    for i in range(k):
        r_tmp = ri_0[i] + (v_prime[i] * delta_t) # r is a vector
        ri.append(r_tmp)
    return ri

def step_3_i(init_p,rho,rho_0,gama,X):
    p = init_p *(math.pow(rho/rho_0,gama)-1) + X
    return p

def calc_delta_p_rho(mj,rho_i,rho_j,p,rij,del_w ,eta):
    p_rho = np.zeros((k,3))
    for i in range(k):
        rho = rho_i[:,1][i] + rho_j[:,1][i]
        p_rho[:,1][i] = mj * 8/math.pow(rho,2) * p * np.dot(rij[i] , del_w[i])/math.pow(math.fabs(rij[i]),2) + math.pow(eta,2)
    return p_rho

def step_4_i(vi_prime , del_p_rho , delta_t):
    vi_n_1 = []
    for i in range(len(v_prime)):
       vi_n = vi_prime[i] - del_p_rho[:,1][i] * delta_t
       vi_n_1.append(vi_n)
    return vi_n_1

def calc_P_W(p,w_i):
    p_w = 0
    for i in range(len(w_i[:,1])):
        p_w += p * w_i[:,1][i]
    return p_w

def calc_r_W(rho,r,w_i):
    r_w = 0
    for i in range(len(w_i[:,1])):
        r_w += rho * r[i] * w_i[:,1][i]
    return r_w

############ End of DEF #################
