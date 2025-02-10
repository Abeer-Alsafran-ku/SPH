import math 
import numpy as np
from sklearn.neighbors import NearestNeighbors
import time
from defs import *
from ic import *
#-----------------------

start_time = time.time()

# the point will chng as time pass #i #j
point = np.array([[0.0, coordinates[0][1], 0.0]])
distance, indxs = distance_ij(point, coordinates,k) # distance and indxs of size k=75
''' the indices of the knn of the first particle
for i in range(k):
    if i == k-1 :
        print(indxs[0][i], indxs[0][k-1]) # last in array 
    else:    
        print(indxs[0][i], indxs[0][i+1])
'''

############ step 1 ####################
r = np.zeros((k,3)) 
for i in range(len(distance[0])):
  r[:, 1][i] = distance[0][i]
rj = calc_rj(r,k)
r_resized = np.resize(r, (n,3)) # resize the array to be compatible with next steps 
rj_resized = np.resize(rj, (n,3)) # resize rj the array to be compatible with next steps 
# calc rij = ri - rj
rij = calc_rij(r,rj_resized,k)
rij_resized = np.resize(rij,(n,3)) # resize rj the array to be compatible with next steps 
# v_in_lam == v0 == vi
vj = calc_vj(vi,i) # to be calc
vij = calc_vij(vi, vj) # vij = vi - vj
vj_resized = np.resize(vj,(n,3))
vij_resized = np.resize(vij, (n,3))
v0 = np.resize(v0,(n,3))
for i in range(k):#here need to specify the length based on k or n??
    v_prime.append(step_1_i(v0,rij_resized,vij_resized.copy(),h,fg))
########### End of step 1 ###############
############### step 2 ##################
'''ri = ri_n + vi . delta t '''
ri = step_2_i(r,v_prime,delta_t)
######## end of step 2 ##################
############### step 3 ##################
'''Pressure P(rho) = P0[(rho/rho0)^gama -1] + X'''
c = 10*Umax
init_p = calc_init_p(rho,c,gama)
w_i= calc_W(ri,h)
rho_star = calc_rho_star(mj , w_i)
rho_i , rho_j = calc_rho_ij(rho_star,indx=0)
p = step_3_i(init_p,rho,rho_i[:,1][0],gama,kai)
######## end of step 3 ##################
############ boundry condition #########
upr = lowr = 0
p = 0
g = 0
rho = 1 
for i in range(len(w_i)):
    upr += vj[i] + w_i[i]
    lowr += w_i[i]

vlcty_of_bndry = 2* va - upr[1]/lowr[1]
wall_pressure_Pw = (calc_P_W(p,w_i) + g*calc_r_W(rho,r,w_i)) / lowr[1] 
print(wall_pressure_Pw)
kai = 0.05*p
wall_dnsty_rho_W = 1
# wall_dnsty_rho_W = rho * math.pow(((wall_pressure_Pw[1]-kai)/p)+1,1/gama) 
############ End boundry condition #####
############### step 4 ##################
'''vi_n_1 = vi_prime - (delta P / rho)i->n+1 * delta_t '''
del_w = calc_del_W(rij,h)
del_p_rho = calc_delta_p_rho(mj,rho_i,rho_j,p,rij,del_w ,eta)
vi_prime = v_prime # copy
vi_n_1 = step_4_i(vi_prime , del_p_rho , delta_t)
############### end of step 4 ###########
############### step 5 ##################
'''vi : the initial voleocity | vi_n_1 : calculated from step 4 '''
ri_n_1 = []
for i in range(len(ri)):
    ri_n= ri[i] + (vi[:,0][0] + vi_n_1[i]) * delta_t/2  
    ri_n_1.append(ri_n)
############### end of step 5 ###########
############### step 6 ##################
'''t : time | delta_t : delta time | n : time step / current step'''
current_step +=1
t += delta_t
############### end of step 6 ###########
end_time = time.time()
time_stamp = end_time - start_time

print(f'''
      Time = {t}\n
      No. of particles = {n}\n
      Pressure = {p}\n
      Velocity = {vi[0]}\n
      Current step = {current_step}\n
      Wall pressure (P W) = {wall_pressure_Pw}\n
      Wall density (rho W) = {wall_dnsty_rho_W}\n  
      Total time of execution = {time_stamp} sec\n
      '''           
      )


