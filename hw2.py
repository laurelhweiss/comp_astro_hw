import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

rc('text', usetex=True)
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})

'''
TYPOS IN DOCUMENTS:
v = u
Eqns 2 & 3 in the HW doc should have n+1 in the superscript not n

useful:
- to compute element_i+1 +/- element_i in array A:
A[1:] +/- A[:-1]
(leaves off last element)
'''

N=200 
gamma = 1.4
q0 = 4
q1 = 0.5
CFL = 0.5
t_choose = 0.245

'''
initialize arrays
x and v arrays defined at edges, so need to be +1 in length
'''
x_array = np.linspace(0, 2, N+1) #position, len=201
v_array = np.zeros(N+1) #velocity, len=201

p_array = np.array([1 if x <= 0.75 else 0.1 for x in x_array[:-1]]) #pressure, len=200
rho_array = np.array([1 if x <= 0.75 else 0.125 for x in x_array[:-1]])#density, len=200
e_array = p_array/((gamma-1)*rho_array) #energy, len=200

#mass arrays stay constant in time
dm_half_array = rho_array*2 / N #mass array #1 (center vals, len=200)
dm_array = np.pad((dm_half_array[:-1] + dm_half_array[1:])/2, 1, 'edge') #mass array #2 (edge vals, len=201)

q_array = np.zeros(N) #viscosity, len=200
cs_array = np.sqrt(gamma*p_array/rho_array)

'''
little function to calculate q values....
slightly easier for me to read. will use in loop below. 
yes, i am 100% aware that I don't have to do this, but
I don't like looking at one giant line
'''
def update_qvals():
  rho_bar = 0.5*(1/rho_array_new + 1/rho_array)
  condition = (v_array_new[1:]-v_array_new[:-1]) / (x_array_new[1:]-x_array_new[:-1])
  #the return of my favorite numpy function: np.where...i love her so much, LOOK AT HER GO 
  q_array_new = np.where(condition < 0, (q0*(v_array_new[1:]-v_array_new[:-1])**2 - q1*(v_array_new[1:]-v_array_new[:-1]))*(cs_array_new/rho_bar), 0)
  return q_array_new

'''
Loop through time...in each step, calculate new arrays using the old/initialized
arrays, then update the old arrays to be the new ones. calculate a dt, add to t,
go to next step.
'''
dt = 0.0001
t=0
while t < t_choose:
  #in this loop, need to pad pressure and viscosity arrays for manipulating in the new v array
  p_array_pad = np.pad(p_array, 1, 'edge')
  q_array_pad = np.pad(q_array, 1, 'edge')

  #begin updating
  v_array_new = v_array - dt*(p_array_pad[1:]-p_array_pad[:-1] + q_array_pad[1:]-q_array_pad[:-1])/ dm_array
  x_array_new = x_array + dt*v_array_new
  
  rho_array_new = dm_half_array / (x_array_new[1:] - x_array_new[:-1])
  cs_array_new = np.sqrt(gamma*p_array/rho_array) #sound speed is lagged, calculate based on old values only
  q_array_new = update_qvals()
  e_array_new = e_array - (p_array+q_array)*(1/rho_array_new - 1/rho_array)
  p_array_new = e_array_new*rho_array_new*(gamma-1)

  #update the old arrays to new arrays
  v_array = v_array_new
  x_array = x_array_new
  rho_array = rho_array_new
  cs_arrary = cs_array_new
  q_array = q_array_new
  e_array = e_array_new
  p_array = p_array_new

  dx = x_array[1:]-x_array[:-1]
  dt = np.min(CFL*dx/(cs_array + v_array[:-1]))
  t+=dt

'''
plotting garbage, i know this could be more nicely written but
im spent to be totally honest
'''
fig, axs = plt.subplots(2, 2)
axs[0, 0].plot(x_array[:-1], rho_array)
axs[0, 0].set_xlabel('x', fontsize=18)
axs[0, 0].set_ylabel(r'$\rho$', fontsize=18)
axs[0, 1].plot(x_array, v_array)
axs[0, 1].set_xlabel('x', fontsize=18)
axs[0, 1].set_ylabel('v', fontsize=18)
axs[1, 0].plot(x_array[:-1], p_array)
axs[1, 0].set_xlabel('x', fontsize=18)
axs[1, 0].set_ylabel('P', fontsize=18)
axs[1, 1].plot(x_array[:-1], e_array)
axs[1, 1].set_xlabel('x', fontsize=18)
axs[1, 1].set_ylabel('Energy', fontsize=18)

axs[0, 0].minorticks_on()
axs[0, 0].tick_params(axis='both', which='both', direction='in', labelsize=14)
axs[0, 1].minorticks_on()
axs[0, 1].tick_params(axis='both', which='both', direction='in', labelsize=14)
axs[1, 0].minorticks_on()
axs[1, 0].tick_params(axis='both', which='both', direction='in', labelsize=14)
axs[1, 1].minorticks_on()
axs[1, 1].tick_params(axis='both', which='both', direction='in', labelsize=14)

plt.tight_layout()
plt.show()


