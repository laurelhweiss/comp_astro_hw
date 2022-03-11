import numpy as np
import time
from threading import Timer
from scipy import random
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from mycolorpy import colorlist as mcp
import matplotlib.patches as mpatches

from hyperion.model import AnalyticalYSOModel
from hyperion.util.constants import rsun, lsun, au, msun, yr, c
from hyperion.model import ModelOutput
from hyperion.util.constants import pc





###Problem 1###

def fxn(x):
  return x**(-3/2)

def mc_integrate(a, b, N):
  #randomly sample points b/w bounds (I know its not the best generator but just to play with)
  points = [random.uniform(a, b) for i in range(N)]
  total = 0
  #eval at each point
  for p in points:
    total+=fxn(p)
  #divide
  integral = (b-a)/float(N)*total
  return integral

#compute integral via trapezoids (from HW1)
def trapezoid_rule(a, b, N):
  del_x = (b-a)/N
  I = 0
  for i in range(1, N+1):
    xi = a+(i-1)*del_x
    xi1 = a+(i)*del_x
    value = (fxn(xi) + fxn(xi1))/2
    I+=value
  return I*del_x




###Problem 2###

def gen_random(seed):
  a = 1664525
  c = 1013904223
  m = 4294967296
  x = (a*seed + c)%m
  return x

R = 70000
l = 0.4
def random_walk():
  x = 0
  y = 0
  total_dist = 0
  r = 0
  
  start = time.time()
  while r < R:
    #get a random theta
    seed = time.time()
    theta = gen_random(seed)

    #calculate x,y coords of step
    x_new = l*np.cos(theta)
    y_new = l*np.sin(theta)

    #calculate distance traveled in step
    r_step = np.sqrt(x_new**2 + y_new**2)
    #update total distance traveled 
    total_dist += r_step

    #update x,y coords
    x+=x_new
    y+=y_new
    #calculate distance from center
    r = np.sqrt(x**2 + y**2)
  
  elapsed_time_lc=(time.time()-start)
  print('t_run:'+str(elapsed_time_lc))
  print('t_esc:'+str(total_dist/(3*10**10))) #time (s) = distance traveled (cm) / speed of light (cm/s)

def extrapolate():
  #from running above for a few r values
  r = [7, 70, 700, 7000, 70000]
  t_esc = [5e-9, 2e-7, 3.5e-6, 2e-5, 0.0007] #time for photon to escape sun
  t_run = [0.008, 0.35, 6, 45, 903] #time for code to run

  def exp_fxn(x, a, b): #assuming an exponential increase
    return a * np.power(x, b)

  #extrapolate out based on data
  extrap_r = np.logspace(0, np.log10(70e9), base=10) #extrapolate up until full size of sun
  popt_tesc, pcov_tesc = curve_fit(exp_fxn, r, t_esc)
  extrap_tesc = exp_fxn(extrap_r, *popt_tesc)
  popt_trun, pcov_trun = curve_fit(exp_fxn, r, t_run)
  extrap_trun = exp_fxn(extrap_r, *popt_trun)

  #plotting nonsense
  fig = plt.figure(figsize=(8, 7))
  ax = fig.add_subplot(1, 1, 1)
  ax.loglog(r, t_esc, label='t escape', marker='o')
  ax.loglog(r, t_run, label='t run', marker='o')
  ax.loglog(extrap_r, extrap_tesc, ls='--', label='t escape extrapolation')
  ax.loglog(extrap_r, extrap_trun, ls='--', label='t run extrapolation')
  ax.axvline(x=70*10**9, label='R sun', color='k', ls='--')
  ax.set_xlabel('R (cm)')
  ax.set_ylabel('t (s)')
  plt.legend()
  plt.show()

  #print out last value 
  print(np.log10(extrap_tesc[-1]), np.log10(extrap_trun[-1]))
  return 
    
  


###Problem 3###
def run_yso_model():

  '''Initalize the model'''
  m = AnalyticalYSOModel()

  '''Set the stellar parameters'''
  m.star.radius = 2*rsun
  m.star.luminosity = 5*lsun
  m.star.mass = (5)**(1/3.5) * msun
  m.star.temperature = 6200.
  '''OR add stellar spectrum'''
  #m.star.spectrum = (nu, fnu)
  #wav, fnu = np.loadtxt('kt04000g+3.5z-2.0.ascii', unpack=True)
  #nu = c / (wav * 1.e-4)


  '''Add a flared disk'''
  disk = m.add_flared_disk()
  disk.mass = 0.01 * msun
  disk.rmin = 10 * rsun
  disk.rmax = 200 * au
  disk.r_0 = m.star.radius
  disk.h_0 = 0.01 * disk.r_0
  disk.p = -1.0
  disk.beta = 1.25
  disk.dust = 'kmh94_3.1_full.hdf5'

  '''Add an Ulrich envelope'''
  envelope = m.add_ulrich_envelope()
  envelope.rc = disk.rmax
  envelope.mdot = 5.e-6 * msun / yr
  envelope.mass = 0.4*msun
  envelope.rmin = 200 * au
  envelope.rmax = 1e4 * au
  envelope.p = -2
  envelope.dust = 'kmh94_3.1_full.hdf5'

  '''
  Add a bipolar cavity
  cavity = envelope.add_bipolar_cavity()
  cavity.power = 1.5
  cavity.theta_0 = 20
  cavity.r_0 = envelope.rmax
  cavity.rho_0 = 5e4 * 3.32e-24
  cavity.rho_exp = 0.
  cavity.dust = 'kmh94_3.1_full.hdf5'
  '''

  '''Use raytracing to improve s/n of thermal/source emission'''
  m.set_raytracing(True)

  '''Use the modified random walk'''
  m.set_mrw(True, gamma=2.)

  '''Set up grid'''
  m.set_spherical_polar_grid_auto(200, 100, 5) #finegrid: 1000, 500, 50 #default: 200, 100, 5 #coarse: 40 25 1

  '''Set up SED'''
  sed = m.add_peeled_images(sed=True, image=False)
  sed.set_viewing_angles(np.linspace(0., 90., 10), np.repeat(45., 10))
  sed.set_wavelength_range(150, 0.02, 2000.)

  '''Set number of photons'''
  m.set_n_photons(initial=1e3, imaging=1e3,
                raytracing_sources=1e3, raytracing_dust=1e3)

  '''Set number of temperature iterations and convergence criterion'''
  m.set_n_initial_iterations(10)
  m.set_convergence(True, percentile=99.0, absolute=2.0, relative=1.1)

  '''Write out file'''
  m.write('class1_example_1e3_defaultgrid.rtin', overwrite=True) #change
  m.run('class1_example_1e3_defaultgrid.rtout', mpi=True, overwrite=True) #change
  return

def plot_sed():
  fig = plt.figure(figsize=(8, 7))
  ax = fig.add_subplot(1, 1, 1)
  
  photon_files = ['class1_example_1e3_defaultgrid.rtout', 'class1_example_1e4_defaultgrid.rtout', 'class1_example_1e5_defaultgrid.rtout', 'class1_example_1e6_defaultgrid.rtout']
  grid_files = ['class1_example_1e5_coarsegrid.rtout','class1_example_1e5_defaultgrid.rtout', 'class1_example_1e5_finegrid.rtout', ]
  colors = mcp.gen_color(cmap="viridis",n=3)
  n_photons = ['1e3', '1e4', '1e5', '1e6']
  grid_types = ['coarse', 'default', 'fine']
  patches = []

  for f, c, gt in zip(grid_files, colors, grid_types):
    mo = ModelOutput(f)
    sed = mo.get_sed(aperture=-1, distance=300. * pc)
    ax.loglog(sed.wav, sed.val.transpose(), color=c)
    ax.set_xlim(0.03, 2000.)
    ax.set_ylim(2.e-16, 5e-9)
    ax.set_xlabel(r'$\lambda$ [$\mu$m]')
    ax.set_ylabel(r'$\lambda F_\lambda$ [ergs/cm$^2/s$]')

    patch = mpatches.Patch(color=c, label=gt+' grid')
    patches.append(patch)

  plt.legend(handles=patches)
  plt.show()







