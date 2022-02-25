import rebound
import matplotlib
import numpy as np
import matplotlib.pyplot as plt

#HAVE TO INITIALIZE NEW SIM EACH TIME, NOT THE NEATEST
#BUT JUST REDEFINE EVERYTHING FOR EACH PROBLEM
#you can't see differences running forward in time. you have to
#run backwards
def problem_2(integrator, dt, add=False):
  #2a - same as example in documentation
  sim = rebound.Simulation()

  sim.add("Sun")
  sim.add("Jupiter")
  sim.add("Saturn")
  sim.add("Churyumov-Gerasimenko", m=5.03e-18)

  #2c - add in more planets

  if add==True:
    sim.add("Mercury")
    sim.add("Venus")
    sim.add("Earth")
    sim.add("Mars")
    sim.add("Uranus")
    sim.add("Neptune")

  Noutputs = 1000 #easier to see differences in models
  year = 2.*np.pi # One year in units where G=1
  times = np.linspace(0.,-70.*year, Noutputs) #forward in time is what you wanted right?
  x = np.zeros((3,Noutputs))
  y = np.zeros((3,Noutputs))
  z = np.zeros((3,Noutputs))
  
  sim.integrator = integrator #ias15, leapfrog, whfast
  sim.move_to_com()
  ps = sim.particles
  sim.dt = dt
  sim.ri_ias15.min_dt = -1

  #1 = jupiter, 2 = saturn, 3 = comet
  for i,time in enumerate(times):
    sim.integrate(time)
 
    x[0][i] = ps[1].x   
    y[0][i] = ps[1].y
    z[0][i] = ps[1].z
    x[1][i] = ps[2].x
    y[1][i] = ps[2].y
    z[1][i] = ps[2].z
    x[2][i] = ps[3].x
    y[2][i] = ps[3].y
    z[2][i] = ps[3].z

  return x[0], y[0], x[1], y[1], x[2], y[2]

def plot_2c():
  x1, y1, x2, y2, x3, y3 = problem_2('whfast', -1, add=False)
  wp_x1, wp_y1, wp_x2, wp_y2, wp_x3, wp_y3 = problem_2('whfast', -1, add=True)
  fig = plt.figure(figsize=(6.5,6.5))
  ax = plt.subplot(111)
  ax.plot(x1, y1, label='Jupiter orbit')
  ax.plot(x2, y2, label='Saturn orbit')
  ax.plot(x3, y3, label='Comet orbit', ls='--')
  ax.plot(wp_x1, wp_y1, label='Jupiter orbit (all planets)')
  ax.plot(wp_x2, wp_y2, label='Saturn orbit (all planets)')
  ax.plot(wp_x3, wp_y3, label='Comet orbit (all planets)', ls='--')
  
  ax.set_ylabel('Y (AU)')
  ax.set_xlabel('X (AU)')
  ax.set_title('whfast, dt=-1')
  plt.legend()
  plt.show()
  

def plot_2b():
  i_x1, i_y1, i_x2, i_y2, i_x3, i_y3 = problem_2('ias15', -0.01, add=False)
  w_x1, w_y1, w_x2, w_y2, w_x3, w_y3 = problem_2('whfast', -0.01, add=False)
  ws_x1, ws_y1, ws_x2, ws_y2, ws_x3, ws_y3 = problem_2('whfast', -1, add=False)
  l_x1, l_y1, l_x2, l_y2, l_x3, l_y3 = problem_2('leapfrog', -0.01, add=False)
  ls_x1, ls_y1, ls_x2, ls_y2, ls_x3, ls_y3 = problem_2('leapfrog', -1, add=False)
  
  fig = plt.figure(figsize=(6.5,6.5))
  ax = plt.subplot(111)
  ax.plot(i_x1, i_y1, label='Jupiter orbit')
  ax.plot(i_x2, i_y2, label='Saturn orbit')
  ax.plot(i_x3, i_y3, label='Comet orbit (IAS15, dt=-0.01)', ls='--')
  ax.plot(w_x3, w_y3, label='Comet orbit (WHFast, dt=-0.01)', ls='--')
  ax.plot(ws_x3, ws_y3, label='Comet orbit (WHFast, dt=-1)', ls='--')
  #ax.plot(l_x3, l_y3, label='Comet orbit (Leapfrog, dt=-0.01)', ls='--')
  #ax.plot(ls_x3, ls_y3, label='Comet orbit (Leapfrog, dt=-1)', ls='--')
  
  ax.set_ylabel('Y (AU)')
  ax.set_xlabel('X (AU)')
  plt.legend()
  plt.show()

  

#need to initialize this separately bc of sim object
def run_pluto_orbit():
  sim = rebound.Simulation()
  
  sim.add("Sun")
  sim.add("Jupiter")
  sim.add("Saturn")
  sim.add("Mercury")
  sim.add("Venus")
  sim.add("Earth")
  sim.add("Mars")
  sim.add("Uranus")
  sim.add("Neptune")
  sim.add("Pluto")

  Noutputs = 100 
  year = 2.*np.pi 
  times = np.linspace(0.,10.*year, Noutputs)
  x = np.zeros((1,Noutputs))
  y = np.zeros((1,Noutputs))
  z = np.zeros((1,Noutputs))
  
  sim.integrator = "ias15"
  sim.move_to_hel()
  ps = sim.particles

  #pluto = 9
  for i,time in enumerate(times):
    sim.integrate(time)
 
    x[0][i] = ps[9].x   
    y[0][i] = ps[9].y
    z[0][i] = ps[9].z

  return x[0], y[0], z[0]

def problem_3():
  pluto_normal_x, pluto_normal_y, pluto_normal_z = run_pluto_orbit()

  #now run with planet X
  sim = rebound.Simulation()
  
  sim.add("Sun")
  sim.add("Jupiter")
  sim.add("Saturn")
  sim.add("Mercury")
  sim.add("Venus")
  sim.add("Earth")
  sim.add("Mars")
  sim.add("Uranus")
  sim.add("Neptune")
  sim.add("Pluto")
  sim.add(m=250, e=0.4, a=460)

  Noutputs = 100 
  year = 2.*np.pi 
  times = np.linspace(0.,10.*year, Noutputs)
  x = np.zeros((2,Noutputs))
  y = np.zeros((2,Noutputs))
  z = np.zeros((2,Noutputs))
  
  sim.integrator = "ias15"
  sim.move_to_hel()
  ps = sim.particles

  #pluto = 9, planet X = 10
  for i,time in enumerate(times):
    sim.integrate(time)
 
    x[0][i] = ps[9].x   
    y[0][i] = ps[9].y
    z[0][i] = ps[9].z
    x[1][i] = ps[10].x   
    y[1][i] = ps[10].y
    z[1][i] = ps[10].z

  fig = plt.figure(figsize=(6.5,6.5))
  ax = plt.subplot(111)
  ax.plot(pluto_normal_x, pluto_normal_y, label='Pluto orbit (w/o Planet X)')
  ax.plot(x[0], y[0], label='Pluto orbit (w/ Planet X)')
  ax.plot(x[1], y[1], label='Planet X orbit')

  ax.set_ylim([-31.0, -27.0])
  ax.set_xlim([15, 25])
  ax.set_ylabel('Y (AU)')
  ax.set_xlabel('X (AU)')
  ax.set_title('IAS15, dt=-0.01')
  plt.legend()
  plt.show()

  #calc when orbit changed by 10% at year 10...this is how im assuming that:
  #final coords w planet X off by 10% of final coords without planet X
  r_normal = np.sqrt(pluto_normal_x[-1]**2 + pluto_normal_y[-1]**2 + pluto_normal_z[-1]**2)
  r_w_planetx = np.sqrt(x[0][-1]**2 + y[0][-1]**2 + z[0][-1]**2)
  print(abs(r_normal - r_w_planetx)/r_normal)

  

  
  
