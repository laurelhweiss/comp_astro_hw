import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate


###----1-----###

def problem1_parta():
  var = 0.1
  var = np.single(var)
  var = format(var, '.60g')
  print(var)
  
  var = 0.1
  var = np.double(var)
  var = format(var, '.60g')
  print(var)
  return

def problem1_partb():
  e = 1
  while e + 1 > 1:
    e = e/1.7
    print(e)
  return

###----2-----###

#function to integrate
def f1(x):
  return (1/(x**(3/2)))

#compute integral via rectangles
def rectangle_rule(f, N, a, b):
  del_x = (b-a)/N
  I = 0
  for i in range(1, N+1):
    xi = a+(i-1)*del_x
    value = f(xi)*del_x
    I+=value
  return I

#compute integral via trapezoids
def trapezoid_rule(f, N, a, b):
  del_x = (b-a)/N
  I = 0
  for i in range(1, N+1):
    xi = a+(i-1)*del_x
    xi1 = a+(i)*del_x
    value = (f(xi) + f(xi1))/2
    I+=value
  return I*del_x

def problem2_parta():
  I_exact = (2 - 2/np.sqrt(5))
  
  #loop through step sizes to get points to plot, try 10000 steps
  steps = np.arange(10, 10000, 100)
  rect_errors, trap_errors = [], []
  for s in steps:
    I_rect = rectangle_rule(f1, s, 1, 5)
    err_rect = abs(I_rect - I_exact)/I_exact #fractional error
    rect_errors.append(err_rect)
    
    I_trap = trapezoid_rule(f1, s, 1, 5)
    err_trap = abs(I_trap - I_exact)/I_exact #fractional error
    trap_errors.append(err_trap)

  #plotting stuff
  plt.figure(figsize=(9, 3))
  plt.loglog(steps, rect_errors, color='teal', label='Rectangle Rule')
  plt.loglog(steps, trap_errors, color='darkorange', label='Trapezoidal Rule')
  plt.axhline(y=0.001, lw=0.75, color='k', ls='--')
  plt.axhline(y=0.00001, lw=0.75, color='k', ls='--') 
  plt.minorticks_on()
  plt.tick_params(axis='both', which='both', direction='in')
  plt.xlabel('Steps')
  plt.ylabel('Errors')
  plt.legend()
  plt.show()


def problem2_partb():
  #use scipy integrate quad
  I_sci, sci_err = integrate.quad(f1, 1, 5)
  print(I_sci, sci_err)
  '''
  The order of the error of scipy.integrate.quad is 1e-13.
  This function uses the Fortran library QUADPACK,
  which, for definite integrals, performs the integral
  using a Clenshaw-Curtis method, which uses Chebyshev
  polynomials (Google then tells me that this essentially
  means integrating using a change of variable x=cos(y),
  converting the integral to a cosine series). To obtain a
  similar degree of accuracy, I would need step sizes on
  the order of ~4/10**7 for the trapezoid rule, and
  ~4/10**12 for the rectangle rule (extrapolating by eye
  from the graph)
  '''

###----3-----###

#function to integrate
O_m = 0.3
O_lam = 0.7
h = 0.7 #choosing this cosmology
def f2(z):
  return (O_m*(1+z)**3 +(1-O_m-O_lam)*(1+z)**2 + O_lam)**(-0.5)

def problem3_partb():
  #say 10000 steps, use trap rule, should be decently accurate (see above)
  I = trapezoid_rule(f2, 10000, 0, 2)
  co_dist = (3000/h)*I
  print(co_dist)#in Mpc
  I_sci, err = integrate.quad(f2, 0, 2)
  co_dist_sci = (3000/h)*I_sci #check
  print(co_dist_sci)

  co_dists = []
  zs = np.linspace(0, 10, 100) #100 points should be decently smooth
  for z in zs:
    I = trapezoid_rule(f2, 10000, 0, z)
    co_dists.append((3000/h)*I)

  #plotting stuff
  plt.figure(figsize=(9, 3))
  plt.plot(zs, co_dists, color='teal') 
  plt.minorticks_on()
  plt.tick_params(axis='both', which='both', direction='in')
  plt.ylabel('Comoving Distance (Mpc)')
  plt.xlabel('z')
  plt.show()

  


  
  
  
  
  














