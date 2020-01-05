#libraries
import matplotlib.pyplot as plt
import numpy as np
import math
from math import pi
import scipy.optimize
from matplotlib.lines import Line2D

#plotting
fig, ax = plt.subplots(figsize=(60,5))
ax.set_aspect("equal")

# select how many circles you want

N = 2

#change bead properties
                        
Diameter =0.019
magnetic_moment = 3
mass = 0.200
gravity = 9.81

# constants
mu_zero_devided_by_four_pi = 4*pi*1*10**-7
R = Diameter/2

#starting magnetic moment directions
angle = np.linspace(0,-pi/2,N)

# scipy.minimize 
def energy_function(params):
    for i in range(1,N):
        angle[i] = params[i-1]
#sets the condition that provides that circles don't overlap
    x = [0]
    y = [0]
    for i in range(1,N):
        x.append(np.cos(angle[i])*2*R+x[i-1])
        y.append(np.sin(angle[i])*2*R+y[i-1])
    
# calculates x and y components
    x_end = []
    y_end = []
    for i in range(N):
        mx_component = (magnetic_moment * np.cos(angle))
        x_end.append(mx_component)
        my_component = (magnetic_moment * np.sin(angle))
        y_end.append(my_component)
        
# first_part (mu zero * -0.5 / 4*pi*r_cubed)    
    first_part = []
    for a in range(N):
        for b in range(a+1,N):
            first_part.append((-0.5* mu_zero_devided_by_four_pi)/(4*pi*math.sqrt(((x[a]-x[b])**2 + (y[a]-y[b])**2)**3)))
            
    first_part = np.asarray(first_part)

#x and y hat
    x_hat = []
    y_hat = []
    
    for a in range(N):
        for b in range(a+1,N):
            x_hat.append((x[b]-x[a])/(math.sqrt((x[a]-x[b])**2 + (y[a]-y[b])**2)))
            y_hat.append((y[b]-y[a])/(math.sqrt((x[a]-x[b])**2 + (y[a]-y[b])**2)))
            
    x_hat = np.asarray(x_hat)
    y_hat = np.asarray(y_hat)
    
#we are finding the minimum of this function
    test = 0
    pot_e = 0
    for a in range(N):
        for b in range(a+1,N):
            partial_energy =  (first_part[test]*
            (3*(mx_component[b]*x_hat[test]+my_component[b]*y_hat[test])*
            (mx_component[a]*x_hat[test] + my_component[a]*y_hat[test])
            -(mx_component[a]*mx_component[b]+my_component[b]*my_component[a])))
            pot_e = pot_e + partial_energy
            test += 1
#including gravity here
    for a in range(N):
        pot_e = pot_e +mass*gravity*y[a]
        
    return pot_e 
    

#sets initial angles
params_guess = [angle]
#minimization method
minimizer_kwargs = {"method":"Nelder-Mead"}

res = scipy.optimize.basinhopping(energy_function,  params_guess,niter=333, T=1.0, stepsize=0.4,
                                  minimizer_kwargs=minimizer_kwargs)
                                
print(res)
solution = []
solution = res.x
print('+++++++++++++++++++++++++++++++++++++++++++++++++++' )
print('Result with gravity = ' + str(energy_function(res.x)))
solution = np.asarray(solution)

#exctracts angle,x and y value
x = [0]
y = [0]
for i in range(1,N):
    x.append(np.cos(angle[i])*2*R+x[i-1])
    y.append(np.sin(angle[i])*2*R+y[i-1])

#calculates x and y component
x_end = []
y_end = []
for i in range(N):
    mx_component = (magnetic_moment * np.cos(angle))
    x_end.append(mx_component)
    my_component = (magnetic_moment * np.sin(angle))
    y_end.append(my_component)

#Draws beads and magnetic moment directions
Factor = 0.00005*magnetic_moment**1/2/Diameter**1/4
"""Warning, this might not be addaptable for all bead combinations, when 
Changing the initial values, you might need to change the Factor value, so the arrow endings
are good enough to get a clear picture"""
for l in range(N):
    plt.arrow(x[l], y[l],(mx_component[l]*Factor), (my_component[l]*Factor),
              length_includes_head=False,head_width=R*0.7,head_length=R*0.8)

for m in range(N):
    c = plt.Circle((x[m],y[m]), R,facecolor="None", edgecolor='green')
    ax.add_artist(c)
    
    
#  r_hat as well as mu zero* -0.5/4*pi*r_cubed
first_part = []
x_hat = []
y_hat = []
for a in range(N):
    for b in range(a+1,N):
        x_hat.append((x[b]-x[a])/(math.sqrt((x[a]-x[b])**2 + (y[a]-y[b])**2)))
        y_hat.append((y[b]-y[a])/(math.sqrt((x[a]-x[b])**2 + (y[a]-y[b])**2)))
        first_part.append((-0.5* mu_zero_devided_by_four_pi)/(4*pi*math.sqrt(((x[a]-x[b])**2 + (y[a]-y[b])**2)**3)))

x_hat = np.asarray(x_hat)
y_hat = np.asarray(y_hat)
first_part = np.asarray(first_part)

#formula for calculating only the dipole-dipole interactions of the potential energy 
energy = []
test = 0
for a in range(N):
    for b in range(a+1,N):
        energy.append(first_part[test]*
        (3*(mx_component[b]*x_hat[test]+my_component[b]*y_hat[test])*
        (mx_component[a]*x_hat[test] + my_component[a]*y_hat[test])
        -(mx_component[a]*mx_component[b]+my_component[b]*my_component[a])))
        test+=1
        
max_energy = sum(energy)
print('Result without gravity:',max_energy)

#psi parameter value
psi = mu_zero_devided_by_four_pi*magnetic_moment**2/4*pi*mass*gravity*Diameter**4
print('Psi = '+str(psi))

# some more plotting
ax.set_xlabel('Y coordinate, m',fontsize=14)
ax.set_ylabel('X coordinate, m',fontsize=14)

legend_elements = [Line2D([0], [0], marker='o', color='w', label='Experimental results',
                          markerfacecolor='None',markeredgecolor='r', markersize=12),
                   Line2D([0], [0], marker='o', color='w', label='Computational results',
                          markerfacecolor='None',markeredgecolor='g', markersize=12),
                   Line2D([0], [0], marker='None', color='w', label='Number of beads: '+str(N),
                          markerfacecolor='None',markeredgecolor='r', markersize=12),
                   Line2D([0], [0], marker='None', color='w', label='Diameter: '+"%.3f" % Diameter +' m',
                          markerfacecolor='None',markeredgecolor='r', markersize=12)]

ax.legend(loc='center left', bbox_to_anchor=(0.8, 0.5), fancybox=True, shadow=True,
          handles=legend_elements)
plt.grid( color='#111111', linestyle='--')
ax.set_title('Magnetic bead chain in a gravitational field\n',fontsize=16)
plt.xlim([-Diameter,Diameter*1.2*N])
plt.ylim([-Diameter*1.2*N,Diameter])
plt.show()