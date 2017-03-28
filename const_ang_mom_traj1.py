"""" Python script: const_ang_mom_traj1.py 
 Problem M1.4   
 Script to compute constant angular momentum trajectories in spherical 
 coordinates with curvature terms included and to show time development.
 Stars mark time at one day intervals.
 Time differencing by 3rd order Adams-Bashforth method.
"""
import matplotlib.pyplot as plt
import numpy as np
import math as mt

#INPUTS
print 'Initial longitude is zero. Specify latitude and speed when asked.'
try:
    init_lat=int(raw_input('(1) Give an initial latitude in degrees '))
except ValueError:
    print "Not a number"
try:
    u0=int(raw_input('(2) Give a zonal wind in m/s '))
except ValueError:
    print "Not a number"
try:
    v0=int(raw_input('(3) Give a meridional wind in m/s '))
except ValueError:
    print "Not a number"
rad = 6.37e6
omega = 7.292e-5
lat = mt.pi*init_lat/180.
dt = 60
dt12 = dt/12
try:
    runtime=int(raw_input('Specify integration time in days'))
except ValueError:
    print "Not a number"
time = runtime*24*3600
ind = 1
X  = np.array([u0, v0, 0, lat])
xprimn= np.zeros((1,4))

# vector xprimn has components of dX/dt
xprimn[0][0] = (2*omega+ind*X[0]/(rad*mt.cos(X[3])))*mt.sin(X[3])*X[2]
xprimn[0][1] = -(2*omega+ind*X[0]/(rad*mt.cos(X[3])))*mt.sin(X[3])*X[0]
xprimn[0][2] = X[0]/(rad*mt.cos(X[3]))
xprimn[0][3] = X[1]/rad
xprim1 = xprimn
xprim2 = xprim1

plt.figure(1)
plt.xlim([-180, 180])
plt.ylim([-60, 60])
plt.xlabel('longitude')
plt.ylabel('latitude')
plt.title('constant angular momentum trajectory')

for t in range(0, time, dt):
    Xn = X +dt12*(23*xprimn -16*xprim1 +5*xprim2)
    X = Xn
    xprim2 = xprim1
    xprim1 = xprimn
    xprimn[0][0] = (2.*omega+ind*X[0]/(rad*cos(X[3])))*sin(X[3])*X[3]
    xprimn[0][1] = -(2.*omega+ind*X[0]/(rad*cos(X[3])))*sin(X[3])*X[0]
    xprimn[0][2] = X[0]/(rad*cos(X[3]))
    xprimn[0][3] = X[1]/rad
    if X[2] > mt.pi:
        X[2]= -2*mt.pi+X[2]
    
    if X[2] < -mt.pi:
        X[2] = 2*mt.pi+X[2]
    
    v = X[2:3]*180./mt.pi
    p = plot(v[0],v[1],':','EraseMode','none', 'MarkerSize',5)
    set(p,'Xdata',v[0],'Ydata',v[1])
    if mod(t,24*3600)==0
        plot(v[0],v[1],'*','EraseMode','none', 'MarkerSize',8)
    end
    drawnow
    t = t +dt