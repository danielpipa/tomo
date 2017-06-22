# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 16:07:22 2017

@author: ajmello
"""
import Tomo_func
import numpy
import math
import scipy
import matplotlib.pyplot
pi = numpy.pi;

#**************************************************
# INITIALIZATION - WFS DATA
#**************************************************
r0 = 0.2;
L0 = 50.;
f0 = 1/L0;
frac = numpy.array([1., 1.]);   ###This will need to be included on the HH_array, leaving at ones for now####
tel_diam = 8.0; #Telescope Diameter (m)
ang_dist = 30.;  #Angular distance between outer guide stars (arcsec)
WFS_angles = ang_dist*numpy.array([[numpy.sin(0),numpy.cos(0)], [numpy.sin(72*pi/180),numpy.cos(72*pi/180)], [numpy.sin(144*pi/180),numpy.cos(144*pi/180)], [numpy.sin(216*pi/180),numpy.cos(216*pi/180)], [numpy.sin(288*pi/180),numpy.cos(288*pi/180)]]);   #for each WFS, in x,y
WFS_target = numpy.array([0.,0.]);
WFS_size = 41;   #Number of subapertures in diameter (Only odd numbers)


#**************************************************
# INITIALIZATION - TURBULENCE DATA
#**************************************************
altitudes = numpy.array([0.,4000.]);     #layers altitudes (m)
layer_size = WFS_size;
alt_p = altitudes/(tel_diam/layer_size);  #layer altitudes relative to pixel size (pixels)
shifts = numpy.tan(ang_dist*pi/648000)*alt_p;        #maximum shift in each layer
extra_p = numpy.zeros(2);
for layer_i in range(0,2):
    extra_p[layer_i] = 2*math.ceil(shifts[layer_i]);        #extra pixels for layer movement
    if extra_p[layer_i]%2!=0:
        extra_p[layer_i] = extra_p[layer_i]+1;    #must be even
sizes = numpy.array([int(layer_size), int(layer_size+extra_p[1])]);


#**************************************************
# CREATE COVARIANCE MATRIX
#**************************************************
c = math.pow((24./5.*scipy.special.gamma(6./5.)),(5./6.))*scipy.special.gamma(11./6.)/(math.pow(2.,(5./6.))*math.pow(pi,(8./3.)))*math.pow((r0*f0),(-5./3.));
x = numpy.arange(0,layer_size);
y = x;
X_index, Y_index = numpy.meshgrid(x, y)
Indexes = numpy.vstack([numpy.reshape(X_index,len(X_index)**2), numpy.reshape(Y_index,len(Y_index)**2)])
mindist = 1e-6; #Avoid zeros
Cov_M1 = numpy.zeros([layer_size**2,layer_size**2]);
for lin in range(0,layer_size**2):
    for col in range(0,layer_size**2):
        #Distance
        dist = numpy.sqrt(numpy.sum((Indexes[:,col] - Indexes[:,lin])**2));
        dist = dist*tel_diam/layer_size + mindist;  #distance between voxels
        Cov_M1[lin,col] = c*frac[0]*math.pow((2*pi*dist*f0),(5./6.))*scipy.special.kv(5./6.,2*pi*f0*dist);

x = numpy.arange(0,sizes[1]);
y = x;
X_index, Y_index = numpy.meshgrid(x, y)
Indexes = numpy.vstack([numpy.reshape(X_index,len(X_index)**2), numpy.reshape(Y_index,len(Y_index)**2)])
mindist = 1e-6; #Avoid zeros
Cov_M2 = numpy.zeros([sizes[1]**2,sizes[1]**2]);
for lin in range(0,sizes[1]**2):
    for col in range(0,sizes[1]**2):
        #Distance
        dist = numpy.sqrt(numpy.sum((Indexes[:,col] - Indexes[:,lin])**2));
        dist = dist*tel_diam/layer_size + mindist;  #distance between voxels
        Cov_M2[lin,col] = c*frac[1]*math.pow((2*pi*dist*f0),(5./6.))*scipy.special.kv(5./6.,2*pi*f0*dist);


#**************************************************
# REGERATE RANDOM DISPLACEMENTS
#**************************************************

L1 = scipy.linalg.cholesky(Cov_M1, lower=True); # Cholesky decomposition
L2 = scipy.linalg.cholesky(Cov_M2, lower=True); # Cholesky decomposition

layersXY = numpy.hstack([L1.dot(numpy.random.normal(0,1,len(L1)))+1j*L1.dot(numpy.random.normal(0,1,len(L1))),  L2.dot(numpy.random.normal(0,1,len(L2)))+1j*L2.dot(numpy.random.normal(0,1,len(L2)))]);     #Two layers with random slopes

layersXY = layersXY/numpy.sqrt(numpy.var(layersXY));  #Set variance to one ##There must be a better way##

del Cov_M1, Cov_M2, L1, L2

#temp = numpy.reshape(layersXY[0:sizes[0]**2],[sizes[0],sizes[0]])
#matplotlib.pyplot.imshow(temp.real)     #Show a layer


#%%
#**************************************************
# simulate WFS
#**************************************************
WFS = Tomo_func.HH(layersXY,sizes,alt_p,WFS_angles,0);
#add noise
WFS = WFS + 0.01*(numpy.random.normal(0,1,WFS.shape)+1j*numpy.random.normal(0,1,WFS.shape));

#temp = numpy.reshape(WFS[0:WFS_size**2],[WFS_size,WFS_size]);
#matplotlib.pyplot.imshow(temp.real)     #Show a WFS

#Testa Transposto
#retorno = Tomo_func.HH(WFS,sizes,alt_p,WFS_angles,1);
#temp = numpy.reshape(WFS[0:WFS_size**2],[WFS_size,WFS_size]);
#matplotlib.pyplot.imshow(temp.real)     #Show a WFS

def matvec(x):
    return Tomo_func.HH(x, sizes, alt_p, WFS_angles, 0)

def rmatvec(x):
    return Tomo_func.HH(x, sizes, alt_p, WFS_angles, 1)


A = scipy.sparse.linalg.LinearOperator( (8405,3890), matvec=matvec, rmatvec=rmatvec);

recon = scipy.sparse.linalg.lsqr(A,WFS);

layerXY_hat = recon[0];

temp = numpy.reshape(layersXY[0:sizes[0]**2],[sizes[0],sizes[0]])
matplotlib.pyplot.imshow(temp.real)     #Show a layer

temp2 = numpy.reshape(layerXY_hat[0:sizes[0]**2],[sizes[0],sizes[0]])
matplotlib.pyplot.imshow(temp2.real)     #Show a layer
