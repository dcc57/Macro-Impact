
# coding: utf-8

# In[4]:

import scipy.io as sio
from scipy.stats import mode
import healpy as hp
import numpy as np
import itertools as it
import math

nside = 16;
R = 1.731E6; #Radius of Moon
AverageDensity = 3.0E3
AverageVelocity = 8.0E3
frequencymax = 2.0E0
DisplacementThreshold = 3.0E-10; #How much does the detector need to be displaced in order to register the ray?

rho = AverageDensity
v = AverageVelocity
D = DisplacementThreshold

A = hp.nside2pixarea(nside) * R**2 #Area of pixel
T = (A**(1/2))/v; #How close in time do two rays have to be in order to be considered hitting together?
kmax = frequencymax * 2 * np.pi * (v**(-1));
N = hp.nside2npix(nside);
Energy = np.zeros(N);

#Loading Data from Matlab
mat_contents = sio.loadmat('/Users/David/Macro-Impact/Ray-Tracing/Points.mat');
points = mat_contents['output']; #Points has the format (x,y,z,time,energy,vTQ)
length = len(points);

#Assigning each point its HEALPix coordinate
hppoints = hp.pixelfunc.vec2pix(nside, points[:,0],points[:,1],points[:,2])
#Mode = mode(hppoints)[1][0]; #Mode is the maximum number of times a pixel is struck by a ray

PixelTimes = np.zeros((N,len(hppoints))); #Each row contains a list of times at which rays struck the nth pixel
PixelEnergies = np.zeros((N,len(hppoints))); #This is the corresponding list of energies
PixelvTQ = np.zeros((N,len(hppoints))); #Quality Factor (ish) of each ray
PixelDisplacements = np.zeros((N,len(hppoints))); #Maximum Displacements measured in a pixel

#Now we assign PixelTimes and PixelEnergies their values
for i in range(0,N - 1):
    counter = 0;
    for j in range(0,length - 1):
        if hppoints[j] == i:
            PixelTimes[i][counter] = points[j,3];
            PixelEnergies[i][counter] = points[j,4];
            PixelvTQ[i][counter] = points[j,5];
            counter = counter + 1;
#New but somehow slower way
#for i in range(0,N - 1):            
#    PixelTimes[i][0 : sum(hppoints == i)] = points[hppoints == i,3]
#    PixelEnergies[i][0 : sum(hppoints == i)] = points[hppoints == i,4]

#New improved code by Prof. Copi
#Now we add together those energies whose times are close together
coeff = (N**(1/2))*((2*np.pi*kmax*v*(rho**(1/2)))**(-1))
for i in range(0,N - 1):
    nz = np.nonzero(PixelTimes[i])[0]
    if len(nz) < 2 : continue
    ind = np.array(list(it.combinations(range(len(nz)),2)))
    d = np.abs(PixelTimes[i,ind[:,0]] - PixelTimes[i,ind[:,1]])
    for ii in np.where(d < T)[0]:
        PixelDisplacements[i,ind[ii,0]] += coeff*(((PixelEnergies[i,ind[ii,1]])**(1/2))*((PixelvTQ[i,ind[ii,1]])**(-1))*(((2*np.pi*((PixelvTQ[i,ind[ii,1]])**(-1)))**(1/2)) * math.erf(((kmax/2)*PixelvTQ[i,ind[ii,1]])**(1/2)) - 2 * (kmax**(1/2)) * math.exp(-(kmax/2)*PixelvTQ[i,ind[ii,1]])))
#Old, slow code, by me.
#for i in range(0,N - 1):
#    for j in range(0,Mode - 2):
#        for k in range(j + 1,Mode - 1):
#            if abs(PixelTimes[i][j] - PixelTimes[i][k]) < T:
#                PixelEnergies[i][j] = PixelEnergies[i][j]+PixelEnergies[i][k];

#Now we find the maximum energy incident on each pixel
Displacement = np.amax(PixelDisplacements, axis=1);

#Lastly, we talley up the number of pixels with energye exceeding the threshold to obtain the fraction of the moon for which the impact is detectable
F = len(np.where(Displacement>D)[0])
print(F/N)

    
            


# In[478]:

np.mean(Displacement)


# In[472]:

((kmax/2)*PixelvTQ[1,1])


# In[414]:

kmax


# In[436]:

math.erf(1)


# In[391]:




# In[386]:

np.pi


# In[ ]:



