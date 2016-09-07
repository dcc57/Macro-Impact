#%reset
import sys
import scipy.io as sio
from scipy.stats import mode
import healpy as hp
import numpy as np
import itertools as it
import math

nside = int(sys.argv[1]) #EG 4
R = float(sys.argv[2]); #Radius of Moon EG 1.737E6
vsurf = float(sys.argv[3]); #Velocity of Rays at Surface EG 1.0E3
frequencymax = float(sys.argv[4]); #Maximum Detectable Frequency EG 2.0E1

A = hp.nside2pixarea(nside) * R**2 #Area of pixel
T = ((A/np.pi)**(1/2))/vsurf; #How close in time do two rays have to be in order to be considered hitting together?
N = hp.nside2npix(nside);

#Loading Data from Matlab
mat_contents = sio.loadmat('/Users/David/Macro-Impact/Data/PointsTestMat'+sys.argv[5]+'.mat');
points = mat_contents['outputmat']; #Points has the format (x,y,z,time,density,REF,ATT,kappa,velocity)
num_contents = sio.loadmat('/Users/David/Macro-Impact/Data/PointsTestNum'+sys.argv[5]+'.mat');
numbers = num_contents['outputnum']; #Points has the format (M,N,L,D)
length = len(points);

hppoints = hp.pixelfunc.vec2pix(nside, points[:,0],points[:,1],points[:,2]) #Assigning each point its HEALPix coordinate
PixelTimes = np.zeros((N,len(hppoints))); #Each row contains a list of times at which rays struck the nth pixel
PreDisplacementPrime = np.zeros((N,len(hppoints))); #Intermediate array to store ray initial condition dependent information
PixelDisplacements = np.zeros((N,len(hppoints))); #Maximum Displacements measured in a pixel

m = numbers[0][0]
n = numbers[0][1]
L = numbers[0][2]
D = numbers[0][3]

#Now we assign PixelTimes and PixelEnergies their values
for i in range(0,N):
    counter = 0;
    for j in range(0,length):
        if hppoints[j] == i:
            PixelTimes[i][counter] = points[j,3];
            PreDisplacementPrime[i][counter] = (points[j,5] * points[j,7] * points[j,4])/(points[j,8]**3)
            counter = counter + 1;
#New improved code by Prof. Copi
#Now we add together those energies whose times are close together
for i in range(0,N):
    nz = np.nonzero(PixelTimes[i])[0]
    if len(nz) < 2 : continue
    ind = np.array(list(it.combinations(range(len(nz)),2)))
    d = np.abs(PixelTimes[i,ind[:,0]] - PixelTimes[i,ind[:,1]])
    for ii in np.where((d < T) & (d > 0))[0]:
        PixelDisplacements[i,ind[ii,0]] += PreDisplacementPrime[i,ind[ii,1]]

#Now we find the maximum displacement factor on each pixel
PreDisplacement = np.amax(PixelDisplacements, axis=1)**(1/2) * ((9 * L) / (2 * np.pi * n * m * A))**(1/2); #Multiply by fmax**(1/2) * sigmaX * (vX**2) * p0**(-1) to get displacement

#We save the max displacements in each pixel (up to position independent info i.e. fmax**(1/2) * sigmaX * (vX**2) * p0**(-1))
sio.savemat('PreDisplacement'+sys.argv[5]+'.mat', {'PreDisplacement':PreDisplacement})