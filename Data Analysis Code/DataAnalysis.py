#%reset
import sys
import scipy.io as sio
from scipy.stats import mode
import healpy as hp
import numpy as np
import itertools as it
import math



nside = int(sys.argv[0]);
R = 1.731E6; #Radius of Moon
sigmaX = 4.0E-11 #Cross-sectional area of macro
vX = 2.5E5 #Velocity of Macro
p0 = 1.0E8 #Hookean Pressure
frequencymax = 2.0E0
DisplacementThreshold = 3.0E-10; #How much does the detector need to be displaced in order to register the ray?
vsurf = 1.0E3; #Velocity of Rays at Surface

A = hp.nside2pixarea(nside) * R**2 #Area of pixel
T = ((A/np.pi)**(1/2))/vsurf; #How close in time do two rays have to be in order to be considered hitting together?
N = hp.nside2npix(nside);

#Loading Data from Matlab
mat_contents = sio.loadmat('/Users/David/Macro-Impact/Data/PointsTestMat1.mat');
points = mat_contents['outputmat']; #Points has the format (x,y,z,time,density,REF,ATT,kappa,velocity)
num_contents = sio.loadmat('/Users/David/Macro-Impact/Data/PointsTestNum1.mat');
numbers = num_contents['outputnum']; #Points has the format (M,N,L,D)
length = len(points);

#Assigning each point its HEALPix coordinate
hppoints = hp.pixelfunc.vec2pix(nside, points[:,0],points[:,1],points[:,2])
#Mode = mode(hppoints)[1][0]; #Mode is the maximum number of times a pixel is struck by a ray

PixelTimes = np.zeros((N,len(hppoints))); #Each row contains a list of times at which rays struck the nth pixel
Density = np.zeros((N,len(hppoints))); #This is the corresponding list of energies
PixelREF = np.zeros((N,len(hppoints))); #Refraction/Reflection Coefficient of each ray
#PixelATT = np.zeros((N,len(hppoints))); #Quality Factor (ish) of each ray
Velocity = np.zeros((N,len(hppoints))); #Velocity at source
kmax = np.zeros((N,len(hppoints))); #Max wavenumber at source
Xi = np.zeros((N,len(hppoints))); #Xi at source
PixelDisplacements = np.zeros((N,len(hppoints))); #Maximum Displacements measured in a pixel

m = numbers[0][0]
n = numbers[0][1]
L = numbers[0][2]
D = numbers[0][3]

#Now we assign PixelTimes and PixelEnergies their values
for i in range(0,N - 1):
    counter = 0;
    for j in range(0,length - 1):
        if hppoints[j] == i:
            PixelTimes[i][counter] = points[j,3];
            Density[i][counter] = points[j,4];
            PixelREF[i][counter] = points[j,5];
            #PixelATT[i][counter] = points[j,6];
            Velocity[i][counter] = points[j,8];
            kmax[i][counter] = frequencymax * 2 * np.pi * (points[j,8]**(-1));
            Xi[i][counter] = points[j,7] * points[j,4] * (points[j,8]**(-2)) * (9 * sigmaX * (vX**2) * ((frequencymax * 2 * np.pi) **2))/(16 * np.pi * p0**2);
            counter = counter + 1;
#id = np.where(PixelATT[:] > 0)
#ATT = np.zeros(N)
#ATT[:] = np.mean(PixelATT[:][id])


#New but somehow slower way
#for i in range(0,N - 1):            
#    PixelTimes[i][0 : sum(hppoints == i)] = points[hppoints == i,3]
#    PixelEnergies[i][0 : sum(hppoints == i)] = points[hppoints == i,4]

#New improved code by Prof. Copi
#Now we add together those energies whose times are close together
for i in range(0,N - 1):
    nz = np.nonzero(PixelTimes[i])[0]
    if len(nz) < 2 : continue
    ind = np.array(list(it.combinations(range(len(nz)),2)))
    d = np.abs(PixelTimes[i,ind[:,0]] - PixelTimes[i,ind[:,1]])
    for ii in np.where((d < T) & (d > 0))[0]:
        #PixelDisplacements[i,ind[ii,0]] += ((2 / (Density[i,ind[ii,1]]*(Velocity[i,ind[ii,1]]**2) * A))**(1/2))*((PixelREF[i,ind[ii,1]]*Xi[i,ind[ii,1]]*(Density[i,ind[ii,1]]*sigmaX*(vX**2) * L / (m*n))/ATT[i,ind[ii,1]])**(1/2))*(kmax[i,ind[ii,1]]**(-1))*math.erf(((PixelATT[i,ind[ii,1]]*kmax[i,ind[ii,1]])/2)**(1/2))
        PixelDisplacements[i,ind[ii,0]] += (4*PixelREF[i,ind[ii,1]]*Xi[i,ind[ii,1]]*sigmaX*(vX**2)*L)/(m*n*np.pi*A*kmax[i,ind[ii,1]]*Velocity[i,ind[ii,1]]**2)
#Old, slow code, by me.
#for i in range(0,N - 1):
#    for j in range(0,Mode - 2):
#        for k in range(j + 1,Mode - 1):
#            if abs(PixelTimes[i][j] - PixelTimes[i][k]) < T:
#                PixelEnergies[i][j] = PixelEnergies[i][j]+PixelEnergies[i][k];

#Now we find the maximum energy incident on each pixel
Displacement = (np.amax(PixelDisplacements, axis=1))**(1/2);

#Lastly, we talley up the number of pixels with energye exceeding the threshold to obtain the fraction of the moon for which the impact is detectable
F = len(np.where(Displacement>DisplacementThreshold)[0])
print(F/N)
