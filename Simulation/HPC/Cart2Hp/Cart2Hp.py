#%reset
import sys
import scipy.io as sio
from scipy.stats import mode
import healpy as hp
import numpy as np
import itertools as it
import math

#Loading Data from Matlab
mat_contents = sio.loadmat('/Users/David/Macro-Impact/Simulation/Current/DataGeneration/Mat'+ sys.argv[1] +'.mat');
points = mat_contents['outputmat'];

nside = int(sys.argv[2]);

#Assigning each point its HEALPix coordinate
hppoints = hp.pixelfunc.vec2pix(nside, points[:,0],points[:,1],points[:,2])

sio.savemat('/Users/David/Macro-Impact/Simulation/Current/Cart2Hp/hppoints' + sys.argv[1] + '.mat', mdict = {'hppoints':hppoints})