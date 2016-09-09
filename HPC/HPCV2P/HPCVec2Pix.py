#%reset
import sys
import scipy.io as sio
from scipy.stats import mode
import healpy as hp
import numpy as np
import itertools as it
import math

mat_contents = sio.loadmat('PointsTestMat'+sys.argv[2]+'.mat')
points = mat_contents['outputmat']
print(points)
hppoints = hp.pixelfunc.vec2pix(int(sys.argv[1]), points[:,0],points[:,1],points[:,2])

sio.savemat('hppoints'+sys.argv[2]+'.mat', {'hppoints':hppoints})