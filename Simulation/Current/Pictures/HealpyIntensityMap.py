import sys
import scipy.io as sio
from scipy.stats import mode
import healpy as hp
import numpy as np
import itertools as it
import math

#Loading Data from Matlab
mat_contents = sio.loadmat('/Users/David/Macro-Impact/Simulation/Current/DataAnalysis/Mat'+ sys.argv[1] +'.mat');
points = mat_contents['outputmat'];

nside = int(sys.argv[2]);