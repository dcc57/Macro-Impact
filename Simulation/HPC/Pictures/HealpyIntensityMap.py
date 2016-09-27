import sys
import scipy.io as sio
from scipy.stats import mode
import healpy as hp
import numpy as np
import itertools as it
import math

#Loading Data from Matlab
mat_contents = sio.loadmat('/Users/David/Macro-Impact/Simulation/Current/Results/Disp'+ sys.argv[1] +'.mat');
Disp = mat_contents['Disp'];

nside = int(sys.argv[2]);

Intensities = hp.read_map(Disp)

hp.mollview(Intensities, coord['G','E'])