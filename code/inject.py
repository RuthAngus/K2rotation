import numpy as np
from recover import inj
import sys
from K2misc import load_K2_data

N = int(sys.argv[1])  # the number of injections
ifs = np.random.uniform(10e-6, 280e-6, N)  # injected freqs
a_s = 10**np.random.uniform(-4, -3, N)  # injection amps
fname = "data/ktwo201121245-c01_lpd-lc.fits"
x, y, basis = load_K2_data(fname)

inj(fname, ifs, a_s)
