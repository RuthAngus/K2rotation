import numpy as np
from recover import inj
import sys
from K2misc import load_K2_data

N = int(sys.argv[1])  # the number of injections
rotation = False

ifs = np.random.uniform(10e-6, 270e-6, N)  # injected freqs (astero)
if rotation:
    ifs = 1./(np.random.uniform(1, 30, N) * 24 * 3600)  # injected freqs (rot)
a_s = 10**np.random.uniform(-7, -2, N)  # injection amps
fname = "data/ktwo201121245-c01_lpd-lc.fits"
x, y, basis = load_K2_data(fname)

inj(fname, ifs, a_s, rotation=rotation)
