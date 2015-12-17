import numpy as np
from recover import inj
import sys
from K2misc import load_K2_data

N = int(sys.argv[1])  # the number of injections
rotation = False

ifs = np.random.uniform(10e-6, 260e-6, N)  # injected freqs (astero)
if rotation:
    ifs = 1./(np.random.uniform(1, 30, N) * 24 * 3600)  # injected freqs (rot)
a_s = 10**np.random.uniform(-7, -3, N)  # injection amps
print(a_s[:10])

fname = "data/ktwo201121245-c01_lpd-lc.fits"
x, y, basis, med = load_K2_data(fname)
# the light curve has been median normalised, so what used to be units of flux
# is now units of flux / median
a_s = 10**np.random.uniform(0, 3, N) * 1e-6 # convert to ppm
print(a_s[:10])

inj(fname, ifs, a_s, rotation=rotation)
