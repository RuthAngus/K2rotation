import numpy as np
from recover import inj
import sys

N = int(sys.argv[1])  # the number of injections
ifs = np.random.uniform(10e-6, 300e-6, N) * 1e-6  # injected freqs
a_s = 10**np.random.uniform(-4, -3, N)  # injection amps
fname = "data/ktwo201121245-c01-lpd-lc.fits"

inj(fname, ifs, a_s)
