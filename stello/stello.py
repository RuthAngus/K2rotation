import numpy as np
import matplotlib.pyplot as plt
import wget
import h5py

epids = np.genfromtxt("ktwo_c1_APO-RGs_llc.dat.epic.list", dtype=str).T

for id in epids:
    a = "%s00000" % id[:4]
    b = "%s000" % id[4:6]
    c = "ktwo%s-c01_lpd-lc.fits" % id
    wget.download("http://bbq.dfm.io/ketu/lightcurves/c1/%s/%s/%s" % (a, b, c))
