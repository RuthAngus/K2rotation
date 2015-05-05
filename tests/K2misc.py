import numpy as np
import matplotlib.pyplot as plt
import fitsio
import h5py
import glob
from K2pgram import K2pgram

def load_K2_data(fname):
    fname = str(int(fname))
    # load K2 light curve
    data = fitsio.read("../data/c1/ktwo%s-c01_lpd-lc.fits" % fname)
    aps = fitsio.read("../data/c1/ktwo%s-c01_lpd-lc.fits" % fname, 2)
    y = data["flux"][:, np.argmin(aps["cdpp6"])]
    x = data["time"]
    q = data["quality"]
    l = np.isfinite(y) * np.isfinite(x) * (q==0)
    y, x = y[l], x[l]
    y /= np.median(y)
    y -= 1
    x *= 24*3600

    # load basis
    with h5py.File("../data/c1.h5", "r") as f:
        basis = f["basis"][:150, l]

    return x, y, basis

def peak_detect(x, y):
    peaks = np.array([i for i in range(1, len(x)-1) if y[i-1] < y[i] and
                     y[i+1] < y[i]])
    l = y[peaks] == max(y[peaks])
    mx, my = x[peaks][l][0], y[peaks][l][0]
    return mx, my

def load_basis():
    # load basis
    with h5py.File("../data/c1.h5", "r") as f:
        basis = f["basis"][:150, l]
    return basis
