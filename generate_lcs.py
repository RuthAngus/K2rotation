import numpy as np
import matplotlib.pyplot as plt
import glob
import pyfits
from mklc import mklc
import fitsio

fname = 201300080

# load K2 light curve
data = fitsio.read("data/ktwo%s-c01_lpd-lc.fits" % fname)
aps = fitsio.read("data/ktwo%s-c01_lpd-lc.fits" % fname, 2)
raw_y = data["flux"][:, np.argmin(aps["cdpp6"])]
raw_x = data["time"]
q = data["quality"]
l = np.isfinite(raw_y) * np.isfinite(raw_x) * (q==0)
raw_y, raw_x = raw_y[l], raw_x[l]
raw_y /= np.median(raw_y)
raw_y -= 1

periods = np.exp(np.random.uniform(np.log(1), np.log(70), 100))
for i, p in enumerate(periods):

    # generate injection
    nspot = 200
    incl = np.pi*5./12.
    amp = 1.
    tau = 30.5
    res0, res1 = mklc(raw_x, nspot, incl, amp, tau, p)
    time = res1[0]
    flux = res1[2]/np.median(res1[2]) - 1
    np.savetxt("injections/%s_lc.txt" % i, np.transpose((time, flux)))

np.savetxt("truth.txt", np.transpose((range(len(periods)), periods)))
