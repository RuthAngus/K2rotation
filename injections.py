import numpy as np
import matplotlib.pyplot as plt
import glob
import pyfits
from mklc import mklc
import scipy.signal as sps
import fitsio
import h5py
from K2pgram import K2pgram
from colours import plot_colours
cols = plot_colours()
from gatspy.periodic import LombScargle

plotpar = {'axes.labelsize': 10,
           'text.fontsize': 10,
           'legend.fontsize': 10,
           'xtick.labelsize': 8,
           'ytick.labelsize': 8,
           'text.usetex': True}
plt.rcParams.update(plotpar)

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

# generate injection
nspot = 200
incl = np.pi*5./12.
amp = 1.
tau = 30.5
p = 10.0,
res0, res1 = mklc(raw_x)
time = res1[0]
flux = res1[2]/np.median(res1[2]) - 1

# add lcs together
a = .01
flux *= a
y = flux + raw_y

SN = np.var(flux) / np.var(raw_y)
print "signal = ", np.var(flux)
print "noise = ", np.var(raw_y)
print "S/N = ", SN

# load basis
with h5py.File("data/c1.h5", "r") as f:
    basis = f["basis"][:10, l]

# construct arrays
AT = np.concatenate((basis, np.ones((3, len(y)))), axis=0)
ATA = np.dot(AT, AT.T)

# calculate K2pgram
print "calculating K2pgram..."
fs = np.arange(1./20, 1./2, .001)
_K2pgram = K2pgram(raw_x, y, fs, AT, ATA)

# calculate periodogram of raw light curve
print "calculating raw pgram..."
y = np.array([i.astype("float64") for i in y])
raw_x = np.array([i.astype("float64") for i in raw_x])
pgram = sps.lombscargle(raw_x, y, 2*np.pi*fs)

# plot
plt.clf()
plt.subplot(3, 1, 1)
plt.plot(raw_x, raw_y, "k.")
plt.plot(raw_x, flux, color=cols.green)
plt.subplot(3, 1, 2)
plt.plot(raw_x, y, "k.")
plt.subplot(3, 1, 3)
plt.axvline(1./10, color=".7", linestyle="--")
plt.plot(fs, _K2pgram/max(_K2pgram), color=cols.blue,
         label="$\mathrm{K2pgram$}")
plt.plot(fs, pgram/max(pgram), color=cols.pink, label="$\mathrm{raw pgram}$")
plt.savefig("%s_result" % fname)
