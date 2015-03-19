import numpy as np
import matplotlib.pyplot as plt
from gatspy.periodic import LombScargle
import fitsio

plotpar = {'axes.labelsize': 20,
           'text.fontsize': 35,
           'legend.fontsize': 20,
           'xtick.labelsize': 15,
           'ytick.labelsize': 15,
           'text.usetex': True}
plt.rcParams.update(plotpar)

eid = "201211472"
fname = "../data/c1/ktwo%s-c01_lpd-lc.fits" % eid

# load raw data
data = fitsio.read(fname)
aps = fitsio.read(fname, 2)
y = data["flux"][:, np.argmin(aps["cdpp6"])]
x = data["time"]
q = data["quality"]
l = np.isfinite(y) * np.isfinite(x) * (q==0)
y, x = y[l], x[l]
y /= np.median(y)
y -= 1
x *= 24*3600  # convert to seconds

fs = np.arange(.1, 300, 4e-2) * 1e-6  # astero

# plot raw data
plt.clf()
model = LombScargle().fit(x, y, np.ones_like(y)*1e-5)
period = 1. / fs
raw_pgram = model.periodogram(period)
plt.plot(fs*1e6, raw_pgram, "k")
plt.xlabel("$\\nu\mathrm{~(}\mu\mathrm{Hz)}$")
plt.ylabel("$\mathrm{Power}$")
plt.xlim(10, 280)
plt.ylim(0, .015)
plt.savefig("../documents/raw_%s.pdf" % eid)

# load andrew's lcs
x, y, _ = np.genfromtxt("../data/c1/ep201211472.csv", delimiter=",").T
x *= 24*3600
model = LombScargle().fit(x, y, np.ones_like(y)*1e-5)
ps = 1. / fs
pgram = model.periodogram(ps)
plt.clf()
plt.plot(fs*1e6, pgram, "k")
plt.xlabel("$\\nu\mathrm{~(}\mu\mathrm{Hz)}$")
plt.ylabel("$\mathrm{Power}$")
plt.xlim(10, 280)
plt.ylim(0, .015)
plt.savefig("vbg_%s.pdf" % eid)
