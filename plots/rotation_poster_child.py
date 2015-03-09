# This script creates the rotation_poster_child.pdf figure

import numpy as np
import matplotlib.pyplot as plt
import emcee
from gatspy.periodic import LombScargle

plotpar = {'axes.labelsize': 12,
           'text.fontsize': 10,
           'legend.fontsize': 10,
           'xtick.labelsize': 8,
           'ytick.labelsize': 8,
           'text.usetex': True}
plt.rcParams.update(plotpar)

def max_peak_detect(x, y):
    peaks = np.array([i for i in range(1, len(x)-1) if y[i-1] < y[i] and
                     y[i+1] < y[i]])
    l = y[peaks] == max(y[peaks])
    return x[peaks][l][0], y[peaks][l][0]

# load poster child light curve
x, y, _ = np.genfromtxt("../data/ep201317002.csv", delimiter=",").T
y = y/np.median(y) - 1

# compute acf
dt = 0.02043359821692  # time resolution of K2 data (from header)
acf = emcee.autocorr.function(y)
lags = np.arange(len(acf)) * dt

# compute LS periodogram
model = LombScargle().fit(x, y, np.ones_like(y)*1e-5)
fmax = max(lags)/100.
fs = np.linspace(1e-6, fmax, 1000)
ps = 1./fs
pgram = model.periodogram(ps)

plt.clf()
plt.subplot(3, 1, 1)
l = x < 2016
plt.plot(x[l], y[l], "k")
plt.plot(x[~l], y[~l], "k")
plt.xlabel("$\mathrm{BJD-2454833~(days)}$")
plt.ylabel("$\mathrm{Normalized~Flux}$")
plt.xlim(min(x), max(x))
plt.ylim(-.02, .02)

plt.subplot(3, 1, 2)
plt.plot(lags, acf, "k")
plt.ylabel("$\mathrm{Autocorrelation}$")
plt.xlabel("$\mathrm{Time~(days)}$")
mx, my = max_peak_detect(lags, acf)
plt.axvline(mx, color=".5", linestyle="--", label="$P_{rot}=%.2f$" % mx)
plt.legend()
plt.xlim(min(lags), max(lags))

plt.subplot(3, 1, 3)
plt.plot(fs, pgram, "k")
plt.xlabel("$\mathrm{Frequency~(days}^{-1}\mathrm{)}$")
plt.ylabel("$\mathrm{Power}$")
mx, my = max_peak_detect(fs, pgram)
plt.xlim(min(fs), fmax)
px = 1./mx
plt.axvline(mx, color=".5", linestyle="--", label="$P_{rot}=%.2f$" % px)
plt.legend()
plt.subplots_adjust(hspace=.4)
plt.savefig("../documents/rotation_poster_child.pdf")
