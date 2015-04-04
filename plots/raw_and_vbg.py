import numpy as np
import matplotlib.pyplot as plt
from gatspy.periodic import LombScargle
import fitsio

plotpar = {'axes.labelsize': 18,
           'text.fontsize': 18,
           'legend.fontsize': 18,
           'xtick.labelsize': 18,
           'ytick.labelsize': 18,
           'text.usetex': True}
plt.rcParams.update(plotpar)

# eid = "201545182"
eid = "201183188"
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
fig = plt.figure()
ax = fig.add_subplot(111)
ax1 = fig.add_subplot(211)
model = LombScargle().fit(x, y, np.ones_like(y)*1e-5)
period = 1. / fs
raw_pgram = model.periodogram(period)
ax1.plot(fs[::3]*1e6, raw_pgram[::3], "k")
ax.set_title("$\mathrm{EPIC~%s}$" % eid)
ax1.set_xlim(10, 280)
ax1.set_ylim(0, .015)

ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='w', top='off', bottom='off', left='off',
               right='off')

# load andrew's lcs
ax2 = fig.add_subplot(212)
x, y, _ = np.genfromtxt("../data/c1/ep%s.csv" % eid, delimiter=",").T
x *= 24*3600
model = LombScargle().fit(x, y, np.ones_like(y)*1e-5)
ps = 1. / fs
pgram = model.periodogram(ps)
ax2.plot(fs[::3]*1e6, pgram[::3], "k")
ax.set_xlabel("$\\nu\mathrm{~(}\mu\mathrm{Hz)}$")
ax2.set_xlim(10, 280)
ax2.set_ylim(0, .015)
fig.subplots_adjust(hspace=0)
fig.text(0.04, 0.5, "$\mathrm{Power}$", ha="center", va="center",
         rotation="vertical")
plt.savefig("../documents/rawvbg_%s.pdf" % eid)
