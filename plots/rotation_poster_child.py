# This script creates the rotation_poster_child.pdf figure

import numpy as np
import matplotlib.pyplot as plt
import emcee
from gatspy.periodic import LombScargle
import scipy.signal as sps
import subprocess

def acf_peak_detect(x, y):
    x, y = x[50:], y[50:]
    peaks = np.array([i for i in range(1, len(x)-1) if y[i-1] < y[i] and
                     y[i+1] < y[i]])
    if len(peaks):  # make sure that a peak is detected
        l = y[peaks] == max(y[peaks])
        mx, my = x[peaks][l][0], y[peaks][l][0]
        if my > 0:
            return mx, my
        else:
            return 0, 0
    else:
        return 0, 0

def max_peak_detect(x, y):
    peaks = np.array([i for i in range(1, len(x)-1) if y[i-1] < y[i] and
                     y[i+1] < y[i]])
    if len(peaks):  # make sure that a peak is detected
        l = y[peaks] == max(y[peaks])
        return x[peaks][l][0], y[peaks][l][0]
    else:
        return 0, 0

# load poster child light curve
# x, y, _ = np.genfromtxt("../data/ep201317002.csv", delimiter=",").T

def p_child_plot(x, y, eid):

    plotpar = {'axes.labelsize': 15,
               'text.fontsize': 15,
               'legend.fontsize': 15,
               'xtick.labelsize': 12,
               'ytick.labelsize': 12,
               'text.usetex': True}
    plt.rcParams.update(plotpar)

    y = y/np.median(y) - 1

    # SUBTRACT LINEAR TREND!
    plv = np.polyfit(x, y, 1)
    print plv
    poly = plv[1] + x * plv[0]
    y -= poly

    # compute acf
    dt = 0.02043359821692  # time resolution of K2 data (from header)
    acf = emcee.autocorr.function(y)
    lags = np.arange(len(acf)) * dt

    # smooth acf
    acf = sps.savgol_filter(acf, 15, 1)

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
    plt.title("$\mathrm{EPIC~%s}$" % eid)

#     plt.ylim(-.02, .02)

    plt.subplot(3, 1, 2)
    plt.plot(lags, acf, "k")
    plt.ylabel("$\mathrm{Autocorrelation}$")
    plt.xlabel("$\mathrm{Time~(days)}$")
    acfx, acfy = acf_peak_detect(lags, acf)
    plt.axvline(acfx, color=".5", linestyle="--", label="$P_{rot}=%.2f$"
                % acfx)
    plt.legend()
    plt.xlim(min(lags), max(lags))

    plt.subplot(3, 1, 3)
#     plt.plot(fs, pgram, "k")
#     plt.xlim(min(fs), fmax)
    plt.plot(1./fs, pgram, "k")
    plt.xlim(min(1./fs), max(lags))
    plt.xlabel("$\mathrm{Period~(days)}$")
#     plt.xlabel("$\mathrm{Frequency~(days}^{-1}\mathrm{)}$")
    plt.ylabel("$\mathrm{Power}$")
    l = fs > 1./100
    mx, my = max_peak_detect(fs[l], pgram[l])
    px = 1./mx
#     plt.axvline(mx, color=".5", linestyle="--", label="$P_{rot}=%.2f$" % px)
    plt.axvline(px, color=".5", linestyle="--", label="$P_{max}=%.2f$" % px)
#     plt.ylim(0, .6)
    plt.legend(loc="best")
    plt.subplots_adjust(hspace=.4)
#     plt.savefig("vbg_%s" % eid)
#     plt.savefig("../documents/rotation_poster_child.pdf")
    print "../documents/rotation%s.pdf" % eid
    plt.savefig("../documents/rotation%s.pdf" % eid)
    return acfx, px

if __name__ == "__main__":

#     eid = "201317002"
#     eid = "201129544"
#     eid = "201132518"
    eids = [201129544, 201132518, 201133037, 201133147, 201135311, 201138638,
            201138849, 201142023, 201142127]
#     eids = [201133037, 201132518]
    eids = [201133037]
#     eids = [201142127]
    for eid in eids:
        try:
            x, y, _ = np.genfromtxt("../data/c1/ep%s.csv" % eid, delimiter=",").T
        except:
            print "copying data"
            print "cp /Users/angusr/data/K2/c1lcsr4/ep%s.csv ../data/c1" % eid
            subprocess.call("cp /Users/angusr/data/K2/c1lcsr4/ep%s.csv ../data/c1"
                            % eid, shell=True)
            x, y, _ = np.genfromtxt("../data/c1/ep%s.csv" % eid, delimiter=",").T
        p_child_plot(x, y, eid)
