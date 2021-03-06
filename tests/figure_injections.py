import numpy as np
import matplotlib.pyplot as plt
import glob
import pyfits
from mklc import mklc
import scipy.signal as sps
import fitsio
import h5py
from K2pgram import K2pgram, eval_freq
from gatspy.periodic import LombScargle
import glob
import sys

plotpar = {'axes.labelsize': 20,
           'text.fontsize': 20,
           'legend.fontsize': 20,
           'xtick.labelsize': 20,
           'ytick.labelsize': 20,
           'text.usetex': True}
plt.rcParams.update(plotpar)

def period2amp(period):
    ps = 10**period * 24 * 3600
    fHz = 1/ps
    fuHz = fHz*1e6
    amp = 860 * fuHz**-.71
    return np.log10(amp)

def histo(amps, Ps, namps, npers, allas, allps, fname, flag):
    amps *= 1e4  # convert to ppm
    allas *= 1e4
    nbins = namps
    Ps = np.log10(Ps)
    allps = np.log10(allps)
    amps = np.log10(amps)
    allas = np.log10(allas)
    max_n_per_bin = int(npers/namps)
    my_yedges = np.linspace(min(amps), max(amps), nbins)
    my_xedges = np.linspace(min(Ps), max(Ps), nbins)
    K2_hist, xedges, yedges = np.histogram2d(Ps, amps, bins=nbins,
                                             range=[[min(my_xedges),
                                                    max(my_xedges)],
                                                    [min(my_yedges),
                                                     max(my_yedges)]])
    truth_hist, xedges, yedges = np.histogram2d(allps, allas, bins=nbins,
                                             range=[[min(my_xedges),
                                                    max(my_xedges)],
                                                    [min(my_yedges),
                                                     max(my_yedges)]])

    fig = plt.figure()
    ax = fig.add_subplot(111)
    X, Y = np.meshgrid(xedges, yedges)
    color = K2_hist.T/truth_hist.T * 100
    color[color==np.inf] = 0
    cax = ax.pcolormesh(X, Y, color, cmap="Blues")
    ax.set_ylabel("$\log_{10}\mathrm{Amplitude~(ppm)}$")
    ax.set_xlabel("$\log_{10}\mathrm{Period~(days)}$")
    plt.colorbar(cax, label="$\mathrm{Completeness~(\%)}$")
    plt.subplots_adjust(bottom=.2)
    plt.plot(np.sort(Ps), period2amp(np.sort(Ps)), color="#FF33CC",
             linestyle="--")
    plt.ylim(1, 3)
#     plt.plot(Ps, amps, "k.")
    plt.plot(allps, allas, "r.", markersize=.5)
    plt.savefig("../injections/sine/%s_hist_%s.pdf" % (fname, flag))
    plt.close(fig)
#     print K2_hist.T
    return K2_hist.T/truth_hist.T * 100, xedges, yedges

def make_histogram_plot(flag, namps=30, npers=1000):

    hist_names = glob.glob("../injections/sine/histogram_*_*_%s.h5" % flag)
    truth_names = glob.glob("../injections/sine/truths_*_*_%s.h5" % flag)
#     print hist_names
    with h5py.File(hist_names[0], "r") as f:
        K2_amps = f["K2"][:, 0] * 100  # convert to %
        K2_Ps = f["K2"][:, 1]
        raw_amps = f["raw"][:, 0] * 100  # convert to %
        raw_Ps = f["raw"][:, 1]
    print min(K2_Ps), max(K2_Ps)
    for i in range(1, len(hist_names)):
        with h5py.File(hist_names[i], "r") as f:
            K2_amps = np.concatenate((K2_amps, f["K2"][:, 0] * 100))
            K2_Ps = np.concatenate((K2_Ps, f["K2"][:, 1]))
            raw_amps = np.concatenate((raw_amps, f["raw"][:, 0] * 100))
            raw_Ps = np.concatenate((raw_Ps, f["raw"][:, 1]))
    with h5py.File(truth_names[0], "r") as f:
        allas = f["K2"][:, 0] * 100  # convert to %
        allps = f["K2"][:, 1]
    for i in range(1, len(truth_names)):
        with h5py.File(truth_names[i], "r") as f:
            allas = np.concatenate((allas, f["K2"][:, 0] * 100))
            allps = np.concatenate((allps, f["K2"][:, 1]))
    print min(allps), max(allps)

    if flag == "r":  # cut off short p_rots
        l = K2_Ps > 10**(.3)
        K2_amps, K2_Ps = K2_amps[l], K2_Ps[l]
        l = raw_Ps > 10**(1)
        raw_amps, raw_Ps = raw_amps[l], raw_Ps[l]

#     K2_amps = np.log(K2_amps)

    K2_hist, xedges, yedges = histo(K2_amps, K2_Ps, namps, npers, allas, allps,
                                    "K2", flag)
    raw_hist, _, _ = histo(raw_amps, raw_Ps, namps, npers, allas, allps,
                           "raw", flag)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    X, Y = np.meshgrid(xedges, yedges)
    b = K2_hist-raw_hist
    l = b < 0
    b[l] = 0
#     print b
    cax = ax.pcolormesh(X, Y, b, cmap="Blues")
    ax.set_ylabel("$\log_{10}\mathrm{Amplitude~(ppm)}$")
    ax.set_xlabel("$\log_{10}\mathrm{(Period)~(days)}$")
    plt.colorbar(cax, label="$\mathrm{Completeness~(\%)}$")
#     plt.plot(raw_Ps, raw_amps, "k.")
    plt.savefig("../injections/sine/both_hist_%s" % flag)
    plt.close(fig)

if __name__ == "__main__":
    make_histogram_plot("a")
    make_histogram_plot("r")
