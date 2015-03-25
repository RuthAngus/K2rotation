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

def histo(amps, Ps, namps, npers, fname):
    amps = np.log(amps)
    nbins = namps
    max_n_per_bin = int(npers/namps)
    Ps = np.log(Ps)
    my_yedges = np.linspace(min(amps), max(amps), nbins)
    my_yedges = np.linspace(min(amps), max(amps), nbins)
    my_xedges = np.linspace(min(Ps), max(Ps), nbins)
    K2_hist, xedges, yedges = np.histogram2d(Ps, amps, bins=nbins,
                                             range=[[min(my_xedges),
                                                    max(my_xedges)],
                                                    [min(my_yedges),
                                                     max(my_yedges)]])

    fig = plt.figure()
    ax = fig.add_subplot(111)
    X, Y = np.meshgrid(xedges, yedges)
    cax = ax.pcolormesh(X, Y, K2_hist.T/float(max_n_per_bin)*100, cmap="Blues")
    ax.set_yticks(np.exp(np.arange(-7, -2)))
    ticks = np.arange(Y.min(), Y.max(), 6)
    labels = range(-7, -2)
    plt.xticks(ticks, labels)
    ax.set_ylabel("$\mathrm{Amplitude~(\%)}$")
    ax.set_xlabel("$\ln\mathrm{Period~(days)}$")
    plt.colorbar(cax, label="$\mathrm{Completeness~(\%)}$")
#     plt.plot(Ps, amps, "k.")
    plt.savefig("../injections/sine/%s_hist_%s.pdf" % (fname, flag))
    plt.close(fig)
#     print K2_hist.T
    return K2_hist, xedges, yedges

def make_histogram_plot(flag, namps=11, npers=1000):
    hist_names = glob.glob("../injections/sine/histogram_*_*_%s.h5" % flag)
#     print hist_names
    with h5py.File(hist_names[0], "r") as f:
        K2_amps = f["K2"][:, 0] * 100  # convert to %
        K2_Ps = f["K2"][:, 1]
        raw_amps = f["raw"][:, 0] * 100  # convert to %
        raw_Ps = f["raw"][:, 1]
    for i in range(1, len(hist_names)):
        with h5py.File(hist_names[0], "r") as f:
            K2_amps = np.concatenate((K2_amps, f["K2"][:, 0] * 100))
            K2_Ps = np.concatenate((K2_Ps, f["K2"][:, 1]))
            raw_amps = np.concatenate((raw_amps, f["raw"][:, 0] * 100))
            raw_Ps = np.concatenate((raw_Ps, f["raw"][:, 1]))

    if flag == "r":  # cut off short p_rots
        l = K2_Ps > np.exp(0.5)
        K2_amps, K2_Ps = K2_amps[l], K2_Ps[l]
        l = raw_Ps > np.exp(1)
        raw_amps, raw_Ps = raw_amps[l], raw_Ps[l]

#     K2_amps = np.log(K2_amps)

    K2_hist, xedges, yedges = histo(K2_amps, K2_Ps, namps, npers, "K2")
    raw_hist, _, _ = histo(raw_amps, raw_Ps, namps, npers, "raw")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    X, Y = np.meshgrid(xedges, yedges)
    b = K2_hist.T-raw_hist.T
    l = b < 0
    b[l] = 0
#     print b
    cax = ax.pcolormesh(X, Y, b, cmap="Blues")
    ax.set_ylabel("$\mathrm{Amplitude~(\%)}$")
    ax.set_xlabel("$\ln\mathrm{(Period)~(days)}$")
    plt.colorbar(cax, label="$\mathrm{Completeness~(\%)}$")
#     plt.plot(raw_Ps, raw_amps, "k.")
    plt.savefig("../injections/sine/both_hist_%s" % flag)
    plt.close(fig)

if __name__ == "__main__":
    make_histogram_plot("r")
    make_histogram_plot("a")
