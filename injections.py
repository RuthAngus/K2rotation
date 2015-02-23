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
import glob

plotpar = {'axes.labelsize': 10,
           'text.fontsize': 10,
           'legend.fontsize': 10,
           'xtick.labelsize': 8,
           'ytick.labelsize': 8,
           'text.usetex': True}
plt.rcParams.update(plotpar)

# grid over amplitudes (the K2 pgram step takes the time)
def grid_over_amps(basis, flux, raw_x, raw_y, truth, fs, amps, true_a,
                   plot=False):

    K2P, rawP, K2a, rawa = [], [], [], []
    for i, a in enumerate(amps):

        tf = 1./truth
        print "period = ", truth
        print "amplitude = ", a
        print "frequency = ", tf

        # add lcs together
        fx = flux * a
        y = fx + raw_y
        SN = np.var(fx) / np.var(raw_y)

        # construct arrays
        AT = np.concatenate((basis, np.ones((3, len(y)))), axis=0)
        ATA = np.dot(AT, AT.T)

        # calculate K2pgram
        _K2pgram = K2pgram(raw_x, y, fs, AT, ATA)
        l = _K2pgram==max(_K2pgram)
        print "recovered frequency", fs[l][0]
        s = 0
        if tf-.1*tf < fs[l][0] and fs[l][0] < tf+.1*tf:
            K2P.append(truth)
            K2a.append(a)
            print "success!", "\n"
            s = 1
            K2_period = fs[l][0]

        # calculate periodogram of raw light curve
        y = np.array([_y.astype("float64") for _y in y])
        raw_x = np.array([_raw_x.astype("float64") for _raw_x in raw_x])
        model = LombScargle().fit(raw_x, y, np.ones_like(y)*1e-5)
        period = 1. / fs
        pgram = model.periodogram(period)
        l = pgram==max(pgram)
        if truth-.1*truth < fs[l][0] and fs[l][0] < truth+.1*truth:
            rawP.append(truth)
            rawa.append(a)

        if plot == True:
            plt.clf()
            plt.subplot(3, 1, 1)
            plt.plot(raw_x, raw_y, "k.")
            plt.plot(raw_x, fx, color=cols.green)
            plt.subplot(3, 1, 2)
            plt.plot(raw_x, y, "k.")
            plt.subplot(3, 1, 3)
            plt.axvline(1./truth, color=".7", linestyle="--")
            if s == 1:
                plt.axvline(K2_period, color=cols.orange, linestyle="--")
            plt.plot(fs, _K2pgram/max(_K2pgram), color=cols.blue,
                     label="$\mathrm{K2pgram$}")
            plt.plot(fs, pgram/max(pgram), color=cols.pink,
                     label="$\mathrm{raw pgram}$")
            plt.savefig("injections/sine/%s_%s_result" % (str(int(truth)), i))
    return np.array(K2a), np.array(K2P), np.array(rawa), np.array(rawP)

# add simulated to real light curves and grid over periods
def grid_over_periods(basis, raw_x, raw_y, true_p, fs, true_a, fnames):
    K2_amps, K2_Ps, raw_amps, raw_Ps = [], [], [], []
    for i, fname in enumerate(fnames):
        print fname
        print true_p[i]
        time, flux = np.genfromtxt(fname).T
        K2a, K2P, rawa, rawP = grid_over_amps(basis, flux, raw_x, raw_y,
                                              true_p[i], fs, amps, true_a[i])
        K2_amps.append(K2a)
        raw_amps.append(rawa)
        K2_Ps.append(K2P)
        raw_Ps.append(rawP)
    K2_amps = np.array([j for i in K2_amps for j in i])
    K2_Ps = np.array([j for i in K2_Ps for j in i])
    raw_amps = np.array([j for i in raw_amps for j in i])
    raw_Ps = np.array([j for i in raw_Ps for j in i])

    f = h5py.File("injections/sine/histogram.h5", "w")
    K2data = f.create_dataset("K2", (len(K2_amps), 2))
    K2data[:, 0] = K2_amps
    K2data[:, 1] = K2_Ps
    rawdata = f.create_dataset("raw", (len(raw_amps), 2))
    rawdata[:, 0] = raw_amps
    rawdata[:, 1] = raw_Ps
    f.close()
    return K2_amps, K2_Ps, raw_amps, raw_Ps

def make_histogram_plot(nbins=10):
#     with h5py.File("histogram.h5", "r") as f:
    with h5py.File("injections/sine/histogram.h5", "r") as f:
        K2_amps = f["K2"][:, 0]
        K2_Ps = f["K2"][:, 1]
        raw_amps = f["raw"][:, 0]
        raw_Ps = f["raw"][:, 1]
    my_xedges = np.exp(np.linspace(np.log(min(K2_amps)), np.log(max(K2_amps)),
                       nbins))
    my_yedges = np.exp(np.linspace(np.log(min(K2_Ps)), np.log(max(K2_Ps)),
                       nbins))
    K2_hist, xedges, yedges = np.histogram2d(K2_amps, K2_Ps, bins=5,
                                             range=[[min(my_xedges),
                                                    max(my_xedges)],
                                                    [min(my_yedges),
                                                     max(my_yedges)]])
    raw_hist, xedges, yedges = np.histogram2d(raw_amps, raw_Ps)
    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    my_xedges = ["{:.4f}".format(i) for i in my_xedges]
    my_yedges = ["{:.4f}".format(i) for i in my_yedges]
    print my_yedges
    print my_xedges
    ax.set_xticklabels(my_xedges)
    ax.set_yticklabels(my_yedges)
    ax.xlabel("$\mathrm{Amplitude}")
    ax.ylabel("$\mathrm{Period~(days)}$")
    im = plt.imshow(K2_hist[::-1, ::-1], interpolation='nearest', origin='low',
                    extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                    cmap="Blues")
    plt.savefig("injections/sine/hist")
    plt.close(fig)

if __name__ == "__main__":

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

    # load basis
    with h5py.File("data/c1.h5", "r") as f:
        basis = f["basis"][:150, l]

    # load injections and truths
    sine = True
    if sine == True:
        fnames = glob.glob("injections/sine/????_lc.txt")
        name, true_p = np.genfromtxt("injections/sine/truth.txt").T
        true_a = np.ones_like(true_p)
    else:
        fnames = glob.glob("injections/*_lc.txt")
        name, true_p, true_a = np.genfromtxt("truth.txt").T

    fs = np.arange(1./70, 1.1, .001)
    amps = np.arange(.0001, .001, .0001)

#     # calculate the 2d histogram of completeness over period and amplitude
#     K2_amps, K2_Ps, raw_amps, raw_Ps = grid_over_periods(basis, raw_x,
#                                                          raw_y, true_p, fs,
#                                                           true_a, fnames)
    make_histogram_plot()
