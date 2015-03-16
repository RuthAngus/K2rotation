import numpy as np
import matplotlib.pyplot as plt
import glob
import pyfits
from mklc import mklc
import scipy.signal as sps
import fitsio
import h5py
from K2pgram import K2pgram, eval_freq
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

# calculate the false alarm probability
def fap(x, y, basis, fs, N, plot=False, sig=False):
    amp2s, s2n, _ = K2pgram(x, y, basis, fs)  # 1st pgram
    if sig: power = s2n
    else: power = amp2s
    mf, ms2n = peak_detect(fs, power)  # find peak
    AT = np.concatenate((basis, np.ones((3, len(y)))), axis=0)
    ATA = np.dot(AT, AT.T)
    # compute trends
    _, _, trends = eval_freq(x, y, mf, AT, ATA, compute_trends=True)
    if plot:
        plt.clf()
        plt.plot(1./fs, power, "k")
    peak_heights = []
    for n in range(N):
        detrended_y = y - trends  # remove trends
        detrended_y = np.random.choice(detrended_y, len(y))  # shuffle
        # add trends back in
        amp2s, s2n, _ = K2pgram(x, detrended_y + trends, basis, fs)
        if sig: power = s2n
        else: power = amp2s
        mx, my = peak_detect(fs, power)
        peak_heights.append(my)
        if plot:
            plt.plot(1./fs, power, alpha=.2)
    fap95 = np.percentile(peak_heights, 95)
    fap90 = np.percentile(peak_heights, 90)
    fap85 = np.percentile(peak_heights, 85)
    fap50 = np.percentile(peak_heights, 50)
    if plot:
        plt.axhline(fap95, color=".5")
        plt.savefig("fap")
    print fap95, fap90, fap85, fap50
    return fap95, fap90, fap85, fap50

def peak_detect(x, y):
    peaks = np.array([i for i in range(1, len(x)-1) if y[i-1] < y[i] and
                     y[i+1] < y[i]])
    l = y[peaks] == max(y[peaks])
    mx, my = x[peaks][l][0], y[peaks][l][0]
    return mx, my

# grid over amplitudes (the K2 pgram step takes the time)
def grid_over_amps(basis, flux, raw_x, raw_y, truth, fs, amps, true_a, fap,
                   flag, plot=True, raw=False):

    # find the threshold level
    _, initial_pgram, _ = K2pgram(raw_x, raw_y, basis, fs)
    mx, threshold = peak_detect(fs, initial_pgram)

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

        # calculate K2pgram
        amp2s, s2n, w = K2pgram(raw_x, y, basis, fs)
        pgram = s2n
        best_f, best_pgram = peak_detect(fs, pgram)  # find peaks
        print "recovered frequency", best_f
        s = 0  # success indicator
        if tf-.1*tf < best_f and best_f < tf+.1*tf and best_pgram > threshold:
            K2P.append(truth)
            K2a.append(a)
            print "success!", "\n"
            s = 1

        # calculate periodogram of raw light curve
        y = np.array([_y.astype("float64") for _y in y])
        raw_x = np.array([_raw_x.astype("float64") for _raw_x in raw_x])
        model = LombScargle().fit(raw_x, y, np.ones_like(y)*1e-5)
        period = 1. / fs
        pg = model.periodogram(period)
        best_f2, best_pg2 = peak_detect(fs, pg)
        if tf-.1*tf < best_f2 and best_f2 < tf+.1*tf:
            rawP.append(truth)
            rawa.append(a)

        if plot:
            plt.clf()
            plt.subplot(2, 1, 1)
            plt.plot(raw_x, y, "k.")
            plt.plot(raw_x, fx, color=cols.green)
            plt.subplot(2, 1, 2)
            plt.axvline(1./tf, color=".7", linestyle="--")
            plt.axhline(threshold, color=".7")
            if s == 1:
                plt.axvline(1./best_f, color="k")
            if flag == "r":
                plt.plot(1./fs, pgram, color=cols.blue,
                         label="$\mathrm{K2pgram$}")
            else:
                plt.plot(1./fs, pgram, color=cols.blue,
                         label="$\mathrm{K2pgram$}")
            plt.savefig("../injections/sine/%s_%s_result_%s"
                        % (str(int(truth)).zfill(2), i, flag))
    return np.array(K2a), np.array(K2P), np.array(rawa), np.array(rawP)

# add simulated to real light curves and grid over periods
def grid_over_periods(basis, raw_x, raw_y, true_p, fs, true_a, fap, fnames,
                      flag):
    K2_amps, K2_Ps, raw_amps, raw_Ps = [], [], [], []
    for i, fname in enumerate(fnames):
        print fname
        print true_p[i]
        time, flux = np.genfromtxt(fname).T
        K2a, K2P, rawa, rawP = grid_over_amps(basis, flux, raw_x, raw_y,
                                              true_p[i], fs, amps, true_a[i],
                                              fap, flag)
        K2_amps.append(K2a)
        raw_amps.append(rawa)
        K2_Ps.append(K2P)
        raw_Ps.append(rawP)
    K2_amps = np.array([j for i in K2_amps for j in i])
    K2_Ps = np.array([j for i in K2_Ps for j in i])
    raw_amps = np.array([j for i in raw_amps for j in i])
    raw_Ps = np.array([j for i in raw_Ps for j in i])

    f = h5py.File("../injections/sine/histogram.h5", "w")
    K2data = f.create_dataset("K2", (len(K2_amps), 2))
    K2data[:, 0] = K2_amps
    K2data[:, 1] = K2_Ps
    rawdata = f.create_dataset("raw", (len(raw_amps), 2))
    rawdata[:, 0] = raw_amps
    rawdata[:, 1] = raw_Ps
    f.close()
    return K2_amps, K2_Ps, raw_amps, raw_Ps

def make_histogram_plot(flag, nbins=10):
    with h5py.File("../injections/sine/histogram.h5", "r") as f:
        K2_amps = f["K2"][:, 0] * 100  # convert to %
        K2_Ps = f["K2"][:, 1]
        raw_amps = f["raw"][:, 0] * 100  # convert to %
        raw_Ps = f["raw"][:, 1]
    my_xedges = np.exp(np.linspace(np.log(min(K2_Ps)), np.log(max(K2_Ps)),
                       nbins))
    my_yedges = np.linspace(min(K2_amps), max(K2_amps), nbins)
    my_xedges = np.linspace(min(K2_Ps), max(K2_Ps), nbins)
    K2_hist, xedges, yedges = np.histogram2d(K2_Ps, K2_amps, bins=nbins,
                                             range=[[min(my_xedges),
                                                    max(my_xedges)],
                                                    [min(my_yedges),
                                                     max(my_yedges)]])

    raw_hist, xedges, yedges = np.histogram2d(raw_Ps, raw_amps, bins=nbins,
                                             range=[[min(my_xedges),
                                                    max(my_xedges)],
                                                    [min(my_yedges),
                                                     max(my_yedges)]])
    fig = plt.figure()
    ax = fig.add_subplot(111)
    X, Y = np.meshgrid(xedges, yedges)
    ax.pcolormesh(X, Y, K2_hist.T, cmap="Blues")
    ax.set_ylabel("$\mathrm{Amplitude~(\%)}$")
    ax.set_xlabel("$\mathrm{Period~(days)}$")
#     plt.plot(K2_Ps, K2_amps, "r.")
    plt.savefig("../injections/sine/K2_hist_%s" % flag)
    plt.close(fig)
    print K2_hist.T

    fig = plt.figure()
    ax = fig.add_subplot(111)
    X, Y = np.meshgrid(xedges, yedges)
    ax.pcolormesh(X, Y, raw_hist.T, cmap="Blues")
    ax.set_ylabel("$\mathrm{Amplitude}$")
    ax.set_xlabel("$\mathrm{Period~(days)}$")
#     plt.plot(raw_Ps, raw_amps, "r.")
    plt.savefig("../injections/sine/raw_hist_%s" % flag)
    plt.close(fig)
    print raw_hist.T

    fig = plt.figure()
    ax = fig.add_subplot(111)
    X, Y = np.meshgrid(xedges, yedges)
    b = K2_hist.T-raw_hist.T
    l = b < 0
    b[l] = 0
    print b
    ax.pcolormesh(X, Y, b, cmap="Blues")
    ax.set_ylabel("$\mathrm{Amplitude}$")
    ax.set_xlabel("$\mathrm{Period~(days)}$")
#     plt.plot(raw_Ps, raw_amps, "r.")
    plt.colorbar()
    plt.savefig("../injections/sine/both_hist_%s" % flag)
    plt.close(fig)

if __name__ == "__main__":

#     fname = 201295312
    fname = 201311941

    # load K2 light curve
    data = fitsio.read("../data/c1/ktwo%s-c01_lpd-lc.fits" % fname)
    aps = fitsio.read("../data/c1/ktwo%s-c01_lpd-lc.fits" % fname, 2)
    raw_y = data["flux"][:, np.argmin(aps["cdpp6"])]
    raw_x = data["time"]
    q = data["quality"]
    l = np.isfinite(raw_y) * np.isfinite(raw_x) * (q==0)
    raw_y, raw_x = raw_y[l], raw_x[l]
    raw_y /= np.median(raw_y)
    raw_y -= 1

    # load basis
    with h5py.File("../data/c1.h5", "r") as f:
        basis = f["basis"][:150, l]

    # load injections and truths
    sine = True
    flag = "r"  # r for rotation or a for asteroseismology
    if sine:
        fnames = glob.glob("../injections/sine/????_lc_%s.txt" % flag)
        name, true_p = np.genfromtxt("../injections/sine/truth_%s.txt"
                                     % flag).T
        true_a = np.ones_like(true_p)
    else:
        fnames = glob.glob("../injections/*_lc.txt")
        name, true_p, true_a = np.genfromtxt("truth.txt").T

    if flag == "r":
        fs = np.arange(1./30, 1.1, .001)
    else:
        f1 = 10. * 1e-6 * 24 * 3600
        f2 = 300. * 1e-6 * 24 * 3600
        fs = np.linspace(f1, f2, 1000)

    amps = np.arange(.0001, .001, .0001)

    fap = 2.40516004879e-06  # amp2s 95% fap
    # calculate the 2d histogram of completeness over period and amplitude
#     K2_amps, K2_Ps, raw_amps, raw_Ps = grid_over_periods(basis, raw_x,
#                                                          raw_y, true_p, fs,
#                                                          true_a, fap, fnames,
#                                                          flag)
    make_histogram_plot(flag)
