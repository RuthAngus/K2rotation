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

plotpar = {'axes.labelsize': 20,
           'text.fontsize': 20,
           'legend.fontsize': 20,
           'xtick.labelsize': 20,
           'ytick.labelsize': 20,
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
                   flag, n, plot=False, raw=False):

    # find the threshold level
    _, initial_pgram, _ = K2pgram(raw_x, raw_y, basis, fs)
    mx, threshold = peak_detect(fs, initial_pgram)

    K2P, rawP, K2a, rawa = [], [], [], []
    for i, a in enumerate(amps):

        tf = 1./truth
        print "period = ", truth
#         print "amplitude = ", a
#         print "frequency = ", tf

        # add lcs together
        fx = flux * a
        y = fx + raw_y
        SN = np.var(fx) / np.var(raw_y)

        # calculate K2pgram
        amp2s, s2n, w = K2pgram(raw_x, y, basis, fs)
        pgram = s2n
        best_f, best_pgram = peak_detect(fs, pgram)  # find peaks
#         print "recovered frequency", best_f
        print "recovered period", 1./best_f
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
            # plt.axhline(1.57445766955e-05, color="r")  # fap line
            c = cols.blue
            if s == 1:
                c = cols.pink
            plt.plot(1./fs, pgram, color=c,
                     label="$\mathrm{K2pgram$}")
            plt.savefig("../injections/sine/%s_%s_result_%s"
                        % (str(n).zfill(2), str(i).zfill(2), flag))
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
                                              fap, flag, i)
        K2_amps.append(K2a)
        raw_amps.append(rawa)
        K2_Ps.append(K2P)
        raw_Ps.append(rawP)
    K2_amps = np.array([j for i in K2_amps for j in i])
    K2_Ps = np.array([j for i in K2_Ps for j in i])
    raw_amps = np.array([j for i in raw_amps for j in i])
    raw_Ps = np.array([j for i in raw_Ps for j in i])

    f = h5py.File("../injections/sine/histogram_%s.h5" % flag, "w")
    K2data = f.create_dataset("K2", (len(K2_amps), 2))
    K2data[:, 0] = K2_amps
    K2data[:, 1] = K2_Ps
    rawdata = f.create_dataset("raw", (len(raw_amps), 2))
    rawdata[:, 0] = raw_amps
    rawdata[:, 1] = raw_Ps
    f.close()
    return K2_amps, K2_Ps, raw_amps, raw_Ps

def histo(amps, Ps, namps, npers, fname):
    nbins = namps
    max_n_per_bin = int(npers/namps)
    Ps = np.log(Ps)
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
    ax.set_ylabel("$\mathrm{Amplitude~(\%)}$")
    ax.set_xlabel("$\ln\mathrm{Period~(days)}$")
    plt.colorbar(cax, label="$\mathrm{Completeness~(\%)}$")
#     plt.plot(Ps, amps, "k.")
    plt.savefig("../injections/sine/%s_hist_%s.pdf" % (fname, flag))
    plt.close(fig)
    print K2_hist.T
    return K2_hist, xedges, yedges

def make_histogram_plot(flag, namps=19, npers=1000):
    with h5py.File("../injections/sine/histogram_%s.h5" % flag, "r") as f:
        K2_amps = f["K2"][:, 0] * 100  # convert to %
        K2_Ps = f["K2"][:, 1]
        raw_amps = f["raw"][:, 0] * 100  # convert to %
        raw_Ps = f["raw"][:, 1]

    if flag == "r":  # cut off short p_rots
        l = K2_Ps > np.exp(0.5)
        K2_amps, K2_Ps = K2_amps[l], K2_Ps[l]
        l = raw_Ps > np.exp(1)
        raw_amps, raw_Ps = raw_amps[l], raw_Ps[l]

    K2_hist, xedges, yedges = histo(K2_amps, K2_Ps, namps, npers, "K2")
    raw_hist, _, _ = histo(raw_amps, raw_Ps, namps, npers, "raw")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    X, Y = np.meshgrid(xedges, yedges)
    b = K2_hist.T-raw_hist.T
    l = b < 0
    b[l] = 0
    print b
    cax = ax.pcolormesh(X, Y, b, cmap="Blues")
    ax.set_ylabel("$\mathrm{Amplitude~(\%)}$")
    ax.set_xlabel("$\ln\mathrm{(Period)~(days)}$")
    plt.colorbar(cax, label="$\mathrm{Completeness~(\%)}$")
#     plt.plot(raw_Ps, raw_amps, "k.")
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
#     flag = "a"  # r for rotation or a for asteroseismology
    if sine:
        fnames = glob.glob("../injections/sine/????_lc_%s.txt" % flag)
        name, true_p = np.genfromtxt("../injections/sine/truth_%s.txt"
                                     % flag).T
        true_a = np.ones_like(true_p)
    else:
        fnames = glob.glob("../injections/*_lc.txt")
        name, true_p, true_a = np.genfromtxt("truth.txt").T

    if flag == "r":
#         ps = np.linspace(.4, 30., 1000)
        ps = np.linspace(1., 50., 1000)
        fs = 1./ps
    elif flag == "a":
        fs = np.linspace(1./4., 26., 1000)

#     amps = np.arange(.0, .001, .00005)
    amps = np.exp(np.linspace(np.log(1e-6), np.log(1e-3), 20))

    fap = 2.40516004879e-06  # amp2s 95% fap
    # calculate the 2d histogram of completeness over period and amplitude
#     K2_amps, K2_Ps, raw_amps, raw_Ps = grid_over_periods(basis, raw_x,
#                                                          raw_y, true_p, fs,
#                                                          true_a, fap, fnames,
#                                                          flag)
    make_histogram_plot("r")
    make_histogram_plot("a")
