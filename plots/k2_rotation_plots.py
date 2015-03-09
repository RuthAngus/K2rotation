# This script creates the K2_rotation_poster_child.pdf figure
# it can also plot the FFT of the eigen light curves,
# the 2nd K2pgram (periodogram of a 2nd frequency),
# the data conditioned on the best freqency and the top
# eigen light curves.

import numpy as np
import matplotlib.pyplot as plt
import fitsio
from K2pgram import K2pgram, eval_freq
from rotation_poster_child import max_peak_detect
import h5py
from gatspy.periodic import LombScargle

plotpar = {'axes.labelsize': 12,
           'text.fontsize': 10,
           'legend.fontsize': 10,
           'xtick.labelsize': 8,
           'ytick.labelsize': 8,
           'text.usetex': True}
plt.rcParams.update(plotpar)

def read_data(epid, nbases):
    # read the data
    data = fitsio.read("../data/c1/ktwo%s-c01_lpd-lc.fits" % epid)
    aps = fitsio.read("../data/c1/ktwo%s-c01_lpd-lc.fits" % epid, 2)
    y = data["flux"][:, np.argmin(aps["cdpp6"])]
    x = data["time"]
    q = data["quality"]
    l = np.isfinite(y) * np.isfinite(x) * (q==0)
    y, x = y[l], x[l]
    y /= np.median(y)
    y -= 1

    # load basis
    with h5py.File("../data/c1.h5", "r") as f:
        basis = f["basis"][:nbases, l]
    return x, y, basis

def bases_FFT(x, y, fs, basis):
    ps = 1./fs
    plt.clf()
    cols = np.linspace(.1, .99, len(basis))
    for i in range(len(basis)):
        model = LombScargle().fit(x, basis[i], np.ones_like(y)*1e-5)
        pgram = model.periodogram(ps)
        plt.plot(fs, pgram, color="%s" % cols[i])
    plt.savefig("all_bases")

def K2pgram2(x, y, fs, basis):
    # construct arrays
    AT = np.concatenate((basis, np.ones((5, len(y)))), axis=0)
    ATA = np.dot(AT, AT.T)

    # compute 2nd k2pgram
    s2n = K2pgram2(x, y, mx, fs, AT, ATA)

    plt.clf()
    plt.subplot(2, 1, 1)
    plt.plot(x, y, "k")
    plt.xlim(min(x), max(x))
    plt.xlabel("$\mathrm{BJD-2454833~(days)}$")
    plt.ylabel("$\mathrm{Normalized~Flux}$")

    plt.subplot(2, 1, 2)
    plt.plot(fs, s2n, "k")
    plt.xlabel("$\mathrm{Frequency~(days}^{-1}\mathrm{)}$")
#     plt.ylabel("$\mathrm{S/N}$")
    plt.ylabel("$\mathrm{Power}$")
    plt.subplots_adjust(hspace=.4)
    mx, my = max_peak_detect(fs, s2n)
    plt.axvline(mx, color=".5", linestyle="--", label="$P_{rot}=%.2f$" % mx)
    plt.savefig("../documents/K22_rotation_poster_child.pdf")

# plot the best frequencies
def plot_best(x, y, fs, AT, ATA, mx):
    plt.clf()
    f = mx
    plt.subplot(2, 1, 1)
    s2n, trends = eval_freq(x, y, f, AT, ATA, compute_trends=True)
    plt.plot(x, y-trends, "k")
    plt.subplot(2, 1, 2)
    s2n, trends = eval_freq(x, y, f, AT, ATA, compute_trends=True)
    plt.plot(x, y, x, trends)
    plt.savefig("test1")

    plt.clf()
    f = mx
    plt.subplot(2, 1, 1)
    s2n, trends = eval_freq(x, y, 1./10.52, AT, ATA, compute_trends=True)
    plt.plot(x, y-trends, "k")
    plt.subplot(2, 1, 2)
    s2n, trends = eval_freq(x, y, 1./10.52, AT, ATA, compute_trends=True)
    plt.plot(x, y, x, trends)
    plt.savefig("test2")

def K2_poster_child_plot(x, y, fs, s2n, epid):

    # find highest peak
    mx, my = max_peak_detect(fs, amp2s)
    print my

    x2, y2, _ = np.genfromtxt("/Users/angusr/data/K2/c1lcsr4/ep%s.csv" % epid,
                              delimiter=",").T
    plt.clf()
    plt.subplot(2, 1, 1)
    l = x < 2016
    plt.plot(x[l], y[l], "k")
    plt.plot(x[~l], y[~l], "k")
    plt.xlim(min(x), max(x))
    plt.xlabel("$\mathrm{BJD-2454833~(days)}$")
    plt.ylabel("$\mathrm{Normalized~Flux}$")

    plt.subplot(2, 1, 2)
    plt.plot(fs, s2n, "k")
    plt.xlabel("$\mathrm{Frequency~(days}^{-1}\mathrm{)}$")
    plt.ylabel("$\mathrm{Power}$")
    plt.ylim(0, my)
    plt.subplots_adjust(hspace=.4)
    plt.axvline(mx, color=".5", linestyle="--",
                label="$P_{rot}=%.2f ~\mathrm{days}$" % (1./mx))
    plt.legend()
    plt.savefig("../documents/K2_rotation_%s.pdf" % epid)

# plot the top 5 components
def top_5(x, basis, w):
    sw = np.sort(w)
    l = np.arange(len(w))[w == sw[0]][0]
    print l
    plt.clf()
    plt.subplot(5, 1, 1)
    plt.plot(x, basis[l, :], "k")
    plt.xticks(visible=False)
    plt.yticks(visible=False)
    plt.xlim(x[0], x[-1])
    plt.subplot(5, 1, 2)
    l = np.arange(len(w))[w == sw[1]][0]
    print l
    plt.plot(x, basis[l, :], "k")
    plt.yticks(visible=False)
    plt.xticks(visible=False)
    plt.xlim(x[0], x[-1])
    plt.subplot(5, 1, 3)
    l = np.arange(len(w))[w == sw[2]][0]
    print l
    plt.plot(x, basis[l, :], "k")
    plt.xticks(visible=False)
    plt.yticks(visible=False)
    plt.xlim(x[0], x[-1])
    plt.subplot(5, 1, 4)
    l = np.arange(len(w))[w == sw[3]][0]
    print l
    plt.plot(x, basis[l, :], "k")
    plt.xticks(visible=False)
    plt.yticks(visible=False)
    plt.xlim(x[0], x[-1])
    plt.subplot(5, 1, 5)
    l = np.arange(len(w))[w == sw[4]][0]
    print l
    plt.plot(x, basis[l, :], "k")
    plt.yticks(visible=False)
    plt.xlabel("$\mathrm{BJD-2456000}$")
    plt.subplots_adjust(hspace=0)
    plt.xlim(x[0], x[-1])
    plt.savefig("../documents/%s_top5.pdf" % epid)

if __name__ == "__main__":

    epid = "201317002"
    # epid = "201372313"
    # epid = "201310077"
    # epid = "201315808"
    # epid = "201370145"
    # epid = "201318479"
    # epid = "201310768"

    x, y, basis = read_data(epid, 200)

    # compute K2 pgram
    try:
        fs, s2n = np.genfromtxt("%spgram.txt" % epid).T
        print "periodogram file found"
    except:
        fs = np.linspace(1e-6, .7, 1000)
        np.savetxt("%spgram.txt" % epid, np.transpose((fs, s2n)))
    amp2s, s2n, w  = K2pgram(x, y, basis, fs)

#     K2_poster_child_plot(x, y, fs, amp2s, epid)
    top_5(x, basis, w)
