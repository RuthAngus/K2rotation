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

# create file containing highest amplitude frequency for all elcs
def bases_FFT(eid, nbases):
    x, y, basis = read_data(eid, nbases)
    fs = np.linspace(1e-6, .7, 1000)
    ps = 1./fs
    plt.clf()
    cols = np.linspace(.1, .99, len(basis))
    freqs = []
    for i in range(len(basis)):
        model = LombScargle().fit(x, basis[i], np.ones_like(x)*1e-5)
        pgram = model.periodogram(ps)
        plt.plot(fs, pgram, color="%s" % cols[i])
        freqs.append(fs[pgram==max(pgram)][0])
    plt.savefig("all_bases")
    freqs = np.array(freqs)
    np.savetxt("elc_freqs.txt", np.transpose((np.arange(nbases), freqs)))

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
    mx, my = max_peak_detect(fs, s2n)
    print my

    plt.clf()
    plt.subplot(2, 1, 1)
    l = x < 2016
    plt.plot(x[l], y[l], "k")
    plt.plot(x[~l], y[~l], "k")
    plt.xlim(min(x), max(x))
    plt.xlabel("$\mathrm{BJD-2454833~(days)}$")
    plt.ylabel("$\mathrm{Normalized~Flux}$")

    plt.subplot(2, 1, 2)
    plt.plot(fs, s2n*1e5, "k")
    plt.xlabel("$\mathrm{Frequency~(days}^{-1}\mathrm{)}$")
    plt.ylabel("$\mathrm{S/N~(} \\times 10^5\mathrm{)}$")
    plt.ylim(0, my*1e5)
    plt.subplots_adjust(hspace=.4)
    plt.axvline(mx, color=".5", linestyle="--",
                label="$P_{rot}=%.2f ~\mathrm{days}$" % (1./mx))
    plt.legend()
    print "../documents/K2_rotation_%s.pdf" % epid
    plt.savefig("../documents/K2_rotation_%s.pdf" % epid)
    return mx

def K2_conditioned_plot(fs, epid):

    x, y, basis = read_data(epid, 150)
    amps2, s2n, w = K2pgram(x, y, basis, fs)

    # find highest peak
    mx, my = max_peak_detect(fs, s2n)

    # construct arrays
    AT = np.concatenate((basis, np.ones((3, len(y)))), axis=0)
    ATA = np.dot(AT, AT.T)
    _, _, trends = eval_freq(x, y, mx, AT, ATA, compute_trends=True)

    plt.clf()
    plt.subplot(2, 1, 1)
    l = x < 2016
    plt.plot(x[l], y[l], "k")
    plt.plot(x[l], y[l]-trends[l])
    plt.plot(x[~l], y[~l], "k")
    plt.plot(x[~l], y[~l]-trends[~l])
    plt.xlim(min(x), max(x))
    plt.xlabel("$\mathrm{BJD-2454833~(days)}$")
    plt.ylabel("$\mathrm{Normalized~Flux}$")

    plt.subplot(2, 1, 2)
    plt.plot(fs, s2n*1e5, "k")
    plt.xlabel("$\mathrm{Frequency~(days}^{-1}\mathrm{)}$")
    plt.ylabel(r"$\mathrm{S/N~(} \\times 10^5\mathrm{)}$")
    plt.ylim(0, my*1e5)
    plt.subplots_adjust(hspace=.4, bottom=.2)
    plt.axvline(mx, color=".5", linestyle="--",
                label="$P_{rot}=%.2f ~\mathrm{days}$" % (1./mx))
    plt.legend()
    plt.savefig("K2_%s" % epid)
    return mx

# plot the top 5 components
def top_5(x, basis, w):
    b = 3
    sw = np.sort(w)
    l = np.arange(len(w))[w == sw[0]][0]
    print l
    plt.clf()
    plt.subplot(5, 1, 1)
    plt.plot(x[::b], basis[l, :][::b], "k")
    plt.xticks(visible=False)
    plt.yticks(visible=False)
    plt.xlim(x[0], x[-1])
    plt.subplot(5, 1, 2)
    l = np.arange(len(w))[w == sw[1]][0]
    print l
    plt.plot(x[::b], basis[l, :][::b], "k")
    plt.yticks(visible=False)
    plt.xticks(visible=False)
    plt.xlim(x[0], x[-1])
    plt.subplot(5, 1, 3)
    l = np.arange(len(w))[w == sw[2]][0]
    print l
    plt.plot(x[::b], basis[l, :][::b], "k")
    plt.xticks(visible=False)
    plt.yticks(visible=False)
    plt.xlim(x[0], x[-1])
    plt.subplot(5, 1, 4)
    l = np.arange(len(w))[w == sw[3]][0]
    print l
    plt.plot(x[::b], basis[l, :][::b], "k")
    plt.xticks(visible=False)
    plt.yticks(visible=False)
    plt.xlim(x[0], x[-1])
    plt.subplot(5, 1, 5)
    l = np.arange(len(w))[w == sw[4]][0]
    print l
    plt.plot(x[::b], basis[l, :][::b], "k")
    plt.yticks(visible=False)
    plt.xlabel("$\mathrm{BJD-2454833}$")
    plt.subplots_adjust(hspace=0)
    plt.xlim(x[0], x[-1])
    plt.savefig("../documents/%s_top5.pdf" % epid)

if __name__ == "__main__":

    plotpar = {'axes.labelsize': 16,
               'text.fontsize': 16,
               'legend.fontsize': 16,
               'xtick.labelsize': 14,
               'ytick.labelsize': 16,
               'text.usetex': True}
    plt.rcParams.update(plotpar)

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

    K2_poster_child_plot(x, y, fs, amp2s, epid)
#     top_5(x, basis, w)
#     K2_conditioned_plot(fs, epid)
