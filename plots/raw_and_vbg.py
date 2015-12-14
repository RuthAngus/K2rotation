import numpy as np
import matplotlib.pyplot as plt
from gatspy.periodic import LombScargle
import fitsio


def raw_and_vbg():

    plotpar = {'axes.labelsize': 12,
               'text.fontsize': 18,
               'legend.fontsize': 18,
               'xtick.labelsize': 14,
               'ytick.labelsize': 14,
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
    l = np.isfinite(y) * np.isfinite(x) * (q == 0)
    y, x = y[l], x[l]
    MAD = np.median(y - np.median(y))
    y /= np.median(y)
    y -= 1
    x *= 24*3600  # convert to seconds

    fs = np.arange(.1, 300, 4e-2) * 1e-6  # astero

    # plot raw data
    fig = plt.figure()
    ax = fig.add_subplot(211)
    ax1 = fig.add_subplot(311)
    model = LombScargle().fit(x, y, np.ones_like(y)*1e-5)
    period = 1. / fs
    raw_pgram = model.periodogram(period)
    ax1.plot(fs[::3]*1e6, raw_pgram[::3], "k", label="$\mathrm{Raw}$")
    ax.set_title("$\mathrm{EPIC~%s}$" % eid)
    ax1.set_xlim(10, 280)
    ax1.set_ylim(0, .015)
    plt.ylabel("$\mathrm{Power}$")
#     plt.legend()
#     leg1 = Rectangle((0, 0), 0, 0, alpha=0.0)
#     plt.legend([leg1], "$\mathrm{Raw}$", handlelength=0)
    plt.text(230, .012, "$\mathrm{Raw}$")
    ticks = ax1.get_yticks()
    ax1.set_yticks(ticks[1:-1])
    ax.set_yticklabels(ax.get_yticklabels(), visible=False)
    ax.set_xticklabels(ax.get_xticklabels(), visible=False)

    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(labelcolor='w', top='off', bottom='off', left='off',
                   right='off')

    # load andrew's lcs
    ax2 = fig.add_subplot(312)
    x, y = np.genfromtxt("../data/c1/ep%s.csv" % eid, skip_header=1).T
    x *= 24*3600

    model = LombScargle().fit(x, y, np.ones_like(y)*1e-5)
    ps = 1. / fs
    pgram = model.periodogram(ps)
    ax2.plot(fs[::3]*1e6, pgram[::3], "k", label="$\mathrm{Detrended}$")
    plt.text(200, .010, "$\mathrm{VJ14~Detrended}$")
#     leg1 = Rectangle((0, 0), 0, 0, alpha=0.0)
#     plt.legend([leg1], "$\mathrm{VJ14~Detrended}$", handlelength=0)
    ax2.set_xlim(10, 280)
    ax2.set_ylim(0, .015)
    plt.ylabel("$\mathrm{Power}$")
    ticks = ax2.get_yticks()
    ax2.set_yticks(ticks[1:-1])
#     fig.text(0.04, 0.5, "$\mathrm{Power}$", ha="center", va="top",
#              rotation="vertical")

    # load sip
    fs, s2n = np.genfromtxt("../astero/%sastero_pgram.txt"
                            % str(int(eid))).T
    ax3 = fig.add_subplot(313)
    if MAD == 0.:
        MAD = 1.
    plt.plot(fs[::3], s2n[::3]*10e4/MAD**2, "k", label="$\mathrm{SIP}$")
#     leg1 = Rectangle((0, 0), 0, 0, alpha=0.0)
#     plt.legend([leg1], "$\mathrm{SIP}$", handlelength=0)
    plt.text(230, 2.2, "$\mathrm{SIP}$")
    ax3.set_xlim(10, 280)
    plt.ylabel("$\mathrm{Relative~(S/N)}^2\mathrm{~(} \\times 10^{4}\mathrm{)}$")
    plt.xlabel("$\\nu\mathrm{~(}\mu\mathrm{Hz)}$")
    fig.subplots_adjust(hspace=0, bottom=.1)
#     ticks = ax3.get_yticks()
#     ax3.set_yticks(ticks[1:-1])
    print("saving as ../documents/rawvbg_%s.pdf" % eid)
    print("saving as poster_rawvbg_%s.pdf" % eid)
    plt.savefig("../documents/rawvbg_%s.pdf" % eid)
    plt.savefig("poster_rawvbg_%s" % eid, transparent=True)

    # compute Fourier transform
    sp = np.fft.fft(y)
    freq = np.fft.fftfreq(x.shape[-1])
    fft = sp.real**2 * np.imag**2
    plt.clf()
    plt.plot(freq, fft)
    plt.savefig("fft")

if __name__ == "__main__":
    raw_and_vbg()
