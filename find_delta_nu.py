import numpy as np
import matplotlib.pyplot as plt
import fitsio
import h5py
from K2pgram import K2pgram, eval_freq, K2pgram_basis
from gatspy.periodic import LombScargle
import fitsio
import glob
import emcee
import scipy.interpolate as spi
import scipy.signal as sps

plotpar = {'axes.labelsize': 10,
           'text.fontsize': 8,
           'legend.fontsize': 10,
           'xtick.labelsize': 10,
           'ytick.labelsize': 10,
           'text.usetex': True}
plt.rcParams.update(plotpar)

def smoothing(x, y):
    f = spi.interp1d(x, y)
    xx = np.linspace(min(x), max(x), 1000)
    yy = f(xx)
    window = sps.gaussian(20, 8)
    smoothed = sps.convolve(yy, window/window.sum(), mode='same')
    return xx, smoothed

def peak_detect(x, y):
    peaks = np.array([i for i in range(1, len(x)-1) if y[i-1] < y[i] and
                     y[i+1] < y[i]])
    l = y[peaks] == max(y[peaks])
    return x[peaks], y[peaks], x[peaks][l], y[peaks][l]

def find_modes(fname, eid, nbasis=150, campaign=1, raw=False):
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

    # plot raw data
    if raw == True:
        plt.clf()
        model = LombScargle().fit(x, y, np.ones_like(y)*1e-5)
        period = 1. / fs
        raw_pgram = model.periodogram(period)
        plt.plot(fs, raw_pgram, "k")
        plt.savefig("astero/raw_%spgram" % eid)

    # load basis
    with h5py.File("data/c%s.h5" % campaign, "r") as f:
        basis = f["basis"][:150, l]

    fs = np.arange(10, 300, 4e-2) * 1e-6
    print len(fs)
    amps2, s2n, w = K2pgram(x, y, basis, fs)

    # plot our pgram
    plt.clf()
    fs *= 1e6
    plt.plot(fs, s2n, "k")
    plt.xlabel("$\mathrm{Frequency~(}\mu \mathrm{Hz)}$")
    plt.ylabel("$\mathrm{Power}$")
    plt.savefig("astero/%sastero_pgram" % eid)

    # save pgram
    np.savetxt("astero/%sastero_pgram.txt" % eid, np.transpose((fs, s2n)))

def find_delta_nu(fs, s2n, eid, dnu, sub=1, smooth=False):

    fps = dnu * len(fs)
    df = (fs[1] - fs[0])  # frequency lag in uHz
    pos = np.arange(len(fs)-fps)[::sub]  # the position of each section
    print len(pos)
    acor = np.zeros((fps, len(pos)))
    for i in range(len(pos)):
        acor[:, i] = emcee.autocorr.function(s2n[i*sub:fps+(i*sub)])
    lags = np.arange(fps)*df

    plt.clf()
    plt.subplot(3, 1, 1)
    plt.plot(fs, s2n, "k")
    plt.xlim(min(fs), max(fs))
#     plt.ylim(0, max(s2n[fs > 10]))
    plt.xlabel("$\mathrm{Frequency~(}\mu \mathrm{Hz)}$")
    plt.ylabel("$\mathrm{Power}$")
    plt.yticks(visible=False)

    plt.subplot(3, 1, 2)
    plt.imshow(acor, cmap="gray_r", interpolation="nearest",
               aspect="auto", vmin=0, vmax=.3)
    plt.subplots_adjust(hspace=.35)
    plt.xticks(visible=False)
    plt.yticks(visible=False)
    plt.ylabel("$\Delta \\nu$")
    plt.xlabel("$\\nu_{max}~\mathrm{location}$")

    plt.subplot(3, 1, 3)
    collapsed_acf = np.sum(acor, axis=1)
    if len(lags) != len(collapsed_acf):
        lags = lags[:-1]
    plt.plot(lags, collapsed_acf, "k")
#     plt.xlim(1, max(lags))
    plt.xlabel("$\Delta \\nu~\mathrm{(}\mu\mathrm{Hz)}$")
    plt.ylabel("$\mathrm{Correlation}$")

    # cut off first part of ACF
#     l = lags > 6
#     lags, collapsed_acf = lags[l], collapsed_acf[l]

    # smooth acf
    if smooth == True:
        smoothx, smoothy = smoothing(lags, collapsed_acf)
        plt.plot(smoothx, smoothy)
        x_peaks, y_peaks, mx, my = peak_detect(smoothx, smoothy)
    else:
        x_peaks, y_peaks, mx, my = peak_detect(lags, collapsed_acf)

    plt.axvline(mx, color="k", alpha=.3, linestyle="--",
                label="$%s~\mu\mathrm{Hz}$" % mx[0])
    plt.legend()
#     plt.ylim(min(collapsed_acf), my)
    plt.ylim(-100, 200)
    plt.savefig("astero/%s_dnu" % eid)
    return mx[0], my[0], lags, pos, collapsed_acf, np.sum(acor, axis=0)
