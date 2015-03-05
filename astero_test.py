# this script is for finding the modes of an asteroseismic star

import numpy as np
import matplotlib.pyplot as plt
import fitsio
import h5py
from K2pgram import K2pgram, eval_freq, K2pgram_basis
from colours import plot_colours
cols = plot_colours()
from gatspy.periodic import LombScargle
from params import plot_params
reb = plot_params()
import fitsio
import glob
import emcee

def peak_detect(x, y):
    peaks = np.array([i for i in range(1, len(x)-1) if y[i-1] < y[i] and
                     y[i+1] < y[i]])
    return x[peaks], y[peaks]

def find_modes(fname, eid, raw=False):
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
    with h5py.File("data/c1.h5", "r") as f:
        basis = f["basis"][:150, l]

    # construct arrays
    AT = np.concatenate((basis, np.ones((3, len(y)))), axis=0)
    ATA = np.dot(AT, AT.T)

    fs = np.arange(10, 300, 4e-2) * 1e-6
    print len(fs)
    amps2, s2n = K2pgram(x, y, fs, AT, ATA)

    # plot our pgram
    plt.clf()
    fs *= 1e6
    plt.plot(fs, s2n, "k")
    plt.xlabel("$\mathrm{Frequency~(}\mu \mathrm{Hz)}$")
    plt.ylabel("$\mathrm{Power}$")
    plt.savefig("astero/%sastero_pgram" % eid)

    # save pgram
    np.savetxt("%sastero_pgram.txt" % eid, np.transpose((fs, s2n)))

# plot the target_pixel_file
def plot_tpf(fname, eid):
    tpf = fitsio.read(fname)
    img = tpf["FLUX"][-1]
    plt.clf()
    plt.imshow(img.T, cmap="gray", interpolation="nearest");
    plt.savefig("astero/%stpf" % eid)

def plot_vbg(fname, eid):
    # plot andrew's lc
    x_vbg, y_vbg, _ = np.genfromtxt("data/ep%s.csv" % eid,
                                    delimiter=",").T
    x_vbg *= 24*3600
    y_vbg = y_vbg/np.median(y_vbg) - 1
    model = LombScargle().fit(x_vbg, y_vbg, np.ones_like(y_vbg)*1e-5)
    period = 1. / fs
    pgram = model.periodogram(period)
    plt.clf()
    plt.plot(fs, pgram, "k")
    plt.xlabel("$\mathrm{Frequency~(}\mu \mathrm{Hz)}$")
    plt.ylabel("$\mathrm{Power}$")
    plt.savefig("astero/vbg_%spgram" % eid)

def delta_nu(fs, s2n, eid, width=.1, sub=1000, tint=1000):
    df = (fs[1] - fs[0]) * 1e6  # frequency lag in uHz
    fps = width * len(fs)  # frequencies per section
    pos = np.arange(len(fs)-fps)[::sub]  # the position of each section
    acor = np.zeros((fps, len(pos)))
    for i in range(len(pos)):
        acor[:, i] = emcee.autocorr.function(s2n[i*sub:fps+i*sub])
    lags = np.arange(fps)*df
    plt.clf()
    plt.subplot(3, 1, 1)
    plt.plot(fs, s2n, "k")
    plt.subplot(3, 1, 2)
    plt.imshow(acor, cmap="gray_r", interpolation="nearest", aspect="auto")
#     fig, ax = plt.subplots()
#     ax.imshow(acor, cmap="gray_r", interpolation="nearest", aspect="auto")
#     ticks = np.arange(len(lags))[::tint]
#     ax.set_yticks(ticks)
#     ax.set_yticklabels(ticks*df)
    plt.subplot(3, 1, 3)
    collapsed_acf = np.sum(acor, axis=1)
    print np.shape(lags), np.shape(collapsed_acf)
    if len(lags) == len(collapsed_acf):
        plt.plot(lags, collapsed_acf, "k")
    else:
        plt.plot(lags[:-1], collapsed_acf, "k")
#     plt.ylim(-5, 10)
#     plt.ylim(-5, 10)
    plt.savefig("astero/%s_dnu" % eid)
    return acor, lags

if __name__ == "__main__":

    poster_child = "201372313"
#     find_modes("data/astero/ktwo201372313-c01_lpd-lc.fits", poster_child)
    fs, s2n = np.genfromtxt("%sastero_pgram.txt" % poster_child).T
    print len(fs)
    delta_nu(fs, s2n, poster_child, width=.1, sub=1, tint=100)

#     eid = "6442183"
#     fs, fft = np.genfromtxt("6442183pgram.txt").T
#     delta_nu(fs, fft, eid, width=.1, sub=1000, tint=1000)

#     fnames = glob.glob("data/astero/*lc.fits")

#     for fname in fnames:
#         eid = fname[16:25]
#         print eid
#         find_modes(fname, str(int(eid)))
#         delta_nu(str(int(eid)))
