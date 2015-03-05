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

def delta_nu(fs, s2n, eid, nsections=100):
    print s2n
    l = s2n > 0
    fs, s2n = fs[l], s2n[l]
    df = fs[1] - fs[0]
    fps = len(fs)/nsections  # frequencies per section
    acor = np.zeros((fps, nsections))
    for i in range(nsections):
        acor[:, i] = emcee.autocorr.function(s2n[i*fps:(i+1)*fps])
        lags = np.arange(len(acor[:, i]))*df
#         plt.clf()
#         plt.plot(lags, acor[:, i], "k")
#         plt.show()
#     fig, ax1 = plt.subplots(211)
    fig, ax = plt.subplots()
    ax.imshow(acor, cmap="gray", interpolation="nearest", aspect="auto")
    ax.set_yticklabels(range(len(lags))*lags)
    plt.savefig("astero/%s_dnu" % eid)
    plt.clf()
    return acor, lags

# find delta nu
def delta_nu_wrap(eid, nsections=100):
    fs, s2n = np.genfromtxt("%sastero_pgram.txt" % eid).T
    acor, lags = delta_nu(fs, s2n, eid, nsections=nsections)

if __name__ == "__main__":

    # Giant stars
#     poster_child = "201372313"
#     find_modes("data/astero/ktwo201372313-c01_lpd-lc.fits", poster_child)
#     delta_nu_wrap(poster_child)

    eid = "6442183"
    fs, fft = np.genfromtxt("6442183pgram.txt").T
    delta_nu(fs, fft, eid)

#     fnames = glob.glob("data/astero/*lc.fits")

#     for fname in fnames:
#         eid = fname[16:25]
#         print eid
#         find_modes(fname, str(int(eid)))
#         delta_nu(str(int(eid)))
