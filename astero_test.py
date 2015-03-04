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

def peak_detect(x, y):
    peaks = np.array([i for i in range(1, len(x)-1) if y[i-1] < y[i] and
                     y[i+1] < y[i]])
    return x[peaks], y[peaks]

# Giant stars
epicids = ["201372313", "201157462", "201690697", "201171432", "201919748"]

for eid in epicids[-1:]:

    # plot the target_pixel_file
    try:
        tpf = fitsio.read("data/ktwo%s-c01_lpd-targ.fits.gz" % eid)
        img = tpf["FLUX"][-1]
        plt.clf()
        plt.imshow(img.T, cmap="gray", interpolation="nearest");
        plt.savefig("%stpf" % eid)
    except:
        "File not found"
        print "no target pixel file"

    data = fitsio.read("data/ktwo%s-c01_lpd-lc.fits" % eid)
    aps = fitsio.read("data/ktwo%s-c01_lpd-lc.fits" % eid, 2)
    y = data["flux"][:, np.argmin(aps["cdpp6"])]
    x = data["time"]
    q = data["quality"]
    l = np.isfinite(y) * np.isfinite(x) * (q==0)
    y, x = y[l], x[l]
    y /= np.median(y)
    y -= 1
    x *= 24*3600  # convert to seconds

    # load basis
    with h5py.File("data/c1.h5", "r") as f:
        basis = f["basis"][:150, l]

    # construct arrays
    AT = np.concatenate((basis, np.ones((3, len(y)))), axis=0)
    ATA = np.dot(AT, AT.T)

    fs = np.arange(10, 300, 4e-2) * 1e-6
    print len(fs)
    amps2, s2n = K2pgram(x, y, fs, AT, ATA)

    # plot raw data
    plt.clf()
    model = LombScargle().fit(x, y, np.ones_like(y)*1e-5)
    period = 1. / fs
    raw_pgram = model.periodogram(period)
    plt.plot(fs, raw_pgram, "k")
    plt.savefig("raw_%spgram" % eid)

    # plot andrew's lc
    try:
        x_vbg, y_vbg, _ = np.genfromtxt("data/ep%s.csv" % eid,
                                        delimiter=",").T
        x_vbg *= 24*3600
        y_vbg = y_vbg/np.median(y_vbg) - 1
        model = LombScargle().fit(x_vbg, y_vbg, np.ones_like(y_vbg)*1e-5)
        period = 1. / fs
        pgram = model.periodogram(period)
        plt.clf()
    #     plt.subplot(2, 1, 1)
    #     plt.plot(x_vbg, y_vbg, "k.")
    #     plt.subplot(2, 1, 2)
        plt.plot(fs, pgram, "k")
        plt.xlabel("$\mathrm{Frequency~(}\mu \mathrm{Hz)}$")
        plt.ylabel("$\mathrm{Power}$")
        plt.savefig("vbg_%spgram" % eid)
    except:
        print "No vbg"

    # plot our pgram
    plt.clf()
    fs *= 1e6
    plt.plot(fs, s2n, "k")
    plt.xlabel("$\mathrm{Frequency~(}\mu \mathrm{Hz)}$")
    plt.ylabel("$\mathrm{Power}$")
    plt.savefig("%sastero_pgram" % eid)

#     # calculate the best model
#     x -= x[0]
#     x /= 24*3600
#     a, trends = eval_freq(x, y, 5.89511413, AT, ATA, compute_trends=True)
#     trends = trends - np.median(trends)
#     plt.clf()
#     plt.plot(x, trends, "r")
#     plt.plot(x, y, "k.", markersize=1.5)
#     plt.plot(x, y-trends+.01, "b.", markersize=1)
#     plt.savefig("%s_best" % eid)
#     print eid
