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
import scipy.interpolate as spi
import scipy.signal as sps
from delta_nu import find_modes, delta_nu

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

def campaign0():
    # campaign 0
    eid = "202064472"
#     eid = "202067981"
#     width = .25
#     eid = "202068435"
#     eid = "202065298"
#     eid = "202126997"
#     eid = "202086286"
    width = .22
    find_modes("/Users/angusr/Downloads/ktwo%s-c00_lpd-lc.fits" % eid, eid,
               nbasis=150)
    fs, s2n = np.genfromtxt("astero/%sastero_pgram.txt" % eid).T
    dnu, peak_height = delta_nu(fs, s2n, eid, width, sub=1)
    print dnu, peak_height

def campaign1():
    fnames = glob.glob("data/c1/*lc.fits")
    for fname in fnames:
        eid = fname[12:21]
        print eid
        find_modes(fname, str(int(eid)))
        fs, s2n = np.genfromtxt("astero/%sastero_pgram.txt"
                                % str(int(eid))).T
        width = .22
        dnu, peak_height = delta_nu(fs, s2n, eid, width, sub=1)
        print dnu, peak_height

def campaign1_vbg():
    fnames = glob.glob("data/c1/*lc.fits")
    for fname in fnames:
        eid = fname[12:21]
        print eid
        fname = "/Users/angusr/data/K2/c1lcsr4"
        # load data
        x, y, _ = np.genfromtxt("%s/ep%s.csv" % (fname, str(int(eid))),
                                delimiter=",").T
        y /= np.median(y)
        y -= 1
        x *= 24*3600  # convert to seconds
        # load basis
        with h5py.File("data/c1.h5", "r") as f:
            basis = f["basis"][:150]
        fs = np.arange(10, 300, 4e-2) * 1e-6
        ps = 1./fs
        model = LombScargle().fit(x, y, np.ones_like(y)*1e-5)
        s2n = model.periodogram(ps)
        # save pgram
        plt.clf()
        plt.plot(fs, s2n, "k")
        plt.savefig("astero/%svbg_pgram" % eid)
        np.savetxt("astero/%svbg_pgram.txt" % eid, np.transpose((fs, s2n)))

def kepler_poster_child():
    eid = "6442183"
    fs, fft = np.genfromtxt("6442183pgram.txt").T
    delta_nu(fs*1e6, fft, eid, .1, sub=1000)

if __name__ == "__main__":

#     poster_child = "201372313"
# #     find_modes("data/astero/ktwo201372313-c01_lpd-lc.fits", poster_child)
#     fs, s2n = np.genfromtxt("%sastero_pgram.txt" % poster_child).T
#     print len(fs)
#     delta_nu(fs, s2n, poster_child, width=.1, sub=1, tint=100)

#     kepler_poster_child()
#     campaign0()
    campaign1_vbg()
