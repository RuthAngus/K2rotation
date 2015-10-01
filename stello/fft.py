import numpy as np
import matplotlib.pyplot as plt
import fitsio
from plotstuff import params
reb = params()
from gatspy.periodic import LombScargle
import nufft
from K2pgram import K2pgram
import h5py
from stello import sigma_clipping, load_data

def load_data_new(id):
    print "loading data"
    fname = "ktwo%s-c01_lpd-lc.fits" % id
    data = fitsio.read(fname)
    aps = fitsio.read(fname, 2)
    y = data["flux"][:, np.argmin(aps["cdpp6"])]
    x = data["time"]
    q = data["quality"]
    l = np.isfinite(y) * np.isfinite(x) * (q==0)
    with h5py.File("../data/c1.h5", "r") as f:
        basis = f["basis"][:150, l]
    return x[l], y[l], basis

def load_vbg(id):
    x, y, _ = np.genfromtxt("/Users/angusr/data/K2/c1lcsr4/ep%s.csv"
                            % id, delimiter=",").T
    return x, y

def LS(x, y, fs):
    ps = 1./fs
    model = LombScargle().fit(x, y, np.ones_like(y)*1e-5)
    return model.periodogram(ps)

def FFT(x, y, fs, nufft=False):
    if nufft: return nufft.nufft3(x, y, fs*2*np.pi)
    pgram = np.fft.fft(y)
    plt.clf()
    plt.plot(np.real(pgram[100:-100])**2 + np.imag(pgram[:-100])**2)
    plt.savefig("test")
    assert 0
    return pgram

if __name__ == "__main__":
    epids = np.genfromtxt("ktwo_c1_APO-RGs_llc.dat.epic.list", dtype=str).T

    for id in epids:
        print "\n", id

        # load data
        x, y, basis = load_data(id)
        med = np.median(y)
        y = y/med - 1

        fs = np.arange(10, 270, 4e-2) * 1e-6

        # compute SIP, ls and fft of raw lc
        amps2, s2n, w = K2pgram(x, y, basis, fs)
        ls = LS(x, y, fs)
        fft = FFT(x, y, fs)
        m = fft > 0

        # plot
        plt.clf()
        plt.subplot(3, 1, 1)
        plt.plot(fs, ls, "k")
        plt.subplot(3, 1, 2)
        plt.plot(fs[m], fft[m], "r")
        plt.subplot(3, 1, 3)
        plt.plot(fs, amps2/max(amps2), "b")
#         plt.plot(fs, s2n/max(s2n), "r", alpha=.5)
        plt.savefig("%s_fft" % id)

        # load vbg
        x, y = load_vbg(id)
        x *= 24*3600
        med = np.median(y)
        y = y/med - 1

        # compute ls and fft of vbg
        ls = LS(x, y, fs)
        fft = FFT(x, y, fs)
        m = fft > 0

        # plot
        plt.clf()
        plt.subplot(2, 1, 1)
        plt.plot(fs, ls, "k")
        plt.subplot(2, 1, 2)
        plt.plot(fs[m], fft[m], "r")
        plt.savefig("%s_vbgfft" % id)
