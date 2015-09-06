import numpy as np
import matplotlib.pyplot as plt
import wget
import h5py
from K2pgram import K2pgram
import fitsio
from plotstuff import params
reb = params()
from WD import EPICSIP
from SIP import SIP, eval_freq

def load_data(id):
        print "loading data"
        fname = "ktwo%s-c01_lpd-lc.fits" % id
        data = fitsio.read(fname)
        aps = fitsio.read(fname, 2)
        y = data["flux"][:, np.argmin(aps["cdpp6"])]
        x = data["time"]
        q = data["quality"]
        l = np.isfinite(y) * np.isfinite(x) * (q==0)
        y, x = y[l], x[l]
        MAD = np.median(y - np.median(y))
        x *= 24*3600

        print "load basis"
        with h5py.File("../data/c1.h5", "r") as f:
            basis = f["basis"][:150, l]

        # inject sinusoid at 1000 ppm
        med = np.median(y)
        amp = med / 1e3
        print med
        print amp
        y += amp*np.sin(x*2*np.pi*200*1e-6)

        y = y / med - 1

        print "computing SIP..."
        fs = np.arange(10, 300, 4e-2) * 1e-6
        amp2s, s2n, w = K2pgram(x, y, basis, fs)
        print max(amp2s**.5)

        plt.clf()
        plt.plot(fs[::5], amp2s[::5]**.5*1e6, "k")
        plt.xlabel("$\\nu\mathrm{~(}\mu\mathrm{Hz)}$")
        plt.ylabel("$\mathrm{Relative~(S/N)}^2$")
        plt.title("$\mathrm{EPIC~%s}$" % str(int(id)))
        plt.subplots_adjust(bottom=.15, left=.2)
        plt.savefig("test")
        assert 0

        return x, y, basis

def list_SIP(epids):
    for id in epids:
        print id

        x, y, basis = load_data(id)
        med = np.median(y)
        y = y/med - 1

        plt.clf()
        plt.plot(x, y, "k.")
        plt.savefig("%s_lc" % id)

        print "computing SIP..."
        fs = np.arange(10, 300, 4e-2) * 1e-6
        amps2, s2n, w = K2pgram(x, y, basis, fs)

        if MAD == 0: MAD = 1
        print MAD
#         signal = s2n/MAD**2  # normalising
        data = np.vstack((fs, amps2**.5*med)).T

        print "saving SIP"
        np.savetxt("%s.txt" % id, data)

if __name__ == "__main__":
    epids = np.genfromtxt("ktwo_c1_APO-RGs_llc.dat.epic.list", dtype=str).T
    list_SIP(epids)
    assert 0

    for id in epids:
        fs, amp2s = np.genfromtxt("%s.txt" % id).T
        print "plotting..."
        plt.clf()
        plt.plot(fs[::5], amp2s[::5], "k")
        plt.xlabel("$\\nu\mathrm{~(}\mu\mathrm{Hz)}$")
        plt.ylabel("$\mathrm{Relative~(S/N)}^2$")
        plt.title("$\mathrm{EPIC~%s}$" % str(int(id)))
        plt.subplots_adjust(bottom=.15, left=.2)
        plt.savefig("%s" % str(int(id)))
