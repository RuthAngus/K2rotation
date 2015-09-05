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

def list_SIP(epids):
    for id in epids:
        print id

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

        print "computing SIP..."
        fs = np.arange(10, 300, 4e-2) * 1e-6
        amps2, s2n, w = K2pgram(x, y, basis, fs)

        if MAD == 0: MAD = 1
        print MAD
#         signal = s2n/MAD**2  # normalising
        signal = amps2
        data = np.vstack((fs, amps2)).T

        print "saving SIP"
        np.savetxt("%s.txt" % id, data)

if __name__ == "__main__":
    epids = np.genfromtxt("ktwo_c1_APO-RGs_llc.dat.epic.list", dtype=str).T
    list_SIP(epids)

    for id in epids:
        fs, amp2s = np.genfromtxt("%s.txt" % id).T
        print "plotting..."
        plt.clf()
        plt.plot(fs[::5], amp2s[::5]/max(amp2s), "k")
        plt.xlabel("$\\nu\mathrm{~(}\mu\mathrm{Hz)}$")
        plt.ylabel("$\mathrm{Relative~(S/N)}^2$")
        plt.title("$\mathrm{EPIC~%s}$" % str(int(id)))
        plt.subplots_adjust(bottom=.15, left=.2)
        plt.savefig("%s" % str(int(id)))
        assert 0
