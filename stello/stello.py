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

def EPICSIP(EPIC, C, periods, nELC=150):
    '''
    Given an EPIC ID, campaign number and period array,
    produce a sip.
    '''
    fname = "ktwo%s-c%s_lpd-lc.fits" % (str(int(EPIC)), str(int(C)).zfill(2))
    print fname, "found"
    # load the data
    data = fitsio.read(fname)
    aps = fitsio.read(fname, 2)
    y = data["flux"][:, np.argmin(aps["cdpp6"])]
    x = data["time"]
    q = data["quality"]
    m = np.isfinite(y) * np.isfinite(x) * (q==0)
    y, x = y[m], x[m]
    y = y / np.median(y) - 1
    nELC = 150
    with h5py.File("c%s.h5" % C, "r") as f:
        basis = f["basis"][:nELC, m]
    freqs = 1./periods
    s2n, amp2s, w = SIP(x, y, basis, freqs)
    return x, y, s2n, amp2s, w

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
        signal = s2n/MAD**2  # normalising
#         signal = amps2
        data = np.vstack((fs, signal)).T

        print "saving SIP"
        np.savetxt("%s.txt" % id, data)
        assert 0

if __name__ == "__main__":
    epids = np.genfromtxt("ktwo_c1_APO-RGs_llc.dat.epic.list", dtype=str).T
    list_SIP(epids)

    for id in epids:
        fs, s2n = np.genfromtxt("%s.txt" % id).T
        print s2n
        print "plotting..."
        plt.clf()
#         oom = 6
#         plt.plot(fs[::3], s2n[::3]*10**oom, "k")
        plt.plot(fs[::3], s2n[::3], "k")
#         plt.xlim(min(fs), 280)
        plt.xlabel("$\\nu\mathrm{~(}\mu\mathrm{Hz)}$")
#         plt.ylabel("$\mathrm{Relative~(S/N)}^2 \mathrm{(} \\times 10^%s\mathrm{)}$"
#                    % oom)
        plt.ylabel("$\mathrm{Relative~(S/N)}^2$")
#         ylims = plt.gca().get_ylim()
        plt.title("$\mathrm{EPIC~%s}$" % str(int(id)))
        plt.subplots_adjust(bottom=.15, left=.2)
        plt.savefig("%s" % str(int(id)))
        assert 0
