import numpy as np
import matplotlib.pyplot as plt
from K2pgram import K2pgram, eval_freq
import wget
import h5py
import fitsio
from colours import plot_colours
cols = plot_colours()

plotpar = {'axes.labelsize': 12,
           'legend.fontsize': 10,
           'xtick.labelsize': 12,
           'ytick.labelsize': 12,
           'text.usetex': True}
plt.rcParams.update(plotpar)

def peak_detect(x, y):
    peaks = np.array([i for i in range(1, len(x)-1) if y[i-1] < y[i] and
                     y[i+1] < y[i]])
    l = y[peaks] == max(y[peaks])
    mx, my = x[peaks][l][0], y[peaks][l][0]
    return mx, my

def download(epids):
    base_url = "http://bbq.dfm.io/ketu/lightcurves/c1"
    for i, e in enumerate(epids):
        e = str(int(e))
        print e
        print i, "of ", len(epids)
        url = "%s/%s00000/%s000/ktwo%s-c01_lpd-lc.fits" \
                % (base_url, e[:4], e[4:6], e)
        print url
        wget.download(url)

def demo(epids, demo):
    for fname in epids:
        fname = str(int(fname))
        print fname
        # load K2 light curve
        data = fitsio.read("../data/c1/ktwo%s-c01_lpd-lc.fits" % fname)
        aps = fitsio.read("../data/c1/ktwo%s-c01_lpd-lc.fits" % fname, 2)
        y = data["flux"][:, np.argmin(aps["cdpp6"])]
        x = data["time"]
        q = data["quality"]
        l = np.isfinite(y) * np.isfinite(x) * (q==0)
        y, x = y[l], x[l]
        med = np.median(y)
        y /= np.median(y)
        y -= 1

        # load basis
        with h5py.File("../data/c1.h5", "r") as f:
            basis = f["basis"][:150, l]

        ps = np.linspace(.1, 1, 1000)
        amp2s, s2n, w = K2pgram(x, y, basis, 1./ps)
        mx, my = peak_detect(ps, s2n)

        # construct arrays
        AT = np.concatenate((basis, np.ones((3, len(y)))), axis=0)
        ATA = np.dot(AT, AT.T)
        _, _, trends = eval_freq(x, y, mx, AT, ATA, compute_trends=True)

        plt.clf()
        plt.subplot(2, 1, 1)
        xmin, xmax = 1980, 1995
        plt.xlim(xmin, xmax)
        l = (x < xmax) * (xmin < x)
#         yt = y[l] - trends[l]
#         ypos = plt.gca().get_ylim()
#         pos = max(yt) + .1 #*max(yt)
#         plt.text(1980.5, pos, "$\mathrm{EPIC~%s}$" % fname)
#         plt.text(1980.5, ypos+.1*ypos, "$\mathrm{EPIC~%s}$" % fname)
        plt.title("$\mathrm{EPIC~%s}$" % fname)
        plt.plot(x[l], y[l]-trends[l], "k")
        plt.xlabel("$\mathrm{BJD~(days) - 2454833}$")
        plt.ylabel("$\mathrm{Normalized~Flux}$")
        plt.subplot(2, 1, 2)
        plt.xlabel("$\mathrm{Period~(days)}$")
        plt.ylabel("$\mathrm{S/N}$")
        plt.ylim(0, max(s2n)+.2*max(s2n))
        plt.axvline(mx, color=".5", linestyle="--",
                    label="$%.2f\mathrm{~days}$" % mx)
        plt.legend()
        plt.plot(ps, s2n, "k")
        plt.subplots_adjust(hspace=.25)
        plt.savefig("../documents/%s_%s.pdf" % (demo, fname))

if __name__ == "__main__":
    # download(epids)
    epids = ["201637175", "201665500"]  # planets
    demo([epids[0]], "planet")
    EB = ["201473612"]
    demo(EB, "EB")
    epids = np.genfromtxt("RRlyrae.txt").T
    demo([epids[0]], "RR")
