import numpy as np
import matplotlib.pyplot as plt
import fitsio
import scipy.stats as sps

def astero_example_plots(names, vbg=False):
    for eid in names:
        print eid
        if vbg:
            fs, s2n = np.genfromtxt("../astero/%svbg_pgram.txt"
                                    % str(int(eid))).T
        else:
            # load raw data
            fname = "../data/c1/ktwo%s-c01_lpd-lc.fits" % eid
            data = fitsio.read(fname)
            aps = fitsio.read(fname, 2)
            y = data["flux"][:, np.argmin(aps["cdpp6"])]
            x = data["time"]
            q = data["quality"]
            l = np.isfinite(y) * np.isfinite(x) * (q==0)
            y, x = y[l], x[l]
            MAD = np.median(y - np.median(y))
            print MAD
            fs, s2n = np.genfromtxt("../astero/%sastero_pgram.txt"
                                    % str(int(eid))).T
        plt.clf()
        if MAD == 0: MAD = 1

        signal = s2n/MAD**2
        oom = 6
        print oom, "oom"

#         plt.plot(fs[::3], s2n[::3]/MAD**2, "k")
        plt.plot(fs[::3], s2n[::3]*10**oom, "k")

        plt.xlim(min(fs), 280)
        plt.xlabel("$\\nu\mathrm{~(}\mu\mathrm{Hz)}$")
        plt.ylabel("$\mathrm{Relative~(S/N)}^2 \mathrm{(} \\times 10^%s\mathrm{)}$"
                   % oom)
        ylims = plt.gca().get_ylim()
        plt.title("$\mathrm{EPIC~%s}$" % str(int(eid)))
        plt.subplots_adjust(bottom=.15, left=.2)
#         plt.text(15, ylims[-1]-.05*ylims[-1], "$\mathrm{EPIC~%s}$"
#                  % str(int(eid)))
        if vbg:
            plt.savefig("../documents/vbg%s.pdf" % str(int(eid)))
        else:
            plt.savefig("../documents/%s.pdf" % str(int(eid)))

def make_figs():
    plotpar = {'axes.labelsize': 20,
               'text.fontsize': 20,
               'legend.fontsize': 20,
               'xtick.labelsize': 20,
               'ytick.labelsize': 20,
               'text.usetex': True}
    plt.rcParams.update(plotpar)
    names, dn, nm = np.genfromtxt("../giant_candidates.txt").T
    # hand picked campaign 1 targets
    c1 = [201433687, 201444854, 201545182, 201607835, 201384413, 201272670]
    # paper
    c1 = [201545182, 201183188, 201211472, 201433687, 201444854, 201607835]
    # c0 targets from Lund et al.
    c0 = ["202064472", "202067981", "202068435", "202065298", "202126997",
            "202086286"]
    # paper
    c0 = ["202068435", "202086286"]
#     astero_example_plots(c0)
    astero_example_plots(c1)

if __name__ == "__main__":
    make_figs()
