import numpy as np
import matplotlib.pyplot as plt

def astero_example_plots(names, vbg=False):
    for eid in names:
        print eid
        if vbg:
            fs, s2n = np.genfromtxt("../astero/%svbg_pgram.txt"
                                    % str(int(eid))).T
        else:
            fs, s2n = np.genfromtxt("../astero/%sastero_pgram.txt"
                                    % str(int(eid))).T
        plt.clf()
        plt.plot(fs, s2n*10e4, "k")

        plt.xlim(min(fs), 280)
        plt.xlabel("$\\nu\mathrm{~(}\mu\mathrm{Hz)}$")
        plt.ylabel("$\mathrm{S/N~(} \\times 10^{-4}\mathrm{)}$")
        ylims = plt.gca().get_ylim()
        plt.title("$\mathrm{EPIC~%s}$" % str(int(eid)))
#         plt.text(15, ylims[-1]-.05*ylims[-1], "$\mathrm{EPIC~%s}$"
#                  % str(int(eid)))
        if vbg:
            plt.savefig("../documents/vbg%s.pdf" % str(int(eid)))
        else:
            plt.savefig("../documents/%s.pdf" % str(int(eid)))

if __name__ == "__main__":

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

    astero_example_plots(c0)
