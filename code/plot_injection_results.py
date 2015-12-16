from __future__ import print_function
import numpy as np
import glob
import matplotlib.pyplot as plt

plt.rcParams.update({"axes.labelsize": 20,
                    "text.fontsize": 20,
                    "legend.fontsize": 20,
                    "xtick.labelsize": 15,
                    "ytick.labelsize": 15,
                    "text.usetex": True})

def success_list(rfs, ras, true_fs, true_as, tau):
    """
    returns an array of "successful recoveries"
    rfs: array of recovered frequencies
    true_fs: array of true frequencies
    tau: tolerance fraction
    """
    diff_frac = abs(rfs - true_fs)/true_fs  # frctnl diff btwn true & rcvrd
    m = diff_frac < tau
    return rfs[m], ras[m], true_fs[m], true_as[m]


def histo(rec_f, true_rec_a, all_f, all_a, nbins):
    """
    make a histogram for plotting
    rec_f: the recovered frequencies that are correct
    true_rec_a: the true amplitudes injected that were correctly recovered
    all_f: all injected frequencies
    all_a: all injected amplitudes
    """
    true_rec_a = np.log10(1e4 * true_rec_a**.5)  # convert amp2 to amp and ppm
    all_a = np.log10(1e4 * all_a**.5)
    rec_f, all_f = rec_f * 1e6, all_f * 1e6  # convert to uHz

    # calculate the bin edges and make histograms
    my_yedges = np.linspace(min(all_a), max(all_a), nbins)
    my_xedges = np.linspace(min(all_f), max(all_f), nbins)

    # recovered histogram
    hist, xedges, yedges = np.histogram2d(rec_f, true_rec_a, bins=nbins,
                                          range = [[min(my_xedges),
                                                   max(my_xedges)],
                                                   [min(my_yedges),
                                                    max(my_yedges)]])

    # true histogram
    all_hist, xedges, yedges = np.histogram2d(all_f, all_a, bins=nbins,
                                    range = [[min(my_xedges),
                                             max(my_xedges)],
                                             [min(my_yedges),
                                             max(my_yedges)]])

    # make the actual plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    X, Y = np.meshgrid(xedges, yedges)
    color = hist.T/all_hist.T * 100 # plot the % of recovered in each bin
    print(color)
    cax = ax.pcolormesh(X, Y, color, cmap="Blues")
    ax.set_ylabel("$\log_{10}\mathrm{Amplitude~(ppm)}$")
    ax.set_xlabel("$\\nu~\mathrm{(\\mu Hz)}$")
    plt.plot(all_f, all_a, "m.", ms=13)
    plt.plot(rec_f, true_rec_a, "y.", ms=5)
    plt.colorbar(cax, label="$\mathrm{Completeness~(\%)}$")
    plt.subplots_adjust(bottom=.2, left=.15)
    plt.savefig("hist")
    plt.close(fig)

if __name__ == "__main__":
    # load truths
    true_ids, true_fs, true_as = np.genfromtxt("truths.txt").T

    # load recoveries
    r_files = sorted(glob.glob("recovered_*.txt"))
    ids, rfs, ras = [], [], []
    for file in r_files:
        data = np.genfromtxt(file).T
        ids.append(data[0])
        rfs.append(data[1])
        ras.append(data[2])
    ids = np.array([i for j in ids for i in j])
    rfs = np.array([i for j in rfs for i in j])
    ras = np.array([i for j in ras for i in j])
    true_fs = true_fs[:len(rfs)]
    true_as = true_as[:len(rfs)]
    print(len(true_fs), "injected, ", len(rfs), "recovered", "\n")
    assert len(true_fs) == len(rfs)

    # find the successful recoveries
    rec_f, rec_a, true_rec_f, true_rec_a = \
            success_list(rfs, ras, true_fs, true_as, 1e-4)

    # make a histogram
    print(rec_f, true_rec_a, true_fs, true_as)
    histo(rec_f, true_rec_a, true_fs, true_as, 10)
