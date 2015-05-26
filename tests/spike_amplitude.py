# For the targets in the galactic archaeology dataset, find how the amplitude
# of the 47 uHz spike depends on Kepler magnitude

import numpy as np
import matplotlib.pyplot as plt
import fitsio
import h5py
import glob
from K2pgram import K2pgram, eval_freq
from K2misc import load_K2_data, peak_detect, load_basis
import kplr
import scipy.stats as sst
from matplotlib.ticker import NullFormatter
from gatspy.periodic import LombScargle
from mpl_toolkits.axes_grid1 import ImageGrid

# using kplr, find the kepler magnitude of each target
# and save to file
def load_kepmags(fnames):
    eids, kepmags = [], []
    for i, fname in enumerate(fnames):
        print fname, i, "of", len(fnames)
        eid = fname[15:24]
        star = client.k2_star(eid)
        eids.append(int(eid))
        kepmags.append(star.kmag)
    eids, kepmags = np.array(eids), np.array(kepmags)
    f = h5py.File("kepmags.h5", "w")
    data = f.create_dataset("stars", (len(eids), 2))
    data[:, 0] = eids
    data[:, 1] = kepmags
    f.close()
    return eids, kepmags

# calculate the periodogram of each target and save to file
# takes a list of globbed file names
def find_spikes(fnames):
    for fname in fnames:
        eid = fname[15:24]
        print eid
        x, y, basis = load_K2_data(eid)
        fs = np.arange(40, 55, 1e-1) * 1e-6
        _, pg, _ = K2pgram(x, y, basis, fs)
        f = h5py.File("spikes/spike_%s.h5" % eid, "w")
        data = f.create_dataset("pgram",  (len(fs), 2))
        data[:, 0] = fs
        data[:, 1] = pg
        f.close()

# calculate the periodogram of each vbg light curve and save to file
def find_spikes_vbg(fnames):
    for i, fname in enumerate(fnames):
        eid = fname[15:24]
        print eid, i, "of", len(fnames)
        x, y, _ = np.genfromtxt("../data/c1/ep%s.csv" % eid, delimiter=",").T
        x *= 24*3600
        y /= np.median(y)
        fs = np.arange(40, 55, 1e-1) * 1e-6
        ps = 1./fs
        model = LombScargle().fit(x, y, np.ones_like(y)*1e-5)
        pg = model.periodogram(ps)
        f = h5py.File("spikes/spike_%s_vbg.h5" % eid, "w")
        data = f.create_dataset("pgram",  (len(fs), 2))
        data[:, 0] = fs
        data[:, 1] = pg
        f.close()

#         plt.clf()
#         plt.subplot(2, 1, 1)
#         plt.plot(x, y, "k")
#         plt.subplot(2, 1, 2)
#         plt.plot(fs, pg, "k")
#         plt.savefig("test")
#         assert 0

# calculate the amplitude at 47uHz
def find_value(fnames):
    s2ns, eids = [], []
    for i, fname in enumerate(fnames):
        eid = fname[15:24]
        print eid, i, "of", len(fnames)
        x, y, basis = load_K2_data(eid)
        f = 47.2281
        # construct arrays
        AT = np.concatenate((basis, np.ones((3, len(y)))), axis=0)
        ATA = np.dot(AT, AT.T)
        _, s2n, _ = eval_freq(x, y, f, AT, ATA)
        s2ns.append(s2n)
        eids.append(eid)
#     f = h5py.File("s2ns.h5", "w")
#     data = f.create_dataset("s2n",  (len(s2ns), 2))
#     data[:, 0] = np.array(eids)
#     data[:, 1] = np.array(s2ns)
#     f.close()
    return s2ns

# calculate the positions and heights of periodogram peaks,
# load the kepmags and save to file
def assemble(fnames, vbg=False):
    mxs, mys = [], []
    for i, fname in enumerate(fnames):
        eid = fname[15:24]
        print eid, i, "of", len(fnames)
        if vbg:
            with h5py.File("spikes/spike_%s_vbg.h5" % eid, "r") as f:
                fs = f["pgram"][:, 0]
                pg = f["pgram"][:, 1]
        else:
            try:
                with h5py.File("spikes/spike_%s.h5" % eid, "r") as f:
                    fs = f["pgram"][:, 0]
                    pg = f["pgram"][:, 1]
            except:
                IOError
        try:
            mx, my = peak_detect(fs, pg)
            mxs.append(mx)
            mys.append(my)
        except:
            IndexError  # some nan problems causing these errors
    if vbg:
        f = h5py.File("kepmag_spike_vbg.h5", "w")
    else:
        f = h5py.File("kepmag_spike.h5", "w")
    data = f.create_dataset("spikes",  (len(mxs), 2))
    data[:, 0] = np.array(mxs)
    data[:, 1] = np.array(mys)
    f.close()
    return mxs, mys

# make a plot of bad peak height vs kepler magnitude
def make_plot(fnames, s2ns):
#     with h5py.File("s2ns.h5", "r") as f:
#         eids1 = f["s2n"][:, 0]
#         s2ns = f["s2n"][:, 1]
#     with h5py.File("kepmags.h5", "r") as f:
#         eids2 = f["stars"][:, 0]
#         kepmags = f["stars"][:, 1]
#     print eids1 == eids2

#     eids = []
#     for fname in fnames:
#         eids.append(fname[15:24])
#     eids = np.array(eids)
#     l = s2ns > 1e-4
#     print eids[l]
#     print s2ns
#     plt.clf()
#     plt.hist(s2ns, 100)
#     plt.savefig("bad")

    with h5py.File("kepmag_spike.h5", "r") as f:
        mxs = f["spikes"][:, 0]
        mys = f["spikes"][:, 1]
        kepmags = f["spikes"][:, 1]
    l = mys

def experimental(mxs, mys, paper=True):
    mxs *= 1e6
    s1, s2, s3, s4 = 45, 46.5, 48, 49.5
    l1 = (mxs < s2) * (s1 < mxs)
    l2 = (mxs < s3) * (s2 < mxs)
    l3 = (mxs < s4) * (s3 < mxs)
    nbins = 50
    plt.clf()
    plt.subplot(2, 1, 1)
    if paper == True:
        plt.axvspan(s1+.05, s2, facecolor="#9999FF", edgecolor="w", hatch="\\")
        plt.axvspan(s2, s3-.05, facecolor="w", edgecolor="#FF66FF", hatch="\\")
        plt.axvspan(s3, s4-.05, facecolor="w", edgecolor="#00CCCC", hatch="/")
    else:
        plt.axvspan(s1, s2, facecolor="#6666FF", edgecolor="w", hatch="\\")
        plt.axvspan(s2, s3, facecolor="#FF66FF", edgecolor="w")
        plt.axvspan(s3, s4, facecolor="#00CCCC", edgecolor="w", hatch="/")
    plt.hist(mxs, nbins, color=".3", edgecolor=".3", rwidth=.7)
    plt.ylim(1e0, 1e5)
    plt.yscale("log")
    nbins=50
    plt.xlim(40, 54)
    plt.xlabel("$\\nu~\mathrm{(}\mu\mathrm{Hz)}$")
    plt.ylabel("$\ln\mathrm{(N}_{\mathrm{peaks}})$")
    plt.subplot(2, 1, 2)
    if paper == True:
        plt.hist(np.log(mys[l2]), nbins, histtype="stepfilled",
                color="w", edgecolor="#FF66FF", hatch="\\")
        plt.hist(np.log(mys[l1]), nbins, histtype="stepfilled",
                 color="#9999FF", edgecolor="w", hatch="\\")
        plt.hist(np.log(mys[l3]), nbins, histtype="stepfilled",
                 color="w", edgecolor="#00CCCC", hatch="/")
    else:
        plt.hist(np.log(mys[l2]), nbins, histtype="stepfilled",
                color="#FF66FF", edgecolor="w")
        plt.hist(np.log(mys[l1]), nbins, histtype="stepfilled",
                 color="#6666FF", edgecolor="w", hatch="\\")
        plt.hist(np.log(mys[l3]), nbins, histtype="stepfilled",
                 color="#00CCCC", edgecolor="w", hatch="/")
    plt.ylabel("$\ln\mathrm{(N}_{\mathrm{peaks}})$")
    plt.xlabel("$\mathrm{\ln(Maximum~peak~height~[Relative~(S/N)}^2\mathrm{])}$")
    plt.yscale("log")
    plt.subplots_adjust(hspace=.3)
    plt.savefig("sip_hist.pdf")

def experimental_vbg(mxs, mys, paper=True):
    mxs *= 1e6
    s1, s2, s3, s4 = 45, 46.5, 48, 49.5
    l1 = (mxs < s2) * (s1 < mxs)
    l2 = (mxs < s3) * (s2 < mxs)
    l3 = (mxs < s4) * (s3 < mxs)
    nbins = 50
    plt.clf()
    plt.subplot(2, 1, 1)
    if paper == True:
        plt.axvspan(s1+.05, s2, facecolor="#9999FF", edgecolor="w", hatch="\\")
        plt.axvspan(s2, s3-.05, facecolor="w", edgecolor="#FF66FF", hatch="\\")
        plt.axvspan(s3, s4-.05, facecolor="w", edgecolor="#00CCCC", hatch="/")
    else:
        plt.axvspan(s1, s2, facecolor="#6666FF", edgecolor="w", hatch="\\")
        plt.axvspan(s2, s3, facecolor="#FF66FF", edgecolor="w")
        plt.axvspan(s3, s4, facecolor="#00CCCC", edgecolor="w", hatch="/")
    plt.hist(mxs, nbins, color=".3", edgecolor=".3", rwidth=.7)
    plt.xlabel("$\\nu~\mathrm{(}\mu\mathrm{Hz)}$")
    plt.ylabel("$\mathrm{N}_{\mathrm{peaks}}$")
    plt.xlim(40, 54)
    plt.yscale("log")
    plt.ylim(10e-1, 10e4)
    plt.ylabel("$\ln(\mathrm{N}_{\mathrm{peaks}})$")
    plt.subplot(2, 1, 2)
    if paper == True:
        plt.hist(np.log(mys[l2]), nbins, histtype="stepfilled",
                color="w", edgecolor="#FF66FF", hatch="\\")
        plt.hist(np.log(mys[l1]), nbins, histtype="stepfilled",
                 color="#9999FF", edgecolor="w", hatch="\\")
        plt.hist(np.log(mys[l3]), nbins, histtype="stepfilled",
                 color="w", edgecolor="#00CCCC", hatch="/")
    else:
        plt.hist(np.log(mys[l2]), nbins, histtype="stepfilled",
                color="#FF66FF", edgecolor="w")
        plt.hist(np.log(mys[l1]), nbins, histtype="stepfilled",
                 color="#6666FF", edgecolor="w", hatch="\\")
        plt.hist(np.log(mys[l3]), nbins, histtype="stepfilled",
                 color="#00CCCC", edgecolor="w", hatch="/")
    plt.ylabel("$\mathrm{N}_{\mathrm{peaks}}$")
    plt.xlabel("$\mathrm{\ln(Maximum~peak~height~(S/N)})$")
    plt.yscale("log")
    plt.ylabel("$\ln(\mathrm{N}_{\mathrm{peaks}})$")
    plt.xlabel("$\ln(\mathrm{Maximum~peak~height~[Power]})$")
    plt.subplots_adjust(hspace=.3)
    plt.savefig("vbg_hist.pdf")

def histhist(mxs, mys):
    mxs *= 1e6
    s1, s2, s3, s4 = 45, 46.5, 48, 49.5
    l1 = (mxs < s2) * (s1 < mxs)
    l2 = (mxs < s3) * (s2 < mxs)
    l3 = (mxs < s4) * (s3 < mxs)

    nbins = 20
    xhist = np.histogram(mxs, nbins)
    bins = xhist[1]
    inds = np.digitize(mxs, bins)
    l = inds == 1

    ax1 = plt.subplot2grid((2, nbins), (0, 0), colspan=nbins)
    if paper == True:
        plt.axvspan(s1+.05, s2, facecolor="#9999FF", edgecolor="w", hatch="\\")
        plt.axvspan(s2, s3-.05, facecolor="w", edgecolor="#FF66FF", hatch="\\")
        plt.axvspan(s3, s4-.05, facecolor="w", edgecolor="#00CCCC", hatch="/")
    else:
        plt.axvspan(s1, s2, facecolor="#6666FF", edgecolor="w", hatch="\\")
        plt.axvspan(s2, s3, facecolor="#FF66FF", edgecolor="w")
        plt.axvspan(s3, s4, facecolor="#00CCCC", edgecolor="w", hatch="/")
    plt.hist(mxs, nbins, color=".3", edgecolor=".3", rwidth=.7)
    plt.xlim(40, 54)
    plt.xlabel("$\\nu~\mathrm{(}\mu\mathrm{Hz)}$")
    plt.ylabel("$\mathrm{N}_{\mathrm{peaks}}$")

    plt.hist(mys[l], 20, orientation="horizontal", color=".3", edgecolor=".3",
             rwidth=.7)
    ax2 = plt.subplot2grid((2, nbins), (1, 0), colspan=1)
    plt.setp(plt.gca(), xticklabels=[], xticks=[])
    plt.ylim(-14, -6)
    for i in range(2, nbins+1):
        l = inds == i
        if len(mys[l]):
            plt.hist(np.log(mys[l]), 20, orientation="horizontal", color=".3",
                     edgecolor="w")
        ax2 = plt.subplot2grid((2, nbins), (1, (i-1)), colspan=1)
        plt.setp(plt.gca(), xticklabels=[], yticklabels=[], xticks=[],
                 yticks=[])
        plt.ylim(-14, -6)
    plt.subplots_adjust(wspace=0)
    plt.subplots_adjust(hspace=.3)
    plt.savefig("test")

def histhist_vbg(mxs, mys):
    mxs *= 1e6
    s1, s2, s3, s4 = 45, 46.5, 48, 49.5
    l1 = (mxs < s2) * (s1 < mxs)
    l2 = (mxs < s3) * (s2 < mxs)
    l3 = (mxs < s4) * (s3 < mxs)

    nbins = 10
    xhist = np.histogram(mxs, nbins)
    bins = xhist[1]
    inds = np.digitize(mxs, bins)
    l = inds == 1

    ax1 = plt.subplot2grid((2, nbins), (0, 0), colspan=nbins)
    plt.axvspan(s1, s2, facecolor="b", alpha=.5, edgecolor="w")
    plt.axvspan(s2, s3, facecolor="m", alpha=.5, edgecolor="w")
    plt.axvspan(s3, s4, facecolor="c", alpha=.5, edgecolor="w")
    plt.hist(mxs, nbins, color=".3", edgecolor=".3", rwidth=.7)
    plt.xlim(40, 54)
    plt.xlabel("$\\nu~\mathrm{(}\mu\mathrm{Hz)}$")
    plt.ylabel("$\mathrm{N}_{\mathrm{peaks}}$")

    plt.hist(mys[l], 20, orientation="horizontal", color=".3", edgecolor=".3",
             rwidth=.7)
    ax2 = plt.subplot2grid((2, nbins), (1, 0), colspan=1)
    plt.setp(plt.gca(), xticklabels=[], xticks=[])
    plt.ylim(-10, -2)
#     plt.xlim(0, 4000)
    for i in range(1, nbins):
        l = inds == i
        if len(mys[l]) > 1:
            plt.hist(np.log(mys[l]), 20, orientation="horizontal", color=".3",
                     edgecolor="w")
        print i
        ax2 = plt.subplot2grid((2, nbins), (1, (i)), colspan=1)
        plt.setp(plt.gca(), xticklabels=[], yticklabels=[], xticks=[],
                 yticks=[])
        plt.ylim(-10, -2)
#         plt.xlim(0, 4000)
    plt.subplots_adjust(wspace=0)
    plt.subplots_adjust(hspace=.3)
    plt.savefig("test")

if __name__ == "__main__":

    plotpar = {'axes.labelsize': 15,
                'text.fontsize': 15,
                'legend.fontsize': 15,
                'xtick.labelsize': 15,
                'ytick.labelsize': 15,
                'text.usetex': True}
    plt.rcParams.update(plotpar)

    client = kplr.API()
#     fnames = glob.glob("../../old_K2rotation/data/c1/*.fits") # SIP
    fnames = glob.glob("../data/c1/*.fits") # SIP

#     find_spikes_vbg(fnames)
#     find_spikes(fnames)
#     assemble(fnames)

    # making histogram plots
    with h5py.File("kepmag_spike_vbg.h5", "r") as f:
        mxs = f["spikes"][:, 0]
        mys = f["spikes"][:, 1]
    experimental_vbg(mxs, mys)

    with h5py.File("kepmag_spike.h5", "r") as f:
        mxs = f["spikes"][:, 0]
        mys = f["spikes"][:, 1]
    experimental(mxs, mys)

#     eids, s2ns = find_value(fnames)
#     s2ns = find_value(fnames)
#      make_plot(fnames, np.array(s2ns))
