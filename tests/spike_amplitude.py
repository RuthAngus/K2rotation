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
    for fname in fnames:
        eid = fname
        print eid
        assert 0
        x, y, _ = np.genfromtxt("../data/c1/ep%s.csv" % eid, delimiter=",").T
        x *= 24*3600
        y /= np.median(y)
        basis = load_K2_data(eid)
        fs = np.arange(40, 55, 1e-1) * 1e-6
        ps = 1./fs
        model = LombScargle().fit(x, y, np.ones_like(y)*1e-5)
        pg = model.periodogram(ps)
        f = h5py.File("spikes/spike_%s_vbg.h5" % eid, "w")
        data = f.create_dataset("pgram",  (len(fs), 2))
        data[:, 0] = fs
        data[:, 1] = pg
        f.close()

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
def assemble(fnames, vbg):
    mxs, mys = [], []
    for i, fname in enumerate(fnames):
        if vbg:
            eid = fname[15:24]
            print eid
            assert 0
        else:
            eid = fname[15:24]
        print eid, i, "of", len(fnames)
        if vbg:
            with h5py.File("spikes/spike_%s_vbg.h5" % eid, "r") as f:
                fs = f["pgram"][:, 0]
                pg = f["pgram"][:, 1]
        else:
            with h5py.File("spikes/spike_%s.h5" % eid, "r") as f:
                fs = f["pgram"][:, 0]
                pg = f["pgram"][:, 1]
        mx, my = peak_detect(fs, pg)
        mxs.append(mx)
        mys.append(my)
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

def experimental(mxs, mys):

    l = (mxs < 48*1e-6) * (46 * 1e-6 < mxs)
    l2 = (mxs < 46*1e-6) * (44 * 1e-6 < mxs)
    l3 = (mxs < 50*1e-6) * (48 * 1e-6 < mxs)
    l4 = (mxs < 44*1e-6) * (42 * 1e-6 < mxs)
    l5 = (mxs < 42*1e-6) * (40 * 1e-6 < mxs)
    l6 = (mxs < 52*1e-6) * (50 * 1e-6 < mxs)
    l7 = (mxs < 54*1e-6) * (52 * 1e-6 < mxs)
    plt.clf()
    plt.subplot(3, 1, 1)
    plt.hist(mxs, 50)
    plt.hist(mxs[l], 6)
    plt.hist(mxs[l2], 6)
    plt.hist(mxs[l3], 6)
    plt.subplot(3, 1, 2)
    plt.hist(np.log(mys[l]), 50, color="g", normed=True)
    plt.hist(np.log(mys[l2]), 50, color="r", normed=True)
    plt.hist(np.log(mys[l3]), 50, color="c", normed=True)

    plt.subplot(3, 1, 3)
    h1 = np.histogram(np.log(mys[l]), 50, normed=True)
    h2 = np.histogram(np.log(mys[l2]), 50, normed=True)
    h3 = np.histogram(np.log(mys[l3]), 50, normed=True)
    h4 = np.histogram(np.log(mys[l4]), 50, normed=True)
    h5 = np.histogram(np.log(mys[l5]), 50, normed=True)
    h6 = np.histogram(np.log(mys[l6]), 50, normed=True)
    h7 = np.histogram(np.log(mys[l7]), 50, normed=True)

#     h1 = np.histogram(np.log(mys[l]), 50)
#     h2 = np.histogram(np.log(mys[l2]), 50)
#     h3 = np.histogram(np.log(mys[l3]), 50)
#     h4 = np.histogram(np.log(mys[l4]), 50)
#     h5 = np.histogram(np.log(mys[l5]), 50)
#     h6 = np.histogram(np.log(mys[l6]), 50)
#     h7 = np.histogram(np.log(mys[l7]), 50)
    plt.plot(h1[1][1:], np.cumsum(h1[0]), "g")
    plt.plot(h2[1][1:], np.cumsum(h2[0]))
    plt.plot(h3[1][1:], np.cumsum(h3[0]))
    plt.plot(h4[1][1:], np.cumsum(h4[0]))
    plt.plot(h5[1][1:], np.cumsum(h5[0]))
    plt.plot(h6[1][1:], np.cumsum(h6[0]))
    plt.plot(h7[1][1:], np.cumsum(h7[0]))
    print sst.ks_2samp(h1[0], h2[0])
    print sst.ks_2samp(h1[0], h3[0])
    print sst.ks_2samp(h2[0], h3[0])

    print sst.ks_2samp(h4[0], h3[0])
    print sst.ks_2samp(h5[0], h4[0])
    print sst.ks_2samp(h6[0], h5[0])
    print sst.ks_2samp(h7[0], h6[0])
    plt.show()


if __name__ == "__main__":

    client = kplr.API()
#     fnames = glob.glob("../../old_K2rotation/data/c1/*.fits") # SIP
    fnames = glob.glob("../data/c1/*.csv")  # vbg

    load_kepmags(fnames)
#     find_spikes(fnames)
#     assemble(fnames)

#     with h5py.File("kepmag_spike.h5", "r") as f:
#         mxs = f["spikes"][:, 0]
#         mys = f["spikes"][:, 1]

#     eids, s2ns = find_value(fnames)
#     s2ns = find_value(fnames)
#      make_plot(fnames, np.array(s2ns))
