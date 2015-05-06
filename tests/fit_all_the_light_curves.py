import numpy as np
import matplotlib.pyplot as plt
import h5py
import glob
import fitsio

"""
Fit all the light curves with the basis from c1 and then sample from those
weights

1. assemble_fnames - make a list of all K2 light curve names
2. load_lc - load x and y values of K2 light curve.
3. fit_lc - fit the basis to one light curve and return the weights.
4. fit_all_lcs - load file names, load data and fit all the light curves.
    Return a 2d array of weights.
5. reconstruct_fake_lc - randomly sample 150 weight values from the 2d array
    and make a new fake light curve.
"""

def assemble_fnames():
    fnames = []
    path = "/export/bbq2/dfm/k2/web/lightcurves/c1"
    i1s = glob.glob("%s/*" % path)
    for i1 in i1s:
        i2s = glob.glob("%s/*" % i1)
        for i2 in i2s:
            i3s = glob.glob("%s/*" % i2)
            fnames.append(i3s)
    fnames = [j for i in fnames for j in i]
    return fnames

def load_lc(fname):
    data = fitsio.read(fname)
    aps = fitsio.read(fname, 2)
    y = data["flux"][:, np.argmin(aps["cdpp6"])]
    x = data["time"]
    q = data["quality"]
    l = np.isfinite(y) * np.isfinite(x) * (q==0)
    y, x = y[l], x[l]
    y /= np.median(y)
    y -= 1
    return x, y, l

def fit_lc(x, y, basis):
    # construct arrays
    AT = np.ones((151, len(y)))
    AT[:-1, :] = basis
    ATA = np.dot(AT, AT.T)
    return np.linalg.solve(ATA, np.dot(AT, y))

def fit_all_lcs():

    # load basis
    with h5py.File("../data/c1.h5", "r") as f:
        basis = f["basis"][:150]

    # generate list of all k2 lc filenames
    fnames = assemble_fnames()

    # find the weight vectors for all the k2 light curves
    weights = np.zeros((151, len(fnames)))
    for i, fname in enumerate(fnames):
        x, y, l = load_lc(fname)
        weights[:, i] = fit_lc(x, y, basis.T[l].T)

    # save the weights
    f = h5py.File("all_weights.h5", "w")
    data = f.create_dataset("weights", np.shape(weights))
    data[:, :] = weights[:, :]
    return weights

def reconstruct_fake_lc(n):

    # load the weights
    with h5py.File("all_weights.h5", "r") as f:
        weights = f["weights"][...]

    # load basis
    with h5py.File("../data/c1.h5", "r") as f:
        basis = f["basis"][:150]

    # generate 150 random numbers between 1 and the number of lcs
    # to select the weights.
    chosen_w = np.array([np.random.choice(weights[i]) for i in range(150)])
    fake_lc = np.sum(basis.T * chosen_w, axis=1)
    return fake_lc

if __name__ == "__main__":

#     weights = fit_all_lcs()
    nc1, nb = 21646, 150

    plt.clf()
    for i in range(10):
        print i
        fake_lc = reconstruct_fake_lc(nc1)
        plt.plot(fake_lc)
    plt.savefig("fake_lc")
