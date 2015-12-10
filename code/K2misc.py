import numpy as np
import fitsio
import h5py


def load_K2_data(fname):
    """
    Load the time and flux for a k2 target
    input: fits file name and path
    output: time, flux and basis
    """
    # load K2 light curve
    data = fitsio.read(fname)
    aps = fitsio.read(fname, 2)
    y = data["flux"][:, np.argmin(aps["cdpp6"])]
    x = data["time"]
    q = data["quality"]
    l = np.isfinite(y) * np.isfinite(x) * (q == 0)
    y, x = y[l], x[l]
    y /= np.median(y)
    y -= 1
    x *= 24*3600

    # load basis
    with h5py.File("data/c1.h5", "r") as f:
        basis = f["basis"][:150, l]
    return x, y, basis


def peak_detect(x, y):
    """
    A VERY simple peak detecting algorithm
    in: x and y
    returns the peak positions and peak heights
    """
    peaks = np.array([i for i in range(1, len(x)-1) if y[i-1] < y[i] and
                     y[i+1] < y[i]])
    l = y[peaks] == max(y[peaks])
    mx, my = x[peaks][l][0], y[peaks][l][0]
    return mx, my


def load_basis(l):
    """
    Loads just the basis.
    takes a mask of non-zero inputs
    returns the basis at the right shape for a given target
    """
    with h5py.File("../data/c1.h5", "r") as f:
        basis = f["basis"][:150, l]
    return basis
