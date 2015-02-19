import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import fitsio
import h5py
import time
from scipy.signal import lombscargle

vdb = pd.read_csv("data/ep201859398corcut.csv", index_col=False,
                  names=["time", "flux"], skiprows=1)

fn = "data/ktwo201859398-c01_lpd-lc.fits"
data = fitsio.read(fn)
aps = fitsio.read(fn, 2)
y = data["flux"][:, np.argmin(aps["cdpp6"])]
x = data["time"]
q = data["quality"]
l = np.isfinite(y) * np.isfinite(x) * (q == 0)
y, x = y[l], x[l]
y /= np.median(y)

with h5py.File("data/c1.h5", "r") as f:
    basis = f["basis"][:20, l]

AT = np.concatenate((basis, np.ones((3, len(y)))), axis=0)
ATA = np.dot(AT, AT.T)


def eval_freq(f, compute_trends=False):
    arg = 2*np.pi*f*x
    AT[-2, :] = np.sin(arg)
    AT[-1, :] = np.cos(arg)

    # AT.shape = (153, nt)
    # shape: (151, nt) * (nt, 2) -> (151, 2)
    v = np.dot(AT[:-2, :], AT[-2:, :].T)
    ATA[:-2, -2:] = v
    ATA[-2:, :-2] = v.T

    # AT[-2:, :].shape = (2, nt)
    # (2, nt), (nt, 2)
    ATA[-2:, -2:] = np.dot(AT[-2:, :], AT[-2:, :].T)

    w = np.linalg.solve(ATA, np.dot(AT, y))
    S = np.linalg.inv(ATA)[-2:, -2:]
    s2n = np.dot(w[-2:], np.linalg.solve(S, w[-2:]))

    if compute_trends:
        trends = np.dot(w[:-2], AT[:-2])
        return np.sum(w[-2:]**2), trends
    return s2n


# Plotting the true model.
a, t = eval_freq(0.1, compute_trends=True)
print a
plt.clf()
plt.subplot(2, 1, 1)
plt.plot(x, y, "k.")
plt.plot(x, t, "r")
plt.subplot(2, 1, 2)
plt.plot(x, y-t, "r")
plt.savefig("data")

fs = np.linspace(.05, .5, 10000)
p = lombscargle(np.array(vdb.time), np.array(vdb.flux) - 1.0, 2*np.pi*fs)
A = np.zeros_like(fs)
strt = time.time()
for i, f in enumerate(fs):
    A[i] = eval_freq(f)
print time.time() - strt

plt.clf()
plt.plot(fs, A / A.max())
plt.plot(fs, p / p.max())
plt.yscale("log")
plt.savefig("pgram")

# Plotting the "best" model.
fpeak = fs[A==max(A)]
a, t = eval_freq(fpeak, compute_trends=True)
print a
plt.clf()
plt.subplot(2, 1, 1)
plt.plot(x, y, "k.")
plt.plot(x, t, "r")
plt.subplot(2, 1, 2)
plt.plot(x, y-t, "r")
plt.savefig("best")
