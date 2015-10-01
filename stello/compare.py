import numpy as np
import matplotlib.pyplot as plt
import glob

# load K2pgram without sigma clipping
fnames = glob.glob("?????????_ns.txt")
for i, fname in enumerate(fnames):
    fs, amps = np.genfromtxt(fname).T
    plt.clf()
    plt.plot(fs, amps)
    plt.savefig("%s_ns" % str(i).zfill(3))

# load K2pgram with sigma clipping
fnames = glob.glob("?????????.txt")
for i, fname in enumerate(fnames):
    fs, amps = np.genfromtxt(fname).T
    plt.clf()
    plt.plot(fs, amps)
    plt.savefig("%s" % str(i).zfill(3))
