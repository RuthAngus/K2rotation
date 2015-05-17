import numpy as np
import matplotlib.pyplot as plt
import glob


fname = "/Users/angusr/data/K2/c1lcsr4/ep201129544.csv"
x, y, _ = np.genfromtxt(fname, delimiter=",").T
plt.clf()
plt.plot(x, y)
plt.show()
assert 0

fnames = glob.glob("/Users/angusr/data/K2/c1lcsr4/*")
for fname in fnames:
    x, y, _ = np.genfromtxt(fname, delimiter=",").T
    print fname
    plt.clf()
    plt.plot(x, y)
    plt.show()
