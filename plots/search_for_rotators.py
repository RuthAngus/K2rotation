import numpy as np
import matplotlib.pyplot as plt
import glob

eids = [201129544, 201132518, 201133037, 201133147, 201135311, 201138638,
        201138849, 201142023, 201142127]
fname = "/Users/angusr/data/K2/c1lcsr4/ep%s.csv" % eids[1]
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
