from raw_and_vbg import raw_and_vbg
from astero_examples import astero_example_plots
from k2_rotation_plots import K2_poster_child_plot
from simple_plots import demo
from K2pgram import K2pgram
import numpy as np

# RAW AND VBG
raw_and_vbg()

# GIANT EXAMPLES
c1 = [201545182, 201211472, 201433687, 201444854, 201607835]
astero_example_plots(c1)

# ROTATION EXAMPLE
epid = "201317002"
x, y, basis = read_data(epid, 200)
try:
    fs, s2n = np.genfromtxt("%spgram.txt" % epid).T
    print "periodogram file found"
except:
    fs = np.linspace(1e-6, .7, 1000)
    np.savetxt("%spgram.txt" % epid, np.transpose((fs, s2n)))
amp2s, s2n, w  = K2pgram(x, y, basis, fs)
K2_poster_child_plot(x, y, fs, amp2s, epid)

# EXOPLANET, EB AND RR LYRAE EXAMPLES
epids = ["201637175", "201665500"]  # planets
demo([epids[0]], "planet")
EB = ["201473612"]
demo(EB, "EB")
epids = np.genfromtxt("RRlyrae.txt").T
demo([epids[0]], "RR")
