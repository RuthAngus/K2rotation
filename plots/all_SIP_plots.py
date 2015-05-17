from raw_and_vbg import raw_and_vbg
from astero_examples import make_figs
from k2_rotation_plots import K2_poster_child_plot, read_data
from simple_plots import demo
from K2pgram import K2pgram
import numpy as np

# RAW AND VBG
raw_and_vbg()

# GIANT EXAMPLES
make_figs()

# ROTATION EXAMPLE
epid = "201317002"
x, y, basis = read_data(epid, 200)
fs = np.linspace(1e-6, .7, 1000)
amp2s, s2n, w  = K2pgram(x, y, basis, fs)
K2_poster_child_plot(x, y, fs, amp2s, epid)

# EXOPLANET, EB AND RR LYRAE EXAMPLES
epids = ["201637175", "201665500"]  # planets
demo([epids[0]], "planet")
EB = ["201473612"]
demo(EB, "EB")
epids = np.genfromtxt("RRlyrae.txt").T
demo([epids[0]], "RR")
