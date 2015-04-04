import numpy as np
import matplotlib.pyplot as plt
import glob
from mklc import mklc
import fitsio
import h5py
from K2pgram import K2pgram

def inject_star_spots(raw_x, periods):
    amps = []
    for i, p in enumerate(periods):
        # generate injection
        nspot = 200
        incl = np.pi*5./12.
        amp = 1.
        tau = 30.5
        res0, res1 = mklc(raw_x, nspot, incl, amp, tau, p)
        time = res1[0]
        flux = res1[2]/np.median(res1[2]) - 1
        np.savetxt("injections/%s_lc.txt" % str(i).zfill(4),
                   np.transpose((time, flux)))
        amps.append(res0[-1])
    np.savetxt("injections/truth.txt", np.transpose((range(len(periods)),
               periods, np.array(amps))))

def inject_sine_wave(raw_x, periods, flag):
    for i, p in enumerate(periods):
        print i, "of", len(periods), "period = ", p
        flux = np.sin(raw_x * (1./p)*2*np.pi)
        np.savetxt("../injections/sine/%s_lc_%s.txt" % (str(i).zfill(4), flag),
                   np.transpose((raw_x, flux)))
    np.savetxt("../injections/sine/truth_%s.txt" % flag,
               np.transpose((range(len(periods)), periods)))

if __name__ == "__main__":

    fnames = [201310650, 201311700, 201311941]
    fname = fnames[-1]

    # load K2 light curve
    data = fitsio.read("../data/c1/ktwo%s-c01_lpd-lc.fits" % fname)
    aps = fitsio.read("../data/c1/ktwo%s-c01_lpd-lc.fits" % fname, 2)
    raw_y = data["flux"][:, np.argmin(aps["cdpp6"])]
    raw_x = data["time"]
    q = data["quality"]
    l = np.isfinite(raw_y) * np.isfinite(raw_x) * (q==0)
    raw_y, raw_x = raw_y[l], raw_x[l]
    raw_y /= np.median(raw_y)
    raw_y -= 1

    # load basis
    with h5py.File("../data/c1.h5", "r") as f:
        basis = f["basis"][:150, l]

    periods = 10**(np.linspace(np.log10(.5), np.log10(30), 1000))
    inject_sine_wave(raw_x, periods, "r")
    periods = 10**(np.linspace(np.log10(1./24), np.log10(2), 1000))
    inject_sine_wave(raw_x, periods, "a")
