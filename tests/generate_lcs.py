import numpy as np
import matplotlib.pyplot as plt
import glob
from mklc import mklc
import fitsio
import h5py
from K2pgram import K2pgram
from fit_all_the_light_curves import load_lc

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

    # load example star to get time array
    path = "/export/bbq2/dfm/k2/web/lightcurves/c1/201100000/21000"
    fname = "ktwo201121245-c01_lpd-lc.fits"
    raw_x, y, l = load_lc("%s/%s" % (path, fname))

    # load basis
    with h5py.File("../data/c1.h5", "r") as f:
        basis = f["basis"][:150, l]

    periods = 10**(np.linspace(np.log10(.5), np.log10(30), 1000))
    inject_sine_wave(raw_x, periods, "r")
    periods = 10**(np.linspace(np.log10(1./24), np.log10(5), 1000))
    inject_sine_wave(raw_x, periods, "a")
