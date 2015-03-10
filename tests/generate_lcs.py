import numpy as np
import matplotlib.pyplot as plt
import glob
import pyfits
from mklc import mklc
import fitsio

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
    amps = []
    for i, p in enumerate(periods):
        print i, "of", len(periods), "period = ", p
        flux = np.sin(raw_x * (1./p)*2*np.pi)
        np.savetxt("../injections/sine/%s_lc_%s.txt" % (str(i).zfill(4), flag),
                   np.transpose((raw_x, flux)))
    np.savetxt("../injections/sine/truth_%s.txt" % flag,
               np.transpose((range(len(periods)), periods)))

if __name__ == "__main__":

    fname = 201295312

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

#     periods = np.arange(1., 70., .5)  # rotation periods
    f1 = 10. * 1e-6 * 24 * 3600
    f2 = 300. * 1e-6 * 24 * 3600
    periods = np.linspace(1./f1, 1./f2, 150)  # asteroseismology
    inject_sine_wave(raw_x, periods, "a")
