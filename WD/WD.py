import numpy as np
import matplotlib.pyplot as plt
from EPIC import EPICSIP
from gatspy.periodic import LombScargle
import fitsio
import h5py
from SIP import SIP, eval_freq

def EPICSIP(EPIC, C, periods, nELC=150):
    '''
    Given an EPIC ID, campaign number and period array,
    produce a sip.
    '''
    fname = "ktwo%s-c%s_lpd-lc.fits" % (str(int(EPIC)), str(int(C)).zfill(2))
    print fname, "found"
    # load the data
    data = fitsio.read(fname)
    aps = fitsio.read(fname, 2)
    y = data["flux"][:, np.argmin(aps["cdpp6"])]
    x = data["time"]
    q = data["quality"]
    m = np.isfinite(y) * np.isfinite(x) * (q==0)
    y, x = y[m], x[m]
    y = y / np.median(y) - 1
    m = y < .04
    x, y = x[m], y[m]
    nELC = 150
    with h5py.File("c%s.h5" % C, "r") as f:
        basis = f["basis"][:nELC, m]
    freqs = 1./periods
    s2n, amp2s, w = SIP(x, y, basis, freqs)
    return x, y, s2n, amp2s, w

EPIC, C = "201563164", 1

freqs = np.arange(.5, 24, .001)
periods = 1./freqs
x, y, s2n, amp2s, w = EPICSIP(EPIC, C, periods, nELC=1000)

plt.clf()
plt.subplot(2, 1, 1)
plt.plot(x, y, "k")
plt.ylim(-.1, .1)
plt.subplot(2, 1, 2)
m = s2n == max(s2n)
plt.plot(freqs, s2n/max(s2n), "k", label="period = %.5f" % periods[m])
plt.legend(loc="best")
plt.xlabel("Period (days)")
plt.ylabel("Relative (S/N)^2")
plt.savefig("%s_sip" % EPIC)

np.savetxt("%s_sip" % EPIC, np.transpose((freqs, amp2s)))

freqs = np.arange(4, 7, .0001)
periods = 1./freqs
x, y, s2n, amp2s, w = EPICSIP(EPIC, C, periods)

plt.clf()
m = s2n == max(s2n)
plt.plot(freqs, s2n/max(s2n), "k", label="period = %.5f" % periods[m])
plt.plot(freqs, amp2s/max(amp2s), "r", label="period = %.5f" % periods[m])
plt.legend(loc="best")
plt.xlabel("Period (days)")
plt.ylabel("Relative (S/N)^2")
plt.savefig("%s_sip_zoom" % EPIC)

# t, flux, _ = np.genfromtxt("ep201563164.csv", delimiter=",").T
# plt.clf()
# plt.subplot(2, 1, 1)
# plt.plot(t, flux, "k")
# plt.ylim(.9, 1.1)
# plt.subplot(2, 1, 2)
# model = LombScargle().fit(t, flux, np.ones_like(t)*1e-5)
# pgram = model.periodogram(periods)
# m = pgram == max(pgram)
# # plt.plot(periods, pgram, "k", label="period = %.5f" % periods[m])
# plt.plot(freqs, pgram, "k", label="period = %.5f" % periods[m])
# plt.xlabel("Period (days)")
# plt.legend()
# plt.savefig("%s_vbg" % EPIC)

# import nufft
# t, flux, _ = np.genfromtxt("ep201563164.csv", delimiter=",").T
# m = flux < 1.05
# t, flux = t[m], flux[m]
# flux -= np.mean(flux)
# plt.clf()
# plt.subplot(2, 1, 1)
# plt.plot(t, flux, "k")
# plt.subplot(2, 1, 2)
# pgram = nufft.nufft3(t, flux, freqs*2*np.pi)
# plt.ylim(0, .001)
# m = pgram == max(pgram)
# # plt.plot(periods, pgram, "k", label="period = %.5f" % periods[m])
# plt.plot(freqs, pgram, "k")
# plt.xlabel("Period (days)")
# plt.savefig("%s_vbg_fft" % EPIC)
