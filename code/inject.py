import numpy as np
import matplotlib.pyplot as plt
# import h5py
# from params import colours, params
from K2misc import load_K2_data, peak_detect
import glob
from SIP import SIP, eval_freq


def prewhiten(x, y, f):
    """
    A simple prewhitening function
    takes time and flux for the timeseries and the freq you want to
    remove
    """

    # firstly, find the amplitudes of the sine and cosine functions at
    # the frequency.
    # construct arrays
    AT = np.concatenate((basis, np.ones((3, len(y)))), axis=0)
    ATA = np.dot(AT, AT.T)
    amps2, s2n, w = eval_freq(x, y, f, AT, ATA)
    sin_amp, cos_amp = w[-2:]
    wave = sin_amp * np.sin(2*np.pi*f*x) + cos_amp * np.cos(2*np.pi*f*x)
    return y - wave

# load the data
fnames = glob.glob("data/ktwo*fits")
x, y, basis = load_K2_data(fnames[0])

# compute a SIP
# fs = np.arange(10, 300, 1e-1) * 1e-6
fs = np.arange(10, 300, 1) * 1e-6
s2n, amps2, w = SIP(x, y, basis, fs)

plt.clf()
plt.plot(fs*1e6, amps2)

peak_f, peak_height = peak_detect(fs, amps2)
new_y = prewhiten(x, y, peak_f)
s2n, amps2, w = SIP(x, y, basis, fs)
plt.plot(fs*1e6, amps2)

peak_f, peak_height = peak_detect(fs, amps2)
new_y = prewhiten(x, y, peak_f)
s2n, amps2, w = SIP(x, y, basis, fs)
plt.plot(fs*1e6, amps2)

peak_f, peak_height = peak_detect(fs, amps2)
new_y = prewhiten(x, y, peak_f)
s2n, amps2, w = SIP(x, y, basis, fs)
plt.plot(fs*1e6, amps2)

plt.savefig("sip")
