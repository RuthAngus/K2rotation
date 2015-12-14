import numpy as np
import matplotlib.pyplot as plt
# import h5py
# from params import colours, params
from K2misc import load_K2_data, peak_detect
import glob
from SIP import SIP, eval_freq


def prewhiten(x, y, f, basis):
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
    return y + wave


# inject a sine wave
def inject(x, y, f, a):
    """
    takes x, y, frequency and amplitude
    returns y + wave
    """
    wave = a * np.sin(f*2*np.pi*x)
    return y + wave


def inject_and_recover(fname, fs, ifs, a_s, plot=False):
    """
    fname: "path/to/file/filename" - should be a K2 fits file
    fs: the SIP freq array
    ifs: an array of frequencies to inject
    a_s: an array of amplitudes to inject
    returns the recovered frequencies and amplitudes
    """
    x, y, basis = load_K2_data(fname)
    recovered = np.zeros_like(ifs)  # array of freq of the highest peak
    recovered_amps = np.zeros_like(ifs)
    for i, f in enumerate(ifs):  # loop over frequencies
        for j, a in enumerate(a_s):
            print(j, "of", len(a_s), "amplitudes")
            print(i, "of", len(ifs), "frequencies")
            iy = y + inject(x, y, f, a)
            s2n, amps2, w = SIP(x, iy, basis, fs)
            peak_f, peak_a = peak_detect(fs, amps2)
            recovered[i] = peak_f
            recovered_amps[i] = peak_a
            if plot:
                plt.clf()
                plt.plot(fs, amps2)
                plt.axvline(peak_f, color="r")
                plt.savefig("{0}{1}".format(j, i))
    return recovered, recovered_amps

if __name__ == "__main__":

    # load the data
    fnames = glob.glob("data/ktwo*fits")

    # inject and recover a bunch of sinewaves
    fs = np.arange(10, 300, 1e-1) * 1e-6
    ifs = np.arange(50, 250, 40) * 1e-6  # the injections frequencies
    a_s = np.arange(1e-5, 1e-3, 5e-4)  # the injection amplitudes
    recovered, recovered_amps = inject_and_recover(fnames[0], fs, ifs, a_s,
                                                   plot=True)
    print(ifs)
    print(recovered)
    print(a_s)
    print(recovered_amps)

    plt.clf()
    plt.hist(ifs - recovered)
    plt.savefig("test")
