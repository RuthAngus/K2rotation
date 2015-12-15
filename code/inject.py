import numpy as np
import matplotlib.pyplot as plt
from K2misc import load_K2_data, peak_detect, detect_all_peaks
from SIP import SIP, eval_freq
import sys


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
    return y - wave


# inject a sine wave
def inject(x, y, f, a):
    """
    takes x, y, frequency and amplitude
    returns y + wave
    """
    wave = a * np.sin(f*2*np.pi*x)
    return y + wave


def inj(fname, ifs, a_s):
    """
    inject sine wave into lc
    fname: "path/to/file/filename" - should be a K2 fits file
    ifs: an array of frequencies to inject
    a_s: an array of amplitudes to inject
    fs: the SIP freq array
    returns the recovered frequencies and amplitudes
    """
    truths, ids = [], []
    x, y, basis = load_K2_data(fname)
    id = 0
    for i, f in enumerate(ifs):  # loop over frequencies
        for j, a in enumerate(a_s):  # loop over amplitudes
            print(id+1, "of", len(a_s) * len(ifs))
            iy = y + inject(x, y, f, a)  # inject sinewave
            np.savetxt("injections/{0}.txt".format(str(id).zfill(5)),
                       np.transpose((x, iy)))
            truths.append(f)  # save the truths and the ids
            ids.append(id)
            id += 1
    np.savetxt("truths.txt", np.transpose((np.array(truths), np.array(ids))))


def recover_SIP(template_id, inj_fnames, fs, oa2, start, stop, plot=True,
                subtract_baseline=True):
    """
    Find frequency and amplitude of the highest peak in the SIP
    fname: the name of the target used for injection
    inj_fnames: the names of the injected lc files
    fs: a grid of frequencies
    oamp2s: the original sip, before sine wave injection
    """
    recovered, recovered_amps = [], []  # array of freq of the highest peak
    _, _, basis = load_K2_data(template_id)  # load original lc
    for i, fname in enumerate(inj_fnames[start:stop]):  # loop over injections
        print(i, "of", len(inj_fnames))
        ix, iy = \
            np.genfromtxt("injections/{0}.txt".format(str(fname).zfill(5))).T
        print("computing SIP")
        s2n, amps2, w = SIP(ix, iy, basis, fs)  # compute a sip
        if subtract_baseline:  # subtract the original sip
            amps2 = amps2 - oa2
        peak_f, peak_a = peak_detect(fs, amps2)  # find the highest peak
        print(peak_f)
        recovered.append(peak_f)
        recovered_amps.append(peak_a)
        if plot:
            plt.clf()
            plt.plot(fs, amps2)
            plt.axvline(peak_f, color="r")
            plt.savefig("{0}".format(str(fname).zfill(5)))

    # save the results
    rf, ra = np.array(recovered), np.array(recovered_amps)
    data = np.vstack((inj_fnames, rf, ra))
    np.savetxt("recovered_{0}_{1}.txt".format(start, stop), data.T)
    return rf, ra


def iterative_prewhiten(N):
    # prewhiten 10 times
    x, y, basis = load_K2_data(fname)
    fs = np.arange(10, 300, 1e-1) * 1e-6
    s2n, amp2s, w = SIP(x, y, basis, fs)  # calculate SIP

    # find the top N peaks
    peak_fs, peaks_as = detect_all_peaks(fs, amp2s)
    peakN = np.sort(peaks_as)[-N:]
    peak_f = np.array([peak_fs[peaks_as == i][0] for i in peakN])

    # prewhiten
    for peak in peak_f[::-1]:  # highest peak to lowest
        y = prewhiten(x, y, peak, basis)
    return y

if __name__ == "__main__":

    # load the data and compute initial sip
    fname = "data/ktwo201121245-c01_lpd-lc.fits"
    x, y, basis = load_K2_data(fname)
    fs = np.arange(10, 300, 1e-1) * 1e-6
    s2n, amps2, w = SIP(x, y, basis, fs)

    # injection and recovery parameters
    ifs = np.random.uniform(10e-6, 300e-6, 10)  # the injections frequencies
    a_s = 1**np.random.uniform(-6, -4, 2)  # the injection amplitudes
    injection_fnames = range(len(ifs) * len(a_s))  # names for file saves

    # parallelisation parameters
    start = sys.argv(1)
    stop = sys.argv(2)

    inj(fname, ifs, a_s)
    recovered, recovered_amps = recover_SIP(fname, injection_fnames, fs,
                                            amps2, start, stop)
