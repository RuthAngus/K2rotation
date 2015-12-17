from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from K2misc import load_K2_data, peak_detect, detect_all_peaks
from SIP import SIP, eval_freq
import sys
from multiprocessing import Pool

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


def iterative_prewhiten(N):
    # prewhiten 10 times
    x, y, basis, _ = load_K2_data(fname)
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


# inject a sine wave
def inject(x, y, f, a):
    """
    takes x, y, frequency and amplitude
    returns y + wave
    """
    wave = a * np.sin(f*2*np.pi*x)
    return y + wave


def inj(fname, ifs, a_s, rotation=False):
    """
    inject sine wave into lc
    fname: "path/to/file/filename" - should be a K2 fits file
    ifs: an array of frequencies to inject
    a_s: an array of amplitudes to inject
    fs: the SIP freq array
    if rotation == True then save the output with a _r at the end
    returns the recovered frequencies and amplitudes
    """
    N = len(ifs)
    true_f, true_a = np.zeros_like(ifs), np.zeros_like(a_s)
    x, y, basis, _ = load_K2_data(fname)
    for i, f in enumerate(ifs):  # loop over frequencies
        print(i, "of", N)
        print("injection frequency = ", f)
        iy = y + inject(x, y, f, a_s[i])  # inject sinewave
        if rotation:
            np.savetxt("injections/{0}_r.txt".format(str(i).zfill(5)),
                       np.transpose((x, iy)))
        else:
            np.savetxt("injections/{0}.txt".format(str(i).zfill(5)),
                       np.transpose((x, iy)))
        true_f[i] = f  # save the truths and the ids
        true_a[i] = a_s[i]
    data = np.vstack((np.arange(N), true_f, true_a))
    if rotation:
        np.savetxt("truths_r.txt", data.T)
    else:
        np.savetxt("truths.txt", data.T)

def recover_SIP(template_id, inj_fnames, fs, oa2, start, stop, plot=False,
                subtract_baseline=True, rotation=False):
    """
    Find frequency and amplitude of the highest peak in the SIP
    fname: the name of the target used for injection
    inj_fnames: the names of the injected lc files
    fs: a grid of frequencies
    oamp2s: the original sip, before sine wave injection
    """
    recovered, recovered_amps = [], []  # array of freq of the highest peak
    _, _, basis, _ = load_K2_data(template_id)  # load original lc
    for i, fname in enumerate(inj_fnames[start:stop]):  # loop over injections
        print(i, "of", len(inj_fnames[start:stop]))
        if rotation:
            ix, iy = \
                np.genfromtxt("injections/{0}_r.txt".format(str(fname).zfill(5))).T
        else:
            ix, iy = \
                np.genfromtxt("injections/{0}.txt".format(str(fname).zfill(5))).T
        print("computing SIP")
        s2n, amps2, w = SIP(ix, iy, basis, fs)  # compute a sip
        if subtract_baseline:  # subtract the original sip

            astero_f = 2.131e-4  # this is a hack for normalising the sip.
                                 # specific to this target only!
            peaks_f, peaks_a = detect_all_peaks(fs, amps2)
            find_nearest_ind = lambda arr, val: np.abs(arr-val).argmin()
            astero_a2 = peaks_a[find_nearest_ind(peaks_f, astero_f)]
            opeaks_f, opeaks_a = detect_all_peaks(fs, oa2)
            astero_a1 = opeaks_a[find_nearest_ind(opeaks_f, astero_f)]
            ratio = astero_a1/astero_a2
            amps2 = amps2*ratio - oa2
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
    data = np.vstack((inj_fnames[start:stop], rf, ra))
    if rotation:
        np.savetxt("recovered_{0}_{1}_r.txt".format(str(start).zfill(5),
                   str(stop).zfill(5)), data.T)
    else:
        np.savetxt("recovered_{0}_{1}.txt".format(str(start).zfill(5),
                   str(stop).zfill(5)), data.T)
    return rf, ra


def recover_set(args):
    start, stop, N = args

    # load the data and compute initial sip
    fname = "data/ktwo201121245-c01_lpd-lc.fits"
    x, y, basis, _ = load_K2_data(fname)

    rotation = False
    fs = np.arange(10, 280, 1e-1) * 1e-6
    if rotation:
        fs = 1./(np.linspace(1., 70., 1000) * 24 * 3600)
    s2n, amps2, w = SIP(x, y, basis, fs)

    # parallelisation parameters
#     start = int(sys.argv[1])
#     stop = int(sys.argv[2])
#     N = int(sys.argv[3])

    # recover injections using SIP
    injection_fnames = range(N)  # names for file saves
    recovered, recovered_amps = recover_SIP(fname, injection_fnames, fs,
                                            amps2, start, stop, plot=False,
                                            rotation=rotation)


if __name__ == "__main__":
    N = 2000
    Nper = 100
    starts = np.arange(N/Nper) * Nper
    stops = (np.arange(N/Nper) + 1) * Nper
    Ns = np.ones_like(starts) * N
    arg = np.vstack((starts, stops, Ns)).T

    pool = Pool()
    results = pool.map(recover_set, arg)
