import numpy as np
import matplotlib.pyplot as plt
import glob
import h5py
from SIP import SIP, eval_freq
from gatspy.periodic import LombScargle
import sys
from fit_all_the_light_curves import load_lc, reconstruct_fake_lc
import time

plotpar = {'axes.labelsize': 20,
           'text.fontsize': 20,
           'legend.fontsize': 20,
           'xtick.labelsize': 20,
           'ytick.labelsize': 20,
           'text.usetex': True}
plt.rcParams.update(plotpar)

def peak_detect(x, y):
    peaks = np.array([i for i in range(1, len(x)-1) if y[i-1] < y[i] and
                     y[i+1] < y[i]])
    l = y[peaks] == max(y[peaks])
    mx, my = x[peaks][l][0], y[peaks][l][0]
    return mx, my


# grid over amplitudes (the K2 pgram step takes the time)
def grid_over_amps(basis, flux, raw_x, raw_y, truth, fs, amps, true_a,
                   flag, n, plot=False, raw=False, random_amps=True):

    # find the threshold level
    _, initial_pgram, _ = SIP(raw_x, raw_y, basis, fs)
    mx, threshold = peak_detect(fs, initial_pgram)

    K2P, rawP, K2a, rawa = [], [], [], []
    alla, allp = [], []
    all_results = []
    for i, a in enumerate(amps):
        if random_amps:
            if flag == "r":
                a = 10**(np.random.uniform(np.log10(1e-5), np.log10(1e-3)))
            elif flag == "a":
                a = 10**(np.random.uniform(np.log10(1e-5), np.log10(1e-3)))

        tf = 1./truth
        print("period = ", truth)

        # add lcs together
#         plt.clf()
#         plt.plot(flux * a, "k.")
        noise = np.random.randn(len(flux))*50*13**.5*1e-6
#         print 50*13**.5*1e-6, a
        fx = flux * a + noise
#         plt.plot(fx, "r.")
#         plt.savefig("Test")
#         assert 0
        y = fx + raw_y
        SN = np.var(fx) / np.var(raw_y)

        if flag == "r":
            threshold = .1
        elif flag == "a":
            threshold = .1

#         # Calculate time
#         start = time.time()
#         amp2s, s2n, w = SIP(raw_x, y, basis, fs[:1000])
#         end = time.time()
#         print("SIP time = ", (end-start), "s")
#         print("for", len(y), "data points and", len(fs), "freqs")

        # calculate SIP
        amp2s, s2n, w = SIP(raw_x, y, basis, fs)
        pgram = s2n
        best_f, best_pgram = peak_detect(fs, pgram)  # find peaks
        print("recovered period", 1./best_f)
        s = 0  # success indicator
        alla.append(a)
        allp.append(truth)
        all_results.append(best_f)

        print(tf-threshold*tf,  best_f, tf+threshold*tf)
        if tf-threshold*tf < best_f and best_f < tf+threshold*tf:
            K2P.append(truth)
            K2a.append(a)
            print("success!", "\n")
            s = 1

        # calculate periodogram of raw light curve
        y = np.array([_y.astype("float64") for _y in y])
        raw_x = np.array([_raw_x.astype("float64") for _raw_x in raw_x])

#         # Calculate time
#         start = time.time()
#         model = LombScargle().fit(raw_x, y, np.ones_like(y)*1e-5)
#         end = time.time()
#         print("SIP time = ", (end-start)*1e3, "ms")
#         print("for", len(y), "data points and", len(fs), "freqs")
#         assert 0

        model = LombScargle().fit(raw_x, y, np.ones_like(y)*1e-5)
        period = 1. / fs
        pg = model.periodogram(period)
        best_f2, best_pg2 = peak_detect(fs, pg)

        if tf-threshold*tf < best_f2 and best_f2 < tf+threshold*tf:
            rawP.append(truth)
            rawa.append(a)

        if plot:
            plt.clf()
            plt.subplot(2, 1, 1)
            plt.plot(raw_x, y, "k.")
            plt.plot(raw_x, fx, color="g")
            plt.title("$\mathrm{Amp = %s, P = %.3f}$" % (a, (1./tf)))
            plt.subplot(2, 1, 2)
            plt.axvline(best_f, color="r", linestyle="-")
            plt.axvline(tf, color="k", linestyle="--")
            print("best f = ", best_f)
            print("true f = ", tf)
            print(tf-threshold*tf, tf+threshold*tf)
            c = "b"
            if s == 1:
                c = "m"
            plt.plot(fs, pgram, color=c,
                     label="$\mathrm{SIP$}")
            plt.savefig("../injections/sine/%s_%s_result_%s"
                        % (str(n).zfill(2), str(i).zfill(2), flag))
            # n is the period index, i is the amplitude index
            print("%s_%s_result_%s" % (str(n).zfill(2), str(i).zfill(2),
                                       flag))
#             raw_input('enter')
    return np.array(K2a), np.array(K2P), np.array(rawa), np.array(rawP), \
            np.array(alla), np.array(allp), np.array(all_results)

# add simulated to real light curves and grid over periods
def grid_over_periods(basis, raw_x, raw_y, true_p, fs, true_a, fnames, flag):
    K2_amps, K2_Ps, raw_amps, raw_Ps = [], [], [], []
    ar = []
    allas, allps = [], []
    for i, fname in enumerate(fnames):
        print(fname)
        print(true_p[i])
        time, flux = np.genfromtxt(fname).T
        K2a, K2P, rawa, rawP, alla, allp, all_results = \
                grid_over_amps(basis, flux, raw_x, raw_y, true_p[i], fs,
                               amps, true_a[i], flag, i, plot=False)
        K2_amps.append(K2a)
        raw_amps.append(rawa)
        K2_Ps.append(K2P)
        raw_Ps.append(rawP)
        allas.append(alla)
        allps.append(allp)
        ar.append(all_results)
    K2_amps = np.array([j for i in K2_amps for j in i])
    K2_Ps = np.array([j for i in K2_Ps for j in i])
    raw_amps = np.array([j for i in raw_amps for j in i])
    raw_Ps = np.array([j for i in raw_Ps for j in i])
    allas = np.array([j for i in allas for j in i])
    allps = np.array([j for i in allps for j in i])
    ar = np.array([j for i in ar for j in i])

    f = h5py.File("../injections/sine/histogram_%s_%s_%s.h5" % (start, stop,
                  flag), "w")
    K2data = f.create_dataset("K2", (len(K2_amps), 2))
    K2data[:, 0] = K2_amps
    K2data[:, 1] = K2_Ps
    rawdata = f.create_dataset("raw", (len(raw_amps), 2))
    rawdata[:, 0] = raw_amps
    rawdata[:, 1] = raw_Ps
    f.close()
    f = h5py.File("../injections/sine/truths_%s_%s_%s.h5" % (start, stop,
                  flag), "w")
    K2data = f.create_dataset("K2", (len(allas), 3))
    K2data[:, 0] = allas
    K2data[:, 1] = allps
    K2data[:, 2] = ar
    f.close()
    return K2_amps, K2_Ps, raw_amps, raw_Ps

if __name__ == "__main__":

    # load example star
    path = "/export/bbq2/dfm/k2/web/lightcurves/c1/201100000/21000"
    fname = "ktwo201121245-c01_lpd-lc.fits"
    raw_x, y, l = load_lc("%s/%s" % (path, fname))

    raw_y = reconstruct_fake_lc()[l]

    # load basis
    with h5py.File("../data/c1.h5", "r") as f:
        basis = f["basis"][:150, l]

    # load injections and truths
    sine = True

    flag = str(sys.argv[1])  # r for rotation or a for asteroseismology
    if sine:
        fnames = glob.glob("../injections/sine/????_lc_%s.txt" % flag)
        fnames = np.sort(fnames)
        name, true_p = np.genfromtxt("../injections/sine/truth_%s.txt"
                                     % flag).T
        true_a = np.ones_like(true_p)
    else:
        fnames = glob.glob("../injections/*_lc.txt")
        name, true_p, true_a = np.genfromtxt("truth.txt").T

    # The sip grid
    if flag == "r":
        ps = np.linspace(.4, 50., 1000)
        fs = 1./ps
#         fs = np.linspace(1/50., 1/.4, 1000)
    elif flag == "a":
        fs = np.linspace(2./4., 26., 5000)

    # this is just a place holder, amps are random
    amps = 10**(np.linspace(np.log10(1e-5), np.log10(1e-3), 20))

    # for parallelisation, provide the starting and stopping indices
    start = int(sys.argv[2])
    stop = int(sys.argv[3])
    fnames = fnames[start:stop]
    true_p = true_p[start:stop]
    true_a = true_a[start:stop]

    # calculate the 2d histogram of completeness over period and amplitude
    K2_amps, K2_Ps, raw_amps, raw_Ps = grid_over_periods(basis, raw_x,
                                                         raw_y, true_p, fs,
                                                         true_a, fnames,
                                                         flag)
