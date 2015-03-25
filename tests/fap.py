import numpy as np

# calculate the false alarm probability
def fap(x, y, basis, fs, N, plot=False, sig=False):
    amp2s, s2n, _ = K2pgram(x, y, basis, fs)  # 1st pgram
    if sig: power = s2n
    else: power = amp2s
    mf, ms2n = peak_detect(fs, power)  # find peak
    AT = np.concatenate((basis, np.ones((3, len(y)))), axis=0)
    ATA = np.dot(AT, AT.T)
    # compute trends
    _, _, trends = eval_freq(x, y, mf, AT, ATA, compute_trends=True)
    if plot:
        plt.clf()
        plt.plot(1./fs, power, "k")
    peak_heights = []
    for n in range(N):
        detrended_y = y - trends  # remove trends
        detrended_y = np.random.choice(detrended_y, len(y))  # shuffle
        # add trends back in
        amp2s, s2n, _ = K2pgram(x, detrended_y + trends, basis, fs)
        if sig: power = s2n
        else: power = amp2s
        mx, my = peak_detect(fs, power)
        peak_heights.append(my)
        if plot:
            plt.plot(1./fs, power, alpha=.2)
    fap95 = np.percentile(peak_heights, 95)
    fap90 = np.percentile(peak_heights, 90)
    fap85 = np.percentile(peak_heights, 85)
    fap50 = np.percentile(peak_heights, 50)
    if plot:
        plt.axhline(fap95, color=".5")
        plt.savefig("fap")
#     print fap95, fap90, fap85, fap50
    return fap95, fap90, fap85, fap50
