import numpy as np

# loop over eval_freq to compute a periodogram
def K2pgram(x, y, fs, AT, ATA):
    pgram = np.zeros_like(fs)
    for i, f in enumerate(fs):
        pgram[i] = eval_freq(x, y, f, AT, ATA)
    return pgram

# calculate periodogram by just updating sin and cos parts of the  matrices
def eval_freq(x, y, f, AT, ATA, compute_trends=False):
    arg = 2*np.pi*f*x
    AT[-2, :] = np.sin(arg)
    AT[-1, :] = np.cos(arg)
    v = np.dot(AT[:-2, :], AT[-2:, :].T)
    ATA[:-2, -2:] = v
    ATA[-2:, :-2] = v.T
    ATA[-2:, -2:] = np.dot(AT[-2:, :], AT[-2:, :].T)
    w = np.linalg.solve(ATA, np.dot(AT, y))

    S = np.linalg.inv(ATA)[-2:, -2:]
    s2n = np.dot(w[-2:], np.linalg.solve(S, w[-2:]))

    if compute_trends:
        trends = np.dot(w[:-2], AT[:-2])
        return np.sum(w[-2:]**2), trends
    return s2n
