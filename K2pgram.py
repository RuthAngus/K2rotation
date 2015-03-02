import numpy as np

# loop over eval_freq to compute a periodogram
def K2pgram(x, y, fs, AT, ATA):
    pgram = np.zeros_like(fs)
    amps2 = np.zeros_like(fs)
    for i, f in enumerate(fs):
        pgram[i], amps2[i] = eval_freq(x, y, f, AT, ATA)
    return pgram, amps2

# calculate periodogram by just updating sin and cos parts of the  matrices
def eval_freq(x, y, f, AT, ATA, compute_trends=False):
    arg = 2*np.pi*f*x
    AT[-2, :] = np.sin(arg)
    AT[-1, :] = np.cos(arg)

    # AT.shape = (153, nt)
    # shape: (151, nt) * (nt, 2) -> (151, 2)
    v = np.dot(AT[:-2, :], AT[-2:, :].T)
    ATA[:-2, -2:] = v
    ATA[-2:, :-2] = v.T

    # AT[-2:, :].shape = (2, nt)
    # (2, nt), (nt, 2)
    ATA[-2:, -2:] = np.dot(AT[-2:, :], AT[-2:, :].T)
    w = np.linalg.solve(ATA, np.dot(AT, y))

    S = np.linalg.inv(ATA)[-2:, -2:]
    s2n = np.dot(w[-2:], np.linalg.solve(S, w[-2:]))

    if compute_trends:
        trends = np.dot(w[:-2], AT[:-2])
        return np.sum(w[-2:]**2), trends
    return np.sum(w[-2:]**2), s2n

# loop over eval_freq to compute a periodogram
def K2pgram_basis(x, y, fs, AT, ATA):
    pgram = np.zeros_like(fs)
    for i, f in enumerate(fs):
        pgram[i], s2n = eval_freq_basis(x, y, f, AT, ATA)
    return pgram, s2n

def eval_freq_basis(x, y, f, AT, ATA, compute_trends=False):
    arg = 2*np.pi*f*x
    AT[-2, :] = np.sin(arg)
    AT[-1, :] = np.cos(arg)

    ATA = np.dot(AT, AT.T)
    w = np.linalg.solve(ATA, np.dot(AT, y))

    return np.sum(w[-2:]**2), w

# loop over eval_freq to compute a periodogram.
# f is the frequency you already found
def K2pgram2(x, y, f1, fs, AT, ATA):
    pgram = np.zeros_like(fs)
    for i, f2 in enumerate(fs):
        if f1 != f2:
            pgram[i] = eval_2nd_freq(x, y, f1, f2, AT, ATA)
    return pgram

# calculate periodogram by just updating sin and cos parts of the  matrices
def eval_2nd_freq(x, y, f1, f2, AT, ATA, compute_trends=False):
    arg1 = 2*np.pi*f1*x
    arg2 = 2*np.pi*f2*x
    AT[-4, :] = np.sin(arg1)
    AT[-3, :] = np.cos(arg1)
    AT[-2, :] = np.sin(arg2)
    AT[-1, :] = np.cos(arg2)

    # AT.shape = (155, nt)
    # shape: (151, nt) * (nt, 4) -> (151, 4)
    v = np.dot(AT[:-4, :], AT[-4:, :].T)
    ATA[:-4, -4:] = v
    ATA[-4:, :-4] = v.T

    # AT[-4:, :].shape = (4, nt)
    # (4, nt), (nt, 4)
    ATA[-4:, -4:] = np.dot(AT[-4:, :], AT[-4:, :].T)
    w = np.linalg.solve(ATA, np.dot(AT, y))

    S = np.linalg.inv(ATA)[-4:, -4:]
    s2n = np.dot(w[-4:], np.linalg.solve(S, w[-4:]))

    if compute_trends:
        trends = np.dot(w[:-4], AT[:-4])
        return np.sum(w[-4:]**2), trends
    return s2n
