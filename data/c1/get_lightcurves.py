import numpy as np
import wget


def get_all_C1_lightcurves():
    fnames = np.genfromtxt("C1asciiwget.sh", dtype=str).T
    urls = fnames[2, :]
    for i, url in enumerate(urls):
        print(i, "of", len(urls))
        wget.download(url)


def get_one_C1_lightcurve(id):
    a = "http://archive.stsci.edu/missions/hlsp/k2sff/c01"
    b = "{0}00000".format(id[:4])
    c = id[4:]
    d = "hlsp_k2sff_k2_lightcurve_"
    e = "{0}-c01_kepler_v1_llc-default-aper.txt".format(id)
    url = "{0}/{1}/{2}/{3}{4}".format(a, b, c, d, e)
    wget.download(url)


if __name__ == "__main__":
    get_one_C1_lightcurve("201183188")
