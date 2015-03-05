# a script to download the asteroseismic target data

import numpy as np
import wget
import kplr
client = kplr.API()

data = np.genfromtxt("../../epic.csv", skip_header=1, delimiter=",").T
campaign = data[2]
epicid = data[0]
kepmag = data[24]
j = data[41]
h = data[43]
k = data[45]

l = (campaign==1) * (9.5 < kepmag) * (kepmag < 11) * (j-h > .75) * (h-k > .1)
eid, km = epicid[l], kepmag[l]

np.savetxt("giant_list.txt", np.transpose((eid, km)))

for i, star in enumerate(eid):
    print i, "of", len(eid)
    s = str(int(star))
    base_url = "http://bbq.dfm.io/ketu/lightcurves/c1"
    url = "%s/20%s00000/%s000/ktwo%s-c01_lpd-lc.fits" \
            % (base_url, s[2:4], s[4:6], s)
    print url
    filename = wget.download(url)
