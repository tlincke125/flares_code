from datetime import datetime
import os
import re
from astropy.io import fits

###### DATA
# Get a list of dates active for a specific harpnumber
def get_dates(harpnum, root, sort = False):
    base = os.path.join(root, "magnetogram", "sharp_" + str(harpnum))
    assert os.path.exists(base)
    ret = []
    for i in os.listdir(base):
        if "Br" in i: # Testing for radial coords
            pattern =   re.compile(rf"""hmi\.sharp_cea_720s\.
                        (?P<region>[0-9]+)\.
                        (?P<date>[0-9]+\_[0-9]+)\_TAI\.Br\.fits""", re.VERBOSE)
            match = pattern.match(i)
            ret.append(datetime.strptime(match.group("date"), "%Y%m%d_%H%M%S"))
    return sorted(ret) if sort else ret


def get_data(harpnum, date, root, nantozero = False):
    date_str = date.strftime("%Y%m%d_%H%M%S")

    files = [("magnetogram", "Br"), ("magnetogram", "Bt"), ("magnetogram", "Bp"), ("continuum", "continuum")]
    data = []
    for f1, f2 in files:
        filename = os.path.join(root, f1, f"sharp_{harpnum}", f"hmi.sharp_cea_720s.{harpnum}.{date_str}_TAI.{f2}.fits")
        assert os.path.isfile(filename)
        with fits.open(filename) as hdul:
            #  hdul.verify('fix')
            d = np.rot90(hdul[1].data).copy()
            if nantozero:
                d = np.nan_to_num(d, 0.01)
            data.append(d)
    return data