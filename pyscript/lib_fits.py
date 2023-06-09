import os
import warnings
import numpy as np
import astropy.io.fits as fits
from astropy.wcs import WCS
from astropy.wcs import FITSFixedWarning


def load_fits(infile, data_type):
    """Parameters
    infile: fits file
    data_type: numpy `dtype`
    """
    with fits.open(infile) as hdulist:
        hdu = hdulist[0]
        data = hdu.data
        data = np.array(data, dtype=data_type)
        with warnings.catch_warnings():
            # Ignore a warning on using DATE-OBS in place of MJD-OBS
            warnings.filterwarnings('ignore', message="'datfix' made the change",
                                    category=FITSFixedWarning)
            # Ignore a warning on Deprecated keywords 
            warnings.filterwarnings('ignore', message="RADECSYS= 'ICRS '",
                                    category=FITSFixedWarning)
            wcs = WCS(hdu.header)
    return data, wcs



def np2fits(np_array, outfile, wcs=None):
    """Parameters
    np_array: numpy 2D array
    outfile: fits file
    """
    # Make outfile path.(If outfile contains dir.)
    outdir = os.path.dirname(outfile)
    if outdir != '':
        os.makedirs(outdir, exist_ok=True)

    if wcs is None:
        hdu = fits.PrimaryHDU(data=np_array)
    else:
        header = wcs.to_header()
        hdu = fits.PrimaryHDU(data=np_array, header=header)
    hdu.writeto(outfile, overwrite=True)
