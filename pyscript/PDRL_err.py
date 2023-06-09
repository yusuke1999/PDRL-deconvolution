import os
from tqdm import tqdm
import numpy as np
import astropy.io.fits as fits
import argparse
import lib_fits


def main():
    # Parameters
    parser = argparse.ArgumentParser(description='Calculate the uncertainty of PDRL method by the law of error propagation.')
    parser.add_argument('thresh_img', type=str, help='input counts map file for PDRL')
    parser.add_argument('expmap', type=str, help='input exposure map file for PDRL')
    parser.add_argument('psfs', type=str, help='input psf file for each infile position (.npz)')
    parser.add_argument('im_deconv', type=str, help='input a PDRL result to get the next iteration PDRL uncertainty')
    parser.add_argument('outfile', type=str, help='output file for PDRL error result')
    parser.add_argument('--data_type', type=str, default='float64', help='default:float64. numpy data type (Note:If `killed` is returned, `float32` is useful to lower the memory.)')
    args = parser.parse_args()

    # Loads the thresh.img, expmap and im_deconv.
    thresh_img, wcs = lib_fits.load_fits(infile=args.thresh_img, data_type=args.data_type)
    expmap, _ = lib_fits.load_fits(infile=args.expmap, data_type=args.data_type)
    im_deconv, _ = lib_fits.load_fits(infile=args.im_deconv, data_type=args.data_type)

    # Loads all PSFs as numpy array.
    print('PSF loading in progress.')
    psfs, psfs_bins = load_all_psf(args.psfs, args.data_type)

    # PDRL error method.
    print('Estimate the PDRL error method is start.')
    outdir = os.path.dirname(args.outfile)
    os.makedirs(outdir, exist_ok=True)
    pdrl_err = Position_Dependent_Richardson_Lucy_Err(args.outfile, thresh_img, expmap, psfs, psfs_bins, im_deconv, wcs, args.data_type)
    deconvs_err = pdrl_err.position_dependent_richardson_lucy_err()
    print('Done.')


def load_all_psf(psfs_infile, data_type):
    repro_psfs = np.load(psfs_infile)
    psfs = np.array(repro_psfs['repro_psfs'], dtype=data_type)
    psfs_bins = repro_psfs['psfs_bins'].tolist()
    return psfs, psfs_bins


class Position_Dependent_Richardson_Lucy_Err:
    def __init__(self, outfile, thresh_img, expmap, psfs, psfs_bins, im_deconv, wcs=None, data_type='float64'):
        self.outfile = outfile
        self.thresh_img = thresh_img
        self.expmap = expmap + 1e-15  # Used to avoid 0 divisions
        self.psfs = psfs
        self.psfs_bins = psfs_bins
        self.im_deconv = im_deconv
        self.wcs = wcs
        self.data_type = data_type


    def _convolve_each_kernel(self, image, kernels):
        y_shape, x_shape = image.shape
        *_, k_y_shape, k_x_shape = kernels.shape
        pad_y, pad_x = (k_y_shape // 2, k_x_shape // 2)
        pad_conv = np.zeros([y_shape + pad_y*2, x_shape + pad_x*2], dtype=self.data_type)  # For convolution edge works
        y_range = range(y_shape)
        x_range = range(x_shape)
        for y in tqdm(y_range):
            for x in x_range:
                pad_conv[y:y+k_y_shape, x:x+k_x_shape] += image[y, x] * kernels[y//self.psfs_bins, x//self.psfs_bins]
        conv = pad_conv[pad_y:-pad_y, pad_x:-pad_x]  # Make array size equal to input size
        return conv


    def _convolve_each_kernel_dependent(self, image, kernels):
        y_shape, x_shape = image.shape
        *_, k_y_shape, k_x_shape = kernels.shape
        k_y_half, k_x_half = (k_y_shape // 2, k_x_shape // 2)
        conv = np.zeros_like(image, dtype=self.data_type)
        y_range = range(k_y_half+1, y_shape-k_y_half-1, 1)
        x_range = range(k_x_half+1, x_shape-k_x_half-1, 1)
        for y in tqdm(y_range):
            for x in x_range:
                conv[y, x] = np.sum(image[y-k_y_half:y+k_y_half+1, x-k_x_half:x+k_x_half+1] * kernels[y//self.psfs_bins, x//self.psfs_bins])
        return conv

    
    def position_dependent_richardson_lucy_err(self):
        psfs_sq = self.psfs ** 2
        expmap_sq = self.expmap ** 2

        eps = 1e-15  # Small regularization parameter used to avoid 0 divisions
        
        print('Bottom calculation in progress.')
        conv = self._convolve_each_kernel(self.im_deconv, self.psfs) + eps
        relative_blur = self.thresh_img / expmap_sq / conv**2

        print('Top calculation in progress.')
        im_deconv_err = self.im_deconv*np.sqrt(self._convolve_each_kernel_dependent(relative_blur, psfs_sq))

        # Save the results of deconvolution err
        lib_fits.np2fits(np_array=im_deconv_err, outfile=self.outfile, wcs=self.wcs)
        
        return


if __name__ == '__main__':
    main()
