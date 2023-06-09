import os
import random
from tqdm import tqdm
import numpy as np
from scipy.signal import convolve
import astropy.io.fits as fits
import argparse
import lib_fits


def main():
    # Parameters
    parser = argparse.ArgumentParser(description='This is PDRL method. Results for all iterations are saved.')
    parser.add_argument('thresh_img', type=str, help='input counts map file for PDRL')
    parser.add_argument('expmap', type=str, help='input exposure map file for PDRL')
    parser.add_argument('psfs', type=str, help='input psf file for each infile position (.npz)')
    parser.add_argument('outdir', type=str, help='output dir for each iteration of PDRL results')
    parser.add_argument('--num_iter', type=int, default=200, help='default:200. the number of iterations for PDRL')
    parser.add_argument('--lambda_tv', type=float, default=0.002, help='default:0.002. the TV regularization lambda')
    parser.add_argument('--data_type', type=str, default='float64', help='default:float64. numpy data type (Note:If `killed` is returned, `float32` is useful to lower the memory.)')
    parser.add_argument("--im_deconv_0_flat", action='store_true', help='default:False. initial value for the 0th iteration of the ldrl method (`False` mean the input file)')
    parser.add_argument("--poisson_err", action='store_true', help='default:False. add a Poisson distribution random number according to the input count file for each iteration')
    parser.add_argument("--boundary_px", type=int, default=1, help='default:1. Randomly select a PSF from `boundary_px` pixels near the boundary of the PSFs')
    args = parser.parse_args()

    # Loads the thresh.img and expmap.
    thresh_img, wcs = lib_fits.load_fits(infile=args.thresh_img, data_type=args.data_type)
    expmap, _ = lib_fits.load_fits(infile=args.expmap, data_type=args.data_type)

    # Loads all PSFs as numpy array.
    print('PSF loading in progress.')
    psfs, psfs_bins = load_all_psf(args.psfs, args.data_type)

    # PDRL method.
    np.random.seed(seed=2023)
    print('PDRL method is start.')
    os.makedirs(args.outdir, exist_ok=True)
    pdrl = Position_Dependent_Richardson_Lucy(args.outdir, thresh_img, expmap, psfs, psfs_bins, args.num_iter, args.lambda_tv, wcs, args.data_type, args.im_deconv_0_flat, args.poisson_err, args.boundary_px)
    deconvs = pdrl.position_dependent_richardson_lucy()


def load_all_psf(psfs_infile, data_type):
    repro_psfs = np.load(psfs_infile)
    psfs = np.array(repro_psfs['repro_psfs'], dtype=data_type)
    psfs_bins = repro_psfs['psfs_bins'].tolist()
    return psfs, psfs_bins


class Position_Dependent_Richardson_Lucy:
    def __init__(self, outdir, thresh_img, expmap, psfs, psfs_bins, num_iter=200, lambda_tv=0.002, wcs=None, data_type='float64', im_deconv_0_flat=False, poisson_err=False, boundary_px=1):
        self.outdir = outdir
        self.thresh_img = thresh_img
        self.expmap = expmap + 1e-15  # Used to avoid 0 divisions
        self.psfs = psfs
        self.psfs_bins = psfs_bins
        self.num_iter = num_iter
        self.lambda_tv = lambda_tv
        self.wcs = wcs
        self.data_type = data_type
        self.im_deconv_0_flat = im_deconv_0_flat
        self.poisson_err = poisson_err
        self.boundary_px = boundary_px


    def _convolve_each_kernel(self, image, kernels, kernels_xgrid, kernels_ygrid):
        y_shape, x_shape = image.shape
        *_, k_y_shape, k_x_shape = kernels.shape
        pad_y, pad_x = (k_y_shape // 2, k_x_shape // 2)
        pad_conv = np.zeros([y_shape + pad_y*2, x_shape + pad_x*2], dtype=self.data_type)  # For convolution edge works
        y_range = range(self.boundary_px, y_shape-self.boundary_px)
        x_range = range(self.boundary_px, x_shape-self.boundary_px)
        for y in tqdm(y_range):
            for x in x_range:
                pad_conv[y:y+k_y_shape, x:x+k_x_shape] += image[y, x] * kernels[(y+kernels_ygrid[y, x])//self.psfs_bins, (x+kernels_xgrid[y, x])//self.psfs_bins]
        conv = pad_conv[pad_y:-pad_y, pad_x:-pad_x]  # Make array size equal to input size
        return conv


    def _convolve_each_kernel_dependent(self, image, kernels, kernels_xgrid, kernels_ygrid):
        y_shape, x_shape = image.shape
        *_, k_y_shape, k_x_shape = kernels.shape
        k_y_half, k_x_half = (k_y_shape // 2, k_x_shape // 2)
        conv = np.zeros_like(image, dtype=self.data_type)
        y_range = range(k_y_half+1+self.boundary_px, y_shape-k_y_half-1-self.boundary_px)
        x_range = range(k_x_half+1+self.boundary_px, x_shape-k_x_half-1-self.boundary_px)
        for y in tqdm(y_range):
            for x in x_range:
                conv[y, x] = np.sum(image[y-k_y_half:y+k_y_half+1, x-k_x_half:x+k_x_half+1] * kernels[(y+kernels_ygrid[y, x])//self.psfs_bins, (x+kernels_xgrid[y, x])//self.psfs_bins])
        return conv

    def _generate_poisson_img(self, thresh_img):
        poisson_func = lambda x: np.random.poisson(x)
        poisson_thresh_img = poisson_func(thresh_img)
        return poisson_thresh_img

    def position_dependent_richardson_lucy(self):
        if self.im_deconv_0_flat:
            im_deconv = np.full(self.thresh_img.shape, 0.5, dtype=self.data_type)
        else:
            im_deconv = np.copy(self.thresh_img) / self.expmap

        eps = 1e-15  # Small regularization parameter used to avoid 0 divisions

        for _iter in range(1, self.num_iter+1):
            print(f'iter {_iter} start.')
            print('Bottom calculation in progress.')
            psfs_xgrid = np.random.randint(-self.boundary_px, self.boundary_px+1, size=im_deconv.shape)
            psfs_ygrid = np.random.randint(-self.boundary_px, self.boundary_px+1, size=im_deconv.shape)
            conv = self._convolve_each_kernel(im_deconv, self.psfs, psfs_xgrid, psfs_ygrid) + eps

            if self.poisson_err:
                relative_blur = self._generate_poisson_img(self.thresh_img) / self.expmap / conv
            else:
                relative_blur = self.thresh_img / self.expmap / conv

            # TVregularization
            grad_x = np.gradient(im_deconv, axis=1)
            grad_y = np.gradient(im_deconv, axis=0)
            grad_norm = np.sqrt(grad_x**2+grad_y**2) + eps
            tv_reg = np.gradient(grad_x / grad_norm, axis=1) + np.gradient(grad_y / grad_norm, axis=0)

            print('Top calculation in progress.')
            im_deconv = im_deconv/(1-self.lambda_tv*tv_reg)*self._convolve_each_kernel_dependent(relative_blur, self.psfs, psfs_xgrid, psfs_ygrid)

            # Save the results of deconvolution for each iteration
            outfile = os.path.join(self.outdir, f'iter_{_iter:0>4}.fits')
            lib_fits.np2fits(np_array=im_deconv, outfile=outfile, wcs=self.wcs)
        
        return


if __name__ == '__main__':
    main()
