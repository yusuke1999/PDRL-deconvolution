import os
import shutil
import datetime
from tqdm import tqdm
import numpy as np
from astropy.wcs.utils import pixel_to_skycoord
from astropy.wcs.utils import skycoord_to_pixel
import argparse
import lib_fits


def main():
    # Parameters
    parser = argparse.ArgumentParser(description='Move the center of PSFs to the image center and adjust all PSFs to the same size.')
    parser.add_argument('infile', type=str, help='input counts map file for made psfs (Note: not evt file)')
    parser.add_argument('psfs_dir', type=str, help='input raw psf dir for each location')
    parser.add_argument('outfile', type=str, help='output all psf .npz file')
    parser.add_argument('psfs_bins', type=int, help='spacing of pixels you made each position psf')
    args = parser.parse_args()


    # Get the size of the image from `infile` you made psfs.
    image, wcs_img = lib_fits.load_fits(infile=args.infile, data_type='float64')

    # PSFs coordinates.
    ys = np.arange(args.psfs_bins // 2, image.shape[0], args.psfs_bins)
    xs = np.arange(args.psfs_bins // 2, image.shape[1], args.psfs_bins)

    # Make outfile path (If outfile contains dir).
    outdir = os.path.dirname(args.outfile)
    if outdir != '':
        os.makedirs(outdir, exist_ok=True)

    # Run reprocess PSFs.
    repro_psfs = Reprocess_PSFs(args.psfs_dir, args.outfile, args.psfs_bins, ys, xs, wcs_img)
    repro_psfs.reprocess_psfs()


class Reprocess_PSFs:
    def __init__(self, psf_dir, outfile, psfs_bins, ys, xs, wcs_img):
        self.psf_dir = psf_dir
        self.outfile = outfile
        self.psfs_bins = psfs_bins
        self.ys = ys
        self.xs = xs
        self.wcs_img = wcs_img


    def _shift_psf_center2image_center(self, psf, y, x, wcs_psf):
        point_radec = pixel_to_skycoord(x, y, self.wcs_img)
        center_x, center_y = skycoord_to_pixel(point_radec, wcs_psf)
        center_x, center_y = (int(np.round(center_x)), int(np.round(center_y)))
        y_shift, x_shift = (psf.shape[0]-1-2*center_y, psf.shape[1]-1-2*center_x)
        if x_shift > 0:
            pad_xs = (x_shift, 0)
        else:
            pad_xs = (0, -x_shift)
        if y_shift > 0:
            pad_ys = (y_shift, 0)
        else:
            pad_ys = (0, -y_shift)
        psf = np.pad(psf, [pad_ys, pad_xs], 'constant')
        return psf


    def _pad_psf2same_size(self, psf_max_y_shape, psf_max_x_shape, psf):
        pad_x = (psf_max_x_shape - psf.shape[1]) // 2
        pad_y = (psf_max_y_shape - psf.shape[0]) // 2
        psf = np.pad(psf, [(pad_y,pad_y), (pad_x,pad_x)], 'constant')
        return psf


    def _adjust_psf_center(self, temp_dir):
        psf_y_shapes = np.zeros(self.ys.shape[0]*self.xs.shape[0], dtype='int64')
        psf_x_shapes = np.zeros(self.ys.shape[0]*self.xs.shape[0], dtype='int64')
    
        index = 0
        for y in tqdm(self.ys):
            for x in self.xs:
                psf_infile =  os.path.join(self.psf_dir, f'{y:0>4}_{x:0>4}.psf')
                psf, wcs_psf = lib_fits.load_fits(infile=psf_infile, data_type='float64')

                shifted_psf = self._shift_psf_center2image_center(psf, y, x, wcs_psf)
                psf_y_shapes[index] = shifted_psf.shape[0]
                psf_x_shapes[index] = shifted_psf.shape[1]
                save_file = os.path.join(temp_dir, f'{y:0>4}_{x:0>4}')
                np.save(save_file, shifted_psf)
                index += 1
        psf_max_y_shape = np.max(psf_y_shapes)
        psf_max_x_shape = np.max(psf_x_shapes)

        return psf_max_y_shape, psf_max_x_shape


    def _adjust_psf_size(self, temp_dir, psf_max_y_shape, psf_max_x_shape):
        adjusted_psfs = np.zeros((self.ys.shape[0]+1, self.xs.shape[0]+1, psf_max_y_shape, psf_max_x_shape), dtype='float64')
        edge_y, edge_x = (self.ys[-1], self.xs[-1])
        for y in tqdm(self.ys):
            for x in self.xs:
                centered_psf = np.load(os.path.join(temp_dir, f'{y:0>4}_{x:0>4}.npy'))
                adjusted_psfs[y//self.psfs_bins, x//self.psfs_bins] = self._pad_psf2same_size(psf_max_y_shape, psf_max_x_shape, centered_psf)
                # For edge x psf
                centered_psf = np.load(os.path.join(temp_dir, f'{edge_y:0>4}_{x:0>4}.npy'))
                adjusted_psfs[edge_y//self.psfs_bins+1, x//self.psfs_bins] = self._pad_psf2same_size(psf_max_y_shape, psf_max_x_shape, centered_psf)
            # For edge y psf
            centered_psf = np.load(os.path.join(temp_dir, f'{y:0>4}_{edge_x:0>4}.npy'))
            adjusted_psfs[y//self.psfs_bins, edge_x//self.psfs_bins+1] = self._pad_psf2same_size(psf_max_y_shape, psf_max_x_shape, centered_psf)

        return adjusted_psfs


    def reprocess_psfs(self):
        # Make dir for shifted psf.
        time_stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        temp_dir = f'{time_stamp}_shifted_psf'
        os.makedirs(temp_dir)

        print('Start to shift PSF to center.')
        psf_max_y_shape, psf_max_x_shape = self._adjust_psf_center(temp_dir)

        print(f'Start to all PSFs to the same size as (psf_max_y_shape, psf_max_x_shape)=({psf_max_y_shape}, {psf_max_x_shape}).')
        adjusted_psfs = self._adjust_psf_size(temp_dir, psf_max_y_shape, psf_max_x_shape)
        np.savez_compressed(self.outfile, repro_psfs=adjusted_psfs, psfs_bins=self.psfs_bins)

        # Rremove shifted psf dir.
        shutil.rmtree(temp_dir)

        return


if __name__ == '__main__':
    main()
