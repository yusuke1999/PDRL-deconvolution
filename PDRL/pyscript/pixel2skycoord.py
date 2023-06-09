import os
import warnings
from tqdm import tqdm
import numpy as np
from astropy.wcs.utils import pixel_to_skycoord
import argparse
import lib_fits

def main():
    # Parameters
    parser = argparse.ArgumentParser(description='Get each position skycoord and save .txt file for CIAO `simulate_psf`.')
    parser.add_argument('infile', type=str, help='input counts map file')
    parser.add_argument('outfile', type=str, help='output .txt file containing image_y,image_x,ra,dec')
    parser.add_argument('psf_bins', type=int, help='spacing of pixels to create psf')
    args = parser.parse_args()

    # Load infile.
    data, wcs = lib_fits.load_fits(infile=args.infile, data_type='float64')

    # Make outfile path (If outfile contains dir).
    outdir = os.path.dirname(args.outfile)
    if outdir != '':
        os.makedirs(outdir, exist_ok=True)

    # Export skycoordinates.
    y_shape, x_shape = data.shape
    with open(args.outfile, mode='w', encoding='utf-8') as f_out:
        print('Writting to file in progress.')
        for  y in tqdm(range(args.psf_bins // 2, y_shape, args.psf_bins)):
            for x in range(args.psf_bins // 2, x_shape, args.psf_bins):
                ra = repr(pixel_to_skycoord(x, y, wcs).ra).split(' ')[1]
                dec = repr(pixel_to_skycoord(x, y, wcs).dec).split(' ')[1]
                f_out.write(f'{y:0>4},{x:0>4},{ra},{dec}\n')


if __name__ == '__main__':
    main()
