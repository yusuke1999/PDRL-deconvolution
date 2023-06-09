#!/bin/bash
# Merge the data.
merge_obs 4636/repro/acisf04636_repro_evt2.fits,\
4637/repro/acisf04637_repro_evt2.fits,\
4639/repro/acisf04639_repro_evt2.fits,\
5319/repro/acisf05319_repro_evt2.fits \
merged_4636_4637_4639_5319/ bands=broad binsize=1

# Create a txt file converted from pixel coordinates to radec coordinates.
python ../pyscript/pixel2skycoord.py merged_4636_4637_4639_5319/broad_thresh.img pixel2sky_35bin.txt 35

# Simulation of psf by marx using Obs. ID 4636 as a representative.
bash ../shscript/simulate_psfs.sh pixel2sky_35bin.txt 4636/repro/acisf04636_repro_evt2.fits 2.3 psfs_35bin

# Adjust the psf so that the center of the psf is at the center of the image.
# Save an array with psf for each location.
python ../pyscript/repro_psfs.py merged_4636_4637_4639_5319/broad_thresh.img psfs_35bin repro_psfs_35bin.npz 35

# Run PDRL method.
python ../pyscript/PDRL.py merged_4636_4637_4639_5319/broad_thresh.img merged_4636_4637_4639_5319/broad_thresh.expmap \
       repro_psfs_35bin.npz PDRL_results --lambda_tv 0 --boundary_px 0

# Run the error of PDRL by the law of error propagation. 
python ../pyscript/PDRL_err.py merged_4636_4637_4639_5319/broad_thresh.img merged_4636_4637_4639_5319/broad_thresh.expmap \
       repro_psfs_35bin.npz PDRL_results/iter_0199.fits iter_0200_err.fits
