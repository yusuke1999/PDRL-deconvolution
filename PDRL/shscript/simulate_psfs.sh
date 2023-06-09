function usage {
    cat <<EOF
usage: $(basename ${0}) [-h] infile_txt infile_fits outdir
Get each position PSF from ra, dec .txt file using CIAO \`simulate_psf\`.
positional arguments:
  infile_txt            input .txt file in image_y,image_x,ra,dec
  infile_fits           input evt fits file (Note: not counts map file)
  monoenergy            set psf monoenergy
  outdir                output dir for each position psf
optional arguments:
  -h, --help            show this help message and exit
EOF
}

if [ "${1}" = "-h" ] ||[ "${1}" = "--help" ]; then
  usage
  exit 0
fi

# Parameters
infile_txt=${1}
infile_fits=${2}
monoenergy=${3}
outdir=${4%/}

mkdir -p ${outdir}
while read line
do
  # Assign the strings in the sequence
  # image_y,image_x,ra,dec to the variables respectively.
  image_y="$(echo ${line} | cut -d ',' -f 1)"
  image_x="$(echo ${line} | cut -d ',' -f 2)"
  ra="$(echo ${line} | cut -d ',' -f 3)"
  dec="$(echo ${line} | cut -d ',' -f 4)"

  echo ${image_y},${image_x},${ra},${dec}
  # If you want better statistical psf,
  # we recommend increasing the `flux` parameter.
  # Please refer to https://cxc.cfa.harvard.edu/ciao/ahelp/simulate_psf.html
  # for more information.
  simulate_psf infile=${infile_fits} ra=${ra} dec=${dec} \
  monoenergy=${monoenergy} flux=1e-2 spectrum=none binsize=1 minsize=200 \
  outroot="${outdir}/${image_y}_${image_x}"
done < ${infile_txt}
