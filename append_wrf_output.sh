#!/bin/bash
# This script will append nolnox wrfout* to lnox wrfout*
# and get LNO and LNO2

# Set dir of wrfout* which have been extracted
datadir='/nuist/u/home/yinyan/xin/work/BEHR/data/wrf_profiles/history'
lnox_dir=${datadir}'/lnox/'
nolnox_dir=${datadir}'/nolnox/'
savedir='/nuist/u/home/yinyan/xin/work/BEHR/data/wrf_profiles/history/processed'

# Check same files in lnox_dir and nolnox_dir
# Then append variables and calculate

files="$(comm -12 <(ls -F $lnox_dir) <(ls -F $nolnox_dir))"
echo -e "Append these files:\n$files"
for file in $files
    do
        if [[ -f $lnox_dir$file ]]; then
            ncks -A -h $nolnox_dir$file $lnox_dir$file
            ncap2 -v -h -O -s 'lno2=no2_lnox-no2_nolnox;lno=no_lnox-no_nolnox' $lnox_dir$file $lnox_dir$file.tmpnc
            ncks -A -h $lnox_dir$file.tmpnc $lnox_dir$file
            ncatted -O -h -a description,lno2,o,c,'LNO2 mixing ratio' $lnox_dir$file
            ncatted -O -h -a description,lno,o,c,'LNO mixing ratio' $lnox_dir$file
            mv $lnox_dir$file $savedir$file
        fi
    done

rm $lnox_dir*.tmpnc
