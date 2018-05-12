#!/bin/bash
# Script to be executed on the HPC cluster through the
# slurmrun_wrf_output.sh shell script. Will iterate
# through the list of files passed to it and pull out/
# calculate the necessary quantities to use WRF output
# as BEHR a priori
# Josh Laughner <joshlaugh5@gmail.com> 1 Jul 2015
# modified Xin Zhang <xinzhang1215@gmail.com> May 2018

# Retrieve the operational mode from the environmental variable
mode=$WRFPROCMODE

# Where the various NCO scripts are located
scriptdir=$WRFSCRIPT_DIR
if [[ $scriptdir == '' ]]
then
    echo "Error at $LINENO in read_wrf_output.sh: WRFSCRIPT_DIR not set"
    exit 1
fi

# Keep a list of the files to concatenate
catfiles=''

wrffiles=$(cat $1)

for file in $wrffiles
do
    # Copy the variables needed to keep track of position & time,
    # plus the species mixing ratios that we're interested in.
    # U and V are temporary for studying the effect of wind on a priori
    # COSALPHA AND SINALPHA are needed to convert the winds from grid
    # relative to earth relative. See http://forum.wrfforum.com/viewtopic.php?f=8&t=3225
    # Xin
    # Add all variables which BEHR a priori and wrfpython need.
    # Since rProfile.m calculate 'pres,z,zlev' and convert no2 units,
    # We just need to extract necessary variables directly.
    echo "        Copying variables..."
    BEHR_variables="Times,XLAT,XLONG,XLONG_U,XLAT_U,XLONG_V,XLAT_V,U,V,COSALPHA,SINALPHA,P,PB,PHB,PH,no2,no,"
    my_variables="T,HGT,QVAPOR,QICE,QCLOUD,QGRAUP,QSNOW,QRAIN"
    variables=$BEHR_variables$my_variables
    ncks -A -v $variables $file $file.tmpnc
    
    echo "        Copying attributes..."
    ncks -A -x $file $file.tmpnc
done