#!/bin/bash
# Script to be executed on the HPC cluster through the
# slurmrun_wrf_output.sh shell script. Will iterate
# through the list of files passed to it and pull out/
# calculate the necessary quantities to use WRF output
# as BEHR a priori
# Josh Laughner <joshlaugh5@gmail.com> 1 Jul 2015
# modified Xin Zhang <xinzhang1215@gmail.com> May 2018

# Retrieve the operational mode and kind from the environmental variable
mode=$WRFPROCMODE
kind=$WRFKIND

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
    my_variables="T,HGT,^Q.?"

    if [[ $kind == "lnox" ]]
    then
        variables=$BEHR_variables$my_variables
        echo "Copy $kind kind"
        echo "        Copy variables: $variables"
        ncks -A -h -v $variables $file $file.tmpnc
        ncrename -h -v no2,no2_lnox -v no,no_lnox $file.tmpnc
        ncks -A -h -v no2 $file $file.tmpnc
    else
        variables="no2,no"
        echo "Copy $kind kind"
        echo "        Copy variables: $variables"
        ncks -A -h -v $variables $file $file.tmpnc
        ncrename -h -v no2,no2_nolnox -v no,no_nolnox $file.tmpnc
    fi
    
    echo "        Copying attributes..."
    ncks -A -h -x $file $file.tmpnc
done
