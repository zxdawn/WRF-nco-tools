#!/bin/bash
#
# This script will run the necessary operations to process WRF-Chem output
# for use as BEHR a priori.  There are 3 modes of processing, which should
# be set by modifying the "mode" variable below.
#
#	"hourly" will just extract all profiles within the OMI overpass times 
# of the continental US/CN, with the intention that the user will select the 
# appropriate profile later. 
#	"daily" will average those profiles for each day based on the 
# longitude and UTC time - profiles will be given more weight the closer 
# they are to 1400 local standard time. Those more than an hour off will 
# have a weight of 0.  
#	"monthly" will average with weights as in "daily" but over a month
# rather than a single day.
#
# "d02"(d03) will eclect wrfout_d02*(wrfout_d03*)
#
# This script does not directly perform any of those calculations, instead
# it collects the names of the WRF output files that belong
# to each group to pass to read_wrf_output.sh to do the actual calculation.
# This keeps the structure of this program the same as on the HPC cluster,
# where it needs to be done this way to launch multiple instances of 
# read_wrf_output.sh in parallel.

# Josh Laughner <joshlaugh5@gmail.com> 2 Jul 2015
# modified Xin Zhang <xinzhang1215@gmail.com> May 2018

# Parse command arguments looking for two things: the averaging mode and which set of 
# output quantities to copy/calculate. Credit to 
# http://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
# for the code outline.

while [[ $# > 0 ]]
do
keyin="$1"
# Ensure input is lower case
key=$(echo $keyin | awk '{print tolower($0)}')
    case $key in
        'monthly'|'daily'|'hourly')
        mode=$key
        shift # shift the input arguments left by one
        ;;
        'behr'|'emis'|'avg')
        varsout=$key
        shift
        ;;
        'd01'|'d02'|'d03')
        domain=$key
        shift
        ;;
        'lnox'|'nolnox')
        kind=$key
        shift
        ;;
        'keep'|'del')
        choice=$key
        shift
        ;;
        *) # catch unrecognized arguments
        echo "The argument \"$key\" is not recognized"
        exit 1
        ;;
    esac
done

# Set the defaults - averaging mode will default to "hourly"
# the outputs to "behr" and domain to "d01"
if [[ $mode == '' ]]
then
    mode='hourly'
fi

if [[ $kind == '' ]] || [[ $kind != 'lnox' && $kind != 'nolnox' ]]
then
    echo "Please input WRFKIND: lnox or nolnox"
    exit 0
else
    echo "kind set to $kind"
fi

if [[ $varsout == '' ]]
then
    varsout='behr'
 fi

if [[ $domain == '' ]]
then
    domain='d01'
fi

if [[ $choice == '' ]]
then
    choice='keep'
fi

# export the mode and the kind of wrfout*
# so that the child scripts can access it
export WRFPROCMODE=$mode
export WRFKIND=$kind

# Where the actual scripts (read_wrf_output.sh) are kept.
scriptdir='/nuist/u/home/yinyan/xin/work/BEHR/WRF-nco-tools/'
export WRFSCRIPT_DIR=$scriptdir

# Where new wrfout* files are saved
savedir='/nuist/u/home/yinyan/xin/work/BEHR/data/wrf_profiles/history/'

# Check the mode selection

if [[ $mode != 'daily' && $mode != 'monthly' && $mode != 'hourly' ]]
then
    echo "Input must be 'daily' or 'monthly'"
    exit 1
else
    echo "mode set to $mode"
fi

# Check the domain selection
if [ "$(echo wrfout_$domain*)" == "wrfout_$domain*" ]; then
    echo "wrfout_$domain* not found. Please input correct domain like d02"
fi

# Find all unique dates - we'll need this to iterate over each day
# If we're doing monthly averages, then we just need to get the year and month
dates=''
olddate=''

for file in ./wrfout_$domain*
do
    # Handle wrfout and wrfout_subset files
    dtmp=$(awk -v a="$file" -v b="$domain" 'BEGIN{print index(a,b)}')
    dstart=$((dtmp+3))
    if [[ $mode == 'monthly' ]]
    then
        newdate=${file:$dstart:7}
    else
        newdate=${file:$dstart:10}
    fi

    if [[ $olddate != $newdate ]]
    then
        dates=$(echo $dates $newdate)
    fi
    olddate=$newdate
done


for day in $dates
do
    echo ""
    echo "Files on $day"
    echo ""
    # WRF file names include output time in UTC. We'll look for the output
    # in the range of UTC times when OMI will be passing over domain
    # for this day
    # North America: {18,19,20,21,22}
    # Southeast China: {06,07}
    # If there are no files for this day or month, then it will try to iterate
    # over the wildcard patterns themselves. Since those contain *, we
    # can avoid doing anything in that case by requiring that the file
    # name does not include a *
    filepattern=$(echo wrfout*_"$domain"_${day}_{06,07}*)
    if [[ $filepattern != *'*'* ]]
    then
        echo -e "Select these wrfout* files:\n $filepattern"
        echo "$filepattern" > read_wrf.conf
        # Choose which command to execute based on the command arguments
        if [[ $varsout == 'behr' ]]
        then
            echo "Calling read_wrf_output"
            $scriptdir/read_wrf_output.sh read_wrf.conf

            if [[ $kind == "lnox" ]]
            then
                savedir=${savedir}"lnox/"
            else
                savedir=${savedir}"nolnox/"
            fi
            mkdir -p $savedir

            echo "Save to $savedir"
            mv *.tmpnc $savedir
            for file in $savedir*.tmpnc; do mv "$file" "${file%%.tmpnc}"; done

        elif [[ $varsout == 'emis' ]]
        then
            $scriptdir/read_wrf_emis.sh read_wrf.conf
        elif [[ $varsout == 'avg' ]]
        then
            $scriptdir/avg_wrf_output.sh read_wrf.conf
        else
            echo "Error at $LINENO in slurmrun_wrf_output.sh: \"$varsout\" is not a recognized operation"
            exit 1
        fi
    fi
done

# Whether delete original wrfout* files
if [[ $choice == "del" ]]; then
    echo "        Remove all wrfout* files in current directory"
    rm wrfout_*
fi