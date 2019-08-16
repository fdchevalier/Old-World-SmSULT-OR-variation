#!/bin/bash
# Title: all_exome_pipeline.sh
# Version: 1.3
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2013-12-19
# Modified in: 2017-01-20
# License : GPL v3



#======#
# Aims #
#======#

aim="Start the pipeline for all the libraries and generate stats."



#==========#
# Versions #
#==========#

# v1.3 - 2017-01-20: update due to exome_pipeline.sh update
# v1.2 - 2016-08-13: chromosome length option removed due to update of other script / default location of status folder changed for current folder
# v1.1 - 2016-06-07: script to list software version added / new info messages added / pqsub bug corrected
# v1.0 - 2016-05-07: script renamed (formerly all_pipeline.sh) / options added for passing to exome_pipeline.sh / checking steps added
# v0.3 - 2015-07-26: email address updated / script structure updated / error and warning functions updated / info and test_dep function added / time function added / verbose function added / remove temp file
# v0.2 - 2014-05-23: error and warning functions added / part about algnt_pipeline updated / sex ratio analysis added / checking presence of bed files upgraded
# v0.1 - 2014-01-12: bug correction regarding /tmp/list file name for multisession
# v0.0 - 2013-12-19: creation

version=$(grep -m 1 -i "version" "$0" | cut -d ":" -f 2 | sed "s/ * //g")



#===========#
# Functions #
#===========#

# Usage message
function usage {
    echo -e "
    \e[32m ${0##*/} \e[00m -d|--dir library_directory -r|--ref file -q|--qrecal file -G|--GATK_loc path -P|--picard_loc path -b|--bed file.bed -w|--windows integer -s|--status folder -h|--help

Aim: $aim

Version: $version

Options:

    -d, --dir           path of the directory containing all the library subdirectories [default: data]
    -r, --ref           path to your reference genome. This genome needs to be indexed with bwa, samtools and picard.
    -q, --qrecal        path to the file that contain a list of known variants to perform Q-score recalibration. See GATK documentation for more details. (optional)
    -G, --GATK_loc      path to GATK, symbolic link allowed [default: $HOME/local/bin/GenomeAnalysisTK.jar]
    -P, --picard_loc    path to picard, symbolic link allowed [default: $HOME/local/bin/picard.jar]
    -b, --bed           bed file containing the target region (eg., baits coordinates)
    -w, --windows       integer(s) (bp) corresponding to the window(s) around the target region
                            * 0 needs to be the fisrt integer [default]
                            * list of intergers can be set and must be space separated (eg., 0 10 25 100 150)
    -s, --status        folder where the status files from Grid Engine will be kept [default: $PWD/status]
    -t, --time          time (in day) during the script wait before checking if alignment is done [default: 1]
    -v, --verb          print info message about step ongoing
    -h, --help          this message
    "
}


# Info message
function info {
    if [[ -t 1 ]]
    then
        echo -e "\e[32mInfo:\e[00m $1"
    else
        echo -e "Info: $1"
    fi
}


# Warning message
function warning {
    if [[ -t 1 ]]
    then
        echo -e "\e[33mWarning:\e[00m $1"
    else
        echo -e "Warning: $1"
    fi
}


# Error message
## usage: error "message" exit_code
## exit code optional (no exit allowing downstream steps)
function error {
    if [[ -t 1 ]]
    then
        echo -e "\e[31mError:\e[00m $1"
    else
        echo -e "Error: $1"
    fi

    if [[ -n $2 ]]
    then
        exit $2
    fi
}


# Dependency test
function test_dep {
    which $1 &> /dev/null
    if [[ $? != 0 ]]
    then
        error "Package $1 is needed. Exiting..." 1
    fi
}


# Ended job testing
function job_test {
    sleep $time_sleep

    while [[ $switch == 0 ]]
    do
        if [[ $(qstat -j $1 2>&1 | grep "Following jobs do not exist:") ]]
        then
            switch=1
        else
            sleep $time_sleep
        fi
    done
}



#==============#
# Dependancies #
#==============#

test_dep bwa
test_dep samtools
test_dep java
test_dep R
test_dep bedtools
test_dep awk
test_dep bc

# Grid Engine
test_dep qsub

# Personal scripts
test_dep sftwr-version.sh



#===========#
# Variables #
#===========#

# Options
while [[ $# -gt 0 ]]
do
    case $1 in
        -d|--dir        ) mydir="$2" ; shift 2 ;;
        -r|--ref        ) myref="$2" ; shift 2 ;;
        -q|--recal      ) mysites="$2" ; shift 2 ;;
        -G|--GATK_loc   ) myGATK_loc="$2" ; shift 2 ;;
        -P|--picard_loc ) mypicard_loc="$2" ; shift 2 ;;
        -b|--bed        ) mybed="$2" ; shift 2 ;;
        -w|--win        ) mywindows="$2" ; shift 2
                            while [[ -n "$1" && $(echo "$1"\ | grep -v "^-") ]]
                            do
                                mywindows=(${mywindows[@]} $1)
                                shift
                            done ;;
        -s|--status     ) mystat_folder="$2" ; shift 2 ;;
        -t|--time       ) mytime=$2 ; shift 2 ;;
        -v|--verb       ) verbose=1 ; shift 1 ;;
        -h|--help       ) usage ; exit 0 ;;
        *               ) error "Invalid option: $1\n$(usage)" 1 ;;
    esac
done


#----------------#
# Checking steps #
#----------------#

# Check library directory
if [[ -z "$mydir" ]]
then
    mydir="data"
elif [[ ! -d "$mydir" ]]
then
    error "$mydir is not a directory. Exiting..." 1
fi


# Check reference genome
if [[ -z "$myref" ]]
then
    error "Reference genome is missing."
    flag_error=1
elif [[ ! -s "$myref" ]]
then
    error "${myref##*/} is not a file."
    flag_error=1
fi

# Check mandatory files for the bwa, samtools and GATK
if [[ -n "$myref" && ! (-s "$myref.amb" && -s "$myref.ann" && -s "$myref.bwt" && -s "$myref.pac" && -s "$myref.sa") ]]
then
    error "The reference genome is not indexed for bwa or some files are missing. You must perfom bwa index on the reference genome with options needed (see bwa manual)."
    flag_error=1
fi

if [[ -n "$myref" && ! -s "$myref.fai" ]]
then
    error "The reference genome is not indexed for samtools. You must perfom samtools faidx on the reference genome."
    flag_error=1
fi

if [[ -n "$myref" && ! -s "${myref%.*}.dict" ]]
then
    error "The reference genome is not indexed for GATK. You must use the CreateSequenceDictionary command from picard on the reference genome. See picard documentation for more information."
    flag_error=1
fi


# Q-score recalibration
if [[ -n "$mysites" && ! -s "$mysites" ]]
then
    error "$mysite does not exist or is empty."
    flag_error=1
fi
[[ -n "$mysites" ]] && myq="-q $(printf '%q' "$mysites")"


# Check bed file
if [[ -z "$mybed" ]]
then
    error "Bed file missing."
    flag_error=1
elif [[ ! -s "$mybed" ]]
then
    error "$mybed does not exist or is empty."
    flag_error=1
fi


# Window array
if [[ -z $mywindows ]]
then
    mywindows=0
elif [[ $(echo "$mywindows" | grep -E "[[:alpha:]]|[[:punct:]]") ]]
then
    error "The window option supports only numbers."
    flag_error=1
elif [[ ${mywindows[0]} != 0 ]]
then
    mywindows=(0 ${mywindows[*]})
fi


# Program variables
if [[ -z "$myGATK_loc" ]]
then
    myGATK_loc="$HOME/local/bin/GenomeAnalysisTK.jar"
elif [[ ! -e "$myGATK_loc" ]]
then
    error "$myGATK_loc does not point on a valid file."
    flag_error=1
fi

if [[ -z "$mypicard_loc" ]]
then
    mypicard_loc="$HOME/local/bin/picard.jar"
elif [[ ! -e "$mypicard_loc" ]]
then
    error "$mypicard_loc does not point on a valid file."
    flag_error=1
fi


if [[ -z "$mystat_folder" ]]
then
    mystat_folder="status"
    [[ ! -d "$mystat_folder" ]] && mkdir "$mystat_folder"
elif [[ ! -d "$mystat_folder" ]]
then
    error "$mystat_folder is not a folder."
    flag_error=1
fi


if [[ -z "$mytime" ]]
then
    mytime=1
elif [[ -n "$mytime" && ! $(echo "$mytime" | egrep -x "[[:digit:]]*||[[:digit:]]*.[[:digit:]]*") ]]
then
    error "Time must be numerical. Exiting..."
    flag_error=1
fi


#-------------------#
# Program variables #
#-------------------#

# Grid Engine variables
pqsub="qsub -V -cwd -o $mystat_folder -j y -S /bin/bash"
pqsubR="qsub -V -cwd -o $mystat_folder -j y -r y -S $(which Rscript)"

# Script variables
root_path=$(dirname "$0")
mydate=$(date "+%d%H%M%S%N")
tmp="/tmp/list_$mydate"

# Check for exome pipeline script
if [[ ! -s "$root_path/exome_pipeline.sh" ]]
then
    error "exome_pipeline.sh missing."
    flag_error=1
fi


# Exit for any error encountered 
if [[ $flag_error == 1 ]]
then 
    exit 1
fi



#============#
# Processing #
#============#

# Set bash options to stop script if a command exit with non-zero status
set -e
set -o pipefail


# Clean the dir variable
dir=$(echo "$mydir" | sed "s,/$,,")


# List all the library subdirectories
find $mydir -maxdepth 1 -type d | sed "1d" > "$tmp"


# Start the algnt_pipeline.sh for each subdirectory
[[ -n $verbose ]] && info "Alignment submission..."

while read lib
do
    "$root_path/exome_pipeline.sh" -d "$lib" -r "$myref" $myq -G "$myGATK_loc" -P "$mypicard_loc" -b "$mybed" -w "$mywindows" -s "$mystat_folder" &
done < "$tmp"

[[ -n $verbose ]] && info "Alignment submitted."


# Store software versions used in this pipeline 
[[ -n $verbose ]] && info "Versions of software used in this pipeline stored in sftwr-version.report."
sftwr-version.sh -i bwa samtools "$myGATK_loc" "$mypicard_loc MarkDuplicates" bedtools


# Wait the specified time
[[ -n $verbose ]] && info "Next step will be launch at the soonest in $mytime day(s)."
sleep ${mytime}d


# Check every hour if all the files needed was created to process statistics
while [[ $(find "$mydir" -name *baits-0.bed -type f | wc -l) -lt $(wc -l < "$tmp") ]]
do
    sleep 1h
done

[[ -n $verbose ]] && info "Statistical analysis ongoing..."

if [[ $(find "$mydir" -name *baits-0.bed -type f | wc -l) == $(wc -l < "$tmp") ]]
then
    sizes=$(find "$mydir" -name *baits-0.bed -type f -exec ls -l {} \; | cut -d " " -f 5) # Store size of the files

    if [[ ! $(echo -e "$sizes" | grep -x "0") ]]
    then
        $pqsubR $(which Coverage_along_baits.R)      # The directory must be data
        $pqsubR $(which schisto_sex_RD.R)
        $pqsub $(which lib_stat.sh) -d "$mydir"
    else
        error "Some *baits-0.bed files are empty" 1
    fi
fi


# Remove temporary file
rm $tmp


exit 0
