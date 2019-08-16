#!/bin/bash
# Title: exome_pipeline.sh
# Version: 1.2
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2013-01-12
# Modified in: 2017-01-20
# License : GPL v3



#======#
# Aims #
#======#

aim="Perform alignment, read depth statistics and capture efficiency statistics."



#==========#
# Versions #
#==========#

# v1.2 - 2017-01-20: update due to alignment.sh update / bug correction
# v1.1 - 2016-08-16: chromosome length option removed due to update of other script / default location of status folder changed for current folder
# v1.0 - 2016-05-09: script renamed (formerly algnt_pipeline) / email address added / functions improved
# v0.4 - 2015-01-30: script formatting updated / script name updated / functions updated / bug correction of the $myscript (absolute path was used to look at status file) / check if status files exist
# v0.3 - 2014-05-23: bug correction in the job_test function (error 1 emitted for job ended or job not reachable which occurs during high server load leading to start the following step before the end of the running step) / usage, error and warning functions added / variables updated
# v0.2 - 2013-12-19: usage message updated / variables rearrangement
# v0.1 - 2013-09-03: some improvements
# v0.0 - 2013-01-12: creation

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

    -d, --dir           directory containing the library
    -r, --ref           path to your reference genome. This genome needs to be indexed with bwa, samtools and picard.
    -q, --qrecal        path to the file that contain a list of known variants to perform Q-score recalibration. See GATK documentation for more details. (optional)
    -G, --GATK_loc      path to GATK, symbolic link allowed [default: $HOME/local/bin/GenomeAnalysisTK.jar]
    -P, --picard_loc    path to picard, symbolic link allowed [default: $HOME/local/bin/picard.jar]
    -b, --bed           bed file containing the target region (eg., baits coordinates)
    -w, --windows       integer(s) (bp) corresponding to the window(s) around the target region
                            * 0 needs to be the fisrt integer [default]
                            * list of intergers can be set and must be space separated (eg., 0 10 25 100 150)
    -s, --status        folder where the status files from Grid Engine will be kept [default: $PWD/status]
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
test_dep qrsh



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
        -w|--windows    ) mywindows="$2" ; shift 2
                            while [[ -n "$1" && $(echo "$1"\ | grep -v "^-") ]]
                            do
                                mywindows=(${mywindows[*]} $1)
                                shift
                            done ;;
        -s|--status     ) mystat_folder="$2" ; shift 2 ;;
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
    error "The directory containing libraries is required. Exiting..." 1
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
elif [[ -n "$mysites" ]] 
then
    myq="-q $(printf '%q' "$mysites")"
    myext="_recal"
fi


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
elif [[ $(echo ${mywindows[@]} | grep -E "[[:alpha:]]|[[:punct:]]") ]]
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


# Check fastq file presence
lib_name=$(basename "$mydir")

pe1=$(ls "$mydir/${lib_name}"*_R1.fastq.gz 2> /dev/null)
pe2=$(ls "$mydir/${lib_name}"*_R2.fastq.gz 2> /dev/null)

if [[ ! -f "$pe1" ]]
then
    error "Unable to locate $pe1. Does the file exist?"
    flag_error=1
elif [[ ! -f "$pe2" ]]
then
    error "Unable to locate $pe2. Does the file exist?"
    flag_error=1
fi


# Exit for any error encountered 
if [[ $flag_error == 1 ]]
then 
    exit 1
fi


#-------------------#
# Program variables #
#-------------------#

root_path=$(dirname $0)
myscripts_dir="$root_path/exome_scripts"

pqsub="qsub -V -cwd -o $mystat_folder -j y -S /bin/bash"

# Variables initiation for job_test function
time_sleep=300
i=0
switch=0



#============#
# Processing #
#============#

# Check for scripts directory just in case something happened
if [[ ! -d "$myscripts_dir" ]]
then
    error "$myscripts_dir does not exist. Exiting..." 1
fi


# Alignment and bam file processing
myscript="alignment.sh"
myscript_opt="-d $mydir -r $myref $myq -G $myGATK_loc -P $mypicard_loc"
myscript_path="$myscripts_dir/$myscript"
job_n=$($pqsub "$myscript_path" $myscript_opt | cut -d " " -f 3)

job_test ${job_n}

job_log="$mystat_folder/${myscript}.o${job_n}"
if [[ ! -e "$job_log" ]]
then
    error "Status file does not exist. Errors cannot be checked. Exiting..." 1
elif [[ $(grep -i error "$job_log") ]]
then
    if [[ $(grep "Ignoring SAM validation error:" $job_log &> /dev/null; echo $?) == 1 ]]
    then
        echo -e "File: ${dir}\nJob Number: ${job_n}\nTime: $(date)" >> ERROR_${myscript}_${job_n}.log
        exit 1
    fi
else
    switch=0
    myfile="$mydir/${lib_name}_sorted_realigned_MD${myext}.bam"
fi


# Read depth of the capture region and around
myscript="target_read_depth.sh"
myscript_opt="-i $myfile -b $mybed"
myscript_path="$myscripts_dir/$myscript"
while [[ $i -lt ${#mywindows[@]} ]]
do
    switch=0
    
    job_n=$($pqsub "$myscript_path" $myscript_opt -w ${mywindows[$i]} | cut -d " " -f 3)

    sleep $time_sleep

    job_test $job_n
    
    job_log="$mystat_folder/${myscript}.o${job_n}"
    if [[ ! -e "$job_log" ]]
    then
       error "Status file does not exist. Errors cannot be checked. Exiting..." 1
    elif [[ $(grep -i "error" "$job_log") || $(grep -i "failed" "$job_log") ]]
    then
        echo -e "File: ${dir}\nWindows: ${i}\nJob Number: ${job_n}\nTime: $(date)" >> ERROR_${myscript}_${job_n}.log
        exit 1
    fi

    ((i++))
done

myfile="$mydir/${lib_name}_sorted_realigned_MD${myext}-coverage_baits-0.bed"

if [[ ! -s "$myfile" ]]
then
    error "The file $myfile is empty. Exiting..." 1
fi


# Capture efficiency
myscript="capture_efficiency.sh"
myscript_opt="-i $myfile -b $mybed"
myscript_path="$myscripts_dir/$myscript"
job_n=$($pqsub "$myscript_path" $myscript_opt | cut -d " " -f 3)

job_test $job_n

job_log="$mystat_folder/${myscript}.o${job_n}"
if [[ ! -e "$job_log" ]]
then
    error "Status file does not exist. Errors cannot be checked. Exiting..." 1
elif [[ $(grep -i "error" "$job_log") ]]
then
    echo -e "File: ${dir}\nJob Number: ${job_n}\nTime: $(date)" >> ERROR_${myscript}_${job_n}.log
    exit 1
#else   # Need if subsequent steps
#    switch=0
#    myfile="$mydir/${lib_name}_sorted_realigned_MD${myext}.bam"
fi


exit 0
