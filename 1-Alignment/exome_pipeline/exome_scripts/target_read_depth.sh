#!/bin/bash
# Title: target_read_depth.sh
# Version: 1.1
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2012-09-25
# Modified in: 2016-08-16
# License : GPL v3



#======#
# Aims #
#======#

aim="Read depth on targeted region." 



#==========#
# Versions #
#==========#

# v1.1 - 2016-08-16: remove the chromosome lenght option because useless
# v1.0 - 2016-05-09: name updated (formerly coverageBed.sh) / email address updated / structure updated / functions added / usage message and options added / bug regarding removing file corrected
# v0.2 - 2013-01-04: R analysis integration
# v0.0 - 2012-09-13: creation

version=$(grep -m 1 -i "version" "$0" | cut -d ":" -f 2 | sed "s/ * //g")



#===========#
# Functions #
#===========#

# Usage message
function usage {
    echo -e "
    \e[32m ${0##*/} \e[00m -i|--input file.bam -b|--bed file.bed -w|--window integer -h|--help

Aim: $aim

Version: $version

Options:
    -i, --input     input bam file
    -b, --bed       bed file containing the target region (eg., baits coordinates)
    -w, --window    integer (bp) corresponding to the window around the target region [default: 0]
    -h, --help      this message
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



#==============#
# Dependencies #
#==============#

test_dep bedtools 
test_dep R



#===========#
# Variables #
#===========#

# Options
while [[ $# -gt 0 ]]
do
    case $1 in
        -i|--input    ) mybam="$2" ; shift 2 ;;
        -b|--bed      ) mybed="$2" ; shift 2 ;;
        -w|--window   ) mywindow="$2" ; shift 2 ;;
        -h|--help     ) usage ; exit 0 ;;
        *             ) error "Invalid option: $1\n$(usage)" 1 ;;
    esac
done


# Check mandatory variables
if [[ -z "$mybam" ]]
then
    error "Input file missing. Exiting...\n$(usage)" 1
fi

if [[ -z "$mybed" ]]
then
    error "Bed file missing. Exiting...\n$(usage)" 1
elif [[ ! -s "$mybed" ]]
then
    error "$mybed does not exist or is empty." 1
fi

# Load default values
if [[ -z "$mywindow" ]]
then
    mywindow=0
fi

# Output filenames
output_1="${mybam%.*}-coverage.txt"
output_2="${mybam%.*}-coverage_baits-$mywindow.bed"



#============#
# Processing #
#============#

# Set bash options to stop script if a command exit with non-zero status
set -e
set -o pipefail


#--------------------------#
# Generating coverage file #
#--------------------------#

info "Generating coverage file..."

if [ ! -s "${output_1%.*}.bed" ]
then
    # Compute the coverage for each base
    bedtools genomecov -ibam "$mybam" -dz > "$output_1"

    # Convert the txt file in bed file
    echo -e "#LOCATION\tSTART\tEND\tCOVERAGE" > "${output_1%.*}.bed"
    sed "/^#/d" "$output_1" | awk '{print $1 "\t" $2 "\t" $2+1 "\t" $3}' >> "${output_1%.*}.bed"

    # Remove temporary files
    rm "$output_1"
fi


#-----------------#
# Window analysis #
#-----------------#

info "Window analysis..."

# Intersect the genome coverage with the bait coordinates + window
bedtools window -a "${output_1%.*}.bed" -b "$mybed" -w "$mywindow" -u > "$output_2"


#------------------#
# Statistic report #
#------------------#

info "R analysis..."

# Title
echo -e "$output_2 $mywindow: $(wc -l "$output_2")" >> "${mybam%.*}.log"

# R code computing mean and median
col=$(cut -f 4 "$output_2")
R --slave << EOF >> "${mybam%.*}.log"
mydata <- scan()
$col

summary(mydata)
EOF

echo -e "\n"


exit 0
