#!/bin/bash
# Title: cqpture_efficiency.sh
# Version: 1.0
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2012-09-25
# Modified in: 2016-05-09
# License : GPL v3



#======#
# Aims #
#======#

aim="Produce statistics regarding the capture: number of baits involved, mean and median coverage." 



#==========#
# Versions #
#==========#

# v1.0 - 2016-05-09: email address updated / functions added / usage message and options added
# v0.3 - 2013-01-18: correction of the baits_merged_bed variable to baits_bed
# v0.2 - 2013-01-10: modification of output variable (naming automation)
# v0.1 - XXXX-XX-XX: addition of baits_merged_bed variable
# v0.0 - 2012-09-25: creation

version=$(grep -m 1 -i "version" "$0" | cut -d ":" -f 2 | sed "s/ * //g")



#===========#
# Functions #
#===========#

# Usage message
function usage {
    echo -e "
    \e[32m ${0##*/} \e[00m -i|--input file.bed -b|--bed file.bed -h|--help

Aim: $aim

Version: $version

Options:
    -i, --input     input bed file
    -b, --bed       bed file containing the target region (eg., baits coordinates)
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
test_dep awk
test_dep bc



#===========#
# Variables #
#===========#

# Options
while [[ $# -gt 0 ]]
do
    case $1 in
        -i|--input  ) myfile="$2" ; shift 2 ;;
        -b|--bed    ) mybed="$2" ; shift 2 ;;
        -h|--help   ) usage ; exit 0 ;;
        *           ) error "Invalid option: $1\n$(usage)" 1 ;;
    esac
done


# Check mandatory variables
if [[ -z "$myfile" ]]
then
    error "Input bed file missing. Exiting...\n$(usage)" 1
fi

if [[ -z "$mybed" ]]
then
    error "Target bed file missing. Exiting...\n$(usage)" 1
fi


# Output filename
nb_baits="${myfile%.*}-capture_efficiency_nb.bed"
cov_baits="${myfile%.*}-capture_efficiency_cov.bed"
report_file="${myfile%.*}-capture_efficiency.log"



#============#
# Processing #
#============#

# Set bash options to stop script if a command exit with non-zero status
set -e
set -o pipefail


#------------------------------------------------#
# Proportion of baits really used in the capture #
#------------------------------------------------#
info "Baits proportion running..."

bedtools intersect -a "$mybed" -b "$myfile" -u > "$nb_baits"
my_pc_baits=$(echo "scale=2 ; ($(wc -l < "$nb_baits")*100) / $(wc -l < "$mybed")" | bc -l)


#--------------------------------------#
# Average and median coverage by baits #
#--------------------------------------#
info "Coverage running..."

bedtools coverage -a "$mybed" -b "$myfile" > "$cov_baits"

# Average
myavg=$(sed "/^#/d" "$cov_baits" | cut -f 7 | awk '{ SUM+= $1 } END { print SUM/NR }')

# Median
# source: http://stackoverflow.com/a/6166483
mymed=$(sed "/^#/d" "$cov_baits" | cut -f 7 | sort -n | awk '{ a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }')


#--------#
# Report #
#--------#
info "Generating report..."

echo -e "#Pc_baits\tAvg_cov\tMed_cov" > "$report_file"
echo -e "$my_pc_baits\t$myavg\t$mymed" >> "$report_file"


#-------------#
# Compressing #
#-------------#
info "Compressing files..."

gzip "$nb_baits"
gzip "$cov_baits"


exit 0
