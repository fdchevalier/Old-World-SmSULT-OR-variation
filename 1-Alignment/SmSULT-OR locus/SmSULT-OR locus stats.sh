#!/bin/bash
# Title: SmSULT-OR locus stats.sh
# Version: 0.0
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2019-08-13
# Modified in:
# Licence : GPL v3



#======#
# Aims #
#======#

# Generate read depth and coverage statistics for the SmSULT-OR exons and CDS.



#==========#
# Versions #
#==========#

# v0.0 - 2019-08-13: creation



#===========#
# Functions #
#===========#

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



#===========#
# Variables #
#===========#

mydirname="data"
[[ ! -d "$mydirname" ]] && error "Cannot detect the data directory. Exiting..." 1

mygff="Smp_089320.gff"
[[ ! -s "$mygff" ]] && error "Cannot detect the Smp_089320.gff. Exiting..." 1



#============#
# Processing #
#============#

# Generate the needed GFF
egrep "CDS|UTR" "$mygff" > Smp_089320_exons.gff

# Check if files exist
[[ $(find . -name *coverage.bed) ]] || error "No coverage bed files detectable. Exiting..." 1

# Get exons coverage (start all bedtools jobs)
for i in "$mydirname"/*
do
    bedtools intersect -a "$i"/*coverage.bed -b Smp_089320_exons.gff > $i/Smp_089320_exons.bed &
    mypid="$mypid|^$!$"
done

# Check if all jobs are done
running=true
mypid=$(echo "$mypid" | sed "s/|//")
while $running
do
    if [[ $(egrep "$mypid" <(ps | awk '{print $1}' | tail -n +2)) ]]
    then
        sleep 30s
    else
        running=false
        echo "All jobs done."
    fi
done

# Check if files exist
[[ $(find . -name Smp_089320_exons.bed) ]] || error "No exon bed files detectable. Exiting..." 1

# Run the R script
./exon_stat.R

# Clean up
rm Smp_089320_exons.gff

exit 0
