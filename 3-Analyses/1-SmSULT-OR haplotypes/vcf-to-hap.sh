#!/bin/bash
# Title: vcf-to-hap.sh
# Version: 0.0
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2019-03-25
# Modified in:
# Licence : GPL v3



#======#
# Aims #
#======#

# Generate haplotype sequences in fasta format from a VCF file.



#==========#
# Versions #
#==========#

# v0.0 - 2019-03-25: creation



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

test_dep bgzip
test_dep tabix
test_dep samtools
test_dep bcftools
test_dep revseq



#===========#
# Variables #
#===========#

myvcf="$1"
myoutput="$2"

myref="$HOME/data/sm_genome/sma_v5.0.chr.fa"
[[ ! -s "$myref" ]] && error "Cannot detect $myref. Exiting..." 1
mygff="$HOME/data/sm_Gene_table/v5.07.08.12.chado.raw.gff"
[[ ! -s "$mygff" ]] && error "Cannot detect $mygff. Exiting..." 1

tmp_d=$(mktemp -d)



#============#
# Processing #
#============#

# VCF preparation
if [[ ! $(file -ibL "$myvcf" | grep gzip) ]]
then
    info "VCF preparation..."
    bgzip -c "$myvcf" > $tmp_d/$(basename "$myvcf").gz
    myvcf="$tmp_d/$(basename "$myvcf").gz"

    tabix "$myvcf"
fi

myint=$(grep "Smp_089320" "$mygff" | grep CDS | sort -k 4 | awk '{print $1":"$4"-"$5}' | tr "\n" " ")
mystrand=$(grep "Smp_089320" "$mygff" | grep CDS | cut -f 7 | uniq)

myspl=$(zgrep -m 1 "#C" "$myvcf" | cut -f 10- | tr "\t" "\n" | sort)

# Generate haplotype sequences fro each sample
for s in $myspl
do
    echo "$s"

    for h in 1 2
    do
        myseq=""
        
        for i in $myint
        do
            myseq="$myseq$(samtools faidx "$myref" "$i" | bcftools consensus -s $s -H $h "$myvcf" | sed "1d")"
        done

        if [[ "$mystrand" == "-" ]]
        then
            revseq <(echo $myseq | tr -d "\n") $tmp_d/a &> /dev/null
            sed -i "1d" "$tmp_d/a"
        else
            fold -w 60 <(echo "$myseq" | tr -d "\n") > $tmp_d/a
        fi

        echo ">${s}_${h}" >> "$myoutput"
        cat $tmp_d/a >> "$myoutput"
    done

done

rm -R $tmp_d
