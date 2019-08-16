#!/bin/bash
# Title: hap-processing.sh
# Version: 0.0
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2019-03-25
# Modified in:
# Licence : GPL v3



#======#
# Aims #
#======#

# Process the haplotype sequences to generate fasta file for DnaSp.



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

test_dep clustalo



#===========#
# Variables #
#===========#

myvcf="../../2-VCF processing/Sm.SN.NE.TZ.OM_HR9.pChr6_3M.fb.raw.snp_indel.flt-dp4.cleaned.norm.phased.vcf.gz"
[[ ! -s "$myvcf" ]] && error "Cannot detect the VCF file. Exiting..." 1

output="SmSULT-OR_haplotypes.fas"



#============#
# Processing #
#============#

info "Generating haplotype sequences..."

# Generate haplotypes from Brazilian genotypes
cd "Brazilian haplotypes"
./GT_to_seq.sh "Table S1 - SmSULT-OR genotypes.csv" Smp_089320_cds.fas Sm.BR_hap.fas
cd ..

# Generate haplotypes from the VCF
./vcf-to-hap.sh "$myvcf" "$output"

# Add Brazilian genotypes
cat "Brazilian haplotypes/Sm.BR_hap.fas" >> "$output"

# Add S. rodhaini outgroup
echo ">Sro" >> "$output"
sed -n "2p" "Brazilian haplotypes/SrSULT_cds_rev-comp.fas" | sed "s/$/\n/" >> "$output"

# Align the data for uploading data into DnaSp
info "Aligning haplotypes..."
clustalo -i "$output" -o "${output%.fas}.aln.fas" --output-order=input-order --thread=$(nproc) --force

exit 0
