#!/bin/bash
# Title: vcf-processing.sh
# Version: 0.0
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2019-02-28
# Modified in:
# Licence : GPL v3



#======#
# Aims #
#======#

# Process the VCF from the raw calling to the phased data



#==========#
# Versions #
#==========#

# v0.0 - 2019-02-28: creation



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


# Clean up function for trap command
## Usage: clean_up file1 file2 ...
function clean_up {
    rm -rf $@
    exit 1
}



#==============#
# Dependencies #
#==============#

test_dep wget
test_dep vcftools 
test_dep vcf-subset 
test_dep vt 
test_dep java



#===========#
# Variables #
#===========#

myfilename="Sm.SN.NE.TZ.OM_HR9.pChr6_3M.fb.raw.snp_indel"

mybeagle="$HOME/local/pckg/Beagle/beagle.27Jan18.7e1.jar"
[[ ! -s "$mybeagle" ]] && error "Cannot detect beagle. Exiting..." 1

mygenome="$HOME/data/sm_genome/sma_v5.0.chr.fa"
[[ ! -s "$mygenome" ]] && error "Cannot detect schistosome genome. Exiting..." 1


#============#
# Processing #
#============#

## Command used for calling variants
#freebayes -f ~/data/sm_genome/sma_v5.0.chr.fa \
#    -r Schisto_mansoni.Chr_6:1-3000000  \
#    -b $(find data -name *.bam | tr "\n" " ") \
#    -=
#    --population pop \
#    -q 20            \
#    -m 30            \
#    -@ sm_dbSNP_Smp_089320.vcf.gz > "${myfilename}".vcf

# Download the raw VCF from Zenodo
wget -O "${myfilename}.vcf.gz" "https://zenodo.org/record/2850876/files/Sm.SN.NE.TZ.OM_HR9.pChr6_3M.fb.raw.snp_indel.vcf.gz?download=1"

# Check file
[[ $(md5sum "${myfilename}.vcf.gz" | cut -d " " -f 1) != "57ab9157c1cce8091894a3142f69db56" ]] && error "The downloaded file is incomplete. Exiting..." 1

# Filter on read depth
vcftools --gzvcf "${myfilename}.vcf.gz" --minDP 4 --recode --recode-INFO-all --out "${myfilename}.flt-dp4"
mv "${myfilename}.flt-dp4.recode.vcf" "${myfilename}.flt-dp4.vcf"

# Clean the VCF of invariable ref variants
vcf-subset -a -e "${myfilename}.flt-dp4.vcf" > "${myfilename}.flt-dp4.cleaned.vcf"

# Normalize VCF
vt normalize -o "${myfilename}.flt-dp4.cleaned.norm.vcf" \
    -r "$mygenome" \
    "${myfilename}.flt-dp4.cleaned.vcf"

# Generate report
./func-table-mut.sh -v "${myfilename}.flt-dp4.cleaned.norm.vcf" -r "$mygenome" -g Smp_089320.gff -p pop -o Smp_089320.flt.norm.tsv

# Phase data
java -Xmx16g -jar "$mybeagle" gtgl="${myfilename}.flt-dp4.cleaned.norm.vcf"         out="${myfilename}.flt-dp4.cleaned.norm.exp_gt"
java -Xmx16g -jar "$mybeagle" gt="${myfilename}.flt-dp4.cleaned.norm.exp_gt.vcf.gz" out="${myfilename}.flt-dp4.cleaned.norm.phased"


exit 0
