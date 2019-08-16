#!/bin/bash
# Title: alignment.sh
# Version: 1.1
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2012-08-23
# Modified in: 2017-01-20
# License : GPL v3



#======#
# Aims #
#======#

aim="Alignment of pair-end reads, re-alignement around indels, duplicates marking, Q-score recalibration and alignment statistics."



#==========#
# Versions #
#==========#

# v1.1 - 2017-01-20: bootstrapping base recalibration added / multithread option added
# v1.0 - 2016-05-08: more general options added (genome, GATK, picard, etc.) / checking steps added / functions added / email address updated / structure updated
# v0.3 - 2015-07-26: command lines updates for new version of GATK and Picard
# v0.2 - 2013-09-03: switching from BWA to BWA-MEM / removing temporary files
# v0.1 - 2013-05-22: improvements
# v0.0 - 2012-08-23: creation

version=$(grep -m 1 -i "version" "$0" | cut -d ":" -f 2 | sed "s/ * //g")



#===========#
# Functions #
#===========#

# Usage message
function usage {
    echo -e "
    \e[32m ${0##*/} \e[00m -d|--dir library_directory -r|--ref file -q|--qrecal file -G|--GATK_loc path -P|--picard_loc path -h|--help

Aim: $aim

Version: $version

Options:
    -d, --dir           directory containing the library
    -r, --ref           path to your reference genome. This genome needs to be indexed with bwa, samtools and picard.
    -q, --qrecal        path to the file that contain a list of known variants to perform Q-score recalibration. See GATK documentation for more details. (optional)
    -G, --GATK_loc      path to GATK, symbolic link allowed [default: $HOME/local/bin/GenomeAnalysisTK.jar]
    -P, --picard_loc    path to picard, symbolic link allowed [default: $HOME/local/bin/picard.jar]
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



#==============#
# Dependencies #
#==============#

test_dep bwa
test_dep samtools
# GATK      # tested later in the script
# picard    # tested later in the script



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
        -h|--help       ) usage ; exit 0 ;;
        *               ) error "Invalid option: $1\n$(usage)" 1 ;;
    esac
done


# Check mandatory options
if [[ -z "$mydir" ]]
then
    error "Library directory is missing. Exiting..."
    flag_error=1
elif [[ ! -d "$mydir" ]]
then
    error "$mydir is not a directory."
    flag_error=1
fi

if [[ -z "$myref" ]]
then
    error "Reference genome is missing. Exiting..."
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
if [[ -n "$mysites" && (! -s "$mysites" || -d "$mysites") ]]
then
    error "$mysites does not exist or is empty."
    flag_error=1
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
myGATK="java -Xmx2g -jar $myGATK_loc -R $myref "
mypicard="java -Xmx2g -jar $mypicard_loc "


# Local variables
lib_name=$(basename "$mydir")

pe1=$(ls "$mydir/${lib_name}"*_R1.fastq.gz 2> /dev/null)
pe2=$(ls "$mydir/${lib_name}"*_R2.fastq.gz 2> /dev/null)

if [[ ! -f "$pe1" ]]
then
    error "Unable to locate the R1 fastq file $pe1. Does the file exist?"
    flag_error=1
elif [[ ! -f "$pe2" ]]
then
    error "Unable to locate the R2 fastq file $pe2. Does the file exist?"
    flag_error=1
fi

# Multithread option
mythreads=$(( $(grep -ic processor /proc/cpuinfo) /2 ))
[[ -z $mythreads || $mythreads == 0 ]] && mythreads=1


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


# Some info are reported in the status file
info "Sample: $lib_name"
info "Number of threads: $mythreads"
[[ -n "$mysites" ]] && info "File used for base quality recalibration: $mysites"


#------------#
# Alignement #
#------------#
info "Alignment running...\n"

# Header for BAM files
rg="@RG\tID:$lib_name\tPL:illumina\tLB:$lib_name\tSM:$lib_name"

# Alignment
bwa mem -M -R "$rg" "$myref" "$pe1" "$pe2" > "$mydir/$lib_name.sam"

# BAM sorting and indexing
samtools view -bS -h "$mydir/$lib_name".sam | samtools sort - "$mydir/${lib_name}_sorted"
samtools index "$mydir/${lib_name}_sorted.bam"

# Indel realignment
$myGATK -I "$mydir/${lib_name}_sorted.bam" -T RealignerTargetCreator -o "$mydir/$lib_name.intervals"
$myGATK -I "$mydir/${lib_name}_sorted.bam" -T IndelRealigner -targetIntervals "$mydir/$lib_name.intervals" -o "$mydir/${lib_name}_sorted_realigned.bam"


#--------------------#
# Duplicates marking #
#--------------------#
info "Marking duplicates...\n"

$mypicard MarkDuplicates INPUT="$mydir/${lib_name}_sorted_realigned.bam" OUTPUT="$mydir/${lib_name}_sorted_realigned_MD.bam" METRICS_FILE="$mydir/${lib_name}_sorted_realigned_MD.log" VALIDATION_STRINGENCY=LENIENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=$(ulimit -n)

# BAM sorting and indexing
samtools index "$mydir/${lib_name}_sorted_realigned_MD.bam"


#-----------------------#
# Q-score recalibration #
#-----------------------#
info "Q-score recalibration...\n"
if [[ -n "$mysites" ]]
then
    # Extension setting
    myext="_recal"
    
    # Recalibration table
    $myGATK -I "$mydir/${lib_name}_sorted_realigned_MD.bam" -T BaseRecalibrator -knownSites "$mysites" -nct $mythreads -o "$mydir/${lib_name}_sorted_realigned_MD.grp"

    # Score rewritting
    $myGATK -I "$mydir/${lib_name}_sorted_realigned_MD.bam" -T PrintReads -BQSR "$mydir/${lib_name}_sorted_realigned_MD.grp" -nct $mythreads -o "$mydir/${lib_name}_sorted_realigned_MD${myext}.bam"

    # BAM indexing
    samtools index "$mydir/${lib_name}_sorted_realigned_MD${myext}.bam"
fi


#------------#
# Statistics #
#------------#
info "Generating statistics...\n"

echo -e "File: $mydir/${lib_name}_sorted_realigned_MD${myext}.bam" > "$mydir/${lib_name}_sorted_realigned_MD${myext}.flagstat"
samtools flagstat "$mydir/${lib_name}_sorted_realigned_MD${myext}.bam" >> "$mydir/${lib_name}_sorted_realigned_MD${myext}.flagstat"
echo -e "\n" >> "$mydir/${lib_name}_sorted_realigned_MD${myext}.flagstat"


#--------------#
# Housekeeping #
#--------------#
info "Cleaning the mess...\n"

rm "$mydir/${lib_name}.sam"* 
rm "$mydir/${lib_name}_sorted.ba"*
rm "$mydir/${lib_name}_sorted_realigned.ba"*
[[ -n "$myext" ]] && rm "$mydir/${lib_name}_sorted_realigned_MD."*  # Remove this file only if recal was done


exit 0
