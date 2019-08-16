#!/bin/bash
# Title: fastqc_report.sh
# Version: 1.1
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2013-10-31
# Modified in: 2016-02-11
# Licence : GPL v3



#======#
# Aims #
#======#

aim="Generate pdf report from FastQC analysis. For each library, reports are generated in html, then transformed in pdf. Finally, all reports are concatenated."



#==========#
# Versions #
#==========#

# v1.1 - 2016-02-11: usage message improved / exit number added for the last error message of the variable section / check for pdf files before concatenation
# v1.0 - 2015-09-20: htmldoc replaced by wkhtmltopdf
# v0.3 - 2015-07-25: email address updated / usage message updated / possibility to use symlinks for fastq files added
# v0.2 - 2015-02-10: usage message updated / info, warning, error functions added / check_dep moved to test_dep
# v0.1 - 2014-01-26: usage message moved to functions section / option loop updated / check presence of html document before pdf creation
# v0.0 - 2013-10-31: creation

version=$(grep -m 1 -i version "$0" | cut -d ":" -f 2 | sed "s/ * //g")



#===========#
# Functions #
#===========#

# Usage message
function usage {
    echo -e "
    \e[32m ${0##*/} \e[00m -d|--dir directory -o|--out name[.pdf] -r|--rep -c|--cat -h|--help

Aim: $aim

Version: $version

Options:
    -d, --dir    directory containing libraries (e.g. data)
    -o, --out    name of output file [default: fastqc_report.pdf]
    -r, --rep    generate reports only (no concatenation)
    -c, --cat    concatenate pdfs only (no pdf generatation)
    -h, --help   this message
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
function error {
    if [[ -t 1 ]]
    then
        echo -e "\e[31mError:\e[00m $1"
        exit $2
    else
        echo -e "Error: $1"
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

test_dep fastqc
test_dep wkhtmltopdf


#===========#
# Variables #
#===========#

# Options
while [[ $# -gt 0 ]]
do
    case $1 in
        -d|--dir  ) dir="$2"; shift 2 ;;
        -o|--out  ) output="${2%.pdf}.pdf"; shift 2;;
        -r|--rep  ) myrep=1 ; shift ;;
        -c|--cat  ) mycat=1 ; shift ;;
        -h|--help ) usage; exit 0 ;;
        *         ) error "Invalid option: $1\n$(usage)" 1
    esac
done


# Check the existence of obligatory options
if [[ -z "$dir" ]]
then
    error "Directory is required. Exiting...\n$(usage)" 1 
fi

if [[ -z "$output" ]]
then
    output="fastqc_report.pdf"
fi

if [[ -n $myrep && -n $mycat ]]
then
    error "The -r and -c options cannot be used at the same time. Exiting..." 1
elif [[ -z $myrep && -z $mycat ]]
then
    myrep=1
    mycat=1
fi



#============#
# Processing #
#============#

# FastQC analysis
if [[ -n $myrep ]]
then
    for i in $(find "$dir" -name *.fastq*)
    do

        # Check mime type of the file
        if [[ ! $(file -L -b --mime-type "$i" | grep gzip) && ! $(file -L -b --mime-type "$i" | grep text) ]]
        then
            warning "${i##*/} is not a fastq or fastq.gz"
            continue
        fi


        # Set output directory
        output_dir="$(dirname $i)/fastqc"
        if [[ ! -d "$output_dir" ]]
        then
            mkdir "$output_dir"
        fi

        # FastQC report
        fastqc --extract -q -o "$output_dir" "$i" &> /dev/null


        # PDF report
        fq_name=${i##*/}
        fastqc_dir="$output_dir/${fq_name%*.fastq.gz}_fastqc"
        
        if [[ -f "$fastqc_dir/fastqc_report.html" ]]
        then
            wkhtmltopdf --print-media-type -s Letter "$fastqc_dir/fastqc_report.html" "$fastqc_dir/fastqc_report.pdf" &> /dev/null
        else
            warning -e "The html file is not present in the $fastqc_dir directory."
        fi

    done
fi


# Report concatenation
if [[ -n $mycat ]]
then
    if [[ $(find $dir/ -name fastqc_report.pdf) ]]
    then
        gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile="$output" $(find $dir/ -name fastqc_report.pdf)
    else
        error "No pdf file found. Concatenation aborted." 1
    fi
fi

exit 0
