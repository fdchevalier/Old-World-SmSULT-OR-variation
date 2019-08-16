#!/bin/bash
# Title: lib_stat.sh
# Version: 1.2
# Author: Frédéric CHEVALIER <fcheval@txbiomedgenetics.org>
# Created in: 2013-10-16
# Modified in: 2016-05-13
# Licence : GPL v3



#======#
# Aims #
#======#

aim="Generate table of stat for each libraries using log and flagstat files generated with the exome pipeline"



#==========#
# Versions #
#==========#

# v1.2 - 2016-05-13: update script with the new stat files from the new version of the exome pipline (v1.x)
# v1.1 - 2015-08-17: script updated regarding the new version of samtools that induces some changes in the flagstat file / version number added in usage
# v1.0 - 2015-01-28: new message and test functions added / usage message updated / test for dependencies updated / option loop updated /initial .tab file extension changed for .tsv / new fields in table added (%_read_mapped, Ratio_ZW, sex) / calculation methods for percentage of read mapped
# v0.1 - 2014-01-12: option loop updated / bug correction in getting results
# v0.0 - 2013-10-16: creation

version=$(grep -m 1 -i "version" "$0" | cut -d ":" -f 2 | sed "s/ * //g")



#===========#
# Functions #
#===========#

# Usage message
function usage {
    echo -e "
    \e[32m ${0##*/} \e[00m -d|--dir path_directory -r|--rd file.log -c|--cap-eff file.log -o|--out name[.tsv] -h|--help

Aim: $aim

Version: $version

Options:
    -d, --dir       directory containing library sub-directories (eg., data)
    -r, --rd        read depth log [default: *_sorted_realigned_MD*.log]
    -c, --cap-eff   capture efficiency log [default: *capture_efficiency.log]
    -o, --out       name of output file [default stat.tab]
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


function test_var {
    if [[ -z $(eval echo \$"$1") ]]
    then
        eval "$1"=.
    fi
}



#==============#
# Dependencies #
#==============#

test_dep bc 



#===========#
# Variables #
#===========#

# Options
while [[ $# -gt O ]]
do
    case $1 in
        -d|--dir     ) dir="$2"; shift 2 ;;
        -o|--out     ) output="${2%.tsv}.tsv"; shift 2 ;;
        -r|--rd      ) myrd="$2" ; shift 2 ;;
        -c|--cap-eff ) myce="$2" ; shift 2 ;;
        -h|--help    ) usage; exit 0 ;;
        *            ) error "Invalid option: $1\n$(usage)" 1 ;;
    esac
done


# Check the existence of obligatory options
if [[ -z "$dir" ]]
then
    error "Directory is missing.Exiting...\n$(usage)" 1
fi

if [[ ! -d "$dir" ]]
then
    error "$dir is not a directory. Exiting..." 1
fi

if [[ -z "$myrd" ]]
then
    myrd=*_sorted_realigned_MD*.log
fi

if [[ -z "$myce" ]]
then
    myce=*capture_efficiency.log
fi

if [[ -z "$output" ]]
then
    output="stat.tsv"
fi



#============#
# Processing #
#============#

# Head of the table
echo -e "Library\tMean_bait_region_read_depth\tMedian_bait_region_read_depth\tTotal_reads\tMapped_reads\tPaired_reads\tSingletons\t%_read_mapped\tDuplicates\t%_baits_region\tMean_RD\tMedian_RD\tRatio_ZW\tSex" > "$output"


# Log and flagstat analysis
for i in $(ls "$dir")
do

    dir=$(echo "$dir" | sed "s,/$,,")
    i=$(echo "$i" | sed "s,/$,,")


    # Get log and flagstat files path
    myrd_log=$(ls -1 "$dir/$i"/$myrd 2> /dev/null | tail -1)
    myce_log=$(ls -1 "$dir/$i"/$myce 2> /dev/null | tail -1)
    flag=$(echo "$dir/$i"/*.flagstat)
    sexlog=$(echo "$dir/$i"/*.sex)
    

    # Log file analysis
    if [[ -s "$myrd_log" ]]
    then
        mean_rd=$(sed -n "/baits-0.bed/,/^data/p" "$myrd_log" | sed "s/ * / /g" | sed -n "/^ [0-9]/p" | cut -d " " -f 5 | sed -n "1p") ; test_var mean_rd
        med_rd=$(sed -n "/baits-0.bed/,/^data/p" "$myrd_log" | sed "s/ * / /g" | sed -n "/^ [0-9]/p" | cut -d " " -f 4 | sed -n "1p")  ; test_var med_rd
        
        # Building the first part of the result line
        part1="$mean_rd\t$med_rd"
    else
        part1=".\t."
        warn_msg="read depth log file is"
    fi


    # Flag file analysis
    if [[ -s "$flag" ]]
    then
        tot=$(grep -i "total" "$flag" | cut -d " " -f 1)              ; test_var tot
        dup=$(grep -i "duplicates" "$flag" | cut -d " " -f 1)         ; test_var dup
        map=$(grep -m 1 -i "mapped" "$flag" | cut -d " " -f 1)        ; test_var map
        paired=$(grep -i "properly paired" "$flag" | cut -d " " -f 1) ; test_var paired
        sing=$(grep -i "singleton" "$flag" | cut -d " " -f 1)         ; test_var sing
        pc_map=$(bc <<< "scale=2; $map/$tot*100")                     ; test_var pc_map # percentage mapped

        # Building the second part of the result line
        part2="$tot\t$map\t$paired\t$sing\t$pc_map\t$dup"
    else
        part2=".\t.\t.\t.\t.\t."
        if [[ -n "$warn_msg" ]]
        then
            warn_msg=$(echo "$warn_msg" | sed "s/ is/, flag file are/g")
        else
            warn_msg="flag file is"
        fi
    fi
    
    # Log file analysis
    if [[ -s "$myce_log" ]]
    then
        # Building the third part of the result line
        part3=$(sed -n "2p" "$myce_log")
    else
        part3=".\t.\t."
        warn_msg="capture efficiency log file is"
    fi


    # Sex log analysis
    if [[ -s "$sexlog" ]]
    then
        ratio=$(sed -n "1p" "$sexlog" | sed "s/\"//g" | grep -o ".\...")  ; test_var ratio
        sex=$(sed -n "2p" "$sexlog" | sed "s/\"//g")                      ; test_var sex

        # Building the fourth part of the result line
        part4="$ratio\t$sex"
    else
        part4=".\t."
        if [[ -n "$warn_msg" ]]
        then
            warn_msg=$(echo "$warn_msg" | sed "s/ is\| are/, sex log are/g")
        else
            warn_msg="sex log is"
        fi
    fi


    # Adding line to the table
    echo -e "${i}\t${part1}\t${part2}\t${part3}\t${part4}" >> "$output"


    # Warning message if any
    if [[ -n "$warn_msg" ]]
    then
        warning "$i: $warn_msg missing."
    fi

    unset warn_msg


done


exit 0
