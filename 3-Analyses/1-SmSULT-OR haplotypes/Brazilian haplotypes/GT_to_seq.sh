#!/bin/bash
# Title: GT_to_seq.sh
# Version: 1.0
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2015-11-11
# Modified in: 2019-03-22
# Licence : GPL v3



#======#
# Aims #
#======#

# Generate fasta file of haplotype sequences from a table of Brazilian genotypes of SmSULT-OR.
# TODO: To allow this script to handle any genotype table, extensive modification of the processing section needs to be done.



#==========#
# Versions #
#==========#

# v1.0 - 2019-03-22: generate haplotype sequences instead of a concsensus sequence for each samples / indels handlded / script reformated
# v0.0 - 2015-11-11: creation



#===========#
# Variables #
#===========#

mygt="$1"
myseq="$2"

output="$3"



#============#
# Processing #
#============#

# Sample
for ((i=2 ; i<=$(wc -l < "$mygt") ; i++))
do

    # Get sequence identifier
	id=">$(sed -n "${i}p" "$mygt" | cut -d "," -f 1)"
	
	# Haplotype per sample
	for ((h=0 ; h<=1 ; h++))
	do
		echo "${id}_${h}" >> "$output"
		
		c=$(( $h + 1 ))
		
		# Replace bases of the initial sequence by the mutation
		## source: http://stackoverflow.com/a/24470022
		## Each position mentionned is the position-1 of the mutation (eg., c.103T>C => 102)
		## Make modification first on SNPs then on indels starting from the last to the first (I have not test another way)
		sed -r "s/^(.{102})T/\1$(sed -n "${i}p" "$mygt" | cut -d "," -f 2 | cut -d "/" -f $c)/ ; s/^(.{199})C/\1$(sed -n "${i}p" "$mygt" | cut -d "," -f 3 | cut -d "/" -f $c)/ ; s/^(.{616})G/\1$(sed -n "${i}p" "$mygt" | cut -d "," -f 6 | cut -d "/" -f $c)/  ; s/^(.{642})A/\1$(sed -n "${i}p" "$mygt" | cut -d "," -f 7 | cut -d "/" -f $c)/ ; s/^(.{423})GAA/\1$(sed -n "${i}p" "$mygt" | cut -d "," -f 5 | cut -d "/" -f $c)/ ; s/^(.{213})A/\1$(sed -n "${i}p" "$mygt" | cut -d "," -f 4 | cut -d "/" -f $c)/" <(sed -n "2p" "$myseq") >> "$output"
		
	done
    
    # Echo identifier
    echo $id | sed "s/^>//"
	
done

exit 0
