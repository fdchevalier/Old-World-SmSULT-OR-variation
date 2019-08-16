#!/usr/bin/env Rscript
# Title: exon_stat.R
# Version: 0.0
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2019-08-13
# Modified in:



#==========#
# Comments #
#==========#

# Aim: compute statistics regarding read depths and coverage of SmSULT-OR exons 



#===========#
# Variables #
#===========#

# GFF file
mygff="Smp_089320_exons.gff"

# Output file
myfl  <- "SmSULT-OR_stat.tsv"

# Data
dirs <- dir("data/")



#=================#
# Data processing #
#=================#

mygff <- read.delim(mygff, header=FALSE)

# RNA and CDS info
mRNA.lg <- sum(mygff[,5]-mygff[,4]) + nrow(mygff)
CDS.row <- grepl("CDS", mygff[,3])
CDS.lg  <- sum(mygff[CDS.row,5]-mygff[CDS.row,4]) + sum(CDS.row)

# Final matrix
mytb <- matrix(NA, ncol=7, nrow=length(dirs))

# Information from each library
for (i in dirs) {
    mybed     <- read.delim(paste0("data/",i,"/Smp_089320_exons.bed"), header=FALSE)
    mybed.CDS <- mybed[ rowSums(sapply(mygff[CDS.row, 4:5], function(x) findInterval((mybed[,3]), rev(x), rightmost.closed=TRUE)) == 1) > 0, ]

    mycov  <- nrow(mybed)/mRNA.lg * 100
    mymean <- mean(mybed[,4])
    mymed  <- median(mybed[,4])
    
    mycov.CDS  <- nrow(mybed.CDS)/CDS.lg * 100
    mymean.CDS <- mean(mybed.CDS[,4])
    mymed.CDS  <- median(mybed.CDS[,4])

    mytb[match(i, dirs), ] <- c(i, mymean, mymed, mycov, mymean.CDS, mymed.CDS, mycov.CDS)
}

# Colnames
colnames(mytb) <- c("Sample", "Exon mean read depth", "Exon median read depth", "Exon coverage", "CDS mean read depth", "CDS median read depth", "CDS coverage")

# Output table
write.table(mytb, myfl, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
