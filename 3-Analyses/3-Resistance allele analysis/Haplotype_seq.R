#!/usr/bin/env Rscript
# Title: Haplotype_seq.R
# Version: 1.0
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2018-12-14
# Modified in: 2019-05-08



#==========#
# Comments #
#==========#

# Aim: Determine the longest haplotype for the mutation of interest and draw a network.



#==========#
# Packages #
#==========#

# System
library(magrittr)

# Genomics and genetics
library(rtracklayer)
library(vcfR)
library(seqinr)
library(phangorn)
library(ape)

# Graph library
library(poppr)
library(gridExtra)
library(igraph) # for using layouts with poppr



#===========#
# Functions #
#===========#

source("Haplotype_seq_func.R")



#===========#
# Variables #
#===========#

#---------------------#
# Genes and genotypes #
#---------------------#

# Variant data files
myvcf.f   <- "data/Sm.SN.NE.TZ.OM_HR9.pChr6_3M.fb.raw.snp_indel.flt-dp4.cleaned.norm.phased.vcf.gz"
vcf.uph.f <- "data/Sm.SN.NE.TZ.OM_HR9.pChr6_3M.fb.raw.snp_indel.flt-dp4.cleaned.norm.vcf"

# Annotation file
## source: ftp://ftp.sanger.ac.uk/pub/pathogens/Schistosoma/mansoni/Archive/S.mansoni/genome/Gene_models/ARCHIVE/
mygff.f <- "~/data/sm_Gene_table/v5.07.08.12.chado.raw.gff"

# p.E142del position
mypos <- 1520405

# Filtering parameter
mygq        <- 20
max.missing <- 0.2 # for network analysis


#-------------#
# Populations #
#-------------#

# Population file
pop.ls.f <- "data/pop"

# Reference sample for p.E142del
ref.spl <- "HR9"


# Directory and file names
graph_dir <- "graphs/"
fig_file  <- "Fig. 6 - Haplotypes.pdf"


#=================#
# Data processing #
#=================#

#--------------#
# Loading data #
#--------------#

cat("\nLoading variant data...\n")

# Variant data
myvcf   <- read.vcfR(myvcf.f)
vcf.uph <- read.vcfR(vcf.uph.f)

# Genome annotation
mygff   <- as.data.frame(readGFF(mygff.f))

# Populations
pop.ls <- read.delim("data/pop", header=FALSE)


#---------------#
# Updating data #
#---------------#

# GQ filtering
myvcf.gq <- myvcf
myvcf.gq@gt[,-1][ extract.gt(myvcf, element="GQ") < mygq ] <- NA

# Update the phased data with missing data (because Beagle output missing data as reference allele)
myvcf@gt[,-1][ is.na(extract.gt(vcf.uph)) ] <- NA


#------------------#
# Gene of interest #
#------------------#

# For v5 GFF
myexons <- mygff[ mygff[,3] == "CDS" & grepl("Smp_089320.1", mygff[,9]), c(1,4:5) ]

# Obtain interval
myitv <- matrix(NA, ncol=nrow(myexons)+1, nrow=nrow(getFIX(myvcf)))

myitv[,1] <- getFIX(myvcf)[,1] == unique(myexons[,1])

for (r in 1:nrow(myexons)) {
    myitv[,r+1] <- sapply(getFIX(myvcf)[,2], function(x) as.numeric(x) >= myexons[r,2] & as.numeric(x) <= myexons[r,3])
}

myitv <- myitv[,1] == TRUE & rowSums(myitv[,2:ncol(myitv)]) > 0

# Create table with info from the fix
myvcf.cds <- getFIX(myvcf)[ myitv, ]


#----------------------#
# Mutation of interest #
#----------------------#

# Select lines of interest (including window around the SNPs of interest)
myvcf.fix <- getFIX(myvcf)
myvcf.row <- myvcf.fix[,1] == myvcf.cds[1,1]

myvcf <- myvcf[myvcf.row,]

# p.E142del genotype and samples
p.E142del.row <- grep(mypos, getPOS(myvcf))
p.E142del.ref <- getREF(myvcf)[p.E142del.row]
p.E142del.alt <- getALT(myvcf)[p.E142del.row] %>% strsplit(",") %>% unlist()
for (a in p.E142del.alt) { if (grepl(a, p.E142del.ref)) {p.E142del.alt <- match(a,p.E142del.alt) ; break } }
p.E142del.spl <- colnames(myvcf@gt[,-1])[grep(p.E142del.alt, extract.gt(myvcf)[p.E142del.row, ])]


#-----------------#
# Haplotype block #
#-----------------#

cat("\nIdentifying haplotype block...\n")

# Genotype data
mygt.flt <- extract.gt(myvcf)
mygt.flt <- mygt.flt[ , colnames(mygt.flt) %in% p.E142del.spl ]

# Haplotype matrix
mymat <- data.frame(matrix(NA, nrow=nrow(mygt.flt), ncol=ncol(mygt.flt)))

i <- 1
for (c in 1:ncol(mygt.flt)) {
    mymat[,i] <- t(sapply(strsplit(as.character(mygt.flt[,c]),"|",fixed=TRUE), function(x) c(x[1],x[2])))
    i <- i+1
}

colnames(mymat) <- colnames(mygt.flt)
rownames(mymat) <- rownames(mygt.flt)

# Transoforming genotype in color
mygt.flt.clr <- as.matrix(mymat)

# Select only variable sites
myrow.flt    <- apply(mygt.flt.clr, 1, function(x) all(! is.na(x)) & length(unique(x[!is.na(x)])) > 1)
mygt.flt.clr <- mygt.flt.clr[myrow.flt , ]

# Determining the longest haplotype
myrow     <- grep(mypos, rownames(mygt.flt.clr))
mypos.blk <- matrix(myrow, ncol=2, nrow=2)
myidx     <- 0

for (i in grep(ref.spl, colnames(mygt.flt.clr))) {
    myidx <- myidx + 1

    for(s in grep(ref.spl, colnames(mygt.flt.clr), invert=TRUE)) {
        # Upstream
        for (r in myrow:1) {
            if (mygt.flt.clr[r,i] != mygt.flt.clr[r,s]) {
                #myup <- r+1
                myup <- r
                break
            }
        }
        
        # Downstream
        for (r in myrow:nrow(mygt.flt.clr)) {
            if (mygt.flt.clr[r,i] != mygt.flt.clr[r,s]) {
                #mydw <- r-1
                mydw <- r
                break
            }
        }
        
        if (mypos.blk[1,myidx] > myup && mypos.blk[2,myidx] < mydw) {
            mypos.blk[,myidx] <- c(myup, mydw)
        }

    }
}

# Select the haplotype carrying for p.E142del (for haplotype tree)
p.E142del.hap <- colnames(mygt.flt.clr)[mygt.flt.clr[grep(mypos, rownames(mygt.flt.clr)),] == 1]

# Select the haplotype of reference
ref.hap <- grep(ref.spl, colnames(mygt.flt.clr))[which.max(apply(mypos.blk,2,diff))]

# Select the longest haplotype block
mypos.blk <- mypos.blk[ , which.max(apply(mypos.blk,2,diff)) ]

# Refine table
mygt.flt.clr <- mygt.flt.clr[mypos.blk[1]:mypos.blk[2], ]

# Relative bp position of the GT regarding the deletion
mygt.pos <- rownames(mygt.flt.clr) %>% strsplit(., "_") %>% lapply(., function(x) rev(x)[1]) %>% unlist() %>% as.numeric() - mypos

# Painting and transposing for easy plotting
mygt.flt.clr <- apply(mygt.flt.clr, 1, function(x) {
      if (length(unique(x[!is.na(x)])) == 2) {
          x[ x == x[ref.hap] ] <- "white" #"black"
          x[ x != x[ref.hap] ] <- "black" #"white"
      } else if (length(unique(x[!is.na(x)])) > 2) {
          x[ x == x[ref.hap] ] <- "white" #"black"

          mygt.tmp   <- unique(x[ x != x[ref.hap] ])
          mygray.tmp <- gray.colors(length(mygt.tmp), start=0, end=1/length(mygt.tmp))

          for (i in 1:length(mygt.tmp)) { x[ x == mygt.tmp[i] ] <- mygray.tmp[i] }
      }
      return(x)
})

# Real block position
mypos.blk <- mygt.pos + mypos
mypos.blk <- mypos.blk[c(2,(length(mypos.blk)-1))] # Exclude the first and last that show the haplotype breaks

# Length of the haplotype
cat("Haplotype block size:", diff(mypos.blk), "bp\n")

# Update VCF for each 
myrow.flt <- findInterval(getPOS(myvcf), mypos.blk, rightmost.closed=TRUE) == 1
myvcf.bk <- myvcf
myvcf <- myvcf[ myrow.flt, ]


#-------------------#
# Haplotype network #
#-------------------#

cat("\nPreparing haplotype network...\n")

# Haplotype block
mypos.blk <- mygt.pos + mypos
mypos.blk <- mypos.blk[2:(length(mypos.blk)-1)]

myvcf <- myvcf[ getPOS(myvcf) %in% mypos.blk, ]

# Genotype data (with maximum missing data threshold)
mygt <- extract.gt(myvcf)
mygt <- extract.gt(myvcf[ ! (is.na(mygt) %>% rowSums() / ncol(mygt)) > max.missing , ])

# Haplotype data
mymat <- data.frame(matrix(NA, nrow=nrow(mygt), ncol=ncol(mygt)))

i <- 1
for (c in 1:ncol(mygt)) {
    mymat[,i] <- t(sapply(strsplit(as.character(mygt[,c]),"|",fixed=TRUE), function(x) c(x[1],x[2])))
    i <- i+1
}

colnames(mymat) <- colnames(mygt)
rownames(mymat) <- rownames(mygt)

mygt <- as.matrix(mymat)

# Conversion to genlight
mygl <- new('genlight', t(mygt), n.cores=NULL, ploidy=1)

# Population names
mypop <- gsub(".[1|2]$", "", mygl@ind.names) %>% sapply(., function(x) pop.ls[ match(x, pop.ls[,1]),2]) %>% as.vector()
mypop[ is.na(mypop) ] <- "Caribbean"

# Mutation type data
mymut <- rep("", length(mygl@ind.names))
mymut[mygl@ind.names %in% p.E142del.hap] <- "del"

# Add strata
strata(mygl) <- data.frame(mypop, mymut)
setPop(mygl) <- ~mypop + mymut

# Distance matrix
mydist <- bitwise.dist(mygl, missing_match=F)



#========#
# Graphs #
#========#

# Check graph folder
if (! dir.exists(graph_dir)) {dir.create(graph_dir)}

myrow <- nrow(mygt.flt.clr)
mycln <- ncol(mygt.flt.clr)

#pdf(paste0(title1,".pdf"), width=mycln/15, height=myrow/15+1.5+10)
pdf(paste0(graph_dir,fig_file), width=mycln/15, height=myrow/15+1.5+5, useDingbats=FALSE)

#layout(matrix(1:2, ncol=1), heights=c(2,5))
layout(matrix(1:2, ncol=1), heights=c(myrow/15+1.5,5))

#layout(matrix(1:3, ncol=1), heights=c(2,10,2))
# Save plot to an object using a null PDF device
# http://stackoverflow.com/a/14742001/120898
#pdf(NULL)
#dev.control(displaylist="enable")

#-----------------#
# Haplotype block #
#-----------------#

cat("\nGenerating haplotype block plot...\n")

mygt.pos2 <- mygt.pos
mygt.pos2[ ! mygt.pos %in% mygt.pos[ c(2, which(mygt.pos == 0), length(mygt.pos)-1) ] ] <- NA

mygraph(mygt.flt.clr, grepl(paste0(p.E142del.spl,"[.]",collapse="|"), rownames(mygt.flt.clr), perl=TRUE), mygt.pos2)

# Panel letter
mtext("A", side=3, line=-1, at=line2user(par("mar")[2],2), cex=par("cex")*2, adj=0)

#-------------------#
# Haplotype network #
#-------------------#

cat("\nGenerating haplotype network plot...\n")

# Populations
mygl.pop <- popNames(mygl)

# Colors
myclr <- rep(NA, length(mygl.pop))
names(myclr) <- mygl.pop
myclr[ grep("del", mygl.pop, invert=T) ] <- gray.colors(length(grep("del", mygl.pop, invert=T)), start=1, end=0)
myclr[ mygl.pop == "Caribbean_del" ]     <- "red"
myclr[ mygl.pop == "W-Af_del" ]          <- "orange"
myclr[ mygl.pop == "E-Af_del" ]          <- "yellow"

# Network
## threshold argument can be use to condense representation (e.g., threshold=0.02)
myntw <- poppr.msn(mygl, mydist, palette=myclr, vertex.label=NA, clustering.algorithm="average_neighbor", showplot=FALSE)

# Replace single pie with circle
V(myntw$graph)$shape[ sapply(V(myntw$graph)$pie, length) == 1 ] <- "circle"
V(myntw$graph)$color <- rep(NA, length(V(myntw$graph)$name))
V(myntw$graph)$color[ sapply(V(myntw$graph)$pie, length) == 1 ] <- unlist(V(myntw$graph)$pie.color[ sapply(V(myntw$graph)$pie, length) == 1 ])

# Plot
par(mar=rep(0,4))
set.seed(930734420)
plot_poppr_msn(mygl, myntw, inds="NA", nodelab=NA, wscale=FALSE, scale.leg=FALSE, pop.leg=FALSE, size.leg=FALSE, layfun = layout_with_lgl)

# Legend
myleg <- gsub("E-Af", "East African samples", myntw$population) %>% gsub("W-Af", "West African samples", .) %>% gsub("Caribbean", "Caribbean sample", .) %>% gsub("Oman", "Omani samples", .) %>% gsub("del", " with p.E142del", .) %>% gsub("_", "", .)
legend("bottomright", myleg, pch=21, bty="n", pt.bg=myntw$colors) #, xjust=0.5, yjust=0.5)

# Panel letter
mtext("B", side=3, line=5, at=line2user(par("mar")[2],2), cex=par("cex")*2, adj=0)

# Close PDF device
dev.off()
