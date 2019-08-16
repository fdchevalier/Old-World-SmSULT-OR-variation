#!/usr/bin/env Rscript
# Title: OXA-R_map.R
# Version: 1.0
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2018-12-14
# Modified in: 2019-05-08



#==========#
# Comments #
#==========#

# Aim: Map the OXA-R allele frequency and OXA-R worms from the different populations.



#==========#
# Packages #
#==========#

# Genetics and genomics
library(vcfR)
library(rtracklayer)
library(seqinr)
library(strataG)

# Graphics
library(rgdal)      # Reading and projecting shapefiles
library(plotrix)    # Creating color scales

library(mapplots)   # add.pie function

library(RColorBrewer)



#===========#
# Functions #
#===========#

# Convert geographical coordinates to xy values with a given projection
coord2xy <- function(lon, lat, proj) {

    # lon   longitude vector or matrix of long x lat
    # lat   latitude
    # proj  projection

    # Check variables
    if (is.vector(lon)) {
        if (missing(lon) | ! is.numeric(lon) | ! is.vector(lon)) { stop("lon must be a vector or a matrix of numeric values.") }
        if (missing(lat) | ! is.numeric(lat) | ! is.vector(lat)) { stop("lat must be a vector of numeric values.") }
        if (length(lon) != length (lat)) { stop("lon and lat must have the same length.") }
    } else if (is.matrix(lon) | is.data.frame(lon)) {
        if (! is.numeric(lon) | ncol(lon) != 2) { stop("lon must be a vector or a matrix of numeric values.") }
        if (missing("proj")) { proj <- lat }
    }

    if (! is.character(proj) | length(proj) != 1) { stop("proj must be a character vector of a single value.") }

    # Transform geographical coordinates into graph coordinates
    ## source: https://gis.stackexchange.com/a/45266
    if (is.vector(lon)) {
        mycoord <- data.frame(lon=lon, lat=lat)
    } else {
        mycoord <- as.data.frame(lon)
        colnames(mycoord) <- c("lon", "lat")
    }
    coordinates(mycoord) <- c("lon", "lat")
    proj4string(mycoord) <- CRS("+init=epsg:4326") # WGS 84
    mycoord <- spTransform(mycoord, CRS(proj))

    return(mycoord@coords)
}



#====================#
# Data and variables #
#====================#

#---------------------#
# Genes and genotypes #
#---------------------#

# Variant data files
myvcf.f   <- "data/Sm.SN.NE.TZ.OM_HR9.pChr6_3M.fb.raw.snp_indel.flt-dp4.cleaned.norm.phased.vcf.gz"
vcf.uph.f <- "data/Sm.SN.NE.TZ.OM_HR9.pChr6_3M.fb.raw.snp_indel.flt-dp4.cleaned.norm.vcf"

# Annotation file
## source: ftp://ftp.sanger.ac.uk/pub/pathogens/Schistosoma/mansoni/Archive/S.mansoni/genome/Gene_models/ARCHIVE/
mygff.f <- "~/data/sm_Gene_table/v5.07.08.12.chado.raw.gff"

# Transcript of interest
myRNA <- "Smp_089320.1"


#-----------#
# Mutations #
#-----------#

# Mutation annotations
mymut.ls.f <- "data/Non-syn_mutation_annotation.tsv"
mymut.f    <- "data/Smp_089320.flt.norm.tsv"

# Missing data (if true missing data will be considered as sensistive alleles)
# This will consider any missing data in the resistant type as sensitive (if any sensitive allele is present)
allow.msg <- TRUE

# Add data from South America (Chevalier et al. 2016, https://doi.org/10.1016/j.ijpara.2016.03.006)
if (allow.msg) {
    myalleles.BR <- c(7,0,371)
} else {
    myalleles.BR <- c(7,0,353)
}
myzyg.BR <- c(2,3)


#-------------#
# Populations #
#-------------#

# Population file
pop.ls.f <- "data/pop"

# Population names
mypop    <- c("BR", "W-Af", "E-Af", "Oman")
mypop.nm <- c("South America", "W.Af", "E.Af", "Oman")

# Type of mutations
mytype     <- c("Deleterious Validated", "Deleterious Putative", "Neutral")
mytype.tag <- c("R", "PR", "S")


#-----#
# Map #
#-----#

# Shape file directory for geographical representation
## source: https://www.naturalearthdata.com/downloads/110m-cultural-vectors/ (download all)
## source: https://www.naturalearthdata.com/downloads/110m-physical-vectors/ (download all)
shp.fld     <- "~/data/natural_earth_data"

# Projection
proj <- "+proj=natearth"

# Country code and point coordinates (GPS average of sampling sites)
## source: https://www.dcode.fr/geographic-coordinates-calculation
myct.code <- c("BR",      "SN",    "NE",    "TZ",  "OM")
## comment below to have the center of country by default
myct.ctr  <- matrix(c(
                -41.503889, -16.752778,
                -16.149123,  15.736788,
                  1.281175,  14.334580,
                 32.749469,  -2.514305,
                 54.305932,  17.187315
                   ), ncol=2, byrow=TRUE)

# Pie coordinates
mypie.cd <- matrix(c(
                -20, -16.75,
                -35,  30,
                 70, -16.5,
                 80,  13
                     ), ncol=2, byrow=TRUE)

# Legend coordinates
leg.cd <- c(-170, 0)

# Directory and file names
graph_dir <- "graphs/"
fig_file  <- "Fig. 5 - OXA-R map.pdf"



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

cat("\nLoading geographical data...\n")

# Shape files
countries   <- readOGR(dsn=shp.fld, layer="ne_110m_admin_0_countries")
graticules  <- readOGR(dsn=shp.fld, layer="ne_110m_graticules_10")
bbox        <- readOGR(dsn=shp.fld, layer="ne_110m_wgs84_bounding_box")

# Mutation annotations
mymut.ls <- read.delim(mymut.ls.f, head=FALSE)
mymut    <- read.delim(mymut.f,    head=TRUE)

# Populations
pop.ls <- read.delim("data/pop", header=FALSE)


#------------------#
# Gene of interest #
#------------------#

# For v5 GFF
myexons <- mygff[ mygff[,3] == "CDS" & grepl(myRNA, mygff[,9]), c(1,4:5,7) ]

# Exon interval
myitv <- matrix(NA, ncol=nrow(myexons)+1, nrow=nrow(getFIX(myvcf)))

myitv[,1] <- getFIX(myvcf)[,1] == unique(myexons[,1])

for (r in 1:nrow(myexons)) {
    myitv[,r+1] <- sapply(getFIX(myvcf)[,2], function(x) as.numeric(x) >= myexons[r,2] & as.numeric(x) <= myexons[r,3])
}

myitv <- myitv[,1] == TRUE & rowSums(myitv[,2:ncol(myitv)]) > 0

# Strand orientation
mystrand <- unique(myexons[,4])

if (length(mystrand) != 1) { stop("Unable to determine DNA strand") }


#---------------------------------#
# Missing data in the phased data #
#---------------------------------#

# Update the phased data with missing data (because Beagle output missing data as reference allele)
myvcf@gt[,-1][ is.na(extract.gt(vcf.uph)) ] <- NA


#------------------#
# Gene of interest #
#------------------#

# Select lines of interest
myvcf <- myvcf[myitv,]

# Update my mutation file to get positions
mymut <- merge(mymut[,c(8,1:2)], mymut.ls, by=1, sort=FALSE) %>% unique()

# Update alternative allele sequence from mutation code
mymut[,3] <- sapply(as.character(mymut[,3]), function(x) {
                        y <- strsplit(x, ">") %>% unlist() %>% rev() %>% .[1] %>% strsplit (., "") %>% unlist() %>% comp(., forceToLower=FALSE)
                        if (mystrand == "-") {
                            rev(y) %>% paste(., collapse="")
                        } else {
                            paste(y, collapse="")
                        }
                    })
colnames(mymut)[3] <- "ALT"

# Build haplotype matrix
mygt      <- extract.gt(myvcf)
myvcf.fix <- getFIX(myvcf)
myhap     <- matrix(NA, nrow=nrow(mygt), ncol=ncol(mygt)*2)

i <- 1
for (c in 1:ncol(mygt)) {
    myhap[,i:(i+1)] <- t(sapply(strsplit(as.character(mygt[,c]),"|",fixed=TRUE), function(x) c(x[1],x[2])))
    i <- i+2
}

colnames(myhap) <- unlist(lapply(colnames(mygt), function(x) paste0(x, c(".1", ".2"))))

# Determine allelic genotypes
mymut[,ncol(mymut)+1] <- NA

for (i in 1:nrow(mymut)) {
    mymut[i,ncol(mymut)] <- grep( paste0("^",mymut[i,3],"$"), strsplit(myvcf.fix[myvcf.fix[,2] == mymut[i,2], 5], ",") %>% unlist() )
}


myalleles <- rep(NA, ncol(mygt)*2)
names(myalleles) <- colnames(myhap)

for (s in 1:ncol(myhap)) {

    na.flag <- FALSE

    for (t in 1:length(mytype)) {
        mytype.tmp <- strsplit(mytype[t], " ") %>% unlist()
        mypos      <- unique(mymut[apply(mymut[,4:5], 1, function(x) sum(grepl(paste0(mytype.tmp,collapse="|"), x)) == length(mytype.tmp)), 2])
        mystate    <- sapply(mypos, function(x) grepl(paste0("^",mymut[ mymut[,2] == x, 6 ],"$",collapse="$|^"), myhap[ myvcf.fix[,2] == x , s ]) )

        if ( ! allow.msg & any(is.na(mystate)) & t < length(mytype) ) na.flag <- TRUE
        
        if ( all( ! mystate, na.rm=TRUE) & t == length(mytype) & ! na.flag) {
            myalleles[s] <- mytype.tag[t]
            break
        }

        if ( all(is.na(mystate)) ) {
            next
        } else if ( any(mystate, na.rm=TRUE) ) {
            myalleles[s] <- mytype.tag[t]
            break
        }

    }

}

# Determine worm status (zygocity state for resistance alleles)
myzyg.r <- rep(NA, ncol(mygt))
names(myzyg.r) <- colnames(mygt)

myzyg.pr <- myzyg.r

for (s in 1:length(myzyg.r)) {
    myzyg.tmp <- myalleles[grep(paste0(names(myzyg.r[s]),"."), names(myalleles), fixed=TRUE)]

    if ( all(myzyg.tmp == mytype.tag[3]) | any(is.na(myzyg.tmp)) ) next
    
    if ( all(myzyg.tmp == mytype.tag[1]) ) {
        myzyg.r[s] <- "hmz"
    } else if ( any(myzyg.tmp == mytype.tag[1]) ) {
        myzyg.r[s] <- "htz"
    } 
    
    if ( all(grepl(paste0("^",mytype.tag[1:2],"$", collapse="|"), myzyg.tmp)) ) {
        myzyg.pr[s] <- "hmz"
    } else if ( any(grepl(paste0("^",mytype.tag[1:2],"$", collapse="|"), myzyg.tmp)) ) {
        myzyg.pr[s] <- "htz"
    }
}

# Summarize data per population
myalleles.pop     <- matrix(NA, ncol=length(mypop), nrow=length(mytype))
myalleles.pop[,1] <- myalleles.BR
colnames(myalleles.pop) <- mypop
rownames(myalleles.pop) <- mytype

myzyg.r.pop  <- matrix(NA, ncol=length(mypop), nrow=2)
myzyg.pr.pop <- myzyg.r.pop

# Add the Brazilian data
myzyg.r.pop[,1]  <- myzyg.BR
myzyg.pr.pop[,1] <- myzyg.BR

for (p in 2:length(mypop)) {
    myspl <- paste0(pop.ls[ pop.ls[,2] == mypop[p], 1], collapse="|")

    for (t in 1:length(mytype)) { myalleles.pop[t,p] <- sum( myalleles[ grepl(myspl,names(myalleles)) ] == mytype.tag[t], na.rm=TRUE ) }
    
    myzyg.r.pop[1,p] <- length( grep("hmz", myzyg.r[ grepl(myspl,names(myzyg.r)) ]) )
    myzyg.r.pop[2,p] <- length( grep("htz", myzyg.r[ grepl(myspl,names(myzyg.r)) ]) )

    myzyg.pr.pop[1,p] <- length( grep("hmz", myzyg.pr[ grepl(myspl,names(myzyg.pr)) ]) )
    myzyg.pr.pop[2,p] <- length( grep("htz", myzyg.pr[ grepl(myspl,names(myzyg.pr)) ]) )
}


#-------------------#
# Statistic summary #
#-------------------#

cat("\nGenetic statistics:\n")

# Transition/Transversion ratio
mydb <- vcfR2DNAbin(myvcf)
cat("Transition/Transversion ratio:\n")
TiTvRatio(mydb)

# Allele frequency
cat("Allele frequency of OXA-R alleles:\n")
print(apply(myalleles.pop, 2, function(x) x[1]/sum(x)*100))

mypop.nb <- apply(myalleles.pop, 2, function(x) sum(x)/2)



#========#
# Graphs #
#========#

cat("\nGenerating map...\n")

#------------#
# Parameters #
#------------#

# Projection
countries   <- spTransform(countries,  CRS(proj))
bbox        <- spTransform(bbox,       CRS(proj))
graticules  <- spTransform(graticules, CRS(proj))

# Center of countries
if (! exists("myct.ctr")) {
    myct.ctr  <- matrix(NA, ncol=2, nrow=length(myct.code))

    for (i in 1:length(myct.code)) {
        j <- which(countries$ISO_A2 == myct.code[i])
        myct.ctr[i,] <- countries@polygons[[j]]@labpt
    }
} else {
    myct.ctr <- coord2xy(myct.ctr, proj)
}

# Pie coordinates
mypie.cd <- coord2xy(mypie.cd, proj)

# Lines coordinates
mylines.cd          <- matrix(NA, ncol=4, nrow=4)
mylines.cd[,c(1,3)] <- as.matrix(myct.ctr[c(1:2,4:5),])
mylines.cd[,c(2,4)] <- mypie.cd

## Adjust line start for SN/NE
mylines.cd[2,1] <- mean(myct.ctr[2:3,1]) 
mylines.cd[2,3] <- mean(myct.ctr[2:3,2]) 

# Legend coordinates
leg.cd <- coord2xy(leg.cd[1], leg.cd[2], proj)


# Limit of the map
mymap.cd <- matrix(NA, nrow=length(myct.code), ncol=4)

for (i in myct.code) {
    myct.id <- grep(i, countries$ISO_A2)
    mycd    <- countries@polygons[[ myct.id ]]@Polygons[[1]]@coords
    mymap.cd[match(i,myct.code),] <- c(min(mycd[,1]), max(mycd[,1]), min(mycd[,2]), max(mycd[,2]))
}

mymap.cd <- c(min(mymap.cd[,1]), max(mymap.cd[,2]), min(mymap.cd[,3]), max(mymap.cd[,4]))

## Adjust map coordinates with pies
if (mymap.cd[1] > min(mylines.cd[,1:2])) { mymap.cd[1] <- min(mylines.cd[,1:2]) }
if (mymap.cd[2] < max(mylines.cd[,1:2])) { mymap.cd[2] <- max(mylines.cd[,1:2]) }
if (mymap.cd[3] > min(mylines.cd[,3:4])) { mymap.cd[3] <- min(mylines.cd[,3:4]) }
if (mymap.cd[4] < max(mylines.cd[,3:4])) { mymap.cd[4] <- max(mylines.cd[,3:4]) }

## Add border
mymap.cd <- (abs(mymap.cd * 1.1) + 1e6) * sign(mymap.cd)

# Colors
clr <- rep("#E6E6E6", length(countries$ISO_A2))
clr[ grepl(paste0(myct.code,collapse="|"), countries$ISO_A2) ] <- "#a5a4a4"

clr.line <- "grey50"
clr.pie  <- rev(brewer.pal(3, "YlOrRd"))

myfct   <- diff(mymap.cd[1:2])/diff(mymap.cd[3:4])


#-------------------#
# Ploting world map #
#-------------------#

# Check graph folder
if (! dir.exists(graph_dir)) {dir.create(graph_dir)}

# Start PDF device
mywidth <- 8
pdf(paste0(graph_dir, fig_file), width=mywidth, height=(mywidth/myfct))

# Remove margins
par(mar=rep(0,4))

# Background, graticules and countries
plot(mymap.cd[1:2], mymap.cd[3:4])
plot(bbox,       col="white", border=NA, add=TRUE)
plot(graticules, col="#CCCCCC33", lwd=1, add=TRUE)
plot(countries,  col=clr, border="#AAAAAA", lwd=0.5, add=TRUE)

# Border of the map
plot(bbox, border="grey90", lwd=1, add=TRUE)

# Finer control of pch than points function allows
## source: https://stackoverflow.com/a/2580209
with(as.data.frame(myct.ctr), symbols(x=myct.ctr[,1], y=myct.ctr[,2], circles=rep(1e5,nrow(myct.ctr)), inches=FALSE, ann=FALSE, bg=clr.line, fg=NULL, add=TRUE))

# Specific line between SN and NE
lines(myct.ctr[2:3,], col=clr.line)

for (i in 1:ncol(myalleles.pop)) {
    
    # Line between center of the country and the pie chart
    lines(mylines.cd[i,1:2], mylines.cd[i,3:4], col=clr.line) #, lwd=0.5)
    
    # Pie
    old.par <- par(lwd=0.1) # Changing of lwd for add.pi
    add.pie(myalleles.pop[-2,i], x=mypie.cd[i,1], y=mypie.cd[i,2], col=clr.pie[-2], border=clr.line, labels=NA, clockwise=FALSE, radius=1e6, init.angle=0, lwd=0.5) 
    par(old.par)

    # Hmz/Htz/N data
    text(mypie.cd[i,1], mypie.cd[i,2]-1e6*1.7, paste0("Hmz: ", myzyg.r.pop[1,i], "\nHtz: ", myzyg.r.pop[2,i], "\nN: ", mypop.nb[i]), adj=0.5, cex=0.8)

}

#legend(leg.cd, legend=c("Resistance allele", "Sensitive allele"), fill=clr.pie, bty="n", cex=0.75, yjust=0.5) # If full world map
legend("bottomright", legend=c("Resistance allele", "Sensitive allele"), fill=clr.pie[-2], bty="n", cex=1, yjust=0.5)

dev.off()
