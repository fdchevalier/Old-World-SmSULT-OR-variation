#!/usr/bin/env Rscript
# Title: Selection.R
# Version: 1.0
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2018-12-14
# Modified in: 2019-05-12



#==========#
# Comments #
#==========#

# Aim: simulate natural selection on standing variation and de novo mutation



#===========#
# Functions #
#===========#

source("Selection_func.R")



#===========#
# Variables #
#===========#

# Directory variables
data_dir  <- "data"
graph_dir <- "graphs"
fig_file  <- "Fig. 7 - Simulation.pdf"

# Simulation variables
Ne     <- 6.5e4 # Crellen et al. 2016
nb.gen <- 1.5e3
nb.rep <- 2.5e3
s      <- c(0.01, 0.05, 0.1, 0.2)

s.sim  <- s[3]
p0.std <- 0.15 # Allele frequency observed in East Africa
p0.dnv <- 1/(2*Ne)

# Graph variable
std.clr <- "red"
dnv.clr <- "blue"

# Width of the graph in inch
mysize <- 7



#=================#
# Data processing #
#=================#

# Data folder
if (! dir.exists(data_dir)) {dir.create(data_dir)}

#------------------------------------------------#
# Simulation from standing and de novo variation #
#------------------------------------------------#

# Simulation for a selection coefficient of 0.1

## Standing variation
if (file.exists(paste(data_dir, "standing.Rdata", sep="/"))) {
    cat("Loadng simulation data of standing variation...\n")
    load(paste(data_dir, "standing.Rdata", sep="/"))
} else {
    cat("Running simulation for standing variation...\n")
    standing <- runPopSim2(gen=nb.gen, p=p0.std, Waa=1, Wab=1-s.sim, Wbb=1-s.sim, n=Ne, nPop=nb.rep)
    save(standing, file=paste(data_dir, "standing.Rdata", sep="/"))
}


## de novo variation
if (file.exists(paste(data_dir, "de.novo.Rdata", sep="/"))) {
    cat("Loading simulation data of de novo mutation...\n")
    load(paste(data_dir, "de.novo.Rdata", sep="/"))
} else {
    cat("Running simulation for de novo mutation...\n")
    run <- TRUE
    while (run) {
        de.novo <- runPopSim2(gen=nb.gen, p=p0.dnv, Waa=1, Wab=1-s.sim, Wbb=1-s.sim, n=Ne, nPop=nb.rep)
        # If an allele reaches fixtion, stop the simulation (should arrive at least after 5 simulations)
        if (sum(apply(de.novo, 2, function(x) any(x >=1))) > 0) { run <- FALSE }
    }
    save(de.novo, file=paste(data_dir, "de.novo.Rdata", sep="/"))
}

# Percentage of simulation to reach fixation for standing variation and de novo mutation
std.pc <- sum(apply(standing, 2, function(x) any(x == 1))) / nb.rep * 100
dnv.pc <- sum(apply(de.novo,  2, function(x) any(x == 1))) / nb.rep * 100

# Average generation time to reach fixation
std.tm <- mean(apply(standing, 2, function(x) which(x == 1)[1]), na.rm=TRUE)
dnv.tm <- mean(apply(de.novo,  2, function(x) which(x == 1)[1]), na.rm=TRUE)
cat("  - Average generation time to reach resistance fixation with: \n")
cat("    - standing variation: ", std.tm, "\n")
cat("    - de novo mutation: ",   dnv.tm, "\n")


#--------------------------------------#
# Modelization of resistance evolution #
#--------------------------------------#

cat("Modelization of allele resistance evolution...\n")

pheno.selection <- vector("list", length(s))

for (i in 1:length(s)) {
    fitness <- c(1,1-s[i],1-s[i])
    pheno.selection[[i]] <- diploid.selection(p0.std, fitness, ngen=nb.gen)
}

# Number of generation to reach 10 and 20% respectively
gen.ps.1 <- lapply(pheno.selection, function(x) which(x >= 0.1)[1])
gen.ps.2 <- lapply(pheno.selection, function(x) which(x >= 0.2)[1])



#=========#
# Figures #
#=========#

cat("Generating graph...\n")

# Folder creation for graphic output (if does not exist)
if (! dir.exists(graph_dir)) {dir.create(graph_dir)}

# Start PDF device
pdf(paste(graph_dir,fig_file,sep="/"), height=mysize*2, width=mysize)

layout(matrix(1:2, ncol=1, byrow=TRUE))
par(mar=c(5,4,0.5,0.5)+0.1)


#-------------#
# Simmulation #
#-------------#

# Max generation to reach 0.5 for standing variation
std.x <- max(apply(standing, 2, function(x) which(x >= 0.5)[1]), na.rm=TRUE)

# Max generation to reach 0.5 for de.novo variation
dnv.x <- max(apply(de.novo, 2, function(x) which(x >= 0.5)[1]), na.rm=TRUE)

# Plot
matplot(standing, ylab="p(OXA-R)", xlab="Generation", ylim=c(0,1), type="l", col=std.clr, lty=1)
matplot(de.novo, type="l", col=dnv.clr, lty=1, add=TRUE) 

# Panel letter
mtext("A", side=3, line=-1, at=line2user(par("mar")[2],2), cex=par("cex")*2, adj=0)

# Text
text(std.x, 0.6, paste0(std.pc, "% of simulations reached fixation"), pos=4, col=std.clr)
text(dnv.x, 0.4, paste0(dnv.pc, "% of simulations reached fixation"), pos=2, col=dnv.clr)

# Legend
legend(nb.gen/2, 0.8, legend=c("Standing variation", expression(italic("de novo ") * "variation")), lty=1, col=c(std.clr, dnv.clr), bty="n", xjust=0.5, yjust=0.5, xpd=TRUE)


#--------------#
# Modelisation #
#--------------#

myclr <- gray.colors(length(s), start=0, end=0.8)

plot(pheno.selection[[1]], ylab="Frequency of OXA-R parasite", xlab="Generation", ylim=c(0,0.3), type="l", col=myclr[1], log="x")
nulo <- mapply(function(x,y) lines(x, col=y), pheno.selection[2:length(s)], myclr[2:length(s)])

# Threshold
abline(h=0.1, col="grey", lty=2)
abline(h=0.2, col="grey", lty=2)

# Number of generations
text(gen.ps.1, 0.1, gen.ps.1, adj=c(0,1.5), col=myclr)
text(gen.ps.2, 0.2, gen.ps.2, adj=c(0,1.5), col=myclr)

# Panel letter
mtext("B", side=3, line=-1, at=line2user(par("mar")[2],2), cex=par("cex")*2, adj=0)

# Legend
legend("topleft", legend=rev(s), lty=1, col=rev(myclr), bty="n", xjust=0.5, yjust=0.5, title="Selection coefficients:", xpd=TRUE)

# close PDF device
dev.off()

