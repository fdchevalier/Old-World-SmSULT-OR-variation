#!/usr/bin/env Rscript
# Title:
# Version: 
# Author: Frédéric CHEVALIER
# Created in:
# Modified in:



#==========#
# Comments #
#==========#





#===================#
# Packages required #
#===================#

library(plotrix)	# for ablineclip function
library(TeachingDemos) # for shadowtext function



#===========#
# Functions #
#===========#
#------------------------------------------------#
# Fuction to test if the number is a wholenumber #
#------------------------------------------------#

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol


#-----------------------------------------#
# Function to plot the SNP by chromosomes #
#-----------------------------------------#

matplot.data <- function (data.tab, col, ylim=2000, cex.axis=1, cex=1) { #ylim=10*trunc(log10(mean(data.tab[,col]))*10^trunc(log10(mean(data.tab[,col])))), cex.axis=1, cex=1) {

	# Smoothing curve using runmed
	rmdt <- data.tab[,col]

	# Chromosome parameters use for search and plotting
	chr.names <- c("i.Chr_1", "i.Chr_2", "i.Chr_3", "i.Chr_4", "i.Chr_5", "i.Chr_6", "i.Chr_7", "i.Chr_W", "i.SC")
	chr.names.true <- c("1", "2", "3", "4", "5", "6", "7", "Z", "Unass. scaff.")
	chr.names.all <- unique(as.vector(data.tab[,1]))

	chr.SNP.length <- c(0)
	for (i in 1:length(chr.names)){
		chr.SNP.length <- c(chr.SNP.length, tail(grep(chr.names[i],data.tab[,1]),1))
	}


	# Marge modification (change it if you need to add title or other things around the graph)
	par(mar=c(3,4,4,0.5)+0.1)
#	par(mar=c(4,4,4,0.5)+0.1)	# original
#	par(mar=c(4,4.5,2,0.5)+0.1)	# for cex=2
#	par(mar=c(4,5,2,0.5)+0.1)	# for cex=4
#	par(mar=c(2,2,1,0)+0.1)

	# Empty plot
	plot(1,1, ylim=c(min(rmdt), ylim), xlim=c(1, ceiling(nrow(data.tab)/1000)*1000), type="n", bty="n", axes=FALSE, xlab="Chromosomes", ylab="Read depth", cex.lab=cex)

	# Axes draw
	axis(1, at=chr.SNP.length, labels=FALSE)
	axis(2, cex.axis=cex.axis)

	# X-axis label
	for (i in 1:length(chr.names)) {
#		text(chr.SNP.length[i]+(chr.SNP.length[i+1]-chr.SNP.length[i])/2, -0.1*ylim, labels=chr.names.true[i], cex=cex.axis, xpd=TRUE)
		text(chr.SNP.length[i]+(chr.SNP.length[i+1]-chr.SNP.length[i])/2, -0.075*ylim, labels=chr.names.true[i], cex=cex.axis*0.9, xpd=TRUE)
	}


	# Plot of the SNP frequency on the "chromosomes" and the scaffold
	for (i in 1:(length(chr.names))) {			# Grey Rectangle
		if (is.wholenumber(i/2) == FALSE){
			rect(head(grep(chr.names[i],data.tab[,1]),1), min(trunc(rmdt)), tail(grep(chr.names[i],data.tab[,1]),1), max(ceiling(rmdt)/2), col="grey90", border=NA)
#			rect(head(grep(chr.names[i],data.tab[,1]),1), min(trunc(rmdt)), tail(grep(chr.names[i],data.tab[,1]),1), ylim, col="grey90", border=NA)
		}
	}

	for (i in 1:length(chr.names)) {			# Plotting the data

		print(paste(i, "-", chr.names[i]), quote=FALSE)

		for (j in head(grep(chr.names[i],chr.names.all),2)) {
			if (j == head(grep(chr.names[i],chr.names.all), 1)) {
				length.plot <- seq(head(grep(paste(chr.names.all[j], "$", sep=""), data.tab[,1]), 1), tail(grep(paste(chr.names.all[j], "$", sep=""), data.tab[,1]), 1))
			} else {
				length.plot <- seq(head(grep(paste(chr.names.all[j], "$", sep=""), data.tab[,1]), 1), tail(grep(paste(chr.names.all[tail(grep(chr.names[i],chr.names.all),1)], "$", sep=""), data.tab[,1]), 1))
			}

			if (i == length(chr.names)) {
				color="black"
			} else if (j == head(grep(chr.names[i],chr.names.all),1)) {
				color="red"
			} else {
				color="black"
			}

		matplot(length.plot, rmdt[length.plot], col=color, type="l", add=TRUE)
			}
		}

}  #End of the function



#===========#
# Variables #
#===========#

# Removing of all existing variables from the memory
#rm(list = ls())

# Data file
mydir.list <- list.dirs("data", recursive=F)

myfile.list <- list.files("data", pattern="baits-0.bed", full.names=T, recursive=T)

mytitle.list <- lapply(strsplit(mydir.list, "/"), function(x) x[[length(x)]])

ylim.list <- c(6000,1500,600,1500,1200,1500)



#=========#
# Figures #
#=========#

# Create the graph folder if it do not already exist
if (file.exists("Graphs/") == FALSE) {dir.create("Graphs")}


for (i in 1:length(myfile.list)) {
	myfile <- myfile.list[i]
    mytitle <- as.character(mytitle.list[i])
	print(mytitle)

    mydata <- read.csv(myfile, header=FALSE, sep="\t")

	mydata.mean <- mean(mydata[,4])
	p <- 10^trunc(log10(mydata.mean))
	y.lim <- round(mydata.mean/p)*(p*15)
#	y.lim <- ylim.list[i]

	par(bg="white")
	mycex=1

    png(file=paste("Graphs/",mytitle,".png",sep=""), width=16*72, height=10*72)

    matplot.data(mydata,4, cex.axis=mycex, cex=mycex,  ylim=y.lim)

	title(mytitle, cex=mycex)

    dev.off()

}
