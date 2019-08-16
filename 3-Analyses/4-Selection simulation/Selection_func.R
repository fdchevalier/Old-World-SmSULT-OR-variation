# Title: Selection_func.R
# Version: 1.0
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2018-12-14
# Modified in: 2019-05-12



#================#
# Data functions #
#================#

#------------#
# Simulation #
#------------#

# This section is a modified version of the dev.R file from driftR project.
## source: https://github.com/cjbattey/driftR/blob/fcf11160871b00710dc2ce5df5498e238e398179/dev.R

library(plyr)
library(reshape)
library(magrittr)
library(viridis)
library(parallel)
library(doParallel)

## This function has been extensively modified to allow parallelisation. Migration has been removed and final table simplified.
runPopSim2 <- function(gen=100,p=0.5,Waa=1,Wab=1,Wbb=1,n=100,nPop=2,Uab=0,Uba=0,infinitePop=F){

  registerDoParallel(cores = detectCores()-1)

  allele.freq <- foreach (j=1:nPop, .combine='cbind') %dopar% {
    p.tmp <- p 
      for(i in 1:gen){ 
        p <- (1 - Uab) * p + Uba * (1 - p) #mutation
        q <- 1-p
        if(p>0 && p<1){ #if alleles are not fixed
          w <- p*p*Waa+2*p*q*Wab+q*q*Wbb #population average fitness
          freq.aa <- (p*p*Waa)/w #post-selection genotype frequencies (weighted by relative fitness)
          freq.ab <- (2*p*q*Wab)/w
          if(infinitePop==F){ 
            Naa <- rbinom(1,n,freq.aa)
            if(freq.aa<1){ 
              Nab <- rbinom(1,(n-Naa),(freq.ab/(1-freq.aa)))
            }
            else {
              Nab <- 0
            }
            p <- ((2*Naa)+Nab)/(2*n)
            q <- 1-p
          } 
          else { #no drift (infinite population) conditions
            p <- freq.aa+(freq.ab/2)
            q <- 1-p
          }
        } else { #if alleles are fixed
          if(p<=0){
            p <- 0
            w <- p*p*Waa+2*p*q*Wab+q*q*Wbb
          } else {
            p <- 1
            w <- p*p*Waa+2*p*q*Wab+q*q*Wbb
          }
        }
        p.tmp[i + 1] <- p
      } #end populations loop
    return(p.tmp)
    } #end generations loop
    #summary stats
    names <- c()
    for(i in 1:nPop){names[i]<-paste0("p",i)}
    colnames(allele.freq) <- names
    return(allele.freq)
}


#-----------#
# Selection #
#-----------#

# This deterministic function computes genotype evolution based on initial allele frequency and genotype fitness
## source: Principles of population genetics, Hartl and Clark, 3rd editions, ISBN: 0-87893-306-9, p.219
diploid.selection <- function (p0 = 0.5, w = c(1, 1, 1), ngen = 400) {

    # Adjust generation number to include starting point
    ngen <- ngen+1

    # Genotype vector
    pp <- rep(NA, ngen)

    # Initial allele frequency
    p <- p0
    
    # Initial genotype frequency
    pp[1] <- p^2
    
    # Evolution of allele and genotype frequency
    for (g in 2:ngen) {
        q <- 1-p
        fitness.avg <- p^2*w[1] + 2*p*q*w[2] + q^2*w[3]

        pp[g] <- (p^2*w[1]) / fitness.avg

        p <- (p^2*w[1] + p*q*w[2]) / fitness.avg
    }

    return(pp)
}



#===================#
# Graphic functions #
#===================#

# Line in units
## source: https://stackoverflow.com/a/30835971
line2user <- function(line, side) {
    lh <- par('cin')[2] * par('cex') * par('lheight')
    x_off <- diff(grconvertX(c(0, lh), 'inches', 'npc'))
    y_off <- diff(grconvertY(c(0, lh), 'inches', 'npc'))
    switch(side,
        `1` = grconvertY(-line * y_off, 'npc', 'user'),
        `2` = grconvertX(-line * x_off, 'npc', 'user'),
        `3` = grconvertY(1 + line * y_off, 'npc', 'user'),
        `4` = grconvertX(1 + line * x_off, 'npc', 'user'),
        stop("Side must be 1, 2, 3, or 4", call.=FALSE))
}
