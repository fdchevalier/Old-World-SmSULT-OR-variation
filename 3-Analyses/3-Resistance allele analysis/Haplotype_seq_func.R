#!/usr/bin/env Rscript
# Title: Haplotype_seq_func.R
# Version: 1.0
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2018-12-14
# Modified in: 2019-05-08



#===================#
# Graphic functions #
#===================#

#-----------#
# Converter #
#-----------#

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

#-------------------------#
# Variants to chromosomes #
#-------------------------#

mygraph <- function(x, vec, pos) {
    # x     data frame
    # vec   vector of logicals that contains samples with deletion
    # pos   vector of the GT position in bp

    mygt.tmp <- NULL
    for (i in seq(1,sum(vec),2)) {
        myrow    <- rownames(x)[vec][i:(i+1)]
        mygt.tmp <- rbind(mygt.tmp, x[ myrow,])
        if(i+1 < sum(vec)) { mygt.tmp <- rbind(mygt.tmp, rep(NA, ncol(x)))}
    }

    myrow <- nrow(mygt.tmp)
    mycln <- ncol(mygt.tmp)
    
    # Make squares squared
    myf <- 1/myrow*2

    # Adequate left margin
    mar.l <- round(log10(mycln))
    if (mar.l == 1) { mar.l <- 2 }
    
    # Margins
    par(mar=c(5,4*mar.l,0.5,0.5))

    myrow <- nrow(mygt.tmp)
    mycol <- ncol(mygt.tmp)
    plot(0,0, xlim=c(0,mycln), ylim=c(0.5,myrow), xlab="", ylab="", type="n", axes=F, xaxs="i")
    for (r in 1:myrow) {
        
        if(all(is.na(mygt.tmp[r,]))) { bdr <- "NA" } else { bdr <- par("fg") }

        for (c in 1:mycln) {
            rect(c-0.5, r-0.5, c+0.5, r+0.5, col=mygt.tmp[r,c], xpd=TRUE, border=bdr)
        }
    }

    mylb.u <- rep(NA, length(rownames(mygt.tmp)))

    for (i in rownames(mygt.tmp)) {
        mylb <- strsplit(i, ".", fixed=T) %>% unlist()
        if(length(mylb) == 0) next
        mylb.u[match(i,rownames(mygt.tmp))] <- paste(mylb[1:(length(mylb)-1)], collapse=".")
    }

    mylb.u <- unique(na.omit(mylb.u))

    axis(1, at=which(! is.na(pos)), labels=na.omit(pos), las=2)
    title(xlab="Relative position to p.E142del (bp)") 
    axis(2, at=seq(1.5,myrow,3), pos=1, labels=mylb.u, las=1, tick=FALSE)

}

