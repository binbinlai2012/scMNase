args <- commandArgs( trailingOnly = TRUE )

tb1 <- scan(args[1], list(x=0,y=0))

png(args[2])
par(mar=c(6,6,2,2)+0.01)
plot(tb1$x*10-1000, tb1$y, type='l', col=color1,  lwd = lwidth, bty="n" )
title(xlab="Distance from TSS", ylab="tag density", cex.lab=2)
dev.off()
