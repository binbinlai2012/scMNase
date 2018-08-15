args <- commandArgs( trailingOnly = TRUE )

tb1 <- scan(args[1], list(x=0,y=0))


rangex=c(0, 250)

png(args[2]);
#par(mar=c(5,5,4,2)+0.1);
lwidth=2
plot(tb1$x, tb1$y, type='l',  xlim=rangex, lwd=lwidth )


box(lwd=2)
dev.off()

