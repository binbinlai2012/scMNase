args <- commandArgs( trailingOnly = TRUE )

tb1 <- scan(args[1], list(x=0,y=0))
tb2 <- scan(args[2], list(x=0,y=0))
tb3 <- scan(args[3], list(x=0,y=0))
tb4 <- scan(args[4], list(x=0,y=0))

color1 = "blue"
color2 = "darkgreen"
color3 = "orange"
color4 = "red"
lwidth=2
xrange=c(-500,500)
yrange=c(0,3)
png(args[5])
par(mar=c(6,6,2,2)+0.01)
plot(tb1$x*10-500, tb1$y, type='l', col=color1, xlim=xrange, ylim=yrange,  ann='F', yaxt="n", xaxt="n", lwd = lwidth, bty="n" )
par(new=T)
plot(tb2$x*10-500, tb2$y, type='l', col=color2, xlim=xrange, ylim=yrange,  ann='F', yaxt="n", xaxt="n", lwd = lwidth, bty="n" )
par(new=T)
plot(tb3$x*10-500, tb3$y, type='l', col=color3, xlim=xrange, ylim=yrange,  ann='F', yaxt="n", xaxt="n", lwd = lwidth, bty="n" )
par(new=T)
plot(tb4$x*10-500, tb4$y, type='l', col=color4, xlim=xrange, ylim=yrange,  ann='F', yaxt="n", xaxt="n", lwd = lwidth, bty="n" )
y=c(0, 1, 2, 3)
axis(2, at=y, labels=y, las=2, lwd=lwidth, cex.axis=2 )
x=c(  -500, 0, 500 )

axis(1, at=x, labels=x, las=1, lwd=lwidth, cex.axis=2 )

title(xlab="Distance from TSS", ylab="Normalized subnucl. counts", cex.lab=2)
dev.off()
