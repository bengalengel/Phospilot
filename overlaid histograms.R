# overlaid histograms of data that has been normalized and quantified (class 1 for phospho) in at least one sample


phoshis1
phoshis2
phoshis3
phoshis4
phoshis5
phoshis6

proth1
proth2
proth3
proth4
proth5
proth6

par(mfrow = c(3,2))

plot(phoshis1, col=rgb(0,0,1,1/4), xlim=c(-6,6))
plot(proth1, col=rgb(1,0,0,1/4), xlim=c(-6,6), add=T)  # second

plot(phoshis2, col=rgb(0,0,1,1/4), xlim=c(-6,6))
plot(proth2, col=rgb(1,0,0,1/4), xlim=c(-6,6), add=T)  # second

plot(phoshis3, col=rgb(0,0,1,1/4), xlim=c(-6,6))
plot(proth3, col=rgb(1,0,0,1/4), xlim=c(-6,6), add=T)  # second

plot(phoshis4, col=rgb(0,0,1,1/4), xlim=c(-6,6))
plot(proth4, col=rgb(1,0,0,1/4), xlim=c(-6,6), add=T)  # second

plot(phoshis5, col=rgb(0,0,1,1/4), xlim=c(-6,6))
plot(proth5, col=rgb(1,0,0,1/4), xlim=c(-6,6), add=T)  # second

plot(phoshis6, col=rgb(0,0,1,1/4), xlim=c(-6,6))
plot(proth6, col=rgb(1,0,0,1/4), xlim=c(-6,6), add=T)  # second

par(mfrow = c(1,1))
plot(phoshis6, col=rgb(0,0,1,1/4), xlim=c(-6,6))
plot(proth6, col=rgb(1,0,0,1/4), xlim=c(-6,6), add=T)  # second



p11 <- datanorm$HL16770_1
p22 <- datanorm$HL16770_1
datanorm$HL16770_1
datanorm$HL16770_1
datanorm$HL16770_1
datanorm$HL16770_1

###plot each with lines intersecting the 1st and 3rd quantiles
qqnorm(datanorm$HL16770_1)
qqline(datanorm$HL16770_1)
GRB2pts <- qqnorm(GRB2, plot.it=FALSE)
points(GRB2pts, col=2)
qqline(GRB2, col=2)
SHP2Npts <- qqnorm(SHP2N, plot.it=FALSE)
points(SHP2Npts, col=3)
qqline(SHP2N, col=3)
NCK1pts <- qqnorm(NCK1, plot.it=FALSE)
points(NCK1pts, col=4)
qqline(NCK1, col=4)




##set the range for the qqplots

ylim=range(ABL1,SHP2N,GRB2,NCK1)


###plot each with lines intersecting the 1st and 3rd quantiles
qqnorm(ABL1, ylim=ylim)
qqline(ABL1)
GRB2pts <- qqnorm(GRB2, plot.it=FALSE)
points(GRB2pts, col=2)
qqline(GRB2, col=2)
SHP2Npts <- qqnorm(SHP2N, plot.it=FALSE)
points(SHP2Npts, col=3)
qqline(SHP2N, col=3)
NCK1pts <- qqnorm(NCK1, plot.it=FALSE)
points(NCK1pts, col=4)
qqline(NCK1, col=4)
