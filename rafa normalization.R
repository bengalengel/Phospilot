library(rafalib)
library(affy)
library(affydata)
library(SpikeIn)


data(SpikeIn95)
spms <- pm(SpikeIn95)#perfect match probes
spd = pData(SpikeIn95)#phenotypic data of somekind (this is the spike in data)


mypar(1, 2)
shist(log2(spms[, 2]), unit = 0.1, type = "n", xlab = "log (base 2) intensity", 
      main = "Five techical replicates")
## show the first 10
for (i in 1:10) shist(log2(spms[, i]), unit = 0.1, col = i, add = TRUE, lwd = 2, 
                      lty = i)
## Pick 9 and 10 and make an MA plot (10 is the plot with a different density
i = 10
j = 9  ##example with two samples
siNames <- colnames(spd)
## show probes with expected FC=2
siNames = siNames[which(spd[i, ]/spd[j, ] == 2)]
M = log2(spms[, i]) - log2(spms[, j])
A = (log2(spms[, i]) + log2(spms[, j]))/2
splot(A, M, ylim = c(-1.5, 1.5))
spikeinIndex = which(probeNames(SpikeIn95) %in% siNames)
points(A[spikeinIndex], M[spikeinIndex], ylim = c(-4, 4), bg = 1, pch = 21)


# In this plot the green dots are genes spiked in to be different in the two samples. The rest of the points (black dots) should be at 0 because other than the spiked-in genes these are technical replicates. Notice that without normalization the black dots on the left side of the plot are as high as most of the green dots. If we were to order by fold change, we would obtain many false positives. In the next units we will introduce three normalization procedures that have proven to work well in practice.


# Loess normalization
# 
# In the MA-plot above we see a non-linear bias in the M that changes as function of A. The general idea behind loess normalization is to estimate this bias and remove it. Because the bias is a curve of no obvious parametric form (it is not a line or parabola or a sine function, etc.) we want to fit a curve to the data. Local weighted regression (loess) provides one way to do this. Loess is inspired by Taylor's theorem that in practice means that at any given point, if one looks at a small enough region around that point, the curve looks like a parabola. If you look even closer it looks like a straight line (note that gardeners can make a curved edge with a straight shovel).
# 
# Loess takes advantage of this mathematical property of functions. For each point in your data set a region is defined considered to be small enough to assume the curve approximated by a line in that region and a line is fit with weighted least squares. The weights depend on the distance from the point of interest. The robust version of loess also weights points down that are considered outliers. The following code makes an animation that shows loess at work:
mypar(1,1)
o <- order(A)
a <- A[o]
m <- M[o]
ind <- round(seq(1, length(a), len = 5000))
a <- a[ind]
m <- m[ind]
centers <- seq(min(a), max(a), 0.1)
plot(a, m, ylim = c(-1.5, 1.5), col = "grey")

windowSize <- 1.5
smooth <- rep(NA, length(centers))
# library(animation) saveGIF({ for(i in seq(along=centers)){
# center<-centers[i] ind=which(a>center-windowSize & a<center+windowSize)
# fit<-lm(m~a,subset=ind)
# smooth[i]<-predict(fit,newdata=data.frame(a=center)) if(center<12){
# plot(a,m,ylim=c(-1.5,1.5),col='grey') points(a[ind],m[ind])
# abline(fit,col=3,lty=2,lwd=2) lines(centers[1:i],smooth[1:i],col=2,lwd=2)
# points(centers[i],smooth[i],col=2,pch=16) } } },'loess.gif', interval =
# .15) Final version
plot(a, m, ylim = c(-1.5, 1.5))
lines(centers, smooth, col = 2, lwd = 2)


o <- order(A)
a <- A[o]
m <- M[o]
ind <- round(seq(1, length(a), len = 5000))
a <- a[ind]
m <- m[ind]
fit <- loess(m ~ a)
bias <- predict(fit, newdata = data.frame(a = A))
nM <- M - bias
mypar(1, 1)
splot(A, M, ylim = c(-1.5, 1.5))
points(A[spikeinIndex], M[spikeinIndex], ylim = c(-4, 4), bg = 1, pch = 21)
lines(a, fit$fitted, col = 2, lwd = 2)


splot(A, nM, ylim = c(-1.5, 1.5))
points(A[spikeinIndex], nM[spikeinIndex], ylim = c(-4, 4), bg = 1, pch = 21)
abline(h = 0, col = 2, lwd = 2)

# Note that the bias is removed and now the highest fold changes are almost all spike-ins. Also note that the we can control the size of the intervals in which lines are fit. The smaller we make these intervals the more flexibility we get. This is controlled with the span argument of the loess function. A span of 0.75 means that the closest points are considered until 3/4 of all points are used. Finally, we are fitting parabolas, but for some datasets these can result in over fitting. For example, a few points can force a parabola to "shoot up" very fast. For this reason using lines (the argument is degree=1) is safer.

fit <- loess(m ~ a, degree = 1, span = 1/2)
bias <- predict(fit, newdata = data.frame(a = A))
nM <- M - bias
mypar(1, 1)
splot(A, M, ylim = c(-1.5, 1.5))
points(A[spikeinIndex], M[spikeinIndex], ylim = c(-4, 4), bg = 1, pch = 21)
lines(a, fit$fitted, col = 2, lwd = 2)


splot(A, nM, ylim = c(-1.5, 1.5))
points(A[spikeinIndex], nM[spikeinIndex], ylim = c(-4, 4), bg = 1, pch = 21)
abline(h = 0, col = 2, lwd = 2)

# ******************************************************************************************
# Does loess normalization produce better histograms and boxplots?!!! YES!

# First I need to convert back to normalized log2 values using some algebra
inorm <- 2^((nM+2*A)/2)
jnorm <- 2^-((nM-2*A)/2)

inorm <- as.vector(inorm)
jnorm <- as.vector(jnorm)

test <- cbind(inorm,jnorm)
test2 <- cbind(inorm,jnorm,spms[, i],spms[,j])

test3 <- log2(test2)
boxplot(test3)
summary(test3)


plot.new()
par(mfrow = c(1, 1))
for (i in 1:(ncol(test3))){
  if(i==1) plot(density(test3[, i], na.rm=T), col = i, ylim = c(0,.8))
  else lines(density(test3[, i], na.rm=T), col = i)
}

# Yes it certainly does!! green and orange are the normalized histograms. This certainly also works with quantile norm.

# ******************************************************************************************
library(preprocessCore)
#quantile normalized 
nspms <- normalize.quantiles(spms)

# before norm
M = log2(spms[, i]) - log2(spms[, j])
A = (log2(spms[, i]) + log2(spms[, j]))/2
splot(A, M, ylim = c(-1.5, 1.5))
points(A[spikeinIndex], M[spikeinIndex], bg = 1, pch = 21)

# after norm
M = log2(nspms[, i]) - log2(nspms[, j])
A = (log2(nspms[, i]) + log2(nspms[, j]))/2
splot(A, M, ylim = c(-1.5, 1.5))
points(A[spikeinIndex], M[spikeinIndex], bg = 1, pch = 21)

# Note that the densities are now identical as expected since we forced this to be the case.

pms <- spms
mypar(1, 1)
shist(log2(pms[, 2]), unit = 0.1, type = "n", xlab = "log (base 2) intensity", 
      main = "Five techical replicates")
for (i in 1:5) shist(log2(pms[, i]), unit = 0.1, col = i, add = TRUE, lwd = 2, 
                     lty = i)

qpms <- normalize.quantiles(pms[, 1:5])
shist(log2(qpms[, 2]), unit = 0.1, type = "n", xlab = "log (base 2) intensity", 
      main = "Five techical replicates")
for (i in 1:5) shist(log2(qpms[, i]), unit = 0.1, col = i, add = TRUE, lwd = 2, 
                     lty = i)

# variance stabilizing normalization section


library(rafalib)
N = 10000
e = rexp(N, 1/1000)
b1 = 24
b2 = 20
A1 = 1
A2 = 1.25
sigma = 1
eta = 0.05
y1 = b1 + rnorm(N, 0, sigma) + A1 * e * 2^rnorm(N, 0, eta)
y2 = b2 + rnorm(N, 0, sigma) + A2 * e * 2^rnorm(N, 0, eta)
mypar(1, 1)
maplot(log2(y1), log2(y2), ylim = c(-1, 1), curve.add = FALSE)


# For this type of data, the variance depends on the mean. We seek a transfromation that stabilizies the variance of the estimates of $\theta$ after we subctract the additive background estimate and divide by the estimate of the gain.

ny1 = (y1 - b1)/A1
ny2 = (y2 - b2)/A2
mypar(1, 2)
maplot(ny1, ny2, curve.add = FALSE, ylim = c(-500, 500))
maplot(log2(ny1), log2(ny2), ylim = c(-2, 2), xlim = c(0, 15))


# if we know how the variance depends on the mean, we can compute a variance stabilizing transform: 

# $$ Y \text{ with } \text{E}(Y)=\mu \text{ and } \text{var}(Y) = v(\mu)\ \text{var}{f(Y)} \text{ does not depend on } \mu $$ In the case of the model above we can derive the following transformation
# 
# $$ \text{arsinh}(y) = \log\left(y + \sqrt{y^2+1} \right) $$
#   
# The vsn library implements this apprach. It estimates $\beta$ and $A$ by assuming that most genes don't change, i.e. $\theta$ does not depend on $i$.


library(vsn)
nspms <- exprs(vsn2(spms))
mypar(1,1)
i = 10
j = 9
M = log2(spms[, i]) - log2(spms[, j])
A = (log2(spms[, i]) + log2(spms[, j]))/2
splot(A, M, ylim = c(-1.5, 1.5))
points(A[spikeinIndex], M[spikeinIndex], bg = 1, pch = 21)


M = nspms[, i] - nspms[, j]
A = (nspms[, i] + nspms[, j])/2
splot(A, M, ylim = c(-1.5, 1.5))
points(A[spikeinIndex], M[spikeinIndex], bg = 1, pch = 21)

# We notice that it corrects the bias in a similar way to loess and quantile normalization.

