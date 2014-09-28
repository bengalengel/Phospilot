# some limma example code from Rafa

source("http://bioconductor.org/biocLite.R")
biocLite('SpikeInSubset')
library(SpikeInSubset)

data(rma95)
fac <- factor(rep(1:2, each = 3))##codes the grouping for the ttests


# We can now perform simple t-tests using the rowttests function in the genefilter package:
  
library(genefilter)

rtt <- rowttests(exprs(rma95), fac)#returns t stat, pvalue, and dm (difference in means 1st from the second)


# We will define colors depending on whether the p-value is small, the absolute difference in means is large, and whether the feature is a spike-in value.

mask <- with(rtt, abs(dm) < 0.2 & p.value < 0.01)## low p-value but very small difference in means
spike <- rownames(rma95) %in% colnames(pData(rma95))
cols <- ifelse(mask, "red", ifelse(spike, "dodgerblue", "black"))

# We now plot the results, using the colors defined above. We multiply the dm by -1, because we are interested in the difference from the second group to the first (this is the difference used by lm and the limma package by default). The spike-in genes are in blue, which have mostly small p-value and large difference in means. The red points indicate genes which have small p-values but also small differences in means. We will see how these points change after using limma.


#Volcano plot 
with(rtt, plot(-dm, -log10(p.value), cex = 0.8, pch = 16, xlim = c(-1, 1), ylim = c(0, 
                                                                                    5), xlab = "difference in means", col = cols))
abline(h = 2, v = c(-0.2, 0.2), lty = 2)

# Note that the red genes have mostly low estimates of standard deviation. Limma corrects for this by using empircal bayes to
# correct for what is unlikely to be (if it were sampled more) such a low variance. 

rtt$s <- apply(exprs(rma95), 1, function(row) sqrt(0.5 * (var(row[1:3]) + var(row[4:6]))))
with(rtt, plot(s, -log10(p.value), cex = 0.8, pch = 16, log = "x", xlab = "estimate of standard deviation", 
               col = cols))

#getyourself some limma

biocLite("limma")
library(limma)

fit <- lmFit(rma95, design = model.matrix(~fac))
colnames(coef(fit))

?lmFit

  
fit <- eBayes(fit)
tt <- topTable(fit, coef = 2)
tt

# topTable will return the top genes ranked by whichever value you define. You can also ask topTable to return all the values, sorted by "none". Note that a column automatically is included which gives the adjusted p-values for each gene. By default the method of Benjamini-Hochberg is used, by calling the p.adjust function.
# ?topTable

dim(topTable(fit, coef = 2, number = Inf, sort.by = "none"))

# Here we will compare the previous volcano plot with the limma results. Note that the red points are now all under the line where -log10(p.value) is equal to 2. Also, the blue points which represent real differences have p-values which are even higher than before.

limmares <- data.frame(dm = coef(fit)[, "fac2"], p.value = fit$p.value[, "fac2"])
with(limmares, plot(dm, -log10(p.value), cex = 0.8, pch = 16, col = cols, xlab = "difference in means", 
                    xlim = c(-1, 1), ylim = c(0, 5)))
abline(h = 2, v = c(-0.2, 0.2), lty = 2)

# Finally, we will construct a plot which shows how limma shrinks the variance estimates towards a common value, eliminating false positives which might arise from too-low estimates of variance.

# Here we pick, for each of 40 bins of different variance estimates, a single gene which falls in that bin. We remove bins which do not have any such genes.

n <- 40
qs <- seq(from = 0, to = 0.2, length = n)
idx <- sapply(seq_len(n), function(i) which(as.integer(cut(rtt$s^2, qs)) == 
                                              i)[1])
idx <- idx[!is.na(idx)]

# Now we will plot a line, from the initial estimate of variance for these genes to the estimate after running limma.

par(mar = c(5, 5, 2, 2))
plot(1, 1, xlim = c(0, 0.21), ylim = c(0, 1), type = "n", xlab = "variance estimates", 
     ylab = "", yaxt = "n")
axis(2, at = c(0.1, 0.9), c("before", "after"), las = 2)
segments((rtt$s^2)[idx], rep(0.1, n), fit$s2.post[idx], rep(0.9, n))




# Volcano plot
volcanoplot(fit,coef=2,highlight=2)

# Q-Q plot of moderated t-statistics
qqt(fit$t[,2],df=fit$df.residual+fit$df.prior)
abline(0,1)

y <- rt(50,df=4)
qqt(y,df=4)
abline(0,1)





