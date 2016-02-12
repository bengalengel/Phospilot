#PCA and prince plots of confounded data before and after PC1 regression
require(swamp)
require(limma)
require(sva)

#start with the median normalized and confounded data
str(medianSub.quantiled)
cdata <- na.omit(medianSub.quantiled)
cdata <- as.matrix(cdata)

# set font family within par
op <- par(family = "serif")

#PCA analysis 
x <- t(cdata)#samples are the rows of the column matrix
pc <- prcomp(x, scale = T, center = T) #I am now scaling and centering x
summary(pc)



cols <- as.factor(substr(colnames(cdata), 3, 7))##use 5 digit exp name.
pdf("PCA_confounded_nocorrection.pdf")
plot(pc$x[, 1], pc$x[, 2], col = as.numeric(cols), xlab = "PC1", ylab = "PC2", pch = rep(c(1,1,6,6),3), 
     cex = 1.5, family = "serif")
legend("bottomright", c(levels(cols), "batch1", "batch2") , col = c(seq(along=levels(cols)), 1, 1), lty = c(1,1,1,NA,NA), pch = c(rep(NA,3), 1, 6), pt.cex = 1.5, bty = "n")
dev.off()

#SVD for calculating variance explained;
cx <- sweep(x, 2, colMeans(x), "-")
sv <- svd(cx)
names(sv)
plot(sv$u[, 1], sv$u[, 2], col = as.numeric(cols), main = "SVD", xlab = "U1", ylab = "U2")  
plot(sv$d^2/sum(sv$d^2), xlim = c(1, 12), type = "b", pch = 16, xlab = "principal components", 
     ylab = "variance explained")


##use the swamp package to regress away PC1

# here I will use the data with at least one observation in each biological replicate
medianSub.quantiled.bio <- medianSub.quantiled[rowSums(is.na(medianSub.quantiled[ , 1:2])) < 2 & rowSums(is.na(medianSub.quantiled[ , 3:4])) < 2 & rowSums(is.na(medianSub.quantiled[ , 5:6])) < 2 
     & rowSums(is.na(medianSub.quantiled[ , 7:8])) < 2 & rowSums(is.na(medianSub.quantiled[ , 9:10])) < 2 & rowSums(is.na(medianSub.quantiled[ , 11:12])) < 2,] 

medianSub.quantiled.bio <- as.matrix(medianSub.quantiled.bio)

PC1.regressed <- kill.pc(medianSub.quantiled.bio, 1, imputeknn = T)
x <- t(PC1.regressed)
pc <- prcomp(x, scale = T, center = T) #I am now scaling and centering x

cols <- as.factor(substr(colnames(cdata), 3, 7))##use 5 digit exp name.
pdf("PCA_confounded_PC1_regressed.pdf")
plot(pc$x[, 1], pc$x[, 2], col = as.numeric(cols), xlab = "PC1", ylab = "PC2", pch = rep(c(1,1,6,6),3), 
     cex = 1.5, family = "serif")
legend("bottomleft", c(levels(cols), "batch1", "batch2") , col = c(seq(along=levels(cols)), 1, 1), lty = c(1,1,1,NA,NA), pch = c(rep(NA,3), 1, 6), pt.cex = 1.5, bty = "n")
dev.off()


##use combat estimated and regressed batch effect
com2 <- na.omit(com2)
com2 <- as.matrix(com2)
x <- t(com2)
pc <- prcomp(x, scale = T, center = T) #I am now scaling and centering x

cols <- as.factor(substr(colnames(cdata), 3, 7))##use 5 digit exp name.
pdf("PCA_confounded_combatcorrected.pdf")
plot(pc$x[, 1], pc$x[, 2], col = as.numeric(cols), xlab = "PC1", ylab = "PC2", pch = rep(c(1,1,6,6),3), 
     cex = 1.5, family = "serif")
legend("bottomleft", c(levels(cols), "batch1", "batch2") , col = c(seq(along=levels(cols)), 1, 1), lty = c(1,1,1,NA,NA), pch = c(rep(NA,3), 1, 6), pt.cex = 1.5, bty = "n")
dev.off()



#reset graphics device
par(op)



