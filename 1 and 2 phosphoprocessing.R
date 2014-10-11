# Merge the two files for clustering
library(plyr)

data1 <- read.csv("data1.csv")
data2 <- read.csv("data2.csv")
names(data1)[7] <- "SequenceWindow"
names(data2)[7] <- "SequenceWindow"
names(data1)[8] <- "Multiplicity"
names(data2)[8] <- "Multiplicity"


#Merging the data using the sequence window and multiplicity
combined <- merge(data1,data2, all=T)

#remove non data columns
drops <- c("SequenceWindow","Multiplicity")
combined <- combined[,!(names(combined) %in% drops)]

#boxplots 
boxplot(combined)

#median normalize 
median.subtract <- function(x){ x - median(x, na.rm = TRUE)}##create a wrapper for median subtraction
combined1 <- colwise(median.subtract)(combined) #create median subtracted data but loose the row names here...

# boxplots after median normalization
boxplot(combined1)

#remove all uncommon
combinedCommon <- na.omit(combined1)

#boxplots again
boxplot(combinedCommon)

# PCA based examination of batch effects using swamp package
library(swamp)

data <- as.matrix(combinedCommon)

##### sample annotations (data.frame)
set.seed(50)
o<-data.frame(Factor1=factor(c(rep("A",6),rep("B",6))),
              Numeric1=rnorm(12),row.names=colnames(data))

# PCA analysis
res1<-prince(data,o,top=10,permute=T)
str(res1)
res1$linp#plot p values
res1$linpperm#plot p values for permuted data
prince.plot(prince=res1)

# to plot p values of linear models: lm(principal components ~ sapmle annotations).
# to see if the variation in the data is associated with sample annotations.
res2<-prince.var.plot(data,show.top=12,npermute=10)

# to see how many principal components carry informative variation
## hierarchical clustering analysis
hca.plot(data,o)

# to show a dendrogram with sample annotations below
res3<-hca.test(data,o,dendcut=2,test="fisher")

# to test if the major clusters show differences in sample annotations
## feature associations
res4a<-feature.assoc(data,o$Factor1,method="correlation")

# to calculate correlation between one sample annotation and each feature
res4b<-feature.assoc(data,o$Factor1,method="t.test",g1=res4a$permuted.data)
res4c<-feature.assoc(data,o$Factor1,method="AUC",g1=res4a$permuted.data)

dense.plot(res4a)
# to plot the distribution of correlations in the observed data

# in comparison to permuted data
dense.plot(res4b)
dense.plot(res4c)

res5<-corrected.p(res4a)
# to correct for multiple testing and find out how many features are
# significantly associated to the sample annotation
names(which(res5$padjust<0.05))
names(which(res5$adjust.permute<0.05))
names(which(res5$adjust.rank<0.05))

## associations between sample annotations
res4<-confounding(o,method="fisher")
# to see how biological and technical annotations are inter-related

## adjustment for batch effects
gadj1<-quickadjust.zero(data,o$Factor1)
# to adjust for batches (o$Factor1)
# using median centering of the features for each batch
prince.plot(prince(gadj1,o,top=10))

gadj2<-quickadjust.ref(data,o$Factor1,"B")
# to adjust for batches (o$Factor)
# by adjusting the median of the features to the median of a reference batch.
prince.plot(prince(gadj2$adjusted.data,o,top=10))


gadj3<-kill.pc(data,pc=1)
# to remove one or more principal components (here pc1) from the data
prince.plot(prince(gadj3,o,top=10))




mypar(2,1)
boxplot(gadj3)
boxplot(combinedCommon)
mypar(1,1)





# here we go!
dataZ <- scale(combinedCommon)##Z-scored column wise
dataZ <- scale(gadj3)##Z-scored column wise

# now all data excepting complete cases (note that the sample dendograms look the same)
# hist(dataZ[,6], breaks = 100)

# dendogram using euclidian distance (default) and ward or complete agglomeration
dend.ward<- as.dendrogram(hclust(dist(t(dataZ)),method="ward"))
dend.complete<- as.dendrogram(hclust(dist(t(dataZ))))

ward.o<- order.dendrogram(dend.ward)
complete.o<- order.dendrogram(dend.complete)

plot(dend.complete,ylab="height", main = "Euclidian/Complete")
plot(dend.ward, leaflab = "perpendicular", ylab = "height", main = "Euclidian/Ward")

dev.off()##turns off plot

# Cluster using euclidian distance and ward linkage for both sites(rows) and samples (columns)
# Note that both dendograms are created independently and row Z scores are presented in the heatmap

# row scaled
r <- t(scale(t(combinedCommon)))#transpose to zscale the rows then transpose back to original format

# sample scaled
c <- scale(combinedCommon)


# install heatmap.2 package
install.packages("gplots")
library(gplots)

# Create dendrogram using the data without NAs
feature.dend<- as.dendrogram(hclust(dist(r),method="ward"))
sample.dend<- as.dendrogram(hclust(dist(t(c)),method="ward"))##note that dist caclculates distance between rows by default


##produce the heatmap. Note that the help page has a nice section on identifying subregions by color. Although I will likely have to cut the dendogram to id clusters of interest

heatmap.2(
  r,#row Z scores
  Colv=sample.dend,
  Rowv=feature.dend,
  col=bluered(25),
  scale="none",
  trace="none",
  key.xlab = "Row Z scores", key.ylab=NULL, key.title = "",
  srtCol=45,  ,adjCol = c(1,1),
  margins = c(6,5),
  cexCol=1,
  labRow = NA#remove row labels
)


dev.off()


###########################################################################################################################

##PCA of data from quick R (I will need to use the column standardized data next)

# fit <- princomp(combinedCommon, cor=TRUE)
# summary(fit) # print variance accounted for
# loadings(fit) # pc loadings
# plot(fit,type="lines") # scree plot
# fit$scores # the principal components
# biplot(fit) 
# plot(fit)
# 
# 
# # PCA two using prcomp
# pca.res <- prcomp(combinedCommon, retx=TRUE)
# pca.res
# summary(pca.res)
# 
# 
# # Get principal component vectors using prcomp instead of princomp
# pc <- prcomp(dataZ)
# 
# # First for principal components
# comp <- data.frame(pc$x[,1:4])
# # Plot
# plot(comp, pch=16, col=rgb(0,0,0,0.5))







# Rafa PCA plots!
x <- t(combinedCommon)#samples are the rows of the column matrix
x <- t(gadj3)#samples are the rows of the column matrix

pc <- prcomp(x)#scale = T, center = T) as of now I am not scaling

names(pc)

cols <- as.factor(substr(colnames(combinedCommon), 3, 7))##check me out. use 5 digit exp name.
plot(pc$x[, 1], pc$x[, 2], col=as.numeric(cols), main = "PCA", xlab = "PC1", ylab = "PC2")
legend("bottomright", levels(cols), col = seq(along=levels(cols)), pch = 1, lwd = 2)


summary(pc)

#SVD for calculating variance explained; see Rafa's notes for an explaination
cx <- sweep(x, 2, colMeans(x), "-")
sv <- svd(cx)
names(sv)
plot(sv$u[, 1], sv$u[, 2], col = as.numeric(cols), main = "SVD", xlab = "U1", ylab = "U2")


plot(sv$d^2/sum(sv$d^2), xlim = c(1, 6), type = "b", pch = 16, xlab = "principal components", 
     ylab = "variance explained")


# scaled? analysis...
# x <- t(dataZ)
# pc <- prcomp(x)
# 
# names(pc)
# 
# cols <- as.factor(substr(colnames(combinedCommon), 3, 7))##check me out. use 5 digit exp name.
# plot(pc$x[, 1], pc$x[, 2], col=as.numeric(cols), main = "PCA", xlab = "PC1", ylab = "PC2")
# legend("bottomleft", levels(cols), col = seq(along=levels(cols)), pch = 1)






# MDS with no Zscaling (MDS using euclidian distance is equivalent to PCA first two components)
d <- dist(t(combinedCommon))##note that their is no Z scaling here ?!
mds <- cmdscale(d)
cols <- as.factor(substr(colnames(combinedCommon), 3, 7))##check me out. use 5 digit exp name..
plot(mds, col=as.numeric(cols), lwd = 1.5)
legend("topleft", levels(cols), col = seq(along=levels(cols)), pch = 1)

# MDS with Zscaling I think...
d <- dist(t(dataZ))
mds <- cmdscale(d)
cols <- as.factor(substr(colnames(combinedCommon), 3, 7))##check me out. got this muthatrucka. use 5 digit exp name..
plot(mds, col=as.numeric(cols), lwd = 1.5)
legend("topleft", levels(cols), col = seq(along=levels(cols)), pch = 1)









#################################################################################################################




