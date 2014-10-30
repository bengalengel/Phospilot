# ComBAT TESTS with sparse matrix
require(limma)
require(sva)
require(plyr)
require(swamp)
#make the sparse matrix from real data (all quants of class 1 here) 

data <- multExpanded1[,newnames]

data <- cbind(data,multExpanded1$Intensity)
names(data)[13] <- "Int" #for MA plots


data <- log2(data)

#only peptides found in all samples
#normalize by col median using columnwise (when to normalize??)

median.subtract <- function(x){ x - median(x, na.rm = TRUE)}##create a wrapper for median subtraction
data <- colwise(median.subtract, newnames)(data) #create median subtracted data but loose the row names here...

# row.names(datanorm) <- row.names(data)##add back the row names


# remove exp obs if not observed in each sample 
data2 <- data[rowSums(is.na(data[ , 1:4])) < 4 & rowSums(is.na(data[ , 5:8])) < 4 & rowSums(is.na(data[ , 9:12])) < 4,]    

# remove exp obs if not observed in each batch
data3 <- data[rowSums(is.na(data[ , c("HL18486_1_1", "HL18486_1_2", "HL18862_1_1", "HL18862_1_2", "HL19160_1_1", "HL19160_1_2")])) < 6 
              & rowSums(is.na(data[, c("HL18486_2_1", "HL18486_2_2", "HL18862_2_1", "HL18862_2_2", "HL19160_2_1", "HL19160_2_2")])) < 6,]    

# remove exp obs if not observed twice in each batch to ensure a variance measurement
data4 <- data[rowSums(is.na(data[ , c("HL18486_1_1", "HL18486_1_2", "HL18862_1_1", "HL18862_1_2", "HL19160_1_1", "HL19160_1_2")])) < 5 
              & rowSums(is.na(data[, c("HL18486_2_1", "HL18486_2_2", "HL18862_2_1", "HL18862_2_2", "HL19160_2_1", "HL19160_2_2")])) < 5,] 


boxplot(data2, ylim= c(-4,4))#note the uneven distributions!

# quantile normalize using ties=T for now
quantiled <- normalizeQuantiles(data2,ties = T)
quantiled <- normalizeQuantiles(data3,ties = T)
quantiled <- normalizeQuantiles(data4,ties = T)



boxplot(quantiled, ylim=c(-4,4))# well thats nice looking...


#batch effect identification 


library(swamp)
quantiled <- na.omit(quantiled)
quantiled <- as.matrix(quantiled)

##### sample annotations (data.frame)
set.seed(50)
o1<-data.frame(Factor1=factor(rep(c("A","A","B","B"),3)),
              Numeric1=rnorm(12),row.names=colnames(quantiled))


# PCA analysis
res1<-prince(quantiled,o1,top=10,permute=T)
str(res1)
res1$linp#plot p values
res1$linpperm#plot p values for permuted data
prince.plot(prince=res1)

#There is a batch effect associated with the process date.
# I must combat this

##batch adjustment for Factor 3
com1<-combat(quantiled,o1$Factor1,batchcolumn=1)
##batch adjustment for Factor 3; with covariate
com2<-combat(quantiled,o1$Factor1,batchcolumn=1) #WORKS AFTER ENSURING AT LEAST TWO IN A BATCH

# how did I do?
prince.plot(prince(com1,o1,top=10)) 
#I did well 

# now for the real shabang
cdata <- na.omit(com2)

prince.plot(prince(cdata,o1,top=10)) #huzzah!

# now with missing data removed perform the clustering************************************************************************
dataZ <- scale(cdata)##Z-scored column wise

# now all data excepting complete cases (note that the sample dendograms look the same)
hist(dataZ[,6], breaks = 100)

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
r <- t(scale(t(cdata)))#transpose to zscale the rows then transpose back to original format

# sample scaled
c <- scale(cdata)


# install heatmap.2 package
# install.packages("gplots")
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

# **********************************************************************************************

# now for some DE
dim(com2)

#first on the common subset cdata
##################################################### LIMMA for DE #####################################################
#Biological replication is needed for a valid comparison 

library(limma)


install.packages("statmod")
library(statmod)

# Calculate the correlation between technical replicates?...
# biolrep <- c(1, 1, 2, 2, 3, 3) 
# corfit <- duplicateCorrelation(datanorm, ndups = 1, block = biolrep)


# Produce dataframe from sample means ignoring missing data

HL18486_1 <- rowMeans(cdata[,1:2], na.rm = T)
HL18486_2 <- rowMeans(cdata[,3:4], na.rm = T)
HL18862_1 <- rowMeans(cdata[,5:6], na.rm = T)
HL18862_2 <- rowMeans(cdata[,7:8], na.rm = T)
HL19160_1 <- rowMeans(cdata[,9:10], na.rm = T)
HL19160_2 <- rowMeans(cdata[,11:12], na.rm = T)


# HL18486_1 <- rowMeans(bfdata[,1:2], na.rm = T)
# HL18486_2 <- rowMeans(bfdata[,3:4], na.rm = T)
# HL18862_1 <- rowMeans(bfdata[,5:6], na.rm = T)
# HL18862_2 <- rowMeans(bfdata[,7:8], na.rm = T)
# HL19160_1 <- rowMeans(bfdata[,9:10], na.rm = T)
# HL19160_2 <- rowMeans(bfdata[,11:12], na.rm = T)




# Better 
# 
# HL18486 <- rowMeans(datanorm[,1:4], na.rm = T)
# HL18862 <- rowMeans(datanorm[,5:8], na.rm = T)
# HL19160 <- rowMeans(datanorm[,9:12], na.rm = T)


pilot <- cbind(HL18486_1, HL18486_2, HL18862_1, HL18862_2, HL19160_1, HL19160_2)

boxplot(pilot)

# pilot <- cbind(HL18486, HL18862, HL19160)

pilot2 <- na.omit(pilot)
#note the strange outlier 

boxplot(pilot2)

#Produce the design matrix

fac <- factor(c(1,1,2,2,3,3))##codes the grouping for the ttests
design <- model.matrix(~0 + fac)
dnames <- levels(as.factor(substr(colnames(pilot2), 1, 7))) ##check me out. use 5 digit exp name.
colnames(design) <- dnames

#limma fit using all common for now
fit <- lmFit(pilot2, design)

#Now to make all pairwise comparisons (from Smyth pg 14)
# contrast.matrix <- makeContrasts(HL16778-HL16770, HL16788-HL16778, HL16788-HL16770, levels = design) 

contrast.matrix <- makeContrasts(HL18862-HL18486, HL19160-HL18862, HL19160-HL18486, levels = design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)



#Look at pairwise DE using toptable and the coef parameter to id which genes you are interested in 
topTable(fit2, coef = 1, adjust = "fdr")

results <- decideTests(fit2)

vennDiagram(results) #shazam but I need to remove outliers and the like


# AND NOW FOR ALL THE DATA around 5K
adata <- com2[rowSums(is.na(com2[ , 1:2])) < 2 & rowSums(is.na(com2[ , 3:4])) < 2 & rowSums(is.na(com2[ , 5:6])) < 2 
              & rowSums(is.na(com2[ , 7:8])) < 2 & rowSums(is.na(com2[ , 9:10])) < 2 & rowSums(is.na(com2[ , 11:12])) < 2,]                    
  
# Produce dataframe from sample means ignoring missing data

HL18486_1 <- rowMeans(adata[,1:2], na.rm = T)
HL18486_2 <- rowMeans(adata[,3:4], na.rm = T)
HL18862_1 <- rowMeans(adata[,5:6], na.rm = T)
HL18862_2 <- rowMeans(adata[,7:8], na.rm = T)
HL19160_1 <- rowMeans(adata[,9:10], na.rm = T)
HL19160_2 <- rowMeans(adata[,11:12], na.rm = T)


# HL18486_1 <- rowMeans(bfdata[,1:2], na.rm = T)
# HL18486_2 <- rowMeans(bfdata[,3:4], na.rm = T)
# HL18862_1 <- rowMeans(bfdata[,5:6], na.rm = T)
# HL18862_2 <- rowMeans(bfdata[,7:8], na.rm = T)
# HL19160_1 <- rowMeans(bfdata[,9:10], na.rm = T)
# HL19160_2 <- rowMeans(bfdata[,11:12], na.rm = T)




# Better 
# 
# HL18486 <- rowMeans(datanorm[,1:4], na.rm = T)
# HL18862 <- rowMeans(datanorm[,5:8], na.rm = T)
# HL19160 <- rowMeans(datanorm[,9:12], na.rm = T)


pilot <- cbind(HL18486_1, HL18486_2, HL18862_1, HL18862_2, HL19160_1, HL19160_2)

boxplot(pilot)

# pilot <- cbind(HL18486, HL18862, HL19160)

pilot2 <- na.omit(pilot)
#note the strange outlier 

boxplot(pilot2)

#Produce the design matrix

fac <- factor(c(1,1,2,2,3,3))##codes the grouping for the ttests
design <- model.matrix(~0 + fac)
dnames <- levels(as.factor(substr(colnames(pilot2), 1, 7))) ##check me out. use 5 digit exp name.
colnames(design) <- dnames

#limma fit using all common for now
fit <- lmFit(pilot, design)

#Now to make all pairwise comparisons (from Smyth pg 14)
# contrast.matrix <- makeContrasts(HL16778-HL16770, HL16788-HL16778, HL16788-HL16770, levels = design) 

contrast.matrix <- makeContrasts(HL18862-HL18486, HL19160-HL18862, HL19160-HL18486, levels = design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)



#Look at pairwise DE using toptable and the coef parameter to id which genes you are interested in 
topTable(fit2, coef = 1, adjust = "fdr")

results <- decideTests(fit2)

vennDiagram(results) #shazam but I need to remove outliers and the like














# Batch effect correction on full matrix using combat
quantiled <- as.matrix(quantiled)

require(sva)
# first start with quantiled
#make the batch covariate string matrix which identifes which sample belongs to which batch
batch <- rep(c("first","first", "second","second"),3)

# now make the model matrix of outcome of interest
# mod <- model.matrix(~batch) NOPE



fac <- factor(c(rep("one",4),rep("two",4), rep("three",4)))##codes the grouping for the ttests
fac <- as.character(fac)


fac <- factor(c(1,1,1,1,2,2,2,2,3,3,3,3))##codes the grouping for the ttests
design <- model.matrix(~0 + fac)
dnames <- levels(as.factor(substr(colnames(quantiled), 1, 7))) ##check me out. use 5 digit exp name.
colnames(design) <- dnames

quantiled1 <- na.omit(quantiled)

# prepare for combat
test <- ComBat(quantiled1,batch,design)



