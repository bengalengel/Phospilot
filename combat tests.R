# Normalization and batch correction with ComBAT followed by normalization again
require(limma)
require(sva)
require(plyr)
require(swamp)
require(statmod)

#make the sparse matrix from real data (all quants of class 1 here) 
expCol <- grep("HL(.*)", colnames(multExpanded1))
data <- multExpanded1[,expCol]

# data <- cbind(data,multExpanded1$Intensity)
# names(data)[13] <- "Int" #for individual MA plots assessing bias of ratios with intensity

# add experiment intensities
intCol <- grep("Intensity", colnames(multExpanded1))

#a quick fix..> intCol <- intCol[1:(length(intCol)-3)]

data <- cbind(data,multExpanded1[,intCol])

# add row names with site id and multiplicity designation 
# row.names(data) <- multExpanded$id#note this won't work because of the multiplicity issue
idmult <- paste(multExpanded1$id, multExpanded1$multiplicity, sep="_")
row.names(data) <- idmult
multExpanded1 <- cbind(multExpanded1,idmult)
data <- log2(data)

# Within Individual experiment 'MA' plots with color change. M is log2(H/L) and A is the log2(total intensity from both heavy and light peptides ACROSS ALL EXPERIMENTS (Intensity) and for the individual experiment)
plot(data$Intensity, data$HL18486_1_1, ylim=c(-6,6))
plot(data$Intensity.18486_1_1, data$HL18486_1_1, ylim=c(-6,6))
plot(data$Intensity, data$HL18486_1_2, ylim=c(-6,6))
plot(data$Intensity.18486_1_2, data$HL18486_1_2, ylim=c(-6,6))
plot(data$Intensity, data$HL18486_2_1, ylim=c(-6,6))
plot(data$Intensity.18486_2_1, data$HL18486_2_1, ylim=c(-6,6))
plot(data$Intensity, data$HL18486_2_2, ylim=c(-6,6))
plot(data$Intensity.18486_2_2, data$HL18486_2_2, ylim=c(-6,6))
plot(data$Intensity, data$HL18862_1_1, ylim=c(-6,6))
plot(data$Intensity.18862_1_1, data$HL18862_1_1, ylim=c(-6,6))
plot(data$Intensity, data$HL18862_1_2, ylim=c(-6,6))
plot(data$Intensity.18862_1_2, data$HL18862_1_2, ylim=c(-6,6))
plot(data$Intensity, data$HL18862_2_1, ylim=c(-6,6))
plot(data$Intensity.18862_2_1, data$HL18862_2_1, ylim=c(-6,6))
plot(data$Intensity, data$HL18862_2_2, ylim=c(-6,6))
plot(data$Intensity.18862_2_2, data$HL18862_2_2, ylim=c(-6,6))
plot(data$Intensity, data$HL19160_1_1, ylim=c(-6,6))
plot(data$Intensity.19160_1_1, data$HL19160_1_1, ylim=c(-6,6))
plot(data$Intensity, data$HL19160_1_2, ylim=c(-6,6))
plot(data$Intensity.19160_1_2, data$HL19160_1_2, ylim=c(-6,6))
plot(data$Intensity, data$HL19160_2_1, ylim=c(-6,6))
plot(data$Intensity.19160_2_1, data$HL19160_2_1, ylim=c(-6,6))
plot(data$Intensity, data$HL19160_2_2, ylim=c(-6,6))
plot(data$Intensity.19160_2_2, data$HL19160_2_2, ylim=c(-6,6))

# Note here for the between array definition of MA plots for single channel data from the limma update notes section. For this reason I have also included 'A'

# The definition of the M and A axes for an MA-plot of single channel
# data is changed slightly.  Previously the A-axis was the average of
# all arrays in the dataset - this has been definition since MA-plots
# were introduced for single channel data in April 2003.  Now an
# artificial array is formed by averaging all arrays other than the
# one to be plotted.  Then a mean-difference plot is formed from the
# specified array and the artificial array.  This change ensures the
# specified and artificial arrays are computed from independent data,
# and ensures the MA-plot will reduce to a correct mean-difference
# plot when there are just two arrays in the dataset.



#Between experiment MA plot. There should be no bias in terms of intensity at this scale if there was no bias within experiments with intensity. here are technical replicates of the same bio rep in the same batch
HL184861 <- rowMeans(data[,1:2])
A <- rowMeans(data[,1:12], na.rm=T)
plot(HL184861, data$HL18486_1_1-data$HL18486_1_2)#besides the huge outlier there is no bias
plot(HL184861, data$HL18486_1_1-data$HL18486_1_2, ylim=c(-6,6), xlim=c(-8,8))#here can't see outlier

plot(A, data$HL18486_1_1-data$HL18486_1_2)#besides the huge outlier there is no bias
plot(A, data$HL18486_1_1-data$HL18486_1_2, ylim=c(-6,6), xlim=c(-8,8))#here can't see outlier

# here is the same biological replicate at different batches
HL184861_2 <- rowMeans(data[,c(1,4)])
plot(HL184861_2, data$HL18486_1_1-data$HL18486_2_2)
plot(HL184861_2, data$HL18486_1_1-data$HL18486_2_2, ylim=c(-6,6), xlim=c(-8,8))#here can't see outlier
# Note here that the blob has fattened, but there doesn't appear to be an intensity bias.


plot(A, data$HL18486_1_1-data$HL18486_2_2)
plot(A, data$HL18486_1_1-data$HL18486_2_2, ylim=c(-6,6), xlim=c(-8,8))#here can't see outlier

# here are two different samples within the same batch
change <- rowMeans(data[,c(1,5)])
plot(change, data$HL18486_1_1-data$HL18862_1_2)
plot(change, data$HL18486_1_1-data$HL18862_1_2, ylim=c(-6,6), xlim=c(-8,8))#here can't see outlier
# The blob is actually smaller than the same sample in different batches

plot(A, data$HL18486_1_1-data$HL18862_1_2)
plot(A, data$HL18486_1_1-data$HL18862_1_2, ylim=c(-6,6), xlim=c(-8,8))#here can't see outlier


# Below there seems to be some sort of curvalinear bias if you use the mean of the two ratios compared as opposed to the means
# across all samples

# Two different samples in separate batches. fattest yet.
change <- rowMeans(data[,c(1,8)])
plot(change, data$HL18486_1_1-data$HL18862_2_2)
plot(change, data$HL18486_1_1-data$HL18862_2_2, ylim=c(-6,6), xlim=c(-8,8))#here can't see outlier
# No not really but it is clear that there is no intensity bias

plot(A, data$HL18486_1_1-data$HL18862_2_2)
plot(A, data$HL18486_1_1-data$HL18862_2_2, ylim=c(-6,6), xlim=c(-8,8))#here can't see outlier

# Removing the intensity columns for now
data <- data[,1:12]

# Density plot overlays by color show a difference in the distributions
plot.new()
par(mfrow = c(1, 1))
for (i in 1:(ncol(data)-1)){
  if(i==1) plot(density(data[, i], na.rm=T), col = i, ylim = c(0,.7))
  else lines(density(data[, i], na.rm=T), col = i)
}

# This difference is due to batch
summary(data)
boxplot(data[,1:12])#extreme outlier in HL18486_1_2
boxplot(data[,1:12], ylim=c(-4,4))


# THE APPROACH NEEDS TO BE NORMALIZATION-BATCH EFFECT CORRECTION-RENORMALIZATION

# First normalization rationale: Given that we expect the same distributions across replicates quantile normalization makes sense. It is unclear if this should be preceded by a median normalization because MQ purportedly performs median normalization on the entire quantified peptide list. From the MQ paper supplement:

#   To correct for mixing errors of total protein amount the SILAC ratios determined in the
# previous section are normalized so that the median of logarithmized ratios is at zero. This
# normalization is done in intensity bins, similarly as described in the section on the
# calculation of protein ratios and significance. It is done separately for lysine and arginine
# labeled peptides to compensate for any possible label-specific bias. This peptide ratio
# normalization is done for each LC-MS run separately, allowing for different protein
# mixing ratios in different runs.

# Well if the samples are median normalized for every run the global median should still be zero since the median of combined median normalized vectors will be the same regardless of the length. However, according to Jurgen:

#   "The normalization is done on all peptides, meaning identified + unidentified, including these that do not even have an MS/MS spectrum. the evidence table contains only identified peptides. it means that there is a bias in this sample for identified peptides to have a smaller ratio."

# So I am going to median normalize again to account for this 'ascertainment bias' in my experiments as well

#***********************************************************************************************************************
# #another option is cyclic loess normalization, which may or may not be more robust to values missing in a non-random way.
# loessdata <- normalizeCyclicLoess(data[1:12])
# plot.new()
# par(mfrow = c(1, 1))
# for (i in 1:(ncol(loessdata)-1)){
#   if(i==1) plot(density(loessdata[, i], na.rm=T), col = i, ylim = c(0,.7))
#   else lines(density(loessdata[, i], na.rm=T), col = i)
# }
# loessdata2 <- colwise(median.subtract, names)(as.data.frame(loessdata))#median subtract using function below
# plot.new()
# for (i in 1:(ncol(loessdata2)-1)){
#   if(i==1) plot(density(loessdata2[, i], na.rm=T), col = i, ylim = c(0,.7))
#   else lines(density(loessdata2[, i], na.rm=T), col = i)
# }


#median normalize
names <- colnames(data)[1:12]
median.subtract <- function(x){ x - median(x, na.rm = TRUE)}##create a wrapper for median subtraction
data <- colwise(median.subtract, names)(data) #create median subtracted data but loose intensity and the row names here...

row.names(data) <- idmult##add back the row names

#summaries
summary(data)
boxplot(data[,1:12])#extreme outlier in HL18486_1_2
boxplot(data[,1:12], ylim=c(-4,4))

#remove outlier which is 3X greater than any other datapoint
data[,2][which.max(data[,2])] <- NaN
boxplot(data[,1:12])#extreme outlier in HL18486_1_2 removed

# density plots
plot.new()
par(mfrow = c(1, 1))
for (i in 1:(ncol(data))){
  if(i==1) plot(density(data[, i], na.rm=T), col = i, ylim = c(0,.7))
  else lines(density(data[, i], na.rm=T), col = i)
}

# quantile normalization. from normalize.quantiles {preprocessCore}  
# "This functions will handle missing data (ie NA values), based on the assumption that the data is missing at random."

quantiled <- normalizeQuantiles(data,ties = T)#ties are all assigned the same value for the common quantile
summary(quantiled)
boxplot(data)
boxplot(quantiled)
# density plots all look the same
plot.new()
par(mfrow = c(1, 1))
for (i in 1:(ncol(quantiled))){
  if(i==1) plot(density(quantiled[, i], na.rm=T), col = i, ylim = c(0,.55))
  else lines(density(quantiled[, i], na.rm=T), col = i)
}


quantiled <- cbind(quantiled,log2(multExpanded1[,intCol]))
multExpanded1[,intCol]
# add back logged intensities to see if ther is any new intra-experiment intensity (global that is) dependent bias
# quantiled <- cbind(quantiled,multExpanded1$Intensity)
# names(quantiled)[13] <- "Int" #for individual MA plots assessing bias of ratios with intensity
# quantiled$Int <- log2(quantiled$Int)

# Within Individual experiment 'MA' plots with color change. M is log2(H/L) and A is the log2(total intensity from both heavy and light peptides ACROSS ALL EXPERIMENTS)  There is experiment specific intensity but I don't think multiplicity explicit intensities are reported in the phospho table. They would have to be calcuated from the raw intensities for each experiment and the unnormalized ratios for each multiplicity.

plot(quantiled$Intensity, quantiled$HL18486_1_1, ylim=c(-6,6))
plot(quantiled$Intensity.18486_1_1, quantiled$HL18486_1_1, ylim=c(-6,6))
plot(quantiled$Intensity, quantiled$HL18486_1_2, ylim=c(-6,6))
plot(quantiled$Intensity.18486_1_2, quantiled$HL18486_1_2, ylim=c(-6,6))
plot(quantiled$Intensity, quantiled$HL18486_2_1, ylim=c(-6,6))
plot(quantiled$Intensity.18486_2_1, quantiled$HL18486_2_1, ylim=c(-6,6))
plot(quantiled$Intensity, quantiled$HL18486_2_2, ylim=c(-6,6))
plot(quantiled$Intensity.18486_2_2, quantiled$HL18486_2_2, ylim=c(-6,6))
plot(quantiled$Intensity, quantiled$HL18862_1_1, ylim=c(-6,6))
plot(quantiled$Intensity.18862_1_1, quantiled$HL18862_1_1, ylim=c(-6,6))
plot(quantiled$Intensity, quantiled$HL18862_1_2, ylim=c(-6,6))
plot(quantiled$Intensity.18862_1_2, quantiled$HL18862_1_2, ylim=c(-6,6))
plot(quantiled$Intensity, quantiled$HL18862_2_1, ylim=c(-6,6))
plot(quantiled$Intensity.18862_2_1, quantiled$HL18862_2_1, ylim=c(-6,6))
plot(quantiled$Intensity, quantiled$HL18862_2_2, ylim=c(-6,6))
plot(quantiled$Intensity.18862_2_2, quantiled$HL18862_2_2, ylim=c(-6,6))
plot(quantiled$Intensity, quantiled$HL19160_1_1, ylim=c(-6,6))
plot(quantiled$Intensity.19160_1_1, quantiled$HL19160_1_1, ylim=c(-6,6))
plot(quantiled$Intensity, quantiled$HL19160_1_2, ylim=c(-6,6))
plot(quantiled$Intensity.19160_1_2, quantiled$HL19160_1_2, ylim=c(-6,6))
plot(quantiled$Intensity, quantiled$HL19160_2_1, ylim=c(-6,6))
plot(quantiled$Intensity.19160_2_1, quantiled$HL19160_2_1, ylim=c(-6,6))
plot(quantiled$Intensity, quantiled$HL19160_2_2, ylim=c(-6,6))
plot(quantiled$Intensity.19160_2_2, quantiled$HL19160_2_2, ylim=c(-6,6))




# 
# plot(quantiled$Intensity, quantiled$HL18486_1_1, ylim=c(-6,6))
# plot(quantiled$Int, quantiled$HL18486_1_2, ylim=c(-6,6))
# plot(quantiled$Int, quantiled$HL18486_2_1, ylim=c(-6,6))
# plot(quantiled$Int, quantiled$HL18486_2_2, ylim=c(-6,6))
# plot(quantiled$Int, quantiled$HL18862_1_1, ylim=c(-6,6))
# plot(quantiled$Int, quantiled$HL18862_1_2, ylim=c(-6,6))
# plot(quantiled$Int, quantiled$HL18862_2_1, ylim=c(-6,6))
# plot(quantiled$Int, quantiled$HL18862_2_2, ylim=c(-6,6))
# plot(quantiled$Int, quantiled$HL19160_1_1, ylim=c(-6,6))
# plot(quantiled$Int, quantiled$HL19160_1_2, ylim=c(-6,6))
# plot(quantiled$Int, quantiled$HL19160_2_1, ylim=c(-6,6))
# plot(quantiled$Int, quantiled$HL19160_2_2, ylim=c(-6,6))


# there is no intensity bias between samples (see below). blobs look pretty much the same with no curvilinearity.

#Between experiment MA plot. There should be no bias in terms of intensity at this scale if there was no bias within experiments with intensity. here are technical replicates of the same bio rep in the same batch
HL184861 <- rowMeans(quantiled[,1:2])
# fit <- loess(quantiled$HL18486_1_1-quantiled$HL18486_1_2 ~ HL184861, na.rm=T)
plot(HL184861, quantiled$HL18486_1_1-quantiled$HL18486_1_2)#besides the huge outlier there is no bias
plot(HL184861, quantiled$HL18486_1_1-quantiled$HL18486_1_2, ylim=c(-6,6), xlim=c(-8,8))#here can't see outlier
# lines(HL184861, fit$fitted, col = 2, lwd = 2)

# here is the same biological replicate at different batches
HL184861_2 <- rowMeans(quantiled[,c(1,4)])
plot(HL184861_2, quantiled$HL18486_1_1-quantiled$HL18486_2_2)
plot(HL184861_2, quantiled$HL18486_1_1-quantiled$HL18486_2_2, ylim=c(-6,6), xlim=c(-8,8))#here can't see outlier
# Note here that the blob has fattened, but there doesn't appear to be an intensity bias.


# here are two different samples within the same batch
change <- rowMeans(quantiled[,c(1,5)])
plot(change, quantiled$HL18486_1_1-quantiled$HL18862_1_2)
plot(change, quantiled$HL18486_1_1-quantiled$HL18862_1_2, ylim=c(-6,6), xlim=c(-8,8))#here can't see outlier
# The blob is actually smaller than the same sample in different batches


# Two different samples in separate batches. fattest yet.
change <- rowMeans(quantiled[,c(1,8)])
plot(change, quantiled$HL18486_1_1-quantiled$HL18862_2_2)
plot(change, quantiled$HL18486_1_1-quantiled$HL18862_2_2, ylim=c(-6,6), xlim=c(-8,8))#here can't see outlier
# No not really but it is clear that there is no intensity bias

# END OF NORMALIZATION


# remove exp obs if not observed in each sample 
quantiled2 <- quantiled[rowSums(is.na(quantiled[ , 1:4])) < 4 & rowSums(is.na(quantiled[ , 5:8])) < 4 & rowSums(is.na(quantiled[ , 9:12])) < 4,]    

# remove exp obs if not observed in each batch
quantiled3 <- quantiled[rowSums(is.na(quantiled[ , c("HL18486_1_1", "HL18486_1_2", "HL18862_1_1", "HL18862_1_2", "HL19160_1_1", "HL19160_1_2")])) < 6 
                        & rowSums(is.na(quantiled[, c("HL18486_2_1", "HL18486_2_2", "HL18862_2_1", "HL18862_2_2", "HL19160_2_1", "HL19160_2_2")])) < 6,]    

# remove exp obs if not observed twice in each batch to ensure a variance measurement
quantiled4 <- quantiled[rowSums(is.na(quantiled[ , c("HL18486_1_1", "HL18486_1_2", "HL18862_1_1", "HL18862_1_2", "HL19160_1_1", "HL19160_1_2")])) < 5 
                        & rowSums(is.na(quantiled[, c("HL18486_2_1", "HL18486_2_2", "HL18862_2_1", "HL18862_2_2", "HL19160_2_1", "HL19160_2_2")])) < 5,] 
quantiled5 <- na.omit(quantiled)##common across all


















# # remove exp obs if not observed in each sample 
# data2 <- data[rowSums(is.na(data[ , 1:4])) < 4 & rowSums(is.na(data[ , 5:8])) < 4 & rowSums(is.na(data[ , 9:12])) < 4,]    
# 
# # remove exp obs if not observed in each batch
# data3 <- data[rowSums(is.na(data[ , c("HL18486_1_1", "HL18486_1_2", "HL18862_1_1", "HL18862_1_2", "HL19160_1_1", "HL19160_1_2")])) < 6 
#               & rowSums(is.na(data[, c("HL18486_2_1", "HL18486_2_2", "HL18862_2_1", "HL18862_2_2", "HL19160_2_1", "HL19160_2_2")])) < 6,]    
# 
# # remove exp obs if not observed twice in each batch to ensure a variance measurement
# data4 <- data[rowSums(is.na(data[ , c("HL18486_1_1", "HL18486_1_2", "HL18862_1_1", "HL18862_1_2", "HL19160_1_1", "HL19160_1_2")])) < 5 
#               & rowSums(is.na(data[, c("HL18486_2_1", "HL18486_2_2", "HL18862_2_1", "HL18862_2_2", "HL19160_2_1", "HL19160_2_2")])) < 5,] 


# boxplots 
mypar(1,1)
boxplot(quantiled[,1:12])
boxplot(quantiled2[,1:12])
boxplot(quantiled3[,1:12])
boxplot(quantiled4[,1:12])
boxplot(quantiled5[,1:12])
# boxplot(data2, ylim= c(-4,4))#note the uneven distributions!

# density plots
plot.new()
par(mfrow = c(1, 1))
for (i in 1:(ncol(quantiled))){
  if(i==1) plot(density(quantiled[, i], na.rm=T), col = i, ylim = c(0,.55))
  else lines(density(quantiled[, i], na.rm=T), col = i)
}
plot.new()
par(mfrow = c(1, 1))
for (i in 1:(ncol(quantiled2))){
  if(i==1) plot(density(quantiled2[, i], na.rm=T), col = i, ylim = c(0,.55))
  else lines(density(quantiled2[, i], na.rm=T), col = i)
}
plot.new()
par(mfrow = c(1, 1))
for (i in 1:(ncol(quantiled3))){
  if(i==1) plot(density(quantiled3[, i], na.rm=T), col = i, ylim = c(0,.55))
  else lines(density(quantiled3[, i], na.rm=T), col = i)
}
plot.new()
par(mfrow = c(1, 1))
for (i in 1:(ncol(quantiled4))){
  if(i==1) plot(density(quantiled4[, i], na.rm=T), col = i, ylim = c(0,.55))
  else lines(density(quantiled4[, i], na.rm=T), col = i)
}
plot.new()
par(mfrow = c(1, 1))
for (i in 1:(ncol(quantiled5))){
  if(i==1) plot(density(quantiled5[, i], na.rm=T), col = i, ylim = c(0,.7))
  else lines(density(quantiled5[, i], na.rm=T), col = i)
}
# # quantile normalize using ties=T for now
# quantiled <- normalizeQuantiles(data2,ties = T)
# quantiled <- normalizeQuantiles(data3,ties = T)
# quantiled <- normalizeQuantiles(data4,ties = T)




#batch effect identification and adjustment using swamp/combat*************************


swamp <- as.matrix(quantiled5)
swamp <- swamp[,1:12]
##### sample annotations (data.frame)
set.seed(50)
o1<-data.frame(Factor1=factor(rep(c("A","A","B","B"),3)),
               Numeric1=rnorm(12),row.names=colnames(swamp))


# PCA analysis
res1<-prince(swamp,o1,top=10,permute=T)
str(res1)
a <- res1$linp#plot p values
b <- res1$linpperm#plot p values for permuted data
prince.plot(prince=res1)

#There is a batch effect associated with the process date.
# I must combat this

##batch adjustment using the fully denuded data
com1<-combat(swamp,o1$Factor1,batchcolumn=1)

##batch adjustment using quantiled4
swamp <- as.matrix(quantiled4)
swamp <- swamp[,1:12]
com2<-combat(swamp,o1$Factor1,batchcolumn=1) #WORKS AFTER ENSURING AT LEAST TWO IN A BATCH. How to interpret plots...
# Found 2 batches
# Found 0 covariate(s)
# Found 25032 Missing Data Values
# Standardizing Data across genes
# Fitting L/S model and finding priors
# Finding parametric adjustments
# Adjusting the Data

##batch effect correction using sva combat and 'covariate' matrix

# how did I do?
prince.plot(prince(com1,o1,top=10)) 
#I did well 

# now for the full dataset n=8560
cdata <- na.omit(com2)
prince.plot(prince(cdata,o1,top=10)) #huzzah!

# PCA analysis
res1<-prince(cdata,o1,top=10,permute=T)
#str(res1)
c <- res1$linp#plot p values
d <- res1$linpperm#plot p values for permuted data
out <- rbind(a,b,c,d)
write.table(out, "PC_ba_batch.csv", sep=',', col.names=T, row.names=F) #PCs before and after batch correction.


##batch corrected EDA********************************************************************************
boxplot(com2, cex.axis = 1, cex.names = .5, cex.lab = .5, las=2)#fix the margins later
summary(com2)
# density plots
plot.new()
par(mfrow = c(1, 1))
for (i in 1:(ncol(com2))){
  if(i==1) plot(density(com2[, i], na.rm=T), col = i, ylim = c(0,.75))
  else lines(density(com2[, i], na.rm=T), col = i)
}
# I am not going to normalize again after batch correction because I am not sure if it makes any sense.



# now with missing data removed perform the clustering and heatmaps*******************************************
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

plot.new()##produces a blank canvas

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
  density.info="none",
  key.xlab = "Row Z scores", key.ylab=NULL, key.title = "",
  srtCol=45,  ,adjCol = c(1,1),
  margins = c(6,5),
  cexCol=1,
  labRow = NA#remove row labels
)


plot.new()

#PCA analysis 
# Rafa PCA plots!
x <- t(cdata)#samples are the rows of the column matrix
pc <- prcomp(x)#scale = T, center = T) as of now I am not scaling

names(pc)

cols <- as.factor(substr(colnames(cdata), 3, 7))##check me out. use 5 digit exp name.
plot(pc$x[, 1], pc$x[, 2], col=as.numeric(cols), main = "PCA", xlab = "PC1", ylab = "PC2")
legend("bottomleft", levels(cols), col = seq(along=levels(cols)), pch = 1)


summary(pc)

#SVD for calculating variance explained; see Rafa's notes for an explaination
cx <- sweep(x, 2, colMeans(x), "-")
sv <- svd(cx)
names(sv)
plot(sv$u[, 1], sv$u[, 2], col = as.numeric(cols), main = "SVD", xlab = "U1", ylab = "U2")


plot(sv$d^2/sum(sv$d^2), xlim = c(1, 12), type = "b", pch = 16, xlab = "principal components", 
     ylab = "variance explained")



# **********************************************************************************************

# now for some DE
# from a Jeff Leek post:
#   Several recent questions have focused on removing batch effects from gene
# expression or other high-throughput data as a cleaning step prior to
# performing other analyses. An important point about batch effect correction
# (whether with sva, combat, or any other currently published approach) is
# that a regression analysis is performed and variation is removed from the
# data. So subsequent analyses using a "cleaned" version of the data should
# be performed with caution. In particular, methods use to infer networks or
# to illustrate patterns (MDS/PCA) should be used with caution after
# regressing out batch effects. All currently published batch effect removal
# methods focus on adjusting batch effects for differential expression.

##################################################### LIMMA for DE #####################################################
#Biological replication is needed for a valid comparison 
# com2

# Calculate the correlation between technical replicates?...
# biolrep <- c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6)
# corfit <- duplicateCorrelation(com2, ndups = 1, block = biolrep)

#on the data with no NAs
#corfit2 <- duplicateCorrelation(cdata, ndups = 1, block = biolrep)

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





# AND NOW FOR ALL THE DATA around 5K*****************************************
adata <- com2[rowSums(is.na(com2[ , 1:2])) < 2 & rowSums(is.na(com2[ , 3:4])) < 2 & rowSums(is.na(com2[ , 5:6])) < 2 
              & rowSums(is.na(com2[ , 7:8])) < 2 & rowSums(is.na(com2[ , 9:10])) < 2 & rowSums(is.na(com2[ , 11:12])) < 2,]                    

# Produce dataframe from sample means ignoring missing data

HL18486_1 <- rowMeans(adata[,1:2], na.rm = T)
HL18486_2 <- rowMeans(adata[,3:4], na.rm = T)
HL18862_1 <- rowMeans(adata[,5:6], na.rm = T)
HL18862_2 <- rowMeans(adata[,7:8], na.rm = T)
HL19160_1 <- rowMeans(adata[,9:10], na.rm = T)
HL19160_2 <- rowMeans(adata[,11:12], na.rm = T)


pilot <- cbind(HL18486_1, HL18486_2, HL18862_1, HL18862_2, HL19160_1, HL19160_2)

boxplot(pilot)

#here pilot and pilot2 are the same
pilot2 <- na.omit(pilot) #4,996 class 1 measurements with at least one quant in each biological replicate

boxplot(pilot2)

#Produce the design matrix

fac <- factor(c(1,1,2,2,3,3))##codes the grouping for the ttests
design <- model.matrix(~0 + fac)
dnames <- levels(as.factor(substr(colnames(pilot2), 1, 7))) ##check me out. use 5 digit exp name.
colnames(design) <- dnames

#limma fit using all common for now.
# The philosophy of the approach is as follows. You have to start by fitting a linear model to
# your data which fully models the systematic part of your data. The model is specified by the design
# matrix. Each row of the design matrix corresponds to an array in your experiment and each column
# corresponds to a coefficient that is used to describe the RNA sources in your experiment.
# The main purpose of this step is to estimate the variability in the data, hence the systematic part needs to be modelled so it can be distinguished from random variation.

#Perhaps I can add replication correlation information. As of now sparseness throws an error. Example below:

# 17.3.6 Within-patient correlations
# The study involves multiple cell types from the same patient. Arrays from the same donor are not
# independent, so we need to estimate the within-dinor correlation:
#   > ct <- factor(targets$CellType)
# > design <- model.matrix(~0+ct)
# > colnames(design) <- levels(ct)
# > dupcor <- duplicateCorrelation(y,design,block=targets$Donor)
# > dupcor$consensus.correlation
# [1] 0.134
# As expected, the within-donor correlation is small but positive.

fit <- lmFit(pilot2, design)

# In practice the requirement to have exactly as many coefficients as RNA sources is too restrictive
# in terms of questions you might want to answer. You might be interested in more or fewer comparisons
# between the RNA source. Hence the contrasts step is provided so that you can take the initial
# coefficients and compare them in as many ways as you want to answer any questions you might have,
# regardless of how many or how few these might be.

#Now to make all pairwise comparisons (group2-1, group3-2, group3-1)
contrast.matrix <- makeContrasts(HL18862-HL18486, HL19160-HL18862, HL19160-HL18486, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)

# For statistical analysis and assessing differential expression, limma uses an empirical Bayes method
# to mo derate the standard errors of the estimated log-fold changes. This results in more stable
# inference and improved p ower, esp ecially for exp eriments with small numb ers of arrays
fit2 <- eBayes(fit2)



#Look at pairwise DE using toptable and the coef parameter to id which genes you are interested in 
sig1 <- topTable(fit2, coef = 1, adjust = "BH", n=Inf, sort="p", p=.05)#sorts by adjusted p up to the threshold of .05, which is the default FDR chosen for differential expression ("results" function). This actually seems a conservative way to sort.
sig2 <- topTable(fit2, coef = 2, adjust = "BH", n=Inf, sort="p", p=.05)
sig3 <- topTable(fit2, coef = 3, adjust = "BH", n=Inf, sort="p", p=.05)

# sig1 - 18862-18486
# sig2 - 19160-18862
# sig3 - 19160-18486

c1up  <- sig1[sig1$logFC > 0,]
c1down <- sig1[sig1$logFC < 0,]
c2up <- sig2[sig2$logFC > 0,]
c2down <- sig2[sig2$logFC < 0,]-
c3up <- sig3[sig3$logFC > 0,]
c3down <- sig3[sig3$logFC < 0,]


tt1 <- topTable(fit2, coef = 1, adjust = "BH", n=Inf)#sorts by adjusted p up to the threshold of .
tt2 <- topTable(fit2, coef = 2, adjust = "BH", n=Inf)#sorts by adjusted p up to the threshold of .
tt3 <- topTable(fit2, coef = 3, adjust = "BH", n=Inf)#sorts by adjusted p up to the threshold of .

hist(tt1$P.Value, nc=40, xlab="P values", main = colnames(contrast.matrix)[1])
hist(tt2$P.Value, nc=40, xlab="P values", main = colnames(contrast.matrix)[2])
hist(tt3$P.Value, nc=40, xlab="P values", main = colnames(contrast.matrix)[3])

plot(tt1$logFC,-log10(tt1$P.Value), xlab = colnames(contrast.matrix)[1], pch = 20, ylab = "-log10(P)",xlim = c(-5, 5))
#sites with sig difference in comparison 1
names <- row.names(sig1)
names2 <- row.names(tt1)
index <- which(names2 %in% names)
points(tt1$logFC[index],-log10(tt1$P.Value)[index], col="red3", pch = 20)


plot(tt2$logFC,-log10(tt2$P.Value), xlab = colnames(contrast.matrix)[2], pch = 20, ylab = "-log10(P)",xlim = c(-5, 5))
#sites with sig difference in comparison 1
names <- row.names(sig2)
names2 <- row.names(tt2)
index <- which(names2 %in% names)
points(tt2$logFC[index],-log10(tt2$P.Value)[index], col="red3", pch = 20)


plot(tt3$logFC,-log10(tt3$P.Value), xlab = colnames(contrast.matrix)[3], pch = 20, ylab = "-log10(P)",xlim = c(-5, 5))
#sites with sig difference in comparison 1
names <- row.names(sig3)
names2 <- row.names(tt3)
index <- which(names2 %in% names)
points(tt3$logFC[index],-log10(tt3$P.Value)[index], col="red3", pch = 20)




results <- decideTests(fit2, adjust.method = "BH", method = "separate")#results is a 'TestResults' matrix
#separate compares each sample individually and is the default approach
summary(results)


vennDiagram(results, cex=c(1.2,1,0.7)) #good DE across conditions
vennDiagram(results, cex=c(1.2,1,0.7), include = "up") #good DE across conditions
vennDiagram(results, cex=c(1.2,1,0.7), include = "down") #good DE across conditions
vennDiagram(results, cex=c(1.2,1,0.7), include = c("up", "down")) #good DE across conditions


table("18862-18486" =results[,1],"19160-18862"=results[,2])


volcanoplot(fit2, coef=1, main = colnames(contrast.matrix)[1], highlight = 900)
abline(v=1))
abline(v=-1)
volcanoplot(fit2, coef=2, main = colnames(contrast.matrix)[2])
abline(v=1)
abline(v=-1)
volcanoplot(fit2, coef=3, main = colnames(contrast.matrix)[3])
abline(v=1)
abline(v=-1)

# F statistic distributions and cuts by DE across individuals
# Fvalues <- as.data.frame(Fvalues)
# row.names(Fvalues) <- sites
# boxplot(log10(Fvalues))
# hist(log10(as.matrix(Fvalues)))
# plot(density(log10(as.matrix(Fvalues))))
# plot(density(as.matrix(Fvalues)))
# density(as.matrix(Fvalues))

plot(density(fit2$F))
plot(density(log10(fit2$F)))

# NOTE THE FOLLOWING WHEN COMPARING ALL THE CONTRASTS AT ONCE PG 62 IN LIMMA USER GUIDE - decideTests method global should be used here.

# method="global" is recommended when a set of closely related contrasts are being tested. This
# method simply appends all the tests together into one long vector of tests, i.e., it treats all the tests
# as equivalent regardless of which probe or contrast they relate to. An advantage is that the raw
# p-value cutoff is consistent across all contrasts. For this reason, method="global" is recommended if
# you want to compare the number of DE genes found for different contrasts, for example interpreting
# the number of DE genes as representing the strength of the contrast. However users need to be aware
# that the number of DE genes for any particular contrasts will depend on which other contrasts are
# tested at the same time. Hence one should include only those contrasts which are closely related to
# the question at hand. Unnecessary contrasts should be excluded as these would affect the results for
# the contrasts of interest. Another more theoretical issue is that there is no theorem which proves that
# adjust.method="BH" in combination with method="global" will correctly control the false discovery
# rate for combinations of negatively correlated contrasts, however simulations, experience and some
# theory suggest that the method is safe in practice.

#DE sites by contrast type
#DE sites using 'separate' contrasts
DE <- results[results[,1] != 0 | results[,2] != 0 | results[,3] != 0,]
#sites only DE in exactly one contrast
absDE <- abs(DE)
DE1 <- absDE[rowSums(absDE)==1,]
#sites DE in exactly two contrast
DE2 <- absDE[rowSums(absDE)==2,]
#site DE in all 3 contrasts
DE3 <- results[results[,1] != 0 & results[,2] != 0 & results[,3] != 0,]

#F stats by contrast type

# #Sorting F Values
# test <- topTable(fit2, coef = c(1,2,3), adjust = "BH", n=Inf, sort.by="F", p=.05)#equivalent to below
# test2 <- topTableF(fit2, adjust = "BH", n=Inf, sort.by="F", p=.05)

Fvals <- topTableF(fit2, adjust = "BH", n=Inf, sort.by="F")#all F values
sigFvals <- topTableF(fit2, adjust = "BH", n=Inf, sort.by="F", p=.05)#gives 1355 compared to 1549 DE total for separate comparisons

#subsets by contrast specific DE
FDE1 <- Fvals[which(row.names(DE1)%in%row.names(Fvals)),5:6]
FDE2 <- Fvals[which(row.names(DE2)%in%row.names(Fvals)),5:6]
FDE3 <- Fvals[which(row.names(DE3)%in%row.names(Fvals)),5:6]

#below gives adjusted pvalue
FDE1 <- Fvals[match(row.names(DE1), row.names(Fvals), nomatch = F),5:7]
FDE2 <- Fvals[match(row.names(DE2), row.names(Fvals), nomatch = F),5:7]
FDE3 <- Fvals[match(row.names(DE3), row.names(Fvals), nomatch = F),5:7]

boxplot(FDE1$F,FDE2$F,FDE3$F)
boxplot(log10(FDE1$F),log10(FDE2$F),log10(FDE3$F))
summary(FDE1$F)
summary(FDE2$F)
summary(FDE3$F)

plot(density(log10(FDE1$F)),xlim = c(0,3))
lines(density(log10(FDE2$F)), col = 2)
lines(density(log10(FDE3$F)), col = 3)

hist(log10(FDE1$F))
hist(FDE1$F)



# GSEA of differentially expressed lists across contrasts

#add annotation to multexpanded DF
head(row.names(pilot2))

#add F test values to the table
multExpanded1$globalFsig = ifelse(multExpanded1$idmult %in% row.names(sigFvals),"+","-")


#add DE to table
multExpanded1$SubtoDE = ifelse(multExpanded1$idmult %in% row.names(pilot2),"+","-")
multExpanded1$DEcont1 = ifelse(multExpanded1$idmult %in% row.names(sig1),"+","-")
multExpanded1$DEcont2 = ifelse(multExpanded1$idmult %in% row.names(sig2),"+","-")
multExpanded1$DEcont3 = ifelse(multExpanded1$idmult %in% row.names(sig3),"+","-")

#add DE direction to table
multExpanded1$cont1up = ifelse(multExpanded1$idmult %in% row.names(c1up),"+","-")
multExpanded1$cont1down = ifelse(multExpanded1$idmult %in% row.names(c1down),"+","-")
multExpanded1$cont2up = ifelse(multExpanded1$idmult %in% row.names(c2up),"+","-")
multExpanded1$cont2down = ifelse(multExpanded1$idmult %in% row.names(c2down),"+","-")
multExpanded1$cont3up = ifelse(multExpanded1$idmult %in% row.names(c3up),"+","-")
multExpanded1$cont3down = ifelse(multExpanded1$idmult %in% row.names(c3down),"+","-")


# write output table to perform enrichment analysis in perseus
write.table(multExpanded1,"multExpanded1_withDE.csv",sep=',',col.names=T,row.names=F)









