rm(list=ls(all=TRUE)) #start with empty workspace


# First perform all processing steps using plyr and related tools.
# load required libraries
library(reshape2)
library(stringr)
library(plyr)
source("loadMQ.R")
source("ExpandPhos.R")
source("counts.R")
source("breakdown.R")
source("BatchNorm.R")
source("DiffPhos.R")

# load phospho and protein files with particular variables populated using "loadMQ"
phospho <- load.MQ(directory = "D:/10_9_14/txt/", type = "phospho")
protein <- load.MQ(directory = "D:/10_9_14/txt/", type = "protein")


# load phospho and protein files with particular variables populated using "loadMQ" at home
phospho <- load.MQ(directory = "E:/My Documents/Pilot/10_9_14/txt/", type = "phospho")
protein <- load.MQ(directory = "E:/My Documents/Pilot/10_9_14/txt/", type = "protein")


# remove contaminants and reverse database hits
phospho <- phospho[(phospho$Potential.contaminant != "+" & phospho$Reverse != "+"),]
protein <- protein[(protein$Potential.contaminant != "+" & protein$Reverse != "+"),]

#expand for observation (multiplicity based analysis). All observations quantified in >=1 experiment.
multExpanded <- ExpandPhos(phospho)

# subset phospho to class 1
phospho1 <- phospho[(phospho$Localization.prob >= .75),]

# "only identified by site" hits CAN BE removed because they tend to have lower PEPs (wouldn't pass the FDR TH anyway) and can't be quantified since they are not idd by non-modified peptides. 
# Note there are some high probability proteins here given some proteins are idd by 20+ phosphopeptides.
# eg is A6NKT7 (PEP = 2.23E-70)
protein1 <- protein[(protein$Only.identified.by.site != "+"),]

# Class 1 sites with each source of quantification for that site (singly/doubly/3+) explicitly accounted for 
multExpanded1 <- ExpandPhos(phospho1)

#make tables of basic counts of proteins and phosphopeptides (may have to update when performing the normalization)
phoscount(phospho,phospho1,multExpanded,multExpanded1)
proteincount(protein)

# make breakdown charts of phospho and protein overlap. This function produces barplots of number of sites/obs and proteins per experiment and
# cumulative over experiments. It also has extensive phospo info, including number of phospho per protein, multiplicity breakdown, and venns
breakdown(protein, phospho, multExpanded, cls=F)
breakdown(protein, phospho1, multExpanded1)


#normalization and batch correction (as of now only for class 1 sites). Accepts expanded phospho file and returns "pilot" dataframe. This dataframe is also written out.
pilot <- BatchNorm(multExpanded1)

#DE analysis uding DiffPhos function. Here I perform limma based DE across contrasts and add annotation to the multExpanded1 file. multiple images/venns are returned. 'multExpandedwithDE' is returned and written
multExpanded1_withDE <- DiffPhos(pilot)

#Zia protein workup using the multExpanded1_withDE data frame and 'pilot' dataframe as input. This program outputs a protein normalized matrix 'pilot2' which can be used for other puproses. Note the protein normalization of phospho data word file for the methods section of the paper. 
multExpanded1_withDE <- DiffPhosNorm(multExpanded1_withDE)

#next is GSEA, networkin, motifs to show that this is real functional biology

#next is the nested analysis to assess variation at the different levels of the experimental design to show biologically significant variation in each?

#next is work at the genome level to explicitly show that genetic variation is driving these changes. nonsynSNPs, pQTLs, nonsynSNPs surrounding the phosphosite. Correction for protein length (bootstrapping here). 

#something at the level of the differential network...







#########################################################################################################################
# Data Analysis. Histograms, Dendograms/Clustering/Heatmaps, PCA, DE, QQplots

##log2 the expression values (updated in main table)
multExpanded[newnames] <- log2(multExpanded[,newnames])##log2 transform

##a quick summary
# summary(multExpanded[,newnames])

# Filter so that there is at least one valid value in each sample
data <- multExpanded[,newnames]

data <- cbind(data,multExpanded$Intensity)
names(data)[13] <- "Int"


data <- log2(data)

#replace row names with phospho id later
# row.names(data) <- multExpanded$id#note this won't work because of the multiplicity issue

#note strange outlier in HL18486_1_2 
boxplot(data)


#normalize by col median using columnwise (when to normalize??)

median.subtract <- function(x){ x - median(x, na.rm = TRUE)}##create a wrapper for median subtraction
datanorm <- colwise(median.subtract, newnames)(data) #create median subtracted data but loose the row names here...

# row.names(datanorm) <- row.names(data)##add back the row names


# remove exp obs if not observed in each sample 
data2 <- data[rowSums(is.na(data[ , 1:4])) < 4 & rowSums(is.na(data[ , 5:8])) < 4 & rowSums(is.na(data[ , 9:12])) < 4,]    

# remove exp obs if not observed in each bio replicate (not sure how to automate this for larger datasets)
data3 <- data[rowSums(is.na(data[ , 1:2])) < 2 & rowSums(is.na(data[ , 3:4])) < 2 & rowSums(is.na(data[ , 5:6])) < 2 & rowSums(is.na(data[ , 7:8])) < 2 & rowSums(is.na(data[ , 9:10])) < 2 & rowSums(is.na(data[ , 11:12])) < 2,]                    

# remove exp obs if not observed in every experiment
data4 <- na.omit(data)

quantified <- nrow(data)
quantifiedsample <- nrow(data2)
quantifiedbio <- nrow(data3)
quantifiedall <- nrow(data4)
#overlap summary stats for the  
output2 <- data.frame(quantified, quantifiedsample, quantifiedbio, quantifiedall)
write.table(output2, "phosphosummary2_combined.csv", sep = ",", col.names = T, row.names = F)


# Histograms of data
par(mfrow = c(3,2))
phoshis1 <- hist(datanorm$HL1_A, breaks = 60)
phoshis2 <- hist(datanorm$HL16770_2, breaks = 60)
phoshis3 <- hist(datanorm$HL16778_1, breaks = 60)
phoshis4 <- hist(datanorm$HL16778_2, breaks = 60)
phoshis5 <- hist(datanorm$HL16788_1, breaks = 60)
phoshis6 <- hist(datanorm$HL16788_2, breaks = 60)


################################################ Dendograms and Clustering  ##########################################

#Zscore scale the samples. Note that NAs have an effect on the standardization and the distance metrics so they are removed
# Z score subtracts the mean of each vector from each element, then divides by the sd of the vector
# Z score makes each vector mean = 0 and sd = 1
boxplot(datanorm)

datanorm <- na.omit(datanorm)#now I have limited to observed in all cases.NOTE WAS NORMALIZED USING MORE DATA POINTS

boxplot(datanorm, ylim= c(-4,4))#note the uneven distributions!

# looking at some MA plots - there is no intensity bias!
plot(data2$Int, data2$HL18486_1_2, log="x", ylim=c(-6,6))
> plot(data2$Int, data2$HL18486_1_1, log="x", ylim=c(-6,6))
> plot(data2$Int, data2$HL18486_1_2, log="x", ylim=c(-6,6))
> plot(data2$Int, data2$HL18486_2_1, log="x", ylim=c(-6,6))
> plot(data2$Int, data2$HL18486_2_2, log="x", ylim=c(-6,6))
> plot(data2$Int, data2$HL18862_1_1, log="x", ylim=c(-6,6))
> plot(data2$Int, data2$HL18862_1_2, log="x", ylim=c(-6,6))
> plot(data2$Int, data2$HL18862_2_1, log="x", ylim=c(-6,6))
> plot(data2$Int, data2$HL18862_2_2, log="x", ylim=c(-6,6))
> plot(data2$Int, data2$HL19160_1_1, log="x", ylim=c(-6,6))
> plot(data2$Int, data2$HL19160_1_2, log="x", ylim=c(-6,6))
> plot(data2$Int, data2$HL19160_2_1, log="x", ylim=c(-6,6))
> plot(data2$Int, data2$HL19160_2_2, log="x", ylim=c(-6,6))

# so normalize quantiles on 'data2' median normalized data
quantiled <- normalizeQuantiles(datanorm,ties = T)


#remove pc1 due to confounding batch effect new test on quantiled

library(swamp)

quantiled <- na.omit(quantiled)
quantiled <- as.matrix(quantiled)

##### sample annotations (data.frame)
set.seed(50)
o<-data.frame(Factor1=factor(rep(c("A","A","B","B"),3)),
              Numeric1=rnorm(12),row.names=colnames(quantiled))


# PCA analysis
res1<-prince(quantiled,o,top=10,permute=T)
str(res1)
res1$linp#plot p values
res1$linpperm#plot p values for permuted data
prince.plot(prince=res1)

#remove confounding PC using svd (batch free data - bfdata)
bfdata<-kill.pc(datanorm,pc=1)
# to remove one or more principal components (here pc1) from the data
prince.plot(prince(bfdata,o,top=10))


boxplot(bfdata)#looks a bit better now











#remove pc1 due to confounding batch effect

library(swamp)

datanorm <- as.matrix(datanorm)

##### sample annotations (data.frame)
set.seed(50)
o<-data.frame(Factor1=factor(rep(c("A","A","B","B"),3)),
              Numeric1=rnorm(12),row.names=colnames(data))


# PCA analysis
res1<-prince(datanorm,o,top=10,permute=T)
str(res1)
res1$linp#plot p values
res1$linpperm#plot p values for permuted data
prince.plot(prince=res1)

#remove confounding PC using svd (batch free data - bfdata)
bfdata<-kill.pc(datanorm,pc=1)
# to remove one or more principal components (here pc1) from the data
prince.plot(prince(bfdata,o,top=10))


boxplot(bfdata)#looks a bit better now


##################################################### LIMMA for DE #####################################################
#Biological replication is needed for a valid comparison 

library(limma)


install.packages("statmod")
library(statmod)

# Calculate the correlation between technical replicates?...
# biolrep <- c(1, 1, 2, 2, 3, 3) 
# corfit <- duplicateCorrelation(datanorm, ndups = 1, block = biolrep)


# Produce dataframe from sample means ignoring missing data

HL18486_1 <- rowMeans(datanorm[,1:2], na.rm = T)
HL18486_2 <- rowMeans(datanorm[,3:4], na.rm = T)
HL18862_1 <- rowMeans(datanorm[,5:6], na.rm = T)
HL18862_2 <- rowMeans(datanorm[,7:8], na.rm = T)
HL19160_1 <- rowMeans(datanorm[,9:10], na.rm = T)
HL19160_2 <- rowMeans(datanorm[,11:12], na.rm = T)


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

vennDiagram(results)

















dataZ <- scale(bfdata)##Z-scored column wise

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
r <- t(scale(t(bfdata)))#transpose to zscale the rows then transpose back to original format

# sample scaled
c <- scale(bfdata)


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


###########################################################################################################################

##PCA of data from quick R (I will need to use the column standardized data next)

# fit <- princomp(datanorm, cor=TRUE)
# summary(fit) # print variance accounted for
# loadings(fit) # pc loadings
# plot(fit,type="lines") # scree plot
# fit$scores # the principal components
# biplot(fit) 
# plot(fit)
# 
# 
# # PCA two using prcomp
# pca.res <- prcomp(datanorm, retx=TRUE)
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
x <- t(datanorm)#samples are the rows of the column matrix
pc <- prcomp(x)#scale = T, center = T) as of now I am not scaling

names(pc)

cols <- as.factor(substr(colnames(datanorm), 3, 7))##check me out. use 5 digit exp name.
plot(pc$x[, 1], pc$x[, 2], col=as.numeric(cols), main = "PCA", xlab = "PC1", ylab = "PC2")
legend("bottomleft", levels(cols), col = seq(along=levels(cols)), pch = 1)


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
# cols <- as.factor(substr(colnames(datanorm), 3, 7))##check me out. use 5 digit exp name.
# plot(pc$x[, 1], pc$x[, 2], col=as.numeric(cols), main = "PCA", xlab = "PC1", ylab = "PC2")
# legend("bottomleft", levels(cols), col = seq(along=levels(cols)), pch = 1)






# MDS with no Zscaling (MDS using euclidian distance is equivalent to PCA first two components)
d <- dist(t(datanorm))##note that their is no Z scaling here ?!
mds <- cmdscale(d)
cols <- as.factor(substr(colnames(datanorm), 3, 7))##check me out. use 5 digit exp name..
plot(mds, col=as.numeric(cols), lwd = 1.5)
legend("topleft", levels(cols), col = seq(along=levels(cols)), pch = 1)

# MDS with Zscaling I think...
d <- dist(t(dataZ))
mds <- cmdscale(d)
cols <- as.factor(substr(colnames(datanorm), 3, 7))##check me out. got this muthatrucka. use 5 digit exp name..
plot(mds, col=as.numeric(cols), lwd = 1.5)
legend("topleft", levels(cols), col = seq(along=levels(cols)), pch = 1)









#################################################################################################################




