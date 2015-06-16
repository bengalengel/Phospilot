

####load protein estimates from protein1 and from Zia to assess the suitability of substituting Zia's data with mine.

#is there a good correlation for protein estimates across my data and protein estimates across Zia's data? (>.75 is good in my view)

#Approach for matching Zia's data with my data at the protein group level.

#1) complete match at the majority protein id level


##load Zia protein data
source("loadMQZ.R")
require(plyr)
require(limma)
require(swamp)

directory = "D:/EnsemblDBProteome/iBAQ proteome/"

# load protein files with particular variables populated using "loadMQ"
Ziaproteins <- load.MQZ(directory)#7710 protein groups/81 variables

# remove contaminants and reverse database hits
Ziaproteins <- Ziaproteins[(Ziaproteins$Potential.contaminant != "+" & Ziaproteins$Reverse != "+"),]#7354/78

colnames(Ziaproteins)<- gsub(colnames(Ziaproteins), pattern = "Ratio.H.L.normalized.", replacement = "HL") ##remove redundant information


# Keep only data from lines of interest
vars <- c("id", "Protein.IDs", "Majority.protein.IDs", "Number.of.proteins", "Peptides", 
          "Razor...unique.peptides", "Unique.peptides", "Sequence.coverage....", "Sequence.length", "Sequence.lengths",
          "Mol..weight..kDa.", "PEP", "Peptide.IDs", "Mod..peptide.IDs", "Only.identified.by.site", 
          "Potential.contaminant", "Reverse", "Razor...unique.peptides.18862", "Razor...unique.peptides.18486", 
          "Razor...unique.peptides.19160")

other_data <- Ziaproteins[,vars]
expression <- Ziaproteins[,grep("HL18862|HL18486|HL19160", colnames(Ziaproteins))]
ibaq <- Ziaproteins[,grep("iBAQ.L.18862|iBAQ.L.18486|iBAQ.L.19160", colnames(Ziaproteins), ignore.case = T)]


Ziaproteins <- cbind(other_data,expression,ibaq)

# "only identified by site" hits CAN BE removed because they tend to have lower PEPs (wouldn't pass the FDR TH anyway) and can't be quantified since they are not idd by non-modified peptides. 
Ziaproteins1 <- Ziaproteins[(Ziaproteins$Only.identified.by.site != "+"),]#6686


#some strangeness sometimes there are two extra rows!
#Ziaproteins1 <- Ziaproteins1[2:length(Ziaproteins1)]

#remove proteins if not quantified in at least one sample
expCol <- grep("HL(.*)", colnames(Ziaproteins1))

Ziaproteins1 <- Ziaproteins1[!rowSums(is.na(Ziaproteins1[,expCol])) >= 1,]##removes rows with any nas

#################

#subset 'protein1' similarly
#########
vars <- c("id", "Protein.IDs", "Majority.protein.IDs", "Number.of.proteins", "Peptides", 
          "Razor...unique.peptides", "Unique.peptides", "Sequence.coverage....", "Mol..weight..kDa.", "Sequence.length", "PEP", 
          "Peptide.IDs", "Mod..peptide.IDs", "Phospho..STY..site.IDs", "Only.identified.by.site", "Potential.contaminant", "Reverse")

other_data <- protein1[,vars]

colnames(protein1)<- gsub(colnames(protein1), pattern = "Ratio.H.L.normalized.", replacement = "HL") ##remove redundant information

expression <- protein1[,grep("HL18862*|HL18486*|HL19160*", colnames(protein1))]
ibaq <- protein1[,grep("iBAQ.L.18862*|iBAQ.L.18486*|iBAQ.L.19160*", colnames(protein1), ignore.case = T)]
protein1 <- cbind(other_data, expression, ibaq)

##at least one measurement per experiment
expCol <- grep("HL(.*)", colnames(protein1))
protein1 <- protein1[rowSums(is.na(protein1[,expCol[1:4]])) < 4 & rowSums(is.na(protein1[,expCol[5:8]])) < 4
                     & rowSums(is.na(protein1[,expCol[9:12]])) < 4,]

#############

# merge dataframes, create four separate DFs by protein estimate type, normalize/batch correct (if necessary), 

#size of each df
nrow(protein1)
nrow(Ziaproteins1)

# merge the two dataframes by majority protein idvto perform comparisons
CombinedProt <- merge(Ziaproteins1,protein1, by = "Majority.protein.IDs", suffixes = c("GelPrep", "PhosPrep"))

#subset away ratios from phospo. Normalize and batch correct.
expCol <- grep("HL(.*)_[0-9]",names(CombinedProt))

PhosPrepRatios <- CombinedProt[,expCol]

#log transform, normalize and batch correct
row.names(PhosPrepRatios) <- CombinedProt$Majority.protein.IDs
PhosPrepRatios <- log2(PhosPrepRatios)
boxplot(PhosPrepRatios)

#median then quantile normalize

#median normalize
names <- colnames(PhosPrepRatios)
median.subtract <- function(x){ x - median(x, na.rm = TRUE)}##create a function for median subtraction
MedianNorm <- colwise(median.subtract, names)(PhosPrepRatios) #create median subtracted data but loose intensity and the row names

#add back protien ids
row.names(MedianNorm) <- CombinedProt$Majority.protein.IDs

#summaries
summary(MedianNorm)
boxplot(MedianNorm)

#quantile normalize data being compared
quantiled <- normalizeQuantiles(MedianNorm,ties = T)#ties are all assigned the same value for the common quantile
summary(quantiled)
boxplot(quantiled)
# density plots all look the same now of course
plot.new()
par(mfrow = c(1, 1))
for (i in 1:(ncol(quantiled))){
  if(i==1) plot(density(quantiled[, i], na.rm=T), col = i, ylim = c(0,1.9))
  else lines(density(quantiled[, i], na.rm=T), col = i)
}

##EDA before batch correction reveals as expected a batch effect
############
cdata <- na.omit(quantiled)#559 observations
dataZ <- scale(cdata)##Z-scored column wise

# dendogram using euclidian distance (default) and ward or complete agglomeration
dend.ward<- as.dendrogram(hclust(dist(t(dataZ)),method="ward"))
dend.complete<- as.dendrogram(hclust(dist(t(dataZ))))

ward.o<- order.dendrogram(dend.ward)
complete.o<- order.dendrogram(dend.complete)

plot(dend.complete,ylab="height", main = "Euclidian/Complete")
plot(dend.ward, leaflab = "perpendicular", ylab = "height", main = "Euclidian/Ward")

#PCA analysis 
x <- t(cdata)#samples are the rows of the column matrix
pc <- prcomp(x)#scale = T, center = T) as of now I am not scaling

cols <- as.factor(substr(colnames(cdata), 3, 7))##use 5 digit exp name.
plot(pc$x[, 1], pc$x[, 2], col=as.numeric(cols), main = "PCA", xlab = "PC1", ylab = "PC2")
legend("bottom", levels(cols), col = seq(along=levels(cols)), pch = 1)
#######################

#Batch correction using Leek's 'Combat' via the wrapper package 'swamp'
##########
# remove exp obs if not observed two or more times in each batch to ensure a variance measurement
quantiled4 <- quantiled[rowSums(is.na(quantiled[ , c("HL18486_1_1", "HL18486_1_2", "HL18862_1_1", "HL18862_1_2", "HL19160_1_1", "HL19160_1_2")])) < 5 & rowSums(is.na(quantiled[, c("HL18486_2_1", "HL18486_2_2", "HL18862_2_1", "HL18862_2_2", "HL19160_2_1", "HL19160_2_2")])) < 5,] 


#batch effect identification and adjustment using swamp/combat*************************
##################
swamp <- as.matrix(cdata)
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
##batch adjustment using quantiled4
swamp <- as.matrix(quantiled4)
##### sample annotations (data.frame)
set.seed(50)
o1<-data.frame(Factor1=factor(rep(c("A","A","B","B"),3)),
               Numeric1=rnorm(12),row.names=colnames(swamp))


ccPhosPrepRatios <-combat(swamp,o1$Factor1,batchcolumn=1) #WORKS AFTER ENSURING AT LEAST TWO IN A BATCH.
######
par(mfrow = c(1,1))
boxplot(ccPhosPrepRatios)

# dendrograms and PCA of batch corrected data
#############
# now for the full dataset n=8560
cdata <- na.omit(ccPhosPrepRatios)#also used below
prince.plot(prince(cdata,o1,top=10)) #huzzah!

# PCA analysis and dendrograms on corrected data
x <- t(cdata)#samples are the rows of the column matrix
pc <- prcomp(x)#scale = T, center = T) as of now I am not scaling

cols <- as.factor(substr(colnames(cdata), 3, 7))##use 5 digit exp name.
plot(pc$x[, 1], pc$x[, 2], col=as.numeric(cols), main = "PCA", xlab = "PC1", ylab = "PC2")
legend("bottomleft", levels(cols), col = seq(along=levels(cols)), pch = 1)

#dendrograms
dataZ <- scale(cdata)##Z-scored column wise

# dendogram using euclidian distance (default) and ward or complete agglomeration
dend.ward<- as.dendrogram(hclust(dist(t(dataZ)),method="ward"))
dend.complete<- as.dendrogram(hclust(dist(t(dataZ))))

ward.o<- order.dendrogram(dend.ward)
complete.o<- order.dendrogram(dend.complete)

plot(dend.complete,ylab="height", main = "Euclidian/Complete")
plot(dend.ward, leaflab = "perpendicular", ylab = "height", main = "Euclidian/Ward")
############################

#medians from each sample to compare with Zia data
ccPhosPrepRatios18486 <- apply(ccPhosPrepRatios[ , 1:4],1,function(x) median(x,na.rm =T))
ccPhosPrepRatios18862 <- apply(ccPhosPrepRatios[ , 5:8],1,function(x) median(x,na.rm =T))
ccPhosPrepRatios19160 <- apply(ccPhosPrepRatios[ , 9:12],1,function(x) median(x,na.rm =T))

PhosPrepRatioMedians <- data.frame(ccPhosPrepRatios18486,ccPhosPrepRatios18862,ccPhosPrepRatios19160)
boxplot(PhosPrepRatioMedians)

##these can be combined directly with the normalized ratios from the protein prep data


##NO I SHOULD COMPARE CORRELATION ACROSS BIOLOGICAL REPLICATES WITH CORRELATIONS ACROSS MEDIAN ESTIMATES FROM PROTEIN DATA

#calculate spearman correlations and do the-damn-thing (it's the hyphens that make it funny). The protein workup values agree more closely with one another across samples. This could be due to an increased number of peptides contributing to the protein group concentration estimates in the protein workup vs the phospho. It also could be due to batch effects obscuring the measurements.  
pairs(PhosPrepRatioMedians)
require(gplots)
require(Hmisc)

#using red gray blue color gradient to represent copynumber relative rank
rgbpal <- colorRampPalette(c("Red","Gray","Blue"))


cn.corr <- rcorr(as.matrix(PhosPrepRatioMedians), type = "spearman")
heatmap.2(
  cn.corr$r,
  key = TRUE,
  density.info = "none",#no histogram in key
  key.xlab = "Spearman R", key.ylab=NULL, key.title = "",
  dendrogram = "none",
  trace = "none",
  cellnote = round(cn.corr$r,3),
  notecex = 1.25,
  notecol = "black",
  col=rev(rgbpal(20)),
  cexCol=1.254,
  cexRow=1.254,
  margins = c(10,10)
)

