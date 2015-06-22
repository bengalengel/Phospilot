

####Protein estimates from the gelprep and phosprep are compared. 

#conclusions (numbers 1-3 will be used to justify gelpreps use for manuscript)
1) gel prep provides many more protein quants
2) gel prep provides more sequence coverage and peptides per protein quant (implies greater accuracy)
3) gel prep clusters together and protein prep cluster together. PC2 and PC3 segregate by genetic background.
4) the correlation across samples (from the same sample workup approach) is lower for the gelprep compared to the phosprep and both are lower than the correlations achieved with with the phosphopeptides

Overall, this comparison implies that the protein level estimates from the gelprep approach have additional sample prep level variation that obscures the relevance of the conclusions. It is unlikely that these sample preparation sepcific biases can be completely removed through normalization and batch effect correction. 

Therefore I will run all tests with protein level estimates from PhosPrep and GelPrep separately. Keeping the advantages of depth and protein level estimate accuracy (for Gelprep) in mind vs the greater variance concordance with the PhosPrep numbers in mind.

Future work should be the 


##load Zia protein data
source("loadMQZ.R")
require(plyr)
require(limma)
require(swamp)
require(gplots)
require(Hmisc)
# directory = "D:/EnsemblDBProteome/iBAQ proteome/"
directory = "E:/My Documents/Pilot/EnsemblDBProteome/iBAQ proteome/"


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
nrow(Ziaproteins1)#4211

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
CombinedProt <- merge(Ziaproteins1, protein1, by = "Majority.protein.IDs", suffixes = c("GelPrep", "PhosPrep"))

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


ccPhosPrepRatios <- combat(swamp, o1$Factor1, batchcolumn=1) #WORKS AFTER ENSURING AT LEAST TWO IN A BATCH.
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


# Examine correlation structure amongst protein estimates from phospho data
#####################
#remove proteins that arent observed in each biological replicate
ccPhosPrepRatios <- ccPhosPrepRatios[rowSums(is.na(ccPhosPrepRatios[, 1:2])) < 2 & rowSums(is.na(ccPhosPrepRatios[, 3:4])) < 2 & 
                                       rowSums(is.na(ccPhosPrepRatios[, 5:6])) < 2 & rowSums(is.na(ccPhosPrepRatios[, 7:8])) < 2 & 
                                       rowSums(is.na(ccPhosPrepRatios[, 9:10])) < 2 & rowSums(is.na(ccPhosPrepRatios[, 11:12])) < 2, ]

#normalize again?
ccPhosPrepRatiosNorm <- normalizeQuantiles(ccPhosPrepRatios,ties = T)#ties are all assigned the same value for the common quantile



#medians from each sample and biological replicate to compare with each other and with Zia's data
ccPhosPrepRatios18486_1 <- apply(ccPhosPrepRatios[ , 1:2],1,function(x) median(x,na.rm =T))
ccPhosPrepRatios18486_2 <- apply(ccPhosPrepRatios[ , 3:4],1,function(x) median(x,na.rm =T))
ccPhosPrepRatios18862_1 <- apply(ccPhosPrepRatios[ , 5:6],1,function(x) median(x,na.rm =T))
ccPhosPrepRatios18862_2 <- apply(ccPhosPrepRatios[ , 7:8],1,function(x) median(x,na.rm =T))
ccPhosPrepRatios19160_1 <- apply(ccPhosPrepRatios[ , 9:10],1,function(x) median(x,na.rm =T))
ccPhosPrepRatios19160_2 <- apply(ccPhosPrepRatios[ , 11:12],1,function(x) median(x,na.rm =T))
ccPhosPrepRatios18486 <- apply(ccPhosPrepRatios[ , 1:4],1,function(x) median(x,na.rm =T))
ccPhosPrepRatios18862 <- apply(ccPhosPrepRatios[ , 5:8],1,function(x) median(x,na.rm =T))
ccPhosPrepRatios19160 <- apply(ccPhosPrepRatios[ , 9:12],1,function(x) median(x,na.rm =T))

#medians from each sample and biological replicate to compare with each other and with Zia's data
ccPhosPrepRatios18486_1 <- apply(ccPhosPrepRatiosNorm[ , 1:2],1,function(x) median(x,na.rm =T))
ccPhosPrepRatios18486_2 <- apply(ccPhosPrepRatiosNorm[ , 3:4],1,function(x) median(x,na.rm =T))
ccPhosPrepRatios18862_1 <- apply(ccPhosPrepRatiosNorm[ , 5:6],1,function(x) median(x,na.rm =T))
ccPhosPrepRatios18862_2 <- apply(ccPhosPrepRatiosNorm[ , 7:8],1,function(x) median(x,na.rm =T))
ccPhosPrepRatios19160_1 <- apply(ccPhosPrepRatiosNorm[ , 9:10],1,function(x) median(x,na.rm =T))
ccPhosPrepRatios19160_2 <- apply(ccPhosPrepRatiosNorm[ , 11:12],1,function(x) median(x,na.rm =T))
ccPhosPrepRatios18486 <- apply(ccPhosPrepRatiosNorm[ , 1:4],1,function(x) median(x,na.rm =T))
ccPhosPrepRatios18862 <- apply(ccPhosPrepRatiosNorm[ , 5:8],1,function(x) median(x,na.rm =T))
ccPhosPrepRatios19160 <- apply(ccPhosPrepRatiosNorm[ , 9:12],1,function(x) median(x,na.rm =T))

names <- grep("ccPhos(.*)_[0-9]$", ls(), value = T)
data <- setNames(lapply(names,function(x) get(x)), names)
PhosPrepRatioMedians <- as.data.frame(data)
boxplot(PhosPrepRatioMedians)

##these can be combined directly with the normalized ratios from the protein prep data


##NO I SHOULD COMPARE CORRELATION ACROSS BIOLOGICAL REPLICATES WITH CORRELATIONS ACROSS MEDIAN ESTIMATES FROM PROTEIN DATA

#calculate spearman correlations and do the-damn-thing (it's the hyphens that make it funny). The protein workup values agree more closely with one another across samples. This could be due to an increased number of peptides contributing to the protein group concentration estimates in the protein workup vs the phospho. It also could be due to batch effects obscuring the measurements.  
pairs(PhosPrepRatioMedians)


#using red gray blue color gradient to represent copynumber relative rank
rgbpal <- colorRampPalette(c("Red","Gray","Blue"))


cn.corr <- rcorr(as.matrix(PhosPrepRatioMedians), type = "spearman")
cn.corr <- rcorr(as.matrix(PhosPrepRatioMedians), type = "pearson")
cn.corr <- rcorr(as.matrix(MedianNorm), type = "pearson")#atrocious
cn.corr <- rcorr(as.matrix(quantiled4), type = "pearson")#still  atrocious until batch correction
cn.corr <- rcorr(as.matrix(ccPhosPrepRatios), type = "pearson") #batch correction doubles corr coef for 18486 biological replicate
cn.corr <- rcorr(as.matrix(ccPhosPrepRatiosNorm), type = "pearson") #negligible improvement in correlations after quantile normalizing again


heatmap.2(
  cn.corr$r,
  key = TRUE,
  density.info = "none",#no histogram in key
  key.xlab = "Pearson R", key.ylab=NULL, key.title = "",
  dendrogram = "none",
  trace = "none",
  cellnote = round(cn.corr$r,3),
  notecex = 1.25,
  notecol = "black",
  col=rev(rgbpal(20)),
  cexCol=1.254,
  cexRow=1.254,
  margins = c(15,15)
)
#######################


For phosphoprep derived protein estimates the protein correlation across individuals is ~ .5. How does this compare to the protein prep derived estimates across individuals? If similar why choose protein prep values? If greater this adds rationale for their use. The accuracy of the protein ratio estimates can be inferred by the number of unique peptides as well. A barplot of median number of unique peptides per protein quantification would also convey this message. Bring the point home with the added number of quantifications. In the manuscript add the caveat that protein assignment is uncertain and quantification is always noisy due to unknown modifications affecting estimates. REGRESSION ATTENUATION ISSUES. THE DEGREE OF ATTENUATION IS PROPORTIONAL TO THE INTENSITY OF THE PEPTIDE/NUMBER OF PEAKS. 
###########

# Protein Prep correlations
####################

#retrieve, transform, invert and normalize
#subset away ratios from phospo. Normalize and batch correct.
expCol <- grep("HL18486$|HL18862$|HL19160$",names(CombinedProt))

ProtPrepRatios <- CombinedProt[,expCol]

row.names(ProtPrepRatios) <- CombinedProt$Majority.protein.IDs
ProtPrepRatios <- 1/ProtPrepRatios   #inverse becausebecause the Heavy sample is the standard for the proteomics work
ProtPrepRatios <- log2(ProtPrepRatios)

##change column names to match inversion
colnames(ProtPrepRatios)<- gsub(colnames(ProtPrepRatios), pattern = "HL", replacement = "LH")

#median normalize
names <- colnames(ProtPrepRatios)
median.subtract <- function(x){ x - median(x, na.rm = TRUE)}##create a function for median subtraction
MedianNorm <- colwise(median.subtract, names)(ProtPrepRatios) #create median subtracted data but loose intensity and the row names
row.names(MedianNorm) <- CombinedProt$Majority.protein.IDs

#summaries
summary(MedianNorm)
boxplot(MedianNorm)
nrow(MedianNorm)#ProteinGroup estimate is found across all of the samples

boxplot(MedianNorm)
#differences in distribution shape for sure with HL18486 and HL19160. Something to this 18486 vs the others that is consistent across replication.
par(mfrow = c(1, 1))
for (i in 1:(ncol(MedianNorm))){
  if(i==1) plot(density(MedianNorm[, i], na.rm=T), col = i, ylim = c(0,2))
  else lines(density(MedianNorm[, i], na.rm=T), col = i)
}

#quantile normalize data being compared
quantiled <- normalizeQuantiles(MedianNorm,ties = T)#ties are all assigned the same value for the common quantile
row.names(quantiled) <- CombinedProt$Majority.protein.IDs
summary(quantiled)
boxplot(quantiled)
# density plots all look the same now of course
for (i in 1:(ncol(quantiled))){
  if(i==1) plot(density(quantiled[, i], na.rm=T), col = i, ylim = c(0,1.9))
  else lines(density(quantiled[, i], na.rm=T), col = i)
}

#what is the correlation across these samples?

cn.corr <- rcorr(as.matrix(ProtPrepRatios), type = "pearson")
cn.corr <- rcorr(as.matrix(MedianNorm), type = "pearson")#median transform does not change correlation
cn.corr <- rcorr(as.matrix(quantiled), type = "pearson")

cn.corr <- rcorr(as.matrix(ProtPrepRatios), type = "spearman")
cn.corr <- rcorr(as.matrix(MedianNorm), type = "spearman")#median transform does not change correlation
cn.corr <- rcorr(as.matrix(quantiled), type = "spearman")


heatmap.2(
  cn.corr$r,
  key = TRUE,
  density.info = "none",#no histogram in key
  key.xlab = "Pearson R", key.ylab=NULL, key.title = "",
  dendrogram = "none",
  trace = "none",
  cellnote = round(cn.corr$r,3),
  notecex = 1.25,
  notecol = "black",
  col=rev(rgbpal(20)),
  cexCol=1.254,
  cexRow=1.254,
  margins = c(15,15)
)

#how does this compare with the protein data obtained from phosprep? 

# Combine 'quantile normalized samples above with median phosprep ratios per sample
PhosPrepRatiosbySample <- cbind(ccPhosPrepRatios18486, ccPhosPrepRatios18862, ccPhosPrepRatios19160)#medians of the norm/batch corrected data
PhosPrepRatiosbySample <- as.data.frame(PhosPrepRatiosbySample)
ProtRatiosbySample <- merge(quantiled, PhosPrepRatiosbySample, by = "row.names")
row.names(ProtRatiosbySample) <- ProtRatiosbySample$Row.names
keepers <- which(names(ProtRatiosbySample)!="Row.names")
ProtRatiosbySample <- ProtRatiosbySample[,keepers]
boxplot(ProtRatiosbySample)
summary(ProtRatiosbySample)

cn.corr <- rcorr(as.matrix(ProtRatiosbySample), type = "pearson")

heatmap.2(
  cn.corr$r,
  key = TRUE,
  density.info = "none",#no histogram in key
  key.xlab = "Pearson R", key.ylab=NULL, key.title = "",
  dendrogram = "none",
  trace = "none",
  cellnote = round(cn.corr$r,3),
  notecex = 1.25,
  notecol = "black",
  col=rev(rgbpal(20)),
  cexCol=1.254,
  cexRow=1.254,
  margins = c(15,15)
)


#RESULT IS THAT THE CORRELATIONS ACROSS SAMPLES ARE LOWER FOR THE PROTEIN PREPS THAN THE PHOSPHOPREPS WHEN COMPARING MEDIANS. PERHAPS THIS CAN BE IMPROVED WITH NORMALIZATION?

cn.corr <- rcorr(as.matrix(ccPhosPrepRatios), type = "pearson") #batch correction doubles corr coef for 18486 biological replicate
heatmap.2(
  cn.corr$r,
  key = TRUE,
  density.info = "none",#no histogram in key
  key.xlab = "Pearson R", key.ylab=NULL, key.title = "",
  dendrogram = "none",
  trace = "none",
  cellnote = round(cn.corr$r,3),
  notecex = 1.25,
  notecol = "black",
  col=rev(rgbpal(20)),
  cexCol=1.254,
  cexRow=1.254,
  margins = c(15,15)
)
##############


#HOW DO THE INDIVIDUAL CORRELATIONS COMPARE TO THE MEDIAN OF THE INTER-SAMPLE CORRELATIONS FROM THE PHOSPHO DERIVED PROTEIN PREPS? IT IS HIGHER
##############
protpearson <- cn.corr$r

#assign NA's to values of the matrix that are intra-sample
protpearson[1:4, 1:4] <- NA
protpearson[5:8, 5:8] <- NA
protpearson[9:12, 9:12] <- NA


#the median intra-sample correlation coefficient (.52) is greater than the median inter sample value for the full proteome data! This may be due to batch effects? I am going crazy now.
summary(protpearson[upper.tri(protpearson)])
###########

# This higher correlation exists despite lower sequence coverage and number of uniqe + razor peptides per protein! Note however no batch correction
############

#unique plus razor proteins
PhosPrep <- median(CombinedProt$Razor...unique.peptidesPhosPrep)
GelPrep <- median(CombinedProt$Razor...unique.peptidesGelPrep)
names(PhosPrep) <- "PhosPrep"
names(GelPrep) <- "GelPrep"
barplot(c(PhosPrep,GelPrep), ylab = "#peptides/protein", main = "Median unique and razor peptides/protein")

#sequence coverage
PhosPrep <- median(CombinedProt$Sequence.coverage....PhosPrep)
GelPrep <- median(CombinedProt$Sequence.coverage....GelPrep)
names(PhosPrep) <- "PhosPrep"
names(GelPrep) <- "GelPrep"
barplot(c(PhosPrep,GelPrep), ylab = "#peptides/protein", main = "Median sequence coverage/protein")
################


###### Combine combat corrected phosprepratios with gelprepratios to examine by clustering and pca. The largest variance contributor is predictably experimental design. However pc2 and pc3 segregate the samples by genetic background.

#############
# ProtRatiosbySample compares the median normalized batch corrected data from phosprep with normalized estimates from gelprep
str(ProtRatiosbySample)
names(ProtRatiosbySample) <- c("GelPrep18486","GelPrep18862","GelPrep19160", "PhosPrep18486", "PhosPrep18862", "PhosPrep19160")

# now for the full dataset n=8560
cdata <- na.omit(ProtRatiosbySample)#also used below

# PCA analysis and dendrograms on corrected data
x <- t(cdata)#samples are the rows of the column matrix
pc <- prcomp(x)#scale = T, center = T) as of now I am not scaling

cols <- as.factor(gsub("[a-zA-z]", "", names(cdata)))
pchs <- as.factor(gsub("[0-9]", "", names(cdata)))
plot(pc$x[, 1], pc$x[, 2], col = as.numeric(cols), pch = c(16,17)[as.numeric(pchs)], main = "PCA", xlab = "PC1", ylab = "PC2")
legend("top", levels(cols), col = seq(along=levels(cols)), pch = 19)

cols <- as.factor(gsub("[a-zA-z]", "", names(cdata)))
pchs <- as.factor(gsub("[0-9]", "", names(cdata)))
plot(pc$x[, 2], pc$x[, 3], col = as.numeric(cols), pch = c(16,17)[as.numeric(pchs)], main = "PCA", xlab = "PC2", ylab = "PC3")
legend("top", levels(cols), col = seq(along=levels(cols)), pch = 19)

#dendrograms
dataZ <- scale(cdata)##Z-scored column wise

# dendogram using euclidian distance (default) and ward or complete agglomeration
dend.ward<- as.dendrogram(hclust(dist(t(dataZ)),method="ward"))
dend.complete<- as.dendrogram(hclust(dist(t(dataZ))))

ward.o<- order.dendrogram(dend.ward)
complete.o<- order.dendrogram(dend.complete)

plot(dend.complete,ylab="height", main = "Euclidian/Complete")
plot(dend.ward, leaflab = "perpendicular", ylab = "height", main = "Euclidian/Ward")






####################

# ibaq based quantifications. Do they cluster by individual as well? the copynumber estimates for zia's work correlate much more closely than the h/l ratios'
####################




