##########################################variance components tests###############################################
##bimodal viarance component signature could be caused by artifacts derived from: norm, batch correct, multiplicity, MS acquisition type, and SILAC pair assignments. These tests are meant to examine the presense of such artifacts.

#test 1. MS acquisition and SILAC pair assignment bias. Result is bimodal signature stil present.

##multExpanded1 with all Identification types 'by MS/MS'. That is, MBRs is not used
index1 <- grep("Identification.type.(.*)", colnames(multExpanded1))
#keep only sites acquired by MS/MS across all samples
multExpanded1MSonly <- multExpanded1[rowSums(multExpanded1[,index1] == "By MS/MS") == 12,]

#keep only sites estimated from SILAC pairs without re-quantify
index2 <- which(colnames(multExpanded1) == "Ratio.H.L.iso.count")
multExpanded1NoRQ <- multExpanded1[(multExpanded1[,index2] == 0),]

#need the BatchNorm function to return list of dataframes, pilotMS is list of 8 DFs. quantileds, com2, adata, and pilot
MSonly <- BatchNorm(multExpanded1=multExpanded1MSonly)
NoRQ <- BatchNorm(multExpanded1NoRQ)

#extract list elements quantiled5 and com2 
quantiled5MSO <- MSonly[[7]]#quantile normalized no missing data
com2MSO <- MSonly[[8]]#normalized and batch corrected (nestedvar function will ensure balanced data (nomissingvalues))
quantiled5NoRQ <- NoRQ[[7]]#normalized no missing data
com2NoRQ <- NoRQ[[8]]#normalized and batch corrected (nestedvar function will ensure balanced data (nomissingvalues))

#feed to nestedvar, which fits the nested random effects model and returns the variance components
varcompNormMSO <- NestedVar(ratios=quantiled5MSO)
varcompBatchCorrectMSO <- NestedVar(ratios=com2MSO)
varcompNormNoRQ <- NestedVar(ratios=quantiled5NoRQ)
varcompBatchCorrectNoRQ <- NestedVar(ratios=com2NoRQ)

#Again result is that the bimodal variance component distribution is present with these sources of bias missing.

#test #2. Impact of mult on variance component signature. Variance signature still present.

#subject to 'NestedVar' normalized and batch corrected SILAC ratios for balanced data (com2 from full MultExpanded1)
CorrectedData <- BatchNorm(multExpanded1=multExpanded1)
com2 <- CorrectedData[[8]]

#remove from com2 all sites with multiplicity greater then 1 and feed 'NestedVar'
com2single <- com2[grepl("[0-9]+_1", row.names(com2)),]
varcompSingleMult <- NestedVar(ratios=com2single)


#test #3. Impact of normalization and batch correction.

#feed 'non-normalized' (MQ still has normalized these ratios by intensity) data into model
RawRatios <- CorrectedData[[1]]
RawRatios <- RawRatios[,1:12]
varcompMQnormonly <- NestedVar(ratios=RawRatios)

#feed the same unnormalized data by batch into the model
RawRatiosB1 <- RawRatios[,c(1:2,5:6,9:10)]
RawRatiosB2 <- RawRatios[,c(3:4,7:8,11:12)]

varcompMQnormonlyB1 <- NestedVar(ratios=RawRatiosB1, batch=T)
varcompMQnormonlyB2 <- NestedVar(ratios=RawRatiosB2, batch=T)

#test #4. Is signature there using MQ non-normalized ratio values and separate batches.
#loading unnormalized values
phosphoRaw <- load.MQ2(directory = "D:/10_9_14/txt/", type = "phospho")
phosphoRaw <- load.MQ2(directory = "E:/My Documents/Pilot/10_9_14/txt/", type = "phospho")
phosphoRaw <- phosphoRaw[(phosphoRaw$Potential.contaminant != "+" & phosphoRaw$Reverse != "+"),]

# subset phospho to class 1
phosphoRaw1 <- phosphoRaw[(phosphoRaw$Localization.prob >= .75),]

#for non normalized data
multExpandedRaw1 <- ExpandPhos2(phospho=phosphoRaw1)

#a curious inversion of the outlier. Here is the code from BatchNorm to extract dataframe for nestedvar

expCol <- grep("HL(.*)", colnames(multExpandedRaw1))
rawdata <- multExpandedRaw1[,expCol]

# add row names with site id and multiplicity designation 
# row.names(data) <- multExpanded$id#note this won't work because of the multiplicity issue
idmult <- paste(multExpandedRaw1$id, multExpandedRaw1$multiplicity, sep="_")
row.names(rawdata) <- idmult
rawdata <- log2(rawdata)

#summaries
summary(rawdata)
boxplot(rawdata[,1:12])#extreme outlier in HL18486_1_2
boxplot(rawdata[,1:12], ylim=c(-4,4))

#remove outlier which is 3X smaller than any other datapoint but was 3X greater for the normalized values
rawdata[,2][which.min(rawdata[,2])] <- NaN
boxplot(rawdata[,1:12])#extreme outlier in HL18486_1_2 removed

TotallyRawRatios <- rawdata
#Do unnormalized values have the same signature? Yes
varcompRaw <- NestedVar(ratios=TotallyRawRatios)

#feed the same data by batch into the model
TotallyRawRatiosB1 <- TotallyRawRatios[,c(1:2,5:6,9:10)]
TotallyRawRatiosB2 <- TotallyRawRatios[,c(3:4,7:8,11:12)]

#Do unnormalized bataches have this signature? Yes
varcompRawB1 <- NestedVar(ratios=TotallyRawRatiosB1, batch=T)
varcompRawB2 <- NestedVar(ratios=TotallyRawRatiosB2, batch=T)


#Hypothesis: Variation in labeling efficiency is responsible for bimodal biological replicate variance component distributions. 
# Data consistent with this proposition is observed segregation of proteins into modes of the distribution to greater extent in the biovarcomp vs the indvarcomp  
#x - percent of proteins unique to one mode & y - percent of proteins shared between two modes

#if any phosphosite from a protein is assigned to another mode it is binned shared otherwise it is binned 'uniquehigh' or 'uniquelow'

#subset relevant data to subject to pnvarcomp
keep <- grep(pattern = "pn.*Var|ppMajorityProteinIDs", names(multExpanded1_withDE), value = T)
SubtopnVarcomp <- multExpanded1_withDE[multExpanded1_withDE$ppSubtoVarcomp == "+", names(multExpanded1_withDE) %in% keep] #3483

#subset to proteins with more than one phosphosite/multiply phosphorylated (that is they have a chance to be shared)
multiplyphosphorylated <- unique(SubtopnVarcomp$ppMajorityProteinIDs[duplicated(SubtopnVarcomp$ppMajorityProteinIDs)])

SubtopnVarcomp2 <- SubtopnVarcomp[SubtopnVarcomp$ppMajorityProteinIDs %in% multiplyphosphorylated,]#2831 entries for 632 proteins

#if any of the peptides for a given protein are found in both 'high' and 'low' they are shared, for both 
require(plyr)
ProteinSpread <- ddply(SubtopnVarcomp2, "ppMajorityProteinIDs", function(x){
  ##is protein shared across modes?
  sharedbio <- ifelse(any(x$pnHighBioVar == "+") & any(x$pnLowBioVar == "+"), "+","-")
  sharedind <- ifelse(any(x$pnHighIndVar == "+") & any(x$pnLowIndVar == "+"), "+","-")
  data.frame(shared.bio = sharedbio, shared.ind = sharedind)
})

barplot(table(ProteinSpread$shared.bio), main = "Shared (+) vs segregated (-) multiply phosphorylated proteins. \n Biological variance comp", 
        ylab = "number of proteins")

barplot(table(ProteinSpread$shared.ind), main = "Shared (+) vs segregated (-) multiply phosphorylated proteins. \n Individual variance comp", 
        ylab = "number of proteins")
