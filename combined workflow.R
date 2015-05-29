rm(list=ls(all=TRUE)) #start with empty workspace


# First perform all processing steps using plyr and related tools.
# load required libraries
library(reshape2)
library(stringr)
library(plyr)
library(seqinr)
source("loadMQ.R")
source("ExpandPhos.R")
source("ExpandPhos2.R") #for expanding the non-normalized data
source("counts.R")
source("breakdown.R")
source("BatchNorm.R")
source("DiffPhos.R")
source("DiffPhosProt.R")
source("loadMQ2.R") #non normalized data
source("NestedVar.R")
source("NormProt.R")
source("ProtAssignment2.R")
source("DiffPhosProt.R")
source("AddAnnotation.R")
source("Enrichment.R")
source("PerseusOut.R"

#load, reformat and characterize MQ outputted mass spectrometry data
#######
# load phospho and protein files with particular variables populated using "loadMQ"
phospho <- load.MQ(directory = "D:/EnsemblDBPhospho/PilotPhosphoensemblDB/combined/txt/", type = "phospho")
# protein <- load.MQ(directory = "D:/EnsemblDBPhospho/PilotPhosphoensemblDB/combined/txt/", type = "protein")
protein <- load.MQ(directory = "D:/EnsemblDBPhospho/iBAQ quantification/", type = "protein")#with ibaq

# load phospho and protein files with particular variables populated using "loadMQ" at home
phospho <- load.MQ(directory = "E:/My Documents/Pilot/EnsemblDBPhospho/PilotPhosphoensemblDB/combined/txt/", type = "phospho")
# protein <- load.MQ(directory = "E:/My Documents/Pilot/EnsemblDBPhospho/PilotPhosphoensemblDB/combined/txt/", type = "protein")
protein <- load.MQ(directory = "E:/My Documents/Pilot/EnsemblDBPhospho/iBAQ quantification/", type = "protein")#with ibaq


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
#####################

#normalization and batch correction
#remove an outlier, normalize (median and quantile), and batch correct (combat). Returned are EDA plots and a list of DFs.See BatchNorm for details
######################
CorrectedData <- BatchNorm(multExpanded1=multExpanded1)#class 1 sites
com2 <- CorrectedData[[8]]#normalized/batch corrected (using ComBat) data frame
adata <- CorrectedData[[9]]#normalized/batch corrected data frame with at lesat 1 obs in each bio rep
pilot <- CorrectedData[[10]]#same as above with mean ratios for each bio replicate
#####################


#Differential phosphorylation analysis using DiffPhos function. Below is limma based DE across contrasts with annotation added to the passed multExpanded1 file. multiple images/venns are returned. 'multExpandedwithDE' is returned.
##########################
#!!BEFORE passing a matrix for diffPhos, remove all ids with leading protein mapping to a REV_ hit (32 for ME1). These ids are still present because there are majority protein ids for these peptides that map to a non reverse hit. These are omitted to make protein mapping to Zia's dataset more straightforward and enrichment analysis results comparable between the confounded and non-confounded datasets. That is what does it mean to compare enrichment for a reversed 'protein'?

#add idmult annotation to multexpanded table
idmult <- paste(multExpanded1$id, multExpanded1$multiplicity, sep="_")
multExpanded1 <- cbind(multExpanded1,idmult)
#remove reverse hits and acquire an index to subset corrected data
RevHits <- grep(multExpanded1$Protein, pattern = "REV")
RevHitsidmult <- multExpanded1$idmult[RevHits]
multExpanded1 <- multExpanded1[!grepl(multExpanded1$Protein, pattern = "REV"),] #18238 now 17774

#remove reverse hits from dataframe passed to Diffphos. Right now 'pilot' is used.
pilot <- pilot[!rownames(pilot)%in%RevHitsidmult,]

multExpanded1_withDE <- DiffPhos(pilot, multExpanded1)
#################

#protein normalization and diffphos on non-confounded phosphopeptides.
########################
#Load and normalize protein data.
#loads MQ output from proteomic analysis of 60 LCL lines, subsets to the three of interest, median then quantile normalizes. Returned is a list of 4 DFs: MQoutput heavy,MQ output just data, median normalized, and quantile normalized. 

# Choose directory containing proteomics data to pass to 'NormProt'
CorrectedDataProt <- NormProt(directory = "E:/My Documents/Pilot/EnsemblDBProteome/txt/")#from home
CorrectedDataProt <- NormProt(directory = "D:/EnsemblDBProteome/iBAQ proteome/")#from laptop INCLUDES IBAQ
ProtQuantiled <- CorrectedDataProt[[4]] #Median and quantile normalized inverted (L/H) protein ratios (with MQ normalization as well).
ProteinZia <- CorrectedDataProt[[1]]#Proteins from 60 human LCLs with no contaminants, reverse hits, or non-quantified IDs (6421)

#ProtAssignment matches the two datasets. It returns a DF with the normalized protein L/H values and majority ids appended to the ME DF. It also returns a protein normalized data frame along with EDA plots corresponding to the batch corrected and normalized phospho dataframe that was passed - "phosphonorm".

##read in the proteome fasta file. Here Ensembl cCDS
proteome <- read.fasta( file = "./FASTA//Homo_sapiens.GRCh37.75.pep.all.parsedCCDS.fa", seqtype = "AA", as.string = TRUE)#updataed to local machine. FASTA file used for search

NormalizedResults <- ProtAssignment2(proteinfull = ProteinZia, proteinnorm = ProtQuantiled, multExpanded1_withDE = multExpanded1_withDE, phosphonorm=adata, proteome)#pass com2 perhaps
multExpanded1_withDE <- NormalizedResults[[1]]
ProtNormalized <- NormalizedResults[[2]]#protein subtracted phospho dataframe

#This function performs diffphos analysis on phosphodata (using normalized protein data as a covariate in progress. will add a flag to function call).
#Accepts ME dataframe with confounded data annotation and normalized protein values.
#Returns proteinnormalizedDE annotated ME dataframe and produces diffphos plots derived from protein corrected data.
multExpanded1_withDE <- DiffPhosProt(multExpanded1_withDE, phosphonorm = pilot)#note REV_entries (5) already removed from pilot! See confounded diffphos section
################


#Nested random effects model. 'NestedVar' performs analysis and I would like this to also return an annotated ME DF with each phosphosites annotated for downstream enrichment analysis. One example of a category might be high biological var/low tech and individual. I think a combinatorial categorization would work well here.
###########################
# remove exp obs if not observed at least two times in each sample
com3 <- com2[rowSums(is.na(com2[ , 1:4])) <= 2 & rowSums(is.na(com2[ , 5:8])) <= 2 & rowSums(is.na(com2[ , 9:12])) <= 2,]

#send unbalanced data to NestedVar. Here a nested random effect model is fitted for each phosphopeptide. The peptide model variance components are returned. 
varcomp <- NestedVar(ratios=com3, balanced = F)
varcomp <- as.data.frame(varcomp)
#send to VarComp. Note this is the same dataframe as ProtNormalized. Com2 used next!!
ProtNormalizedVar <- ProtNormalized[rowSums(is.na(ProtNormalized[ , 1:4])) <= 2 & rowSums(is.na(ProtNormalized[ , 5:8])) <= 2 
                                    & rowSums(is.na(ProtNormalized[ , 9:12])) <= 2,]
 
VarcompProt <- NestedVar(ratios=ProtNormalizedVar, balanced = F)#same result as the confounded data. Perhaps can run with protein as a covariate. 
VarcompProt <- as.data.frame(VarcompProt)
# Is this signature unique to phospho? What do the small number of protein estimates show?



#add the new categorizations from the varcompdata to the multexpanded DF for enrichment analyses. SOME DESCREPANCY BETWEEN IDMULT AND ROWNAMES!! 11 MISSING
multExpanded1_withDE$SubtoVarcomp <- ifelse(multExpanded1_withDE$idmult %in% row.names(varcomp), "+", "-")
multExpanded1_withDE$HighIndVar <- ifelse(multExpanded1_withDE$idmult %in% row.names(varcomp[varcomp$high_ind_var=="+",]), "+", "-")
multExpanded1_withDE$LowIndVar <- ifelse(multExpanded1_withDE$idmult %in% row.names(varcomp[varcomp$low_ind_var=="+",]), "+", "-")
multExpanded1_withDE$HighBioVar <- ifelse(multExpanded1_withDE$idmult %in% row.names(varcomp[varcomp$high_bio_var=="+",]), "+", "-")
multExpanded1_withDE$LowBioVar <- ifelse(multExpanded1_withDE$idmult %in% row.names(varcomp[varcomp$low_bio_var=="+",]), "+", "-")

multExpanded1_withDE$ppSubtoVarcomp <- ifelse(multExpanded1_withDE$idmult %in% row.names(VarcompProt), "+", "-")#no descrepancy with pp numbers
multExpanded1_withDE$pnHighIndVar <- ifelse(multExpanded1_withDE$idmult %in% row.names(VarcompProt[VarcompProt$high_ind_var=="+",]), "+", "-")
multExpanded1_withDE$pnLowIndVar <- ifelse(multExpanded1_withDE$idmult %in% row.names(VarcompProt[VarcompProt$low_ind_var=="+",]), "+", "-")
multExpanded1_withDE$pnHighBioVar <- ifelse(multExpanded1_withDE$idmult %in% row.names(VarcompProt[VarcompProt$high_bio_var=="+",]), "+", "-")
multExpanded1_withDE$pnLowBioVar <- ifelse(multExpanded1_withDE$idmult %in% row.names(VarcompProt[VarcompProt$low_bio_var=="+",]), "+", "-")

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
varcompMQnormonly <- NestedVar(ratios=RawRatios)

#feed the same unnormalized data by batch into the model
RawRatiosB1 <- RawRatios[,c(1:2,5:6,9:10)]
RawRatiosB2 <- RawRatios[,c(3:4,7:8,11:12)]
  
varcompMQnormonlyB1 <- NestedVar(ratios=RawRatiosB1, batch=T)
varcompMQnormonlyB2 <- NestedVar(ratios=RawRatiosB2, batch=T)

#test #4. Is signature there using MQ non-normalized ratio values and separate batches.
#loading unnormalized values
phosphoRaw <- load.MQ2(directory = "D:/10_9_14/txt/", type = "phospho")
phosphoRaw <- phosphoRaw[(phosphoRaw$Potential.contaminant != "+" & phosphoRaw$Reverse != "+"),]

# subset phospho to class 1
phosphoRaw1 <- phosphoRaw[(phosphoRaw$Localization.prob >= .75),]

#for non normalized data
multExpandedRaw1 <- ExpandPhos2(phospho=phosphoRaw1)

#a curious inversion of the outlier. Here is the code from BatchNorm to extract dataframe for nestedvar

expCol <- grep("HL(.*)", colnames(multExpandedRaw1))
data <- multExpandedRaw1[,expCol]

# add row names with site id and multiplicity designation 
# row.names(data) <- multExpanded$id#note this won't work because of the multiplicity issue
idmult <- paste(multExpandedRaw1$id, multExpandedRaw1$multiplicity, sep="_")
row.names(data) <- idmult
data <- log2(data)

#summaries
summary(data)
boxplot(data[,1:12])#extreme outlier in HL18486_1_2
boxplot(data[,1:12], ylim=c(-4,4))

#remove outlier which is 3X smaller than any other datapoint but was 3X greater for the normalized values
data[,2][which.min(data[,2])] <- NaN
boxplot(data[,1:12])#extreme outlier in HL18486_1_2 removed

TotallyRawRatios <- data
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


#############################

#Add GOID, Reactome, Entrez, HGNCID, HGNC symbol, and HGNC derived description of each protein gene annotation to multExpanded DF
multExpanded1_withDE <- AddAnnotation(multExpanded1_withDE)


# (Ontology) Enrichment analysis performed on diffphos omnibus F significant and enrichment for each of the four combinations (high/low ind/bio) of variance component estimates. 
enrichment_tables <- Enrichment(multExpanded1_withDE)
#NOTE THE STRANGE REACTOME ISSUE FOR THE CONFOUNDED DATA. annotation added but no corresponding information for that id.

#copynumber estimates using perseus total protein algorithm (can't use protein ruler approach with ensembl database). 
#Here copynumbers across samples and experimental approaches are compared (correlations) and bias for protein expression level is explored within the variance component clusters and diffphos for confounded and non-confounded samples. A resampling based approach may be needed here because of multiple phosphorylation sites/protein requires that I threshold proteins into high/med/low expression (or similar) categories. This function will be self contained but MAY add cn estimates to multexpanded table. 








#next is work at the genome level to explicitly show that genetic variation is driving these changes. nonsynSNPs, pQTLs, nonsynSNPs surrounding the phosphosite, etc


# Absoulte protein concentration estimates (iBAQ or 'protein ruler'), GO, reactome, corum, phosphositeplus, nonsynSNPs, pQTLs, within motif nonsynsnps

#motif based analysis: enrichment of kinase/PBD substrates/targets. motif description across all contrasts

#Perhaps GSEA type approaches within limma


#something at the level of the differential network...





