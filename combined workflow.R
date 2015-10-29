rm(list=ls(all=TRUE)) #start with empty workspace


# load required libraries and scripts -----------------------------------
library(reshape2)
library(stringr)
library(plyr)
library(seqinr)
library(vioplot)
source("loadMQ.R")
source("ExpandPhos.R")
source("ExpandPhos2.R") #for expanding the non-normalized data
source("counts.R")
source("breakdown.R")
source("BatchNorm.R")
source("DiffPhos.R")
source("loadMQ2.R") #non normalized data
source("NestedVar.R")
source("NormProt.R")
source("ProtAssignment.R")
source("ProtAssignment2.R")
source("AddAnnotation.R")
source("Enrichment.R")
source("PerseusOut.R")
source("AddAnnotation.R")


#load, reformat and characterize MQ outputted mass spectrometry  --------
# load phospho and protein files with particular variables populated using "loadMQ" on cluster
phospho <- load.MQ(directory = "/mnt/lustre/home/bengelmann/MQfiles/EnsemblDBPhospho/", type = "phospho")
# protein <- load.MQ(directory = "D:/EnsemblDBPhospho/PilotPhosphoensemblDB/combined/txt/", type = "protein")
protein <- load.MQ(directory = "/mnt/lustre/home/bengelmann/MQfiles/EnsemblDBPhospho/", type = "protein")#file has ibaq values

# load phospho and protein files with particular variables populated using "loadMQ" on laptop
phospho <- load.MQ(directory = "D:/EnsemblDBPhospho/PilotPhosphoensemblDB/combined/txt/", type = "phospho")
# protein <- load.MQ(directory = "D:/EnsemblDBPhospho/PilotPhosphoensemblDB/combined/txt/", type = "protein")
protein <- load.MQ(directory = "D:/EnsemblDBPhospho/iBAQ quantification/", type = "protein")#with ibaq

# load phospho and protein files with particular variables populated using "loadMQ" at home
phospho <- load.MQ(directory = "E:/My Documents/Pilot/EnsemblDBPhospho/PilotPhosphoensemblDB/combined/txt/", type = "phospho")
# protein <- load.MQ(directory = "E:/My Documents/Pilot/EnsemblDBPhospho/PilotPhosphoensemblDB/combined/txt/", type = "protein")
protein <- load.MQ(directory = "E:/My Documents/Pilot/EnsemblDBPhospho/iBAQ quantification/", type = "protein")#with ibaq

# remove contaminants and reverse database hits
phospho <- phospho[(phospho$Potential.contaminant != "+" & phospho$Reverse != "+"),]

#write out phospho for publication
phospho.table <- phospho[, 49:71]
drop <- c("Reverse", "Potential.contaminant", "Positions.within.proteins")
phospho.table <- phospho.table[ , !names(phospho.table) %in% drop]
write.table(phospho.table, file = "./phospho.table.txt", sep = "\t", row.names = F)
rm(phospho.table)

protein <- protein[(protein$Potential.contaminant != "+" & protein$Reverse != "+"),]

# subset phospho to class 1
phospho1 <- phospho[(phospho$Localization.prob >= .75), ]



# subset protein to remove those identifications observed solely by modification site
# "only identified by site" hits CAN BE removed because they tend to have lower PEPs (wouldn't pass the FDR TH anyway) and can't be quantified since they are not idd by non-modified peptides. 
# Note there are some high probability proteins here given some proteins are idd by 20+ phosphopeptides.
# eg is A6NKT7 (PEP = 2.23E-70)
protein1 <- protein[(protein$Only.identified.by.site != "+"), ]

#expand phosphorylation data for observation (multiplicity based analysis). 
# All observations quantified in >=1 experiment.
multExpanded <- ExpandPhos(phospho)

# Class 1 sites with each source of quantification for that site (singly/doubly/3+) explicitly accounted for 
multExpanded1 <- ExpandPhos(phospho1)


# summary statistics of Mass Spec data ------------------------------------

#make tables of basic counts of proteins and phosphopeptides (may have to update when performing the normalization)
phoscount(phospho,phospho1,multExpanded,multExpanded1)
proteincount(protein)

# make breakdown charts of phospho and protein overlap. This function produces barplots of number of sites/obs and proteins per experiment and
# cumulative over experiments. It also has extensive phospo info, including number of phospho per protein, multiplicity breakdown, and venns
breakdown(protein, phospho, multExpanded, cls=F)
breakdown(protein, phospho1, multExpanded1)

# Remove peptides mapping to protein groups with any reverse hits ---------

#Perhaps I should just remove the reverse hits and produce a flag?

#!!BEFORE passing a matrix for diffPhos, remove all ids with leading protein mapping to a REV_ hit (31 for ME1). These ids are still present because there are majority protein ids for these peptides that map to a non reverse hit. These are omitted to make protein mapping to Zia's dataset more straightforward and enrichment analysis results comparable between the confounded and non-confounded datasets. That is, what does it mean to compare enrichment for a reversed 'protein'?

#add idmult annotation to multexpanded table
idmult <- paste(multExpanded1$id, multExpanded1$multiplicity, sep="_")
multExpanded1 <- cbind(multExpanded1,idmult)
#remove reverse hits and acquire an index to subset corrected data
RevHits <- grep(multExpanded1$Protein, pattern = "REV")
RevHitsidmult <- multExpanded1$idmult[RevHits]
multExpanded1 <- multExpanded1[!grepl(multExpanded1$Protein, pattern = "REV"),] #18238 now 17774

# Normalization and batch correction of confounded data -------------------

#remove an outlier, normalize (median and quantile), and batch correct (combat). Returned are EDA plots and a list of DFs.See BatchNorm for details
if(!file.exists("./CorrectedData.rds")) {
CorrectedData <- BatchNorm(multExpanded1=multExpanded1)#class 1 sites
}else{
  CorrectedData <- readRDS("./CorrectedData.rds")
}
com2 <- CorrectedData[[8]]#normalized/batch corrected (using ComBat) data frame
adata <- CorrectedData[[9]]#normalized/batch corrected data frame with at lesat 1 obs in each bio rep
pilot <- CorrectedData[[10]]#same as above with mean ratios for each bio replicate


# protein level assignment and normalization and  using phosprep and gelprep data --------
##read in the proteome fasta file that was used for search. Here Ensembl cCDS. This will be used for protein assignment.
if(!file.exists("./FASTA/Homo_sapiens.GRCh37.75.pep.all.parsedCCDS.fa")){
  source("./FASTA/ensembl FASTA for DB search.R")#produce file for db search
}
proteome <- read.fasta( file = "./FASTA/Homo_sapiens.GRCh37.75.pep.all.parsedCCDS.fa", seqtype = "AA", as.string = TRUE)


#### Add protein level information from PhosPrep workup
#Protein assignment adds protein ids, positions within protein, H/L values and ibaq values to phosphosites for protein level normalization using the "protein groups" file produced from the phospho workup.

#note that this function performs protein level batch correction/normalization and assignment to phosphopeptides.

if(!file.exists("./PhosPrepMatrices.rds")) {
  PhosPrepMatrices <- ProtAssignment(protein, proteome, multExpanded1)
}else{
  PhosPrepMatrices <- readRDS("./PhosPrepMatrices.rds")
}
multExpanded1 <- PhosPrepMatrices[[1]]
PhosPrepquantiledBio <- PhosPrepMatrices[[9]]#med/norm with one measurement in each bio rep
PhosPrepCombat <- PhosPrepMatrices[[10]]#normalized/batch corrected (using ComBat) data frame
PhosPrepCombatBio <- PhosPrepMatrices[[11]]#normalized/batch corrected data frame with at lesat 1 obs in each bio rep
PhosPrepCombatPilot <- PhosPrepMatrices[[12]]#same as above with mean ratios for each bio replicate


#### Add protein level information from GelPrep workup

#For GelPrep protein normalization and assignment of proteins to phosphopeptides is done in two steps. 

# Normprot performs 'traditional' quantile normalization AND calls multiple scripts to perform a pQTL optimization routing by regressing PCs.
# See the PCregression folder for more details

# Choose directory containing proteomics data to pass to 'NormProt'.
#note that this function now requires "genotypes_imputed" files and ensembl_4_30.gz from Jack to run PCregression pQTL optimization anew. 
CorrectedDataProt <- NormProt(directory = "/mnt/lustre/home/bengelmann/MQfiles/EnsemblDBProteome/")#on cluster INCLUDES IBAQ
CorrectedDataProt <- NormProt(directory = "E:/My Documents/Pilot/EnsemblDBProteome/iBAQ proteome/")#from home INCLUDES IBAQ
CorrectedDataProt <- NormProt(directory = "D:/EnsemblDBProteome/iBAQ proteome/")#from laptop INCLUDES IBAQ
ProtQuantiled <- CorrectedDataProt[[4]] #Median and quantile normalized inverted (L/H) protein ratios (with MQ normalization as well).
ProteinZia <- CorrectedDataProt[[1]]#Proteins from 60 human LCLs with no contaminants, reverse hits, or non-quantified IDs (6421)
RegressedCommon <- CorrectedDataProt[[6]]#Gelprep with 13 PCs regressed.

#violin plots of number of peptides/protein and sequence coverage for the proteins present in all three lines (for SF4 A)
GelProtSummary <- ProteinZia[ProteinZia$id %in% row.names(RegressedCommon), c("Razor...unique.peptides", "Sequence.coverage....")]

pdf("peptides_per_prot.pdf", 6, 5)
plot(1, 1, ylim = c(0,300), type = 'n', xaxt = 'n', ylab = "Peptides assigned to protien group", xlab = "Protein Groups \n (n = 3885)", 
     family = "serif")
vioplot(GelProtSummary$Razor...unique.peptides, rectCol="gray", wex = .5, col = "cornflowerblue", add = T)
dev.off()

pdf("sequence_coverage_per_prot.pdf", 6, 5)
plot(1, 1, ylim = c(0,100), type = 'n', xaxt = 'n', ylab = "Protein sequence coverage (%)", xlab = "Protein Groups \n (n = 3885)", 
     family = "serif")
vioplot(GelProtSummary$Sequence.coverage...., rectCol="gray", wex = .5, col = "firebrick3", add = T)
dev.off()
median(GelProtSummary$Razor...unique.peptides)#18
median(GelProtSummary$Sequence.coverage....)#43.8
rm(GelProtSummary)

###!!! Ensure RegressedCommon contains row.names = proteingroup ids. If not delete "./PCregression/PRO_raw.RData" and call NormProt again.


#ProtAssignment2 matches the two datasets. It returns a DF with the PCregressed protein L/H values and majority ids appended to the ME DF. It also returns a protein normalized data frame along with EDA plots corresponding to the batch corrected and normalized phospho dataframe that was passed - "phosphonorm".

#for the moment I am using ME1 as ME with DE. Sorry for this future self, I must move forward. Also note this takes awhile!!!
if(!file.exists("./NormalizedResults.rds")){
NormalizedResults <- ProtAssignment2(proteinfull = ProteinZia, proteinnorm = RegressedCommon, multExpanded1_withDE = multExpanded1, phosphonorm=adata, proteome)
} else {
  NormalizedResults <- readRDS("./NormalizedResults.rds")
}
multExpanded1 <- NormalizedResults[[1]]
ProtNormalized <- NormalizedResults[[2]]#protein subtracted phospho dataframe
GelPrep <- NormalizedResults[[3]]#Protein level

# Differential phosphorylation on confounded and protein normalized data (2 ways) --------

# Differential phosphorylation analysis using DiffPhos function on confounded data. 
# Below is limma based DE across contrasts with annotation added to the passed multExpanded1 file. multiple images/venns are returned. 'multExpandedwithDE' is returned.

# A new combined approach. DiffPhos on confounded and 'corrected' (via covariate or normalization) data.
if(!file.exists("./multExpanded1_withDE.rds")){
  multExpanded1_withDE <- DiffPhos(phosdata = adata, PhosPrep = PhosPrepCombatBio, GelPrep = GelPrep, multExpanded1)
}else{
  multExpanded1_withDE <- readRDS("./multExpanded1_withDE.rds")
}

#still need comparison venns and perhaps discussion/analysis of these results. Power differences with various covariates need to be considered.


# Nested Random Effects modeling ------------------------------------------

# 'NestedVar' performs analysis and I would like this to also return an annotated ME DF 
# with each phosphosites annotated for downstream enrichment analysis. One example of 
# a category might be high biological var/low tech and individual. I think a combinatorial categorization would work well here.

# remove exp obs if not observed at least two times in each sample
com3 <- com2[rowSums(is.na(com2[ , 1:4])) <= 2 & rowSums(is.na(com2[ , 5:8])) <= 2 & rowSums(is.na(com2[ , 9:12])) <= 2,]

# send unbalanced data to NestedVar. Here a nested random effect model is fitted for each phosphopeptide. The peptide model variance components are returned. 

varcomp <- NestedVar(ratios = ProtNormalized, noMissing = F)
varcomp <- as.data.frame(varcomp)

# assign flag to varcomp file according to four categories. high ind/high biological. high ind/low bio. low ind/high bio and low ind/low bio
# cutoff for high/low individual is log10 = -5. cutoff for high/low biological variance is log10 = -6.
varcomp$high_ind_var <- ifelse(log10(varcomp$individual) >= -5, "+", "-")
varcomp$low_ind_var <- ifelse(log10(varcomp$individual) < -5, "+", "-")
varcomp$high_bio_var <- ifelse(log10(varcomp$biorep) >= -6, "+", "-")
varcomp$low_bio_var <- ifelse(log10(varcomp$biorep) < -6, "+", "-")

# Send to VarComp. Note this is the same dataframe as ProtNormalized. Com2 used next!!
ProtNormalizedVar <- ProtNormalized[rowSums(is.na(ProtNormalized[ , 1:4])) <= 2 & rowSums(is.na(ProtNormalized[ , 5:8])) <= 2 
                                    & rowSums(is.na(ProtNormalized[ , 9:12])) <= 2,]#3257
 
VarcompProt <- NestedVar(ratios=ProtNormalizedVar, noMissing = F)#same result as the confounded data. Perhaps can run with protein as a covariate. 
VarcompProt <- as.data.frame(VarcompProt)

#assign flag to VarcompProt file according to four categories. high ind/high biological. high ind/low bio. low ind/high bio and low ind/low bio
#cutoff for high/low individual is log10 = -5. cutoff for high/low biological variance is log10 = -6.
VarcompProt$high_ind_var <- ifelse(log10(VarcompProt$individual) >= -5, "+", "-")
VarcompProt$low_ind_var <- ifelse(log10(VarcompProt$individual) < -5, "+", "-")
VarcompProt$high_bio_var <- ifelse(log10(VarcompProt$biorep) >= -6, "+", "-")
VarcompProt$low_bio_var <- ifelse(log10(VarcompProt$biorep) < -6, "+", "-")

# Is this signature unique to phospho? What do the small number of protein estimates show? They also show this signature

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


#######Enrichment analyses. quite a few of them ----

#Add GOID, Reactome, Entrez, HGNCID, HGNC symbol, and HGNC derived description of each protein gene annotation to multExpanded DF
if(!file.exists("multExpanded1_withDE_annotated.rds")){
multExpanded1_withDE_annotated <- AddAnnotation(multExpanded1_withDE)
}else{
  multExpanded1_withDE_annotated <- readRDS("multExpanded1_withDE_annotated.rds")
}

# add ELM motif instances ----
source("./ELM/elm_processing.R")

#combine with hprd motifs. processed using perseus on 9.19.15
motifs <- read.table("./Perseus/motif_annotations.txt", sep = "\t", header = T, stringsAsFactors = F)
motifs <- motifs[[1]]
multExpanded1_withDE_annotated <- cbind(multExpanded1_withDE_annotated, motifs)
saveRDS(multExpanded1_withDE_annotated, file = "./multExpanded1_withDE_annotated.rds")


#protein level enrichments



#peptide level enrichments



# snp enrichments






# pqtl enrichment








# (Ontology) Enrichment analysis performed on diffphos omnibus F significant and enrichment for each of the four combinations (high/low ind/bio) of variance component estimates. 
enrichment_tables <- Enrichment(multExpanded1_withDE_annotated)
#NOTE THE STRANGE REACTOME ISSUE FOR THE CONFOUNDED DATA. annotation added but no corresponding information for that id.

#copynumber estimates using perseus total protein algorithm (can't use protein ruler approach with ensembl database). 
#Here copynumbers across samples and experimental approaches are compared    (correlations) and bias for protein expression level is explored within the variance component clusters and diffphos for confounded and non-confounded samples. A resampling based approach may be needed here because of multiple phosphorylation sites/protein requires that I threshold proteins into high/med/low expression (or similar) categories. This function will be self contained but MAY add cn estimates to multexpanded table. 


#next is work at the genome level to explicitly show that genetic variation is driving these changes. nonsynSNPs, pQTLs, nonsynSNPs surrounding the phosphosite, etc

# working with "SNPenrich.R" here. Will report a small table with enrichment results as outlined above as well as any resampling results/charts if needed

results <- SNPenrich(multExpanded1_withDE)


# Absoulte protein concentration estimates (iBAQ or 'protein ruler'), GO, reactome, corum, phosphositeplus, nonsynSNPs, pQTLs, within motif nonsynsnps

#motif based analysis: enrichment of kinase/PBD substrates/targets. motif description across all contrasts

#Perhaps GSEA type approaches within limma


#something at the level of the differential network...





