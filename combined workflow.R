rm(list=ls(all=TRUE)) #start with empty workspace


# load required libraries and scripts -----------------------------------
library(reshape2)
library(stringr)
library(plyr)
library(seqinr)
library(vioplot)
source("loadMQ.R")
source("ExpandPhos.R")
source("counts.R")
source("breakdown.R")
source("BatchNorm.R")
source("DiffPhos.R")
source("NestedVar.R")
source("NormProt.R")
source("ProtAssignment.R")
source("ProtAssignment2.R")
source("AddAnnotation.R")
source("Enrichment.R")
source("PerseusOut.R")


#load, reformat and characterize MQ outputted mass spectrometry  --------
# load phospho and protein files with particular variables populated using "loadMQ" on cluster
phospho <- load.MQ(directory = "/mnt/lustre/home/bengelmann/MQfiles/EnsemblDBPhospho/", type = "phospho")
protein <- load.MQ(directory = "/mnt/lustre/home/bengelmann/MQfiles/EnsemblDBPhospho/", type = "protein")#file has ibaq values

# load phospho and protein files with particular variables populated using "loadMQ" on laptop
phospho <- load.MQ(directory = "D:/EnsemblDBPhospho/PilotPhosphoensemblDB/combined/txt/", type = "phospho")
protein <- load.MQ(directory = "D:/EnsemblDBPhospho/iBAQ quantification/", type = "protein")#with ibaq

# load phospho and protein files with particular variables populated using "loadMQ" at home
phospho <- load.MQ(directory = "E:/My Documents/Pilot/EnsemblDBPhospho/PilotPhosphoensemblDB/combined/txt/", type = "phospho")
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
protein1 <- protein[(protein$Only.identified.by.site != "+"), ]

#expand phosphorylation data for observation (multiplicity based analysis). 
# All observations quantified in >=1 experiment.
multExpanded <- ExpandPhos(phospho)

# Class 1 sites with each source of quantification for that site (singly/doubly/3+) explicitly accounted for 
multExpanded1 <- ExpandPhos(phospho1)


# summary statistics of Mass Spec data ------------------------------------

#make tables of basic counts of proteins and phosphopeptides
phoscount(phospho,phospho1,multExpanded,multExpanded1)
proteincount(protein)

# make breakdown charts of phospho and protein overlap. This function produces barplots of number of sites/obs and proteins per experiment and cumulative over experiments. It also has extensive phospo info, including number of phospho per protein, multiplicity breakdown, and venns
breakdown(protein, phospho, multExpanded, cls=F)
breakdown(protein, phospho1, multExpanded1)

# Remove peptides mapping to protein groups with any reverse hits ---------

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
phosdata <- CorrectedData[[9]]#normalized/batch corrected data frame with at lesat 1 obs in each bio rep
pilot <- CorrectedData[[10]]#same as above with mean ratios for each bio replicate
medianSub.quantiled <- CorrectedData[[3]]
medianSub.quantiled.two.per.batch <- CorrectedData[[6]]

# protein level assignment and normalization with phosprep and gelprep data --------
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

# Normprot performs quantile normalization AND calls multiple scripts to perform a pQTL optimization routine by regressing PCs.
# See the PCregression folder for more details

# Choose directory containing proteomics data to pass to 'NormProt'.
#note that this function now requires "genotypes_imputed" files and ensembl_4_30.gz from Jack to run PCregression pQTL optimization anew. 
CorrectedDataProt <- NormProt(directory = "/mnt/lustre/home/bengelmann/MQfiles/EnsemblDBProteome/")#on cluster INCLUDES IBAQ
CorrectedDataProt <- NormProt(directory = "E:/My Documents/Pilot/EnsemblDBProteome/iBAQ proteome/")#from home INCLUDES IBAQ
CorrectedDataProt <- NormProt(directory = "D:/EnsemblDBProteome/iBAQ proteome/")#from laptop INCLUDES IBAQ
ProtQuantiled <- CorrectedDataProt[[4]] #Median and quantile normalized inverted (L/H) protein ratios (with MQ normalization as well).
ProteinZia <- CorrectedDataProt[[1]]#Proteins from 60 human LCLs with no contaminants, reverse hits, or non-quantified IDs (6421)
RegressedCommon <- CorrectedDataProt[[6]]#Gelprep with 13 PCs regressed. 



###########Protein summary stats ---------

#GelProt output
#protein table outputs. Keep columns 244:length and remove unwanted manually. Only keep proteingroups that were PC corrected.
ProteinZia.out <- ProteinZia[ProteinZia$id %in% row.names(RegressedCommon), 244:length(ProteinZia)]
write.table(ProteinZia.out, file = "./GelPrep.table.txt", sep = "\t", row.names = F)
rm(ProteinZia.out)

#PhosProt output
#use protein1 dataframe so that protiens identified only with a site that was modified aren't included
Protein1.out <- protein1[ , 52:length(protein1)]
write.table(Protein1.out, file = "./PhosPrep.table.txt", sep = "\t", row.names = F)
rm(Protein1.out)

#violin plots of number of peptides/protein and sequence coverage for the proteins present in all three lines
PhosProtSummary <- protein1[rowSums(is.na(protein1[ , 1:4])) <= 3 & rowSums(is.na(protein1[ , 5:8])) <= 3 & 
                              rowSums(is.na(protein1[ , 9:12])) <= 3, c("Razor...unique.peptides", "Sequence.coverage....")]
GelProtSummary <- ProteinZia[ProteinZia$id %in% row.names(RegressedCommon), c("Razor...unique.peptides", "Sequence.coverage....")]

#combine into one dataframe...
library(gdata)
names(PhosProtSummary) <- paste0(names(PhosProtSummary), "PhosProt")
names(GelProtSummary) <- paste0(names(GelProtSummary), "GelProt")
proteomes <- cbindX(PhosProtSummary, GelProtSummary)

#make a 'number of peptides and sequence coverage dataframe
PeptideCount <- data.frame( peptide_count = c(PhosProtSummary$Razor...unique.peptidesPhosProt, 
                                              GelProtSummary$Razor...unique.peptidesGelProt),
                            peptide_source = factor( rep( c(1:2), c(dim(PhosProtSummary)[1], dim(GelProtSummary)[1])),
                                                  labels = c("SCX-TiO2", "SDS-PAGE")) )


SequenceCoverage <- data.frame( sequence_coverage = c(PhosProtSummary$Sequence.coverage....PhosProt, 
                                                      GelProtSummary$Sequence.coverage....GelProt),
                                sequence_source = factor( rep( c(1:2), c(dim(PhosProtSummary)[1], dim(GelProtSummary)[1])),
                                                          labels = c("SCX-TiO2", "SDS-PAGE")) )

require(ggplot2)
require(extrafont)


p <- ggplot(PeptideCount,
            aes(x = peptide_source, y = peptide_count ))  +
  geom_violin(alpha = .6) + geom_boxplot(width = .2, alpha = .2, outlier.size = 0) + 
  ylab("Peptides assigned to protein") + xlab("Peptide source") + 
  scale_y_continuous(limits = quantile(GelProtSummary$Razor...unique.peptidesGelProt, c(0,0.95))) + #for clarity
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text=element_text(family="Times New Roman", size=12))   
print(p)
ggsave("peptides_per_group.pdf", p, width = 4, height = 5.5)

  
p <- ggplot(SequenceCoverage,
            aes(x = sequence_source, y = sequence_coverage ))  +
  geom_violin(alpha = .6) + geom_boxplot(width = .25, alpha = .2, outlier.size = 0) + 
  ylab("Protein sequence coverage (%)") + xlab("Peptide source") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text=element_text(family="Times New Roman", size=12))   
print(p)
ggsave("sequence_coverage.pdf", p, width = 4, height = 5.5)


#ProtAssignment2 matches the two datasets. It returns a DF with the PCregressed protein L/H values and majority ids appended to the ME DF. It also returns a protein normalized data frame along with EDA plots corresponding to the batch corrected and normalized phospho dataframe that was passed - "phosphonorm".

#For the moment I am using ME1 as ME with DE.
if(!file.exists("./NormalizedResults.rds")){
NormalizedResults <- ProtAssignment2(proteinfull = ProteinZia, proteinnorm = RegressedCommon, multExpanded1_withDE = multExpanded1, phosphonorm=phosdata, proteome)
} else {
  NormalizedResults <- readRDS("./NormalizedResults.rds")
}
multExpanded1 <- NormalizedResults[[1]]
ProtNormalized <- NormalizedResults[[2]]#protein subtracted phospho dataframe
GelPrep <- NormalizedResults[[3]]#Protein level

# Differential phosphorylation on confounded and protein normalized (run with protein as covariate) data --------

# Below is limma based DE across contrasts with annotation added to the passed multExpanded1 file. multiple images/venns are returned.
if(!file.exists("./multExpanded1_withDE.rds")){
  multExpanded1_withDE <- DiffPhos(phosdata = phosdata, PhosPrep = PhosPrepCombatBio, GelPrep = GelPrep, multExpanded1)
}else{
  multExpanded1_withDE <- readRDS("./multExpanded1_withDE.rds")
}

# Nested Random Effects modeling ------------------------------------------

#A bayesian approach with a naive prior is used to estimate variance components within each hierarchical level of the experimental design:

# Confounded data 
# remove exp obs if not observed at least two times in each sample
confounded.batch.corrected <- com2[rowSums(is.na(com2[ , 1:4])) <= 2 & rowSums(is.na(com2[ , 5:8])) <= 2 & rowSums(is.na(com2[ , 9:12])) <= 2,]

# The peptide model variance components are returned. 
mcmcVarcomp.confounded <- NestedVar(ratios = confounded.batch.corrected, noMissing = TRUE)
colnames(mcmcVarcomp.confounded) <- c("individual","culture","residual")

# Confounded without batch effect correction
confounded.batch.raw <- medianSub.quantiled
mcmcVarcomp.confounded.batchfit <- NestedVar(ratios = confounded.batch.raw, noMissing = TRUE, BatchCorrected = FALSE)
colnames(mcmcVarcomp.confounded.batchfit) <- c("individual","culture","residual")


#Protein as a covariate data.
colnames(GelPrep) <- c("HL18862", "HL18486", "HL19160")
PhosProtGel <- merge(phosdata, GelPrep, by = "row.names", 
                     suffixes = c("_peptide", "_GelPrep") ) #3257 observations
rownames(PhosProtGel) <- PhosProtGel$Row.names
PhosProtGel <- PhosProtGel[ , -1]
PhosProtGel <- as.matrix(PhosProtGel)
mcmcVarcomp.proteinCov <- NestedVar(PhosProtGel, includeProteinCovariate = TRUE)
colnames(mcmcVarcomp.proteinCov) <- c("individual","culture","residual")

#Protein Normalized data. 
colnames(ProtNormalized) <- colnames(confounded.batch.corrected)
mcmcVarcomp.proteinNorm <- NestedVar(ProtNormalized, includeProteinCovariate = FALSE, noMissing = TRUE)
colnames(mcmcVarcomp.proteinNorm) <- c("individual","culture","residual")

#Gelprot protein as a covariate with non-be corrected phospho data
colnames(GelPrep) <- c("HL18862", "HL18486", "HL19160")
PhosProtGelBatch <- merge(medianSub.quantiled, GelPrep, by = "row.names", 
                     suffixes = c("_peptide", "_GelPrep") ) #3257 observations after omitting missing values. 
rownames(PhosProtGelBatch) <- PhosProtGelBatch$Row.names
PhosProtGelBatch <- PhosProtGelBatch[ , -1]
PhosProtGelBatch <- as.matrix(PhosProtGelBatch)
mcmcVarcomp.proteinCov.Batch <- NestedVar(PhosProtGelBatch, includeProteinCovariate = TRUE, BatchCorrected = FALSE)
colnames(mcmcVarcomp.proteinCov.Batch) <- c("individual","culture","residual")

#PhosPrep protein as a covariate with non be corrected phospho dta and protein estimates
PhosProtPhosBatch <- merge(medianSub.quantiled, PhosPrepquantiledBio, by = "row.names", 
                  suffixes = c("_peptide", "_PhosPrep") ) #1308 observations
rownames(PhosProtPhosBatch) <- PhosProtPhosBatch$Row.names
PhosProtPhosBatch <- PhosProtPhosBatch[ , -1]
PhosProtPhosBatch <- as.matrix(PhosProtPhosBatch)
mcmcVarcomp.PhosProteinCov.Batch <- NestedVar(PhosProtPhosBatch, includeProteinCovariate = TRUE, BatchCorrected = FALSE, PhosPrep = TRUE)
colnames(mcmcVarcomp.PhosProteinCov.Batch) <- c("individual","culture","residual")


##------ NRE fit Results ------

#absolute and standardized variance component plots for protein covariate and for confounded.

#combine confounded and protein as a covariate results in to a list of data frames for common processing.

results <- list(Confounded = mcmcVarcomp.confounded, #confounded after combat batch effect correction
                Confounded.Batchfit = mcmcVarcomp.confounded.batchfit, #med/quantile norm confounded with batch as a covariate
                ProteinCorrected = mcmcVarcomp.proteinCov, # pQTL regressed gelprot protein dataframe as a covariate with combat corrected phospho
                ProteinNormalized = mcmcVarcomp.proteinNorm, # pQTL regressed gelprot protein dataframe subtracted from combat corrected phospho
                ProteinCovariateBatch = mcmcVarcomp.proteinCov.Batch, # pQTL regressed gelprot protein dataframe and batch as covariates with median/quantile normalized phospho
                ProteinCovariateBatchPhos = mcmcVarcomp.PhosProteinCov.Batch
)

#test of confounded vs confounded with phosprep genetic variance central tendency
# relative.results <- lapply(names(results), function(x){
#   res <- results[[x]]
#   res/rowSums(res)
# })
# 
# names(relative.results) <- names(results)
# 
# #performs mann-whitney test
# wilcox.test(results$Confounded.Batchfit[ , 1], results$ProteinCovariateBatchPhos[ , 1], alternative="two.sided")$p.value
# 
# #performs mann-whitney test
# wilcox.test(relative.results$Confounded.Batchfit[ , 1], relative.results$ProteinCovariateBatchPhos[ , 1], alternative="two.sided")$p.value


#absolute values of the variance components boxplots
require(ggplot2)

for (ii_result in 1:6) {
  res <- data.frame( results[[ii_result]] )
  p <- ggplot(data.frame(var_estimate = do.call(c, res),
                          var_source = factor( rep( c(1:3), each = dim(res)[1]),
                                               labels = c("individual", "culture", "technical")) ),
               aes(x = var_source, y = log10(var_estimate)) )  +
          geom_violin(alpha = .6) + geom_boxplot(width = .2, alpha = .2) +
          xlab("Source of variation") + ylab(expression(log[10](Variance~Component))) +  
          ggtitle(names(results[ii_result])) +
          theme_bw() +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
          theme(text=element_text(family="Times New Roman", size=12)) 
  ggsave(paste0(names(results[ii_result]), ".pdf"), p, width=5.5, height=4)
}
  

##protein confounded batch effect corrected with axis common to non-confounded (y 7.5, -5.0)
res <- data.frame( results$Confounded )
p <- ggplot(data.frame(var_estimate = do.call(c, res),
                       var_source = factor( rep( c(1:3), each = dim(res)[1]),
                                            labels = c("individual", "culture", "technical")) ),
            aes(x = var_source, y = log10(var_estimate)) )  +
  geom_violin(alpha = .6) + geom_boxplot(width = .2, alpha = .2) +
  ylim(-5,7.5) + 
  xlab("Source of variation") + ylab(expression(log[10](Variance~Component))) +  
  ggtitle(names(results[1])) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text=element_text(family="Times New Roman", size=12)) 
ggsave(paste0("ConfoundedCommonScale", ".pdf"), p, width=5.5, height=4)


##protein corrected batch effect corrected with common axis (y 7.5, -5.0)
res <- data.frame( results$ProteinCorrected )
p <- ggplot(data.frame(var_estimate = do.call(c, res),
                       var_source = factor( rep( c(1:3), each = dim(res)[1]),
                                            labels = c("individual", "culture", "technical")) ),
            aes(x = var_source, y = log10(var_estimate)) )  +
  geom_violin(alpha = .6) + geom_boxplot(width = .2, alpha = .2) +
  ylim(-5,7.5) + 
  xlab("Source of variation") + ylab(expression(log[10](Variance~Component))) +  
  ggtitle(names(results[3])) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text=element_text(family="Times New Roman", size=12)) 
ggsave(paste0("ProteinCorrectedCommonScale", ".pdf"), p, width=5.5, height=4)



#relative values for the variance components boxplots
for (ii_result in 1:6) {
  res <- data.frame( results[[ii_result]] )
  varprop <- res/rowSums(res)
  p <- ggplot(data.frame(var_proportion = do.call(c, varprop),
                         var_source = factor( rep( c(1:3), each = dim(res)[1]),
                                              labels = c("individual", "culture", "technical")) ),
              aes(x = var_source, y = var_proportion) )  +
    geom_violin(alpha = .6) + geom_boxplot(width = .05, alpha = .2) +
    xlab("Source of variation") + ylab("Fraction of total variation") +  
    ggtitle(names(results[ii_result])) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(text=element_text(family="Times New Roman", size=12)) 
  ggsave(paste0(names(results[ii_result]), "standardized.pdf"), p, width=5.5, height=4)
}


#embed the fonts to that I don't have continuous issues with Illustrator on different devices...

#So that I can work with the embedded fonts
# install.packages("extrafont")
# library(extrafont)
# font_import() #only once
# loadfonts(device = "pdf")       #Register fonts for pdf output device. only once
# fonts()    #show available fonts
# Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.18/bin/gswin64c.exe")
# embed_fonts("Confounded.pdf", outfile="Confounded_embed.pdf")
# embed_fonts("ProteinCorrected.pdf", outfile = "ProteinCorrected_embed.pdf")
# embed_fonts("ProteinCorrectedstandardized.pdf", outfile="ProteinCorrectedstandardized_embed.pdf")
# embed_fonts("Confoundedstandardized.pdf", outfile = "Confoundedstandardized_embed.pdf")


#######Enrichment analyses ----

#Add GOID, Reactome, Entrez, HGNCID, HGNC symbol, and HGNC derived description of each protein gene annotation to multExpanded DF
if(!file.exists("multExpanded1_withDE_annotated.rds")){
multExpanded1_withDE_annotated <- AddAnnotation(multExpanded1_withDE)
}else{
  multExpanded1_withDE_annotated <- readRDS("multExpanded1_withDE_annotated.rds")
}

# add ELM motif instances
source("./ELM/elm_processing.R")

#combine with hprd motifs. processed using perseus on 9.19.15
motifs <- read.table("./Perseus/motif_annotations.txt", sep = "\t", header = T, stringsAsFactors = F)
motifs <- motifs[[1]]
multExpanded1_withDE_annotated <- cbind(multExpanded1_withDE_annotated, motifs)
saveRDS(multExpanded1_withDE_annotated, file = "./multExpanded1_withDE_annotated.rds")

# (Ontology) Enrichment analysis 
enrichment_tables <- Enrichment(multExpanded1_withDE_annotated)

#protein enrichments
# ProtLevelEnrichment.R

#peptide enrichments
# PeptideLevelEnrichment.R

# snp enrichments
# SNPenrich.R

# pqtl enrichment
# pqtlEnrich.R






