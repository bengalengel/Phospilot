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
phosdata <- CorrectedData[[9]]#normalized/batch corrected data frame with at lesat 1 obs in each bio rep
pilot <- CorrectedData[[10]]#same as above with mean ratios for each bio replicate
medianSub.quantiled <- CorrectedData[[3]]

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




#protein table output. Keep columns 244:length and remove unwanted manually. Only keep proteingroups that were PC corrected.
ProteinZia.out <- ProteinZia[ProteinZia$id %in% row.names(RegressedCommon), 244:length(ProteinZia)]
write.table(ProteinZia.out, file = "./GelPrep.table.txt", sep = "\t", row.names = F)
rm(ProteinZia.out)


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

# CONFOUNDED data 
# remove exp obs if not observed at least two times in each sample
confounded.batch.corrected <- com2[rowSums(is.na(com2[ , 1:4])) <= 2 & rowSums(is.na(com2[ , 5:8])) <= 2 & rowSums(is.na(com2[ , 9:12])) <= 2,]


# The peptide model variance components are returned. 
mcmcVarcomp.confounded <- NestedVar(ratios = confounded.batch.corrected, noMissing = TRUE)
colnames(mcmcVarcomp.confounded) <- c("individual","culture","residual")


#confounded without batch effect correction
# remove exp obs if not observed at least two times in each sample
confounded.batch.raw <- medianSub.quantiled


# The peptide model variance components are returned. 
mcmcVarcomp.confounded.batchfit <- NestedVar(ratios = confounded.batch.raw, noMissing = TRUE, BatchCorrected = FALSE)
colnames(mcmcVarcomp.confounded.batchfit) <- c("individual","culture","residual")


## incorporating protien data ----

# PROTEIN as a covariate data. It is better to run with protein as a covariate that use normalized data.
colnames(GelPrep) <- c("HL18862", "HL18486", "HL19160")
PhosProtGel <- merge(phosdata, GelPrep, by = "row.names", 
                     suffixes = c("_peptide", "_GelPrep") ) #3257 observations
rownames(PhosProtGel) <- PhosProtGel$Row.names
PhosProtGel <- PhosProtGel[ , -1]
PhosProtGel <- as.matrix(PhosProtGel)

#missing values are removed
mcmcVarcomp.proteinCov <- NestedVar(PhosProtGel, includeProteinCovariate = TRUE)
colnames(mcmcVarcomp.proteinCov) <- c("individual","culture","residual")

# PROTEIN NORMALIZED data. Do we see the same switch in culture vs technical?
colnames(ProtNormalized) <- colnames(confounded.batch.corrected)
mcmcVarcomp.proteinNorm <- NestedVar(ProtNormalized, includeProteinCovariate = FALSE, noMissing = TRUE)
colnames(mcmcVarcomp.proteinNorm) <- c("individual","culture","residual")



#GELPROT PROTEIN AS A COVARIATE WITH NON BE CORRECTED PHOSPHO DATA
colnames(GelPrep) <- c("HL18862", "HL18486", "HL19160")
PhosProtGelBatch <- merge(medianSub.quantiled, GelPrep, by = "row.names", 
                     suffixes = c("_peptide", "_GelPrep") ) #3257 observations after omitting missing values. suffix only used if names are different
rownames(PhosProtGelBatch) <- PhosProtGelBatch$Row.names
PhosProtGelBatch <- PhosProtGelBatch[ , -1]
PhosProtGelBatch <- as.matrix(PhosProtGelBatch)

#Estimate variance componenets
mcmcVarcomp.proteinCov.Batch <- NestedVar(PhosProtGelBatch, includeProteinCovariate = TRUE, BatchCorrected = FALSE)
colnames(mcmcVarcomp.proteinCov.Batch) <- c("individual","culture","residual")


#PhosPrep PROTEIN AS A COVARIATE WITH NON BE CORRECTED PHOSPHO DATA
PhosProtPhosBatch <- merge(medianSub.quantiled, PhosPrepCombatBio, by = "row.names", 
                  suffixes = c("_peptide", "_PhosPrep") ) #1308 observations
rownames(PhosProtPhosBatch) <- PhosProtPhosBatch$Row.names
PhosProtPhosBatch <- PhosProtPhosBatch[ , -1]
PhosProtPhosBatch <- as.matrix(PhosProtPhosBatch)


#Estimate variance componenets include protein estimates and batch as a covariate
mcmcVarcomp.PhosProteinCov.Batch <- NestedVar(PhosProtPhosBatch, includeProteinCovariate = TRUE, BatchCorrected = FALSE, PhosPrep = TRUE)
colnames(mcmcVarcomp.PhosProteinCov.Batch) <- c("individual","culture","residual")



##------ NRE fit Results ------

#absolute and standardized variance component plots for protein covariate and for confounded. (protein as a covariate makes biorep the largest contributor as opposed to the protein normalized and confounded data)


#combine confounded and protein as a covariate results in to a list of data frames for common processing.

results <- list(Confounded = mcmcVarcomp.confounded, #confounded after combat batch effect correction
                Confounded.Batchfit = mcmcVarcomp.confounded.batchfit, #med/quantile norm confounded with batch as a covariate
                ProteinCorrected = mcmcVarcomp.proteinCov, # pQTL regressed gelprot protein dataframe subtracted from combat corrected phospho
                ProteinNormalized = mcmcVarcomp.proteinNorm, # pQTL regressed gelprot protein dataframe as a covariate with combat corrected phospho
                ProteinCovariateBatch = mcmcVarcomp.proteinCov.Batch, # pQTL regressed gelprot protein dataframe and batch as covariates with median/quantile normalized phospho
                ProteinCovariateBatchPhos = mcmcVarcomp.PhosProteinCov.Batch #you get the idea
)


#summaries. Note that these are complete cases
lapply(results, summary)
lapply(results, dim)

lapply(names(results), function(x){
  res <- results[[x]]
  varprop <- res/rowSums(res)
  summary(varprop)
})
  


#So that I can work with the embedded fonts
# install.packages("extrafont")
library(extrafont)
# font_import() #only once
# loadfonts(device = "pdf")       #Register fonts for pdf output device. only once
fonts()    #show available fonts



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
  ggsave(paste0(names(results[ii_result]), ".pdf"), p, width=11, height=8.5)
}
  

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
  ggsave(paste0(names(results[ii_result]), "standardized.pdf"), p, width=11, height=8.5)
}



#embed the fonts to that I don't have continuous issues with Illustrator on different devices! (must install and point device to ghostscript freeware)
# For Windows - in each session
# Adjust the path to match your installation of Ghostscript
Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.18/bin/gswin64c.exe")
embed_fonts("Confounded.pdf", outfile="Confounded_embed.pdf")
embed_fonts("ProteinCorrected.pdf", outfile = "ProteinCorrected_embed.pdf")
embed_fonts("ProteinCorrectedstandardized.pdf", outfile="ProteinCorrectedstandardized_embed.pdf")
embed_fonts("Confoundedstandardized.pdf", outfile = "Confoundedstandardized_embed.pdf")



# Plot the variance component distributions
colnames(mcmcVarcomp) <- c("individual","biorep","residual")
boxplot(log10(mcmcVarcomp), ylab = "log10 variance component")
summary(mcmcVarcomp)
head(mcmcVarcomp)

# Histograms of log10 variance 
for (i in 1:ncol(mcmcVarcomp) ) {
  plot( density(log10(mcmcVarcomp[ ,i]), na.rm = T), xlab = "log10 variance", 
        main = paste(colnames(mcmcVarcomp)[i], "variance") )
}

# Scatter plots of log10 variance components
plot(log10(mcmcVarcomp[,1]),log10(mcmcVarcomp[,3]), 
     main = "log10 variance", xlab = colnames(mcmcVarcomp)[1], ylab = colnames(mcmcVarcomp)[3])
plot(log10(mcmcVarcomp[,1]),log10(mcmcVarcomp[,2]), 
     main = "log10 variance", xlab = colnames(mcmcVarcomp)[1], ylab = colnames(mcmcVarcomp)[2])
plot(log10(mcmcVarcomp[,2]),log10(mcmcVarcomp[,3]), 
     main = "log10 variance", xlab = colnames(mcmcVarcomp)[2], ylab = colnames(mcmcVarcomp)[3])


# Here we'd like to identify phosphpeptides with little or no variability 
# at the individual level and at the biological replicate level. To do so, 
# we standardized the values of the variance components for each phosphopeptides 
# with respect to its sum of variance components. The standardized variance 
# components are the proportion of the total variation in each phosphopeptides
# attributed to individuals, biological replicates, and technical replicates. 

# Boxplots of the standardized VCs confirm our observations from the raw VC values. 
# Proportion of variability attributed to biological replicates is the smallest, 
# followed by technical replicates, with individaul samples contributing the largest 
# portion of variabilty in expression levels. 
par(mfrow = c(1,1))
varprop <- mcmcVarcomp/rowSums(mcmcVarcomp)
labs = c("individual","biorep","tech")
boxplot((varprop), axes = F)
axis(1, at = c(1, 2, 3), labels = labs, col = "white"); axis(2)


# Heatmap representation of the standardized VCs. 
varprop <- na.omit(varprop)
require(gplots)
require(RColorBrewer)
colnames(varprop) = c("individual","bio","tech")
heatmap.2(as.matrix(varprop),
          col=brewer.pal(9,"YlGnBu"),
          Colv=F,
          labRow="",
          trace="none",
          srtCol=45,  ,adjCol = c(1,1),
          margins = c(6,5),
          cexCol=1.5,
          key.xlab = "Standardized VC", key.ylab=NULL, key.title = "",
)

return(mcmcVarcomp)
}









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





