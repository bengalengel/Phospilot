#here I am going to add SNP based annotation to multExpaned table for categorical enrichment test

# SNPenrich <- function(multExpanded1_withDE){
require(seqinr)

####create variant data frame and subset to those present in any of the three samples including the standard --------
#load snpeff_final dataset dataset (see snpeff folder readme file for construction).
SNPeffFinal <- read.table("E:/My Documents/Pilot/snpeff_final.txt", sep = "\t", header = T, stringsAsFactors = F, quote = "")
SNPeffFinal <- read.table("D:/snpeff_final.txt", sep = "\t", header = T, stringsAsFactors = F, quote = "")

#subset to 4 samples of interest: 18486, 18862, 19160, and the 19238 standard
sampleNames <- grep("18486|18862|19160|19238", names(SNPeffFinal), value = T)
variables <- names(SNPeffFinal)[1:13]
SNPeffFinal <- SNPeffFinal[,c(variables,sampleNames)]

#a quick view of the structure of the dataframe
head(SNPeffFinal[,14:length(SNPeffFinal)])
str(SNPeffFinal)

#There is a 1:1 mapping between transcripts and peptides. That is ensemble considers only one ORF per transcript. Genes, transcripts, peptides and snps are obviously duplicated but the combinations thereof are unique. See below for two examples of repeated snps at 1) different positions of resulting isoforms and 2) the same position of different isoforms. Note 2) forms a phosphorylatable residue.
which(duplicated(SNPeffFinal$snp))[1:3]#first is ""
SNPeffFinal[SNPeffFinal$snp == SNPeffFinal$snp[653],]
SNPeffFinal[SNPeffFinal$snp == SNPeffFinal$snp[654],]

#subset to those variants present in at least one of the four lines
hapTypes <- c("0|1", "1|0", "1|1")
#get logical index for subsetting
index <- apply(SNPeffFinal[,sampleNames], 1, function(x){
  any(hapTypes %in% x)})
SNPeffFinal <- SNPeffFinal[index,]

# How many snps are only represented in the standard? 19238
sampleNames <- sampleNames[!grepl("NA19238",sampleNames)]
index <- apply(SNPeffFinal[,sampleNames], 1, function(x){
  any(hapTypes %in% x)})

#Because the standard line cannot contribute to the observed variation, variants unique to it are removed. (note if standard is homozygous positive for variant the peptide itself cannot be observed. See section on effect size estimates)
SNPeffFinal <- SNPeffFinal[index,]#43234

#roughly 18K coding variants in at least 1 line (including standard)
length(unique(SNPeffFinal$snp))#17947





# Add sequence annotation information -------------------------------------

##add sequence information to each entry from the fasta file. Then add a flag column indicating proximity to nearest S/T/Y. Then a flag column indicating yes/no for +- 7 from S/T/Y. The same should be done with phosphorylated site and phosphosite annotation information.

# Just for the phosphosites observed; Are there snps in the vacinity of the phosphosite? The background for enrichmentwill be subtodiffphos + all snps.

#for each ENSPID find the squence contained within the FASTA file. Executed with an sapply call. 
hits <- sapply(SNPeffFinal$peptide, function(x) {
  hit <- grep(x,names(proteome))
  if(length(hit)==0){
    "no match"}else{
      return(hit)
    }
}
)
SNPeffFinal$ProteomeIndex <- hits


#   To be used for motif proximity enrichment analysis, that is for each protein/hit, is the variant within 7 residues of a S/T/Y? 

#getsequence from 'hits' index and position from snpeff table. Return distance to nearest S/T/Y residue. 'NA' results from either no Tyr in sequence or protein mapped to snp is not searched in database. 

# Annotation: Distance to nearest phosphorylatable residue ----------------


#for this work only missense variants are considered
SNPeffFinalMSonly <- SNPeffFinal[SNPeffFinal$effect == "missense_variant" | SNPeffFinal$effect == "missense_variant&splice_region_variant",]#41,388 nonsynonymous snps

#Dist to closest Tyrosine  
TyrDist <- function(ProtSeq, VariantPosition){
  #if ProtSeq index is valid 
  if(ProtSeq != "no match"){
    #retrieve sequence and Tyr index
    seq <- getSequence(proteome[[as.numeric(ProtSeq)]])
    Tyr <- grep("Y", seq)
    #if sequence contains tyrosines
    if(length(Tyr) > 0){
      #retrive variant position
      Variant <- as.integer(gsub("[^0-9]+", "", VariantPosition))#removes any non-numeric elements
      #return minimum distance
      min(abs(Variant - Tyr))}else{
        NA
      }
  }else
  {NA}
}

#Dist to closest Serine  
SerDist <- function(ProtSeq, VariantPosition){
  #if ProtSeq index is valid 
  if(ProtSeq != "no match"){
    #retrieve sequence and Syr index
    seq <- getSequence(proteome[[as.numeric(ProtSeq)]])
    Ser <- grep("S", seq)
    #if sequence contains Serines
    if(length(Ser) > 0){
      #retrive variant position
      Variant <- as.integer(gsub("[^0-9]+", "", VariantPosition))#removes any non-numeric elements
      #return minimum distance
      min(abs(Variant - Ser))}else{
        NA
      }
  }else
  {NA}
}

#Dist to closest Threonine  
ThrDist <- function(ProtSeq, VariantPosition){
  #if ProtSeq index is valid 
  if(ProtSeq != "no match"){
    #retrieve sequence and Thr index
    seq <- getSequence(proteome[[as.numeric(ProtSeq)]])
    Thr <- grep("T", seq)
    #if sequence contains Serines
    if(length(Thr) > 0){
      #retrive variant position
      Variant <- as.integer(gsub("[^0-9]+", "", VariantPosition))#removes any non-numeric elements
      #return minimum distance
      min(abs(Variant - Thr))}else{
        NA
      }
  }else
  {NA}
}

#apply minimum distance functions
SNPeffFinalMSonly$NearestTyrDist <- mapply(TyrDist, SNPeffFinalMSonly$ProteomeIndex, SNPeffFinalMSonly$aa)
SNPeffFinalMSonly$NearestSerDist <- mapply(SerDist, SNPeffFinalMSonly$ProteomeIndex, SNPeffFinalMSonly$aa)
SNPeffFinalMSonly$NearestThrDist <- mapply(ThrDist, SNPeffFinalMSonly$ProteomeIndex, SNPeffFinalMSonly$aa)

#####################min distance to nearest observed phosphorylation site ------

#for each protein group/site combination return comma delimited character string with the minimum distance between all snps and the phosphorylated site position


DistToPhos <- function(ProteinGroup, ProteinGroupPosition){
  if(any(unlist(strsplit(as.character(ProteinGroup), ";")) %in% SNPeffFinalMSonly$peptide)){
    proteins <- unlist(strsplit(as.character(ProteinGroup), ";"))
    positions <- unlist(strsplit(as.character(ProteinGroupPosition), ";"))
    MinDist <- c()
    for(i in seq_along(proteins)){
      #if the protein in the protein group has a snp
      if(proteins[i] %in% SNPeffFinalMSonly$peptide){
        #find the positions of snps within this protein
        VariantPositions <- SNPeffFinalMSonly[SNPeffFinalMSonly$peptide == proteins[i], "aa"]
        VariantPositions <- as.integer(gsub("[^0-9]+", "", VariantPositions))
        #calculate the minimum distance to the phosphosite amongst all variants
        dist <- min(abs(VariantPositions - as.numeric(positions[i])))#position and proteins should have same index
        MinDist <- c(MinDist,dist)
      }else{
        MinDist <- c(MinDist,NA)
      }
    }
    MinDist <- paste(MinDist, collapse = ";")
    MinDist
  }else{
    NA
  }
}

#apply closest snp to phosphorylation site function. UPDATED TO USE ONLY NS VARIANTS
multExpanded1_withDE$ClosestSNPtoSite <- mapply(DistToPhos, multExpanded1_withDE$Proteins, multExpanded1_withDE$Positions.within.proteins)

#apply closest snp to phosphorylation site function. UPDATED TO USE ONLY NS VARIANTS
multExpanded1_withDE$ClosestSNPtoSite <- mapply(DistToPhos, multExpanded1_withDE$Proteins, multExpanded1_withDE$Positions.within.proteins)



# calculate minimum of all the distances for each protein group
multExpanded1_withDE$ClosestSNPtoSiteMin <- sapply(multExpanded1_withDE$ClosestSNPtoSite, function(x){
  distances <- as.numeric(unlist(strsplit(x, ";")))
  if(!all(is.na(distances))){
    min(distances, na.rm = T)
  }else
  {NA}
})  

#Hypothesis: Biological phosphosite variance correlates positively with the distance to observed phosphorylation site   

#is there a positive correlation between bio Varcomp and closest SNP?
holder <- multExpanded1_withDE[,c("idmult", "ClosestSNPtoSiteMin")]

VarcompDist <- merge(varcomp, holder, by.x = "row.names", by.y = "idmult")

index <- !is.na(VarcompDist$ClosestSNPtoSiteMin)
VarcompDist <- VarcompDist[index,]#length of 2427

#There may be a relationship between the distance between phosphosite and the closest SNP and individual variance component magnitude. (per [0-5] AA window interval there may be more sites with disproportionate representation in the 'high' individual variance component section. Note that all of these sites belong to a protein group with at least one member having a nonsyn snp. Therefore these proteins are overrepresented in the diffexp subset. This graphic effectively tests for overrepresentation IN ADDITION TO OVERREPRESENTATION IN HIGH VARCOMP DUE TO HAVING A SNP? If there is an enrichment it seems very small.

plot(log10(VarcompDist$ClosestSNPtoSiteMin), log10(VarcompDist$individual), xlim = c(-.1, 4.0), ylab = "log10(Individual Variance Component)", xlab = "log10(AA Distance between phosphosite and closest SNP)")
plot(log10(VarcompDist$ClosestSNPtoSiteMin), log10(VarcompDist$biorep), ylab = "log10(Biological Variance Component)", xlab = "log10(AA Distance between phosphosite and closest SNP)")
plot(log10(VarcompDist$ClosestSNPtoSiteMin), log10(VarcompDist$residual), ylab = "log10(Individual Variance Component)", xlab = "log10(AA Distance between phosphosite and closest SNP)")

#showing the 0 data point
plot(VarcompDist$ClosestSNPtoSiteMin, log10(VarcompDist$individual), xlim = c(0,10), ylab = "log10(Individual Variance Component)", xlab = "AA Distance between phosphosite and closest SNP")
plot(VarcompDist$ClosestSNPtoSiteMin, log10(VarcompDist$biorep), xlim = c(0,10), ylab = "log10(Individual Variance Component)", xlab = "AA Distance between phosphosite and closest SNP")
plot(VarcompDist$ClosestSNPtoSiteMin, log10(VarcompDist$residual), xlim = c(0,10), ylab = "log10(Individual Variance Component)", xlab = "AA Distance between phosphosite and closest SNP")

#out to 100
plot(VarcompDist$ClosestSNPtoSiteMin, log10(VarcompDist$individual), xlim = c(0,100), ylab = "log10(Individual Variance Component)", xlab = "AA Distance between phosphosite and closest SNP")

#For the protein normalized data



#min distance to nearest annotated phosphorylation or other site modification site







#####adding SNP presence/absence annotation to multExpanded1_withDE --------------------
#is there a match to any nonsyn snp for any protiens assigned to phosphopeptide 
AnyMatchProtSNPeff <- function(queryproteins){
  #This function looks for a match between queryproteins and the proteins containing a non-synonymous snp
  #queryproteins are proteins assigned to the given phosphopeptide
  if(any(unlist(strsplit(as.character(queryproteins), ";")) %in% SNPeffFinal$peptide)){
    "+"
  }else{
    "-"
  }
}
#apply anymatch function to all 'leading proteins' and 'majority protein ids'
multExpanded1_withDE$NsSnpPositive <- mapply(AnyMatchProtSNPeff, multExpanded1_withDE$Leading.proteins)
multExpanded1_withDE$ppNsSnpPositive <- mapply(AnyMatchProtSNPeff, multExpanded1_withDE$ppMajorityProteinIDs)

#how many nonsyn snps per protein group assigned to a phosphopeptide?
NumMatchesProtSNPeff <- function(queryproteins){
  #This function counts the matches between queryproteins and the proteins containing a non-synonymous snp
  #queryproteins are proteins assigned to the peptide
  sum(unlist(strsplit(as.character(queryproteins), ";")) %in% SNPeffFinal$peptide)
}

#apply count function to all 'leading proteins' and 'majority protein ids'
multExpanded1_withDE$NsSnpCount <- mapply(NumMatchesProtSNPeff, multExpanded1_withDE$Leading.proteins)
multExpanded1_withDE$ppNsSnpCount <- mapply(NumMatchesProtSNPeff, multExpanded1_withDE$ppMajorityProteinIDs)

# Less than 1% of the phosphopeptides are mapped to a protein group that (collectively) contains > 1 snp #THIS NEEDS TO BE CORRECTED FOR UNIQUENESS
table(multExpanded1_withDE$NsSnpCount)

#For proteins mapped to peptides using protein prep data the picture is more complicated because more information is being used (more protein ids/peptide). That is the same snp is present in multiple isoforms, which cannot be disambiguated geven the shotgun level information. Proteotypic peptides are need for this.
table(multExpanded1_withDE$ppNsSnpCount)


###Mumblings ---------
# first list; #s 3-4 and some more
# 3) For a given snp found in the phosphoproteomics data, which lines have it?
# 4) For each line that has the snp, what is its genotype?

# Questions about variability should be addressed in the differential phosphorylation analysis via enrichment analysis. Perhaps snps within a domain are enriched, within unstructured regions are enriched, or sites near to a phosphorylatable residue are *especially* enriched for variability. The effect of snps would have to be controlled for somehow by comparing if phosphorylation sites within these regions are intrinsically  

# This information will be used to ask questions about observation and effect size of genetically induced differences in phosphorylation at steady state. 
# 
# hypothesis 1) Missing phosphosites where an individual does not have phosphorylatable residue in one condition but does in another. heterozygote/homozygote cases...OK getting closer
# 
# If it is not observed. must plot the normalized by genotype results for instances where there is differential exprssion within the the DE subset?

######### Effect Size estimate attempt  ------------------
# teaser only 350 possible snps and I only found 1 that was identified, let alone quantified
# I want to estimate effect size for instances where the snp removes a phosphorylation site that was identified in another condition. sites normalized relative to homozygote null. I need the following:
#   1) snp that causes mutation of phosphorylatable residue to non-phosphorylatable residue
#   2) standard line is heterozygous positive or homozygous negative for the mutation
#   3) affected lines only heterozygous. 
#   4) site quantified in homozygous null and heterozygous state. (If still present in heterozygous positive it must be due to one of the other majority protein ids)



#subset snpeff st conditions 1 and 2 are met

#snp modifies a phosphorylatable residue. Regexp based subset
index <- grep("p.Ser|p.Thr|p.Tyr", SNPeffFinal$aa)
SNPeffFinalPhos <- SNPeffFinal[index,]#8200 obs

#standard line is at least heterozygous or homo negative (SILAC standard must be present to quantify reference peptide)
SNPeffFinalPhos <- SNPeffFinalPhos[SNPeffFinalPhos$NA19238 == "0|0" | SNPeffFinalPhos$NA19238 == "1|0"
                                   | SNPeffFinalPhos$NA19238 == "0|1",]#n = 6750

#homozygous null and heterozygous represented in experimental lines.
sampleNames <- sampleNames[!grepl("NA19238",sampleNames)]
hetero <- c("0|1", "1|0")

index <- apply(SNPeffFinalPhos[,sampleNames], 1, function(x){
  any(hetero %in% x) && any("1|1" %in% x)
}
)
SNPeffFinalPhos <- SNPeffFinalPhos[index,]#782....

#how many unique snps? This is a really small number...
length(unique(SNPeffFinalPhos$snp))#349

#which phosphopeptides match these SNPs? add these ids to snpeff dataframe. After this add the normalized numbers
#match ESNPID and *any* 'position.within.proteins'. This seems like a job for plyr or dplyr. two logical tests using mapply for now. 1st choose snp to assign to protein groups assigned to phosphopeptide.

#if ENSPID assigned to snp common with any member of protein group AND any of *Positions within proteins* match snp position, return that snp. If multiple snps match ensure uniqueness. If multiple unique snps match the same peptide then you are in Valhalla and rejoice. But really if that happens assign the peptide to both snps.
SNPprotein <- SNPeffFinalPhos$peptide
SNPproteinposition <- gsub("[a-zA-Z\\.]*" , replacement = "", SNPeffFinalPhos$aa)#deletions give underscore followed by second position
SNPproteinposition <- gsub("_[0-9]*" , replacement = "", SNPproteinposition)#putative position of phosphopeptide

#this will be applied mapply style to get an index to subset snpeff data frame
PepSnpMatch <- function(ProteinGroup, ProteinGroupPosition){
  #things this function does
  if(any(unlist(strsplit(as.character(ProteinGroup), ";")) %in% SNPprotein)
     && any(unlist(strsplit(as.character(ProteinGroupPosition), ";")) %in% SNPproteinposition)){
    #what are the matching indices and are they equal? return the matches
    protmatch <- which(SNPprotein %in% unlist(strsplit(as.character(ProteinGroup), ";")))
    positionmatch <- which(SNPproteinposition %in% unlist(strsplit(as.character(ProteinGroupPosition), ";")))
    index <- intersect(protmatch,positionmatch)
    if(length(index) > 0){
      return(index)
    }else{
      NA
    }
  }else{#can add 'NAP' to see how many no matching proteins
    NA
  }
}

#only one match with leading proteins!
newindex <- mapply(PepSnpMatch,multExpanded1_withDE$Leading.proteins,multExpanded1_withDE$Positions.within.proteins)
which(!is.na(newindex))
newindex[7649]

#only one match with proteins!
newindex <- mapply(PepSnpMatch,multExpanded1_withDE$Proteins,multExpanded1_withDE$Positions.within.proteins)
which(!is.na(newindex))
newindex[7649]

##### code check confirms only 1 match ------
##OK so how many peptides/positions matched but didn't have a matching index? (this could happen when a phosphopeptide positions for one of the proteins in the group matches a protein and a position in the snpeff file but the position is in another protein)
PepSnpMatch2 <- function(ProteinGroup, ProteinGroupPosition){
  if(any(unlist(strsplit(as.character(ProteinGroup), ";")) %in% SNPprotein)
     && any(unlist(strsplit(as.character(ProteinGroupPosition), ";")) %in% SNPproteinposition)){
    "+"
  }else{
    "-"
  }
}

#more hits.
newindex <- mapply(PepSnpMatch2,multExpanded1_withDE$Leading.proteins,multExpanded1_withDE$Positions.within.proteins)
table(newindex)##75 hits

newindex <- mapply(PepSnpMatch2,multExpanded1_withDE$Proteins,multExpanded1_withDE$Positions.within.proteins)
table(newindex)##78 hits

#subset a few and check to see if they actually match (my function wasn't working test)
index <- which(newindex=="+")
checkDF <- multExpanded1_withDE[index,c("Proteins","Positions.within.proteins")]
head(checkDF)
protmatch <- which(SNPprotein %in% unlist(strsplit(as.character(checkDF[1,1]), ";")))
positionmatch <- which(SNPproteinposition %in% unlist(strsplit(as.character(checkDF[1,2]), ";")))
protmatch
positionmatch
intersect(protmatch,positionmatch)
#my script was functional









##############Enrichment tests ------
# 1)  Test for enrichment in diffphos. background is all sites subject to DiffPhos. Foreground is omnibus F significance. Category is 'with snp' or without snp at the phosphopeptide level. Contingency matrix is of the form:

#                            in category  not in category  
#                    DErow
#                    NotDErow
# 'in' category is any majority/leading protein(s) assigned to this phosphosite has a snp.  

subtoDE <- multExpanded1_withDE[multExpanded1_withDE$SubtoDE == "+",] #4738
subtoDEpn <- multExpanded1_withDE[multExpanded1_withDE$SubtoDEpn == "+",] #3488

#confounded analysis
row1 <- c(nrow(subtoDE[subtoDE$globalFsig == "+" & subtoDE$NsSnpPositive == "+",]), 
          nrow(subtoDE[subtoDE$globalFsig == "+" & subtoDE$NsSnpPositive == "-",]))

row2 <- c(nrow(subtoDE[subtoDE$globalFsig == "-" & subtoDE$NsSnpPositive == "+",]), 
          nrow(subtoDE[subtoDE$globalFsig == "-" & subtoDE$NsSnpPositive == "-",]))

#FEtest
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g")
result$p.value
4.492123e-05 


#protnormalized analysis using Zia's data (much less significant and speaks to penetrance at the post translational level)
row1 <- c(nrow(subtoDEpn[subtoDEpn$globalFsigpn == "+" & subtoDEpn$ppNsSnpPositive == "+",]), 
          nrow(subtoDEpn[subtoDEpn$globalFsigpn == "+" & subtoDEpn$ppNsSnpPositive == "-",]))

row2 <- c(nrow(subtoDEpn[subtoDEpn$globalFsigpn == "-" & subtoDEpn$ppNsSnpPositive == "+",]), 
          nrow(subtoDEpn[subtoDEpn$globalFsigpn == "-" & subtoDEpn$ppNsSnpPositive == "-",]))

#FEtest
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g")
result$p.value
0.001188539

#here is the relative proportion difference
apply(contmatrix,1,function(x) x[1]/sum(x))

row1      row2 
0.4036872 0.3530026 

# 2) TEST FOR ENRICHMENT IN QUADRANTS
################################
#                                  in category  not in category  
#                    INQUADRANT
#                    NotINQUADRANT
# 'in' category is any majority/leading protein(s) assigned to this phosphosite has a snp.  

subtoVC <- multExpanded1_withDE[multExpanded1_withDE$SubtoVarcomp == "+",] #6360
subtoVCpn <- multExpanded1_withDE[multExpanded1_withDE$ppSubtoVarcomp == "+",] #3485

#CONFOUNDED ANALYSIS
# a) High individual hish biological RESULT: NOT SIG
row1 <- c(nrow(subtoVC[subtoVC$HighIndVar == "+" & subtoVC$HighBioVar == "+" & subtoVC$NsSnpPositive == "+",]), 
          nrow(subtoVC[subtoVC$HighIndVar == "+" & subtoVC$HighBioVar == "+" & subtoVC$NsSnpPositive == "-",]))

row2 <- c(nrow(subtoVC[!(subtoVC$HighIndVar == "+" & subtoVC$HighBioVar == "+") & subtoVC$NsSnpPositive == "+",]), 
          nrow(subtoVC[!(subtoVC$HighIndVar == "+" & subtoVC$HighBioVar == "+") & subtoVC$NsSnpPositive == "-",]))

#FEtest
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g")
result$p.value
0.168596 

# b) high individual and low biological variance. RESULT: NOT SIG ENRICHED

row1 <- c(nrow(subtoVC[subtoVC$HighIndVar == "+" & subtoVC$LowBioVar == "+" & subtoVC$NsSnpPositive == "+",]), 
          nrow(subtoVC[subtoVC$HighIndVar == "+" & subtoVC$LowBioVar == "+" & subtoVC$NsSnpPositive == "-",]))

row2 <- c(nrow(subtoVC[!(subtoVC$HighIndVar == "+" & subtoVC$LowBioVar == "+") & subtoVC$NsSnpPositive == "+",]), 
          nrow(subtoVC[!(subtoVC$HighIndVar == "+" & subtoVC$LowBioVar == "+") & subtoVC$NsSnpPositive == "-",]))

#FEtest
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g")
result$p.value
0.3198735 

# c) low individual and low biological variance. RESULT: NOT SIG ENRICHED

row1 <- c(nrow(subtoVC[subtoVC$LowIndVar == "+" & subtoVC$LowBioVar == "+" & subtoVC$NsSnpPositive == "+",]), 
          nrow(subtoVC[subtoVC$LowIndVar == "+" & subtoVC$LowBioVar == "+" & subtoVC$NsSnpPositive == "-",]))

row2 <- c(nrow(subtoVC[!(subtoVC$LowIndVar == "+" & subtoVC$LowBioVar == "+") & subtoVC$NsSnpPositive == "+",]), 
          nrow(subtoVC[!(subtoVC$LowIndVar == "+" & subtoVC$LowBioVar == "+") & subtoVC$NsSnpPositive == "-",]))

#FEtest
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g")
result$p.value
0.5806627

# d) low individual and high biological variance. RESULT: NOT SIG ENRICHED, BUT SIGNIFICANTLY DEPLETED

row1 <- c(nrow(subtoVC[subtoVC$LowIndVar == "+" & subtoVC$HighBioVar == "+" & subtoVC$NsSnpPositive == "+",]), 
          nrow(subtoVC[subtoVC$LowIndVar == "+" & subtoVC$HighBioVar == "+" & subtoVC$NsSnpPositive == "-",]))

row2 <- c(nrow(subtoVC[!(subtoVC$LowIndVar == "+" & subtoVC$HighBioVar == "+") & subtoVC$NsSnpPositive == "+",]), 
          nrow(subtoVC[!(subtoVC$LowIndVar == "+" & subtoVC$HighBioVar == "+") & subtoVC$NsSnpPositive == "-",]))

#FEtest
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g")
result$p.value
0.9655049 
result <- fisher.test(contmatrix, alternative = "l")
result$p.value
0.0400571

# d) low individual and high biological variance. RESULT: NOT SIG ENRICHED, BUT SIGNIFICANTLY DEPLETED

row1 <- c(nrow(subtoVCpn[subtoVCpn$LowIndVar == "+" & subtoVC$HighBioVar == "+" & subtoVC$NsSnpPositive == "+",]), 
          nrow(subtoVC[subtoVC$LowIndVar == "+" & subtoVC$HighBioVar == "+" & subtoVC$NsSnpPositive == "-",]))

row2 <- c(nrow(subtoVC[!(subtoVC$LowIndVar == "+" & subtoVC$HighBioVar == "+") & subtoVC$NsSnpPositive == "+",]), 
          nrow(subtoVC[!(subtoVC$LowIndVar == "+" & subtoVC$HighBioVar == "+") & subtoVC$NsSnpPositive == "-",]))

#FEtest
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g")
result$p.value
0.9655049 
result <- fisher.test(contmatrix, alternative = "l")
result$p.value
0.0400571



# e) just high individual variance RESULT:significantly enriched (barely. fails two sided)

row1 <- c(nrow(subtoVC[subtoVC$HighIndVar == "+"  & subtoVC$NsSnpPositive == "+",]), 
          nrow(subtoVC[subtoVC$HighIndVar == "+"  & subtoVC$NsSnpPositive == "-",]))

row2 <- c(nrow(subtoVC[subtoVC$HighIndVar == "-"  & subtoVC$NsSnpPositive == "+",]), 
          nrow(subtoVC[subtoVC$HighIndVar == "-"  & subtoVC$NsSnpPositive == "-",]))

#FEtest
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g")
result$p.value
0.0494056
result <- fisher.test(contmatrix, alternative = "l")
result$p.value
0.9565258

# f) just low individual variance RESULT: significantly depleted (since this is the mirror image of just high)

row1 <- c(nrow(subtoVC[subtoVC$LowIndVar == "+"  & subtoVC$NsSnpPositive == "+",]), 
          nrow(subtoVC[subtoVC$LowIndVar == "+"  & subtoVC$NsSnpPositive == "-",]))

row2 <- c(nrow(subtoVC[subtoVC$LowIndVar == "-"  & subtoVC$NsSnpPositive == "+",]), 
          nrow(subtoVC[subtoVC$LowIndVar == "-"  & subtoVC$NsSnpPositive == "-",]))

#FEtest
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g")
result$p.value
0.9565258
result <- fisher.test(contmatrix, alternative = "l")
result$p.value
0.0494056


# g) just high bio variance RESULT: not sig enriched or depleted

row1 <- c(nrow(subtoVC[subtoVC$HighBioVar == "+"  & subtoVC$NsSnpPositive == "+",]), 
          nrow(subtoVC[subtoVC$HighBioVar == "+"  & subtoVC$NsSnpPositive == "-",]))

row2 <- c(nrow(subtoVC[subtoVC$HighBioVar == "-"  & subtoVC$NsSnpPositive == "+",]), 
          nrow(subtoVC[subtoVC$HighBioVar == "-"  & subtoVC$NsSnpPositive == "-",]))

#FEtest
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g")
result$p.value
0.6622874
result <- fisher.test(contmatrix, alternative = "l")
result$p.value
0.3579017


# h) just low bio variance RESULT: not sig enriched or depleted (again note that this is the mirror image)

row1 <- c(nrow(subtoVC[subtoVC$LowBioVar == "+"  & subtoVC$NsSnpPositive == "+",]), 
          nrow(subtoVC[subtoVC$LowBioVar == "+"  & subtoVC$NsSnpPositive == "-",]))

row2 <- c(nrow(subtoVC[subtoVC$LowBioVar == "-"  & subtoVC$NsSnpPositive == "+",]), 
          nrow(subtoVC[subtoVC$LowBioVar == "-"  & subtoVC$NsSnpPositive == "-",]))

#FEtest
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g")
result$p.value
0.3579017
result <- fisher.test(contmatrix, alternative = "l")
result$p.value
0.6622874




#sidebar comparing enrichment of DE in varcomp sanity check
###################
# a) high individual and low biological variance enriched in DE. RESULT: SUPER SIG (phew)

# first must subset to those subjected to DE so that "-" rows make sense.
subtoVC <- subtoVC[subtoVC$SubtoDE=="+",]#4732

row1 <- c(nrow(subtoVC[subtoVC$HighIndVar == "+" & subtoVC$LowBioVar == "+" & subtoVC$globalFsig == "+",]), 
          nrow(subtoVC[subtoVC$HighIndVar == "+" & subtoVC$LowBioVar == "+" & subtoVC$globalFsig == "-",]))

row2 <- c(nrow(subtoVC[!(subtoVC$HighIndVar == "+" & subtoVC$LowBioVar == "+") & subtoVC$globalFsig == "+",]), 
          nrow(subtoVC[!(subtoVC$HighIndVar == "+" & subtoVC$LowBioVar == "+") & subtoVC$globalFsig == "-",]))

#FEtest
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g")
result$p.value
4.720298e-53 

# b) high individual and high biological variance enriched in DE. RESULT: SIG but less so due to  (phew)

# first must subset to those subjected to DE so that "-" rows make sense.

row1 <- c(nrow(subtoVC[subtoVC$HighIndVar == "+" & subtoVC$HighBioVar == "+" & subtoVC$globalFsig == "+",]), 
          nrow(subtoVC[subtoVC$HighIndVar == "+" & subtoVC$HighBioVar == "+" & subtoVC$globalFsig == "-",]))

row2 <- c(nrow(subtoVC[!(subtoVC$HighIndVar == "+" & subtoVC$HighBioVar == "+") & subtoVC$globalFsig == "+",]), 
          nrow(subtoVC[!(subtoVC$HighIndVar == "+" & subtoVC$HighBioVar == "+") & subtoVC$globalFsig == "-",]))

#FEtest
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g")
result$p.value
1.344748e-07 

# c) high individual variance enriched in DE? RESULT: SUPER SIG as it should be (phew)

# first must subset to those subjected to DE so that "-" rows make sense.

row1 <- c(nrow(subtoVC[subtoVC$HighIndVar == "+"  & subtoVC$globalFsig == "+",]), 
          nrow(subtoVC[subtoVC$HighIndVar == "+"  & subtoVC$globalFsig == "-",]))

row2 <- c(nrow(subtoVC[subtoVC$HighIndVar == "-"  & subtoVC$globalFsig == "+",]), 
          nrow(subtoVC[subtoVC$HighIndVar == "-"  & subtoVC$globalFsig == "-",]))

#FEtest
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g")
result$p.value
row1 
1.023127e-174 
> contmatrix
[,1] [,2]
row1 1329 2342
row2    1 1060
#This confirms that what I am doing is consistent. Only 1 Fsig site is designated low ind variance. 

# d) for completeness low individual variance enriched in DE? RESULT: 
# first must subset to those subjected to DE so that "-" rows make sense.

row1 <- c(nrow(subtoVC[subtoVC$LowIndVar == "+"  & subtoVC$globalFsig == "+",]), 
          nrow(subtoVC[subtoVC$LowIndVar == "+"  & subtoVC$globalFsig == "-",]))

row2 <- c(nrow(subtoVC[subtoVC$LowIndVar == "-"  & subtoVC$globalFsig == "+",]), 
          nrow(subtoVC[subtoVC$LowIndVar == "-"  & subtoVC$globalFsig == "-",]))

#FEtest
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g") #less is significant as well as two sided of course
result$p.value


# e) for completeness low individual variance low biological variance enriched in DE? RESULT: 

row1 <- c(nrow(subtoVC[subtoVC$LowIndVar == "+" & subtoVC$LowBioVar == "+" & subtoVC$globalFsig == "+",]), 
          nrow(subtoVC[subtoVC$LowIndVar == "+" & subtoVC$LowBioVar == "+" & subtoVC$globalFsig == "-",]))

row2 <- c(nrow(subtoVC[!(subtoVC$LowIndVar == "+" & subtoVC$LowBioVar == "+") & subtoVC$globalFsig == "+",]), 
          nrow(subtoVC[!(subtoVC$LowIndVar == "+" & subtoVC$LowBioVar == "+") & subtoVC$globalFsig == "-",]))

#FEtest
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g") 
result$p.value
1
#again it makes sense that the one miss from limma is low bio and technical. I bet more of these guys would have been picked up if limma had not been used. This would make a good part of a course.
contmatrix
[,1] [,2]
row1    1  244
row2 1329 3158

# f) confirmt that there is no enrichment within the lowind and high bio
row1 <- c(nrow(subtoVC[subtoVC$LowIndVar == "+" & subtoVC$HighBioVar == "+" & subtoVC$globalFsig == "+",]), 
          nrow(subtoVC[subtoVC$LowIndVar == "+" & subtoVC$HighBioVar == "+" & subtoVC$globalFsig == "-",]))

row2 <- c(nrow(subtoVC[!(subtoVC$LowIndVar == "+" & subtoVC$HighBioVar == "+") & subtoVC$globalFsig == "+",]), 
          nrow(subtoVC[!(subtoVC$LowIndVar == "+" & subtoVC$HighBioVar == "+") & subtoVC$globalFsig == "-",]))

#FEtest
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g") 
result$p.value

contmatrix
row1    0  816
row2 1330 2586


###################





#protnormalized analysis using Zia's data
row1 <- c(nrow(subtoDEpn[subtoDEpn$globalFsig == "+" & subtoDEpn$NsSnpPositive == "+",]), 
          nrow(subtoDEpn[subtoDEpn$globalFsig == "+" & subtoDEpn$NsSnpPositive == "-",]))

row2 <- c(nrow(subtoDEpn[subtoDEpn$globalFsig == "-" & subtoDEpn$NsSnpPositive == "+",]), 
          nrow(subtoDEpn[subtoDEpn$globalFsig == "-" & subtoDEpn$NsSnpPositive == "-",]))

#FEtest
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g")
result$p.value





###########


# pqtl import and enrichment tests.
#######
#import pqtl table from Battle to use as well. Here I will need to match to protein names
pqtl <- read.table("C:/Users/Brett/Dropbox/Postdoc-Gilad/Yannick SNP Phos/1260793_DatafileS1_pQTLs.csv", sep = ",", header = T, 
                   stringsAsFactors = F)
pqtl <- read.table("E:/My Documents/Dropbox/Postdoc-Gilad/Yannick SNP Phos/1260793_DatafileS1_pQTLs.csv", sep = ",", header = T, 
                   stringsAsFactors = F)

#must expand this table to include all ENSPids associated with an ENSGid using the proper database. What ensembl version did Zia use to map ids?
#########

# What is the 








#######

anysnp <- multExpanded1_withDE[multExpanded1_withDE$Ind18486_SNPiso == "+" | multExpanded1_withDE$Ind18862_SNPiso == "+" | multExpanded1_withDE$Ind19160_SNPiso == "+",]
#make background as data frame and add snp annotation
background <- as.data.frame(background2)#background2 is the subtoDE df but with 'rev' entries removed.
colnames(background) <- "Protein"
background$anysnp <- ifelse(background$Protein %in% anysnp$Protein, "+","-")
#add the site ids as a column in the data frame
background$idmult <- subtoDE$idmult
#add the DE in any contrast information at the site level!
GlobalFs <- subtoDE[subtoDE$globalFsig=="+",]#1355 obs
background$DEany <- ifelse(background$idmult %in% GlobalFs$idmult, "+", "-")

#create the first row of the contingency matrix
row1 <- c(nrow(background[background$DEany=="+" & background$anysnp == "+",]), nrow(background[background$DEany == "+" & background$anysnp == "-",])) 
# 530 and 824 without specific isoform annotation (appropriate given this is bottom up MS) #note that there is a blank protein which may be causing the slight counting issue


#second row
row2 <- c(nrow(background[background$DEany=="-" & background$anysnp == "+",]), nrow(background[background$DEany == "-" & background$anysnp == "-",])) 
# [1] 1072 2565

#FEtest
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g")
result$p.value

#p-value = 7.376e-11  





















#subset to those variants 

##adding annotation to multExpanded1_withDE. The goal is to use a categorical enrichment test to see if phosphoproteins with nonsyn-snps unique to at least one of the lines are more likely to be over-represented in DE than those that do not have non-syn snps. 

# The next goal is to identify those sites with missing phosphorylated residues or regions surrounding the residues to get a sense of effect size after normalizing by protein.

#import yannick files
yannick18486 <- read.table("C:/Users/Brett/Desktop/pilot_18486_allNS.csv", sep = ',', header = T)
yannick18862 <- read.table("C:/Users/Brett/Desktop/pilot_18862_allNS.csv", sep = ',', header = T)
yannick19160 <- read.table("C:/Users/Brett/Desktop/pilot_19160_allNS.csv", sep = ',', header = T)
pqtl <- read.table("E:/My Documents/Dropbox/Postdoc-Gilad/Yannick SNP Phos/1260793_DatafileS1_pQTLs.csv", sep = ',', header = T)

#add a +/- based on UNIPROT.ID annotation, also add if isoform is present. Did I really add the annotation even if there is an
# isoform? NO see below
multExpanded1_withDE$Ind18486_SNP = ifelse(multExpanded1_withDE$Protein %in% yannick18486$UNIPROT.ID, "+","-")
multExpanded1_withDE$Ind18862_SNP = ifelse(multExpanded1_withDE$Protein %in% yannick18862$UNIPROT.ID, "+","-")
multExpanded1_withDE$Ind19160_SNP = ifelse(multExpanded1_withDE$Protein %in% yannick19160$UNIPROT.ID, "+","-")

# test <- as.character(multExpanded1_withDE$Protein)
# isos <- multExpanded1_withDE[which(nchar(test)>6),]
# any(isos$Ind18486_SNP=="+")
# any(isos$Ind18486_SNPiso=="+")



#subsetting to first six characters to ensure isoforms get annoatated. Ignoring isoform maps from phosphopeptides for the moment
multExpanded1_withDE$Ind18486_SNPiso = ifelse(substr(multExpanded1_withDE$Protein,1,6) %in% yannick18486$UNIPROT.ID, "+","-")
multExpanded1_withDE$Ind18862_SNPiso = ifelse(substr(multExpanded1_withDE$Protein,1,6) %in% yannick18862$UNIPROT.ID, "+","-")
multExpanded1_withDE$Ind19160_SNPiso = ifelse(substr(multExpanded1_withDE$Protein,1,6) %in% yannick19160$UNIPROT.ID, "+","-")


##subset to those phosphoobservations subject to DE (the background set of possibilities)
subtoDE <- multExpanded1_withDE[multExpanded1_withDE$SubtoDE == "+",]
subtoDE <- subtoDE[!grepl("REV",subtoDE$Protein),]#now removed REVerse hits


##breakdown to the level of protein for categorical enrichment analysis
background <- subtoDE$Protein

#note that 5 of these proteins have a REVerse designation as the leading protein.
library(plyr)
count(grep("REV",subtoDE$Protein))
#they will be removed from this analysis for now, but I think I can go back and assign them to a different protein within the protein groups file later. This should also be done within the major file as well. I will ask about this in MQ help board.
background2 <- background[!grepl("REV", background)]
# > any(grep("REV",background2))
# [1] FALSE
# > any(grep("REV",background))
# [1] TRUE

#also note one protein is missing!
#missingprotein <- write.csv(background2,"background2.csv")


##**********************************using protien counts as background (not correcting for multiply phosphorylated proteins)
BGtable <- as.matrix(table(background2))#the parenthetical should be a factor vector of ids passed to this function
#remove entries with 0
BGtable <- BGtable[BGtable!=0,,drop=F]#1991 proteins

##add the snp annotation
BGtable <- as.data.frame(BGtable)
#now use the isoform information
anysnp <- multExpanded1_withDE[multExpanded1_withDE$Ind18486_SNPiso == "+" | multExpanded1_withDE$Ind18862_SNPiso == "+" | multExpanded1_withDE$Ind19160_SNPiso == "+",]
BGtable$anysnp <- ifelse(row.names(BGtable)%in%anysnp$Protein,"+","-")

# 
# anysnp <- multExpanded1_withDE[multExpanded1_withDE$Ind18486_SNP == "+" | multExpanded1_withDE$Ind18862_SNP == "+" |                                  multExpanded1_withDE$Ind19160_SNP == "+",]
# BGtable$anysnp <- ifelse(row.names(BGtable)%in%anysnp$Protein,"+","-")

#somehow the first entry has no protein but is found with a SNP!

#also note the isoform issues. Correct annotation to map to any isoform which is a uniprot identifier followed by a dash.

# Now add the DE in any contrast information (already performed below)
BGtable$DEany <- ifelse(row.names(BGtable)%in%row.names(DEtable),"+","-")

#first row is DEany protein with and without a SNP
row1 <- c(nrow(BGtable[BGtable$DEany=="+" & BGtable$anysnp == "+",]), nrow(BGtable[BGtable$DEany == "+" & BGtable$anysnp == "-",])) 
#61 and 208 without isoforms
#96 and 173 with isoforms

#second row
row2 <- c(nrow(BGtable[BGtable$DEany=="-" & BGtable$anysnp == "+",]), nrow(BGtable[BGtable$DEany == "-" & BGtable$anysnp == "-",])) 
#290 and 1432 without isoforms
#481 and 1241 with isoforms

# This is going to be close...
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g")
result$p.value

#p-value = 0.006206
#two sided p = 0.01131


##****************************using full counts as background************************************8
anysnp <- multExpanded1_withDE[multExpanded1_withDE$Ind18486_SNPiso == "+" | multExpanded1_withDE$Ind18862_SNPiso == "+" | multExpanded1_withDE$Ind19160_SNPiso == "+",]
#make background as data frame and add snp annotation
background <- as.data.frame(background2)#background2 is the subtoDE df but with 'rev' entries removed.
colnames(background) <- "Protein"
background$anysnp <- ifelse(background$Protein %in% anysnp$Protein, "+","-")
#add the site ids as a column in the data frame
background$idmult <- subtoDE$idmult
#add the DE in any contrast information at the site level!
GlobalFs <- subtoDE[subtoDE$globalFsig=="+",]#1355 obs
background$DEany <- ifelse(background$idmult %in% GlobalFs$idmult, "+", "-")

#create the first row of the contingency matrix
row1 <- c(nrow(background[background$DEany=="+" & background$anysnp == "+",]), nrow(background[background$DEany == "+" & background$anysnp == "-",])) 
# 530 and 824 without specific isoform annotation (appropriate given this is bottom up MS) #note that there is a blank protein which may be causing the slight counting issue


#second row
row2 <- c(nrow(background[background$DEany=="-" & background$anysnp == "+",]), nrow(background[background$DEany == "-" & background$anysnp == "-",])) 
# [1] 1072 2565

#FEtest
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g")
result$p.value

#p-value = 7.376e-11

##****************************using full counts as background and Yannicks full dataset***********8
#input Yannick's full dataset
#import yannick files
yannick18486full <- read.table("C:/Users/Brett/Desktop/pilot_18486_prot.csv",sep = ',', header = T)
yannick18862full <- read.table("C:/Users/Brett/Desktop/pilot_18862_prot.csv",sep = ',', header = T)
yannick19160full <- read.table("C:/Users/Brett/Desktop/pilot_19160_prot.csv",sep = ',', header = T)


#subsetting to first six characters to ensure isoforms get annoatated. Ignoring isoform maps from phosphopeptides for the moment
multExpanded1_withDE$Ind18486_SNPisofull = ifelse(substr(multExpanded1_withDE$Protein,1,6) %in% yannick18486full$UNIPROT.ID, "+","-")
multExpanded1_withDE$Ind18862_SNPisofull = ifelse(substr(multExpanded1_withDE$Protein,1,6) %in% yannick18862full$UNIPROT.ID, "+","-")
multExpanded1_withDE$Ind19160_SNPisofull = ifelse(substr(multExpanded1_withDE$Protein,1,6) %in% yannick19160full$UNIPROT.ID, "+","-")

#add the phosphomotif information



anysnp <- multExpanded1_withDE[multExpanded1_withDE$Ind18486_SNPisofull == "+" | multExpanded1_withDE$Ind18862_SNPisofull == "+" | multExpanded1_withDE$Ind19160_SNPisofull == "+",]
#make background as data frame and add snp annotation
background <- as.data.frame(background2)#background2 is the subtoDE df but with 'rev' entries removed.
colnames(background) <- "Protein"
background$anysnp <- ifelse(background$Protein %in% anysnp$Protein, "+","-")
#add the site ids as a column in the data frame
background$idmult <- subtoDE$idmult
#add the DE in any contrast information at the site level!
GlobalFs <- subtoDE[subtoDE$globalFsig=="+",]#1355 obs
background$DEany <- ifelse(background$idmult %in% GlobalFs$idmult, "+", "-")

#create the first row of the contingency matrix
row1 <- c(nrow(background[background$DEany=="+" & background$anysnp == "+",]), nrow(background[background$DEany == "+" & background$anysnp == "-",])) 
# 530 and 824 without specific isoform annotation (appropriate given this is bottom up MS) #note that there is a blank protein which may be causing the slight counting issue
#now row1  558  796

#second row
row2 <- c(nrow(background[background$DEany=="-" & background$anysnp == "+",]), nrow(background[background$DEany == "-" & background$anysnp == "-",])) 
# [1] 1072 2565
#now row2 1164 2473

#FEtest
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g")
result$p.value

#p-value = 9.664017e-10 

#some summary numbers at the site level
# contmatrix
snp  nosnp
DE     row1  558  796
not DE row2 1164 2473
# Total sites that map to a protein with a non-syn snp and those that don't (first and second numbers respectively)
# > colSums(contmatrix)
# [1] 1722 3269
#34% of the phosphopeptides map to a protein that has a nonsyn snp

#contrast matrix for proteins
96 173
481 1241
protmatrix <- matrix(c(96,173,481,1241),nrow=2,ncol=2)
colSums(protmatrix)
[1]  269 1722
#13% of the proteins have a snp (yet 34% of the phosphopeptides map to a snp protein. proteins with snp tend to be longer and therefore have more phosphorylations)

##now for the tests of regions near the phosphosite***********************************************************888
#Are motif altering snps enriched for Diffphos?

anysnp <- multExpanded1_withDE[multExpanded1_withDE$Ind18486_SNPisofull == "+" | multExpanded1_withDE$Ind18862_SNPisofull == "+" | multExpanded1_withDE$Ind19160_SNPisofull == "+",]
#make background as data frame and add snp annotation
##breakdown to the level of protein for categorical enrichment analysis
background <- subtoDE$Protein
background <- as.data.frame(background)#background2 is the subtoDE df but with 'rev' entries removed.
colnames(background) <- "Protein"
background$anysnp <- ifelse(background$Protein %in% anysnp$Protein, "+","-")
#add the site ids as a column in the data frame
background$idmult <- subtoDE$idmult
#add the DE in any contrast information at the site level!
GlobalFs <- subtoDE[subtoDE$globalFsig=="+",]#1355 obs
background$DEany <- ifelse(background$idmult %in% GlobalFs$idmult, "+", "-")

#add the phosphovariant information to the anysnp dataframe, then subset this dataframe if variant present in any line
yannick18486motifvar <- yannick18486full[yannick18486full$Phosphomotif.positive..variant.is.within.7.Amino.acids.N.or.C.terminal.to.central.S.T.Y.=="yes",]##871 observations
yannick18862motifvar <- yannick18862full[yannick18862full$Phosphomotif.positive..variant.is.within.7.Amino.acids.N.or.C.terminal.to.central.S.T.Y.=="yes",]#487 observations
yannick19160motifvar <- yannick19160full[yannick19160full$Phosphomotif.positive..variant.is.within.7.Amino.acids.N.or.C.terminal.to.central.S.T.Y.=="yes",]#372 observations
###
anysnp$Ind18486_motifproximal = ifelse(substr(anysnp$Protein,1,6) %in% yannick18486motifvar$UNIPROT.ID, "+","-")
anysnp$Ind18862_motifproximal = ifelse(substr(anysnp$Protein,1,6) %in% yannick18862motifvar$UNIPROT.ID, "+","-")
anysnp$Ind19160_motifproximal = ifelse(substr(anysnp$Protein,1,6) %in% yannick19160motifvar$UNIPROT.ID, "+","-")
#now subset to get all the 'motif'variants
motifvars <- anysnp[anysnp$Ind18486_motifproximal == "+" | anysnp$Ind18862_motifproximal == "+" | anysnp$Ind19160_motifproximal == "+",]
#add PROTEIN LEVEL proximal annotation to 'background'
background$nearmotif <- ifelse(background$Protein %in% motifvars$Protein, "+", "-")#isoforms will match but weren't important for mapping
## only 586 sites are annotated as near the motif out of 5000 mapping to 151 unique proteins
## add pqtl information
background$pqtl <- ifelse(substr(background$Protein,1,6) %in% pqtl$UNIPROT.ID, "+","-")## got 168/278!!

#enrichment of motif proximal variants at the protein level
#create the first row of the contingency matrix
row1 <- c(nrow(background[background$DEany=="+" & background$nearmotif == "+",]), nrow(background[background$DEany == "+" & background$nearmotif == "-",])) 
#second row
row2 <- c(nrow(background[background$DEany=="-" & background$nearmotif == "+",]), nrow(background[background$DEany == "-" & background$nearmotif == "-",])) 

# contmatrix
# [,1] [,2]
# row1  181 1173
# row2  405 3232

#FEtest
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g")
result$p.value
##just barely significant .01743108 

##Note that this assessment asks the question "are proteins with motif variants more likely than those without motif variants to be differentially phosphorylated?"

##"this is different than: "are phosphosites with a modification surrounding the residue more likely to be DiffPhos than those without a proximal modification?"

# The answers to these two questions speaks to the nature of PTM crosstalk.

#enrichment of pqtls in diff phos!!

#create the first row of the contingency matrix
row1 <- c(nrow(background[background$DEany=="+" & background$pqtl == "+",]), nrow(background[background$DEany == "+" & background$pqtl == "-",])) 
#second row
row2 <- c(nrow(background[background$DEany=="-" & background$pqtl == "+",]), nrow(background[background$DEany == "-" & background$pqtl == "-",])) 

# [,1] [,2]
# row1   50 1304
# row2  118 3519

#FEtest
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g")
result$p.value
## Not significant! 0.2422782 this is just "is a protein that has a pqtl (that may not be different across individuals) enriched in DiffPhos"

##I need "are pqtls - where there is a difference in alleles between individuals in this subset - enriched in DE subset?"




#******************************************now for the protein normalized data using all snps and full counts**************************
##subset to those phosphoobservations subject to DE (the background set of possibilities)
subtoDEpn <- multExpanded1_withDE[multExpanded1_withDE$SubtoDEpn == "+",]#3603
subtoDEpn <- subtoDEpn[!grepl("REV",subtoDEpn$Protein),]#now removed REVerse hits and still 3603
##breakdown to the level of protein for categorical enrichment analysis
background <- subtoDEpn$Protein
#add the snp annotation for this analysis
anysnp <- multExpanded1_withDE[multExpanded1_withDE$Ind18486_SNPisofull == "+" | multExpanded1_withDE$Ind18862_SNPisofull == "+" | multExpanded1_withDE$Ind19160_SNPisofull == "+",]
#make background as data frame and add snp annotation
background <- as.data.frame(background)#background2 is the subtoDE df but with 'rev' entries removed.
colnames(background) <- "Protein"
background$anysnp <- ifelse(background$Protein %in% anysnp$Protein, "+","-")
#add the site ids as a column in the data frame
background$idmult <- subtoDEpn$idmult
#add the DE in any contrast information at the site level!
GlobalFs <- subtoDEpn[subtoDEpn$globalFsigpn=="+",]#972 obs with this truncated set
background$DEany <- ifelse(background$idmult %in% GlobalFs$idmult, "+", "-")

#add the variant proximal annotation neglecting isoform information
anysnp$Ind18486_motifproximal = ifelse(substr(anysnp$Protein,1,6) %in% yannick18486motifvar$UNIPROT.ID, "+","-")
anysnp$Ind18862_motifproximal = ifelse(substr(anysnp$Protein,1,6) %in% yannick18862motifvar$UNIPROT.ID, "+","-")
anysnp$Ind19160_motifproximal = ifelse(substr(anysnp$Protein,1,6) %in% yannick19160motifvar$UNIPROT.ID, "+","-")
#now subset to get all the 'motif'variants
motifvars <- anysnp[anysnp$Ind18486_motifproximal == "+" | anysnp$Ind18862_motifproximal == "+" | anysnp$Ind19160_motifproximal == "+",]

#add PROTEIN LEVEL proximity annotation to 'background'
background$nearmotif <- ifelse(background$Protein %in% motifvars$Protein, "+", "-")#isoforms will match but weren't important for mapping


#create the first row of the contingency matrix to test enrichment of snp containing proteins in DE subset.
row1 <- c(nrow(background[background$DEany=="+" & background$anysnp == "+",]), nrow(background[background$DEany == "+" & background$anysnp == "-",])) 

row2 <- c(nrow(background[background$DEany=="-" & background$anysnp == "+",]), nrow(background[background$DEany == "-" & background$anysnp == "-",])) 


#FEtest
contmatrix <- rbind(row1,row2)
# [,1] [,2]
# row1  386  586
# row2  895 1736

result <- fisher.test(contmatrix, alternative = "g")
result$p.value
## 0.0009148777

##now for the tests of regions near the phosphosite
#create the first row of the contingency matrix
row1 <- c(nrow(background[background$DEany=="+" & background$nearmotif == "+",]), nrow(background[background$DEany == "+" & background$nearmotif == "-",])) 
#second row
row2 <- c(nrow(background[background$DEany=="-" & background$nearmotif == "+",]), nrow(background[background$DEany == "-" & background$nearmotif == "-",])) 
#FEtest
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g")
result$p.value
##not significant this time.
## 0.4805825 

# > contmatrix
# [,1] [,2]
# row1  122  850
# row2  327 2304





