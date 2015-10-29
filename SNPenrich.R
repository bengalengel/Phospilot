#here I am going to add SNP based annotation to multExpaned table for categorical enrichment test

# SNPenrich <- function(multExpanded1_withDE){
require(seqinr)
require(iterators)
require(foreach)
require(doParallel)
require(stringr)

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
SNPeffFinal <- SNPeffFinal[index,]
nrow(SNPeffFinal)#51427

#roughly 21K coding variants in at least 1 line (not including standard)
length(unique(SNPeffFinal$snp))#21147

#
length(unique(SNPeffFinal$gene))#9208
length(unique(SNPeffFinal$peptide))#25705





# Add sequence annotation information -------------------------------------

##add sequence information to each entry from the fasta file. Then add a flag column indicating proximity to nearest S/T/Y. Then a flag column indicating yes/no for +- 7 from S/T/Y. The same should be done with phosphorylated site and phosphosite annotation information.

# Just for the phosphosites observed; Are there snps in the vicinity of the phosphosite? The background for enrichmentwill be subtodiffphos + all snps.

#for each ENSPID find the squence contained within the FASTA file. Executed with an sapply call. 

proteome <- read.fasta( file = "./FASTA/Homo_sapiens.GRCh37.75.pep.all.parsedCCDS.fa", seqtype = "AA", as.string = TRUE)


cl <- makeCluster(5)#I have 8 cores but had a crash when using all 8
registerDoParallel(cl)
SNPeffFinal$ProteomeIndex <- foreach(i = 1:length(SNPeffFinal$peptide), .combine = "c") %dopar% {
  peptide <- SNPeffFinal$peptide[i]
  hit <- grep(peptide, names(proteome))
  if(length(hit)==0){
    "no match" 
    } else {
      return(hit)
    }
}
stopCluster(cl)



# Annotation: Distance to nearest phosphorylatable residue (on hold) ----------------

#Can be used for motif proximity enrichment analysis, that is for each protein/hit, is the variant within 7 residues of a S/T/Y? This can provide insight into potential longer range interactions, in cis.

#getsequence from 'hits' index and position from snpeff table. Return distance to nearest S/T/Y residue. 'NA' results from either no Tyr in sequence or protein mapped to snp is not searched in database. 

# 
# #for this work only missense variants are considered
# SNPeffFinalMSonly <- SNPeffFinal[SNPeffFinal$effect == "missense_variant" | SNPeffFinal$effect == "missense_variant&splice_region_variant",]#41,388 nonsynonymous snps
# 
# #Dist to closest Tyrosine  
# TyrDist <- function(ProtSeq, VariantPosition){
#   #if ProtSeq index is valid 
#   if(ProtSeq != "no match"){
#     #retrieve sequence and Tyr index
#     seq <- seqinr::getSequence(proteome[[as.numeric(ProtSeq)]])
#     Tyr <- grep("Y", seq)
#     #if sequence contains tyrosines
#     if(length(Tyr) > 0){
#       #retrive variant position
#       Variant <- as.integer(gsub("[^0-9]+", "", VariantPosition))#removes any non-numeric elements
#       #return minimum distance
#       min(abs(Variant - Tyr))}else{
#         NA
#       }
#   }else
#   {NA}
# }
# 
# #Dist to closest Serine  
# SerDist <- function(ProtSeq, VariantPosition){
#   #if ProtSeq index is valid 
#   if(ProtSeq != "no match"){
#     #retrieve sequence and Syr index
#     seq <- seqinr::getSequence(proteome[[as.numeric(ProtSeq)]])
#     Ser <- grep("S", seq)
#     #if sequence contains Serines
#     if(length(Ser) > 0){
#       #retrive variant position
#       Variant <- as.integer(gsub("[^0-9]+", "", VariantPosition))#removes any non-numeric elements
#       #return minimum distance
#       min(abs(Variant - Ser))}else{
#         NA
#       }
#   }else
#   {NA}
# }
# 
# #Dist to closest Threonine  
# ThrDist <- function(ProtSeq, VariantPosition){
#   #if ProtSeq index is valid 
#   if(ProtSeq != "no match"){
#     #retrieve sequence and Thr index
#     seq <- seqinr::getSequence(proteome[[as.numeric(ProtSeq)]])
#     Thr <- grep("T", seq)
#     #if sequence contains Serines
#     if(length(Thr) > 0){
#       #retrive variant position
#       Variant <- as.integer(gsub("[^0-9]+", "", VariantPosition))#removes any non-numeric elements
#       #return minimum distance
#       min(abs(Variant - Thr))}else{
#         NA
#       }
#   }else
#   {NA}
# }
# 
# #apply minimum distance functions
# SNPeffFinalMSonly$NearestTyrDist <- mapply(TyrDist, SNPeffFinalMSonly$ProteomeIndex, SNPeffFinalMSonly$aa)
# SNPeffFinalMSonly$NearestSerDist <- mapply(SerDist, SNPeffFinalMSonly$ProteomeIndex, SNPeffFinalMSonly$aa)
# SNPeffFinalMSonly$NearestThrDist <- mapply(ThrDist, SNPeffFinalMSonly$ProteomeIndex, SNPeffFinalMSonly$aa)

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

#apply closest snp to phosphorylation site function. UPDATED TO USE ONLY MISSENSE VARIANTS
multExpanded1_withDE_annotated$ClosestSNPtoSite <- mapply(DistToPhos, multExpanded1_withDE$Proteins, multExpanded1_withDE$Positions.within.proteins)

#apply closest snp to phosphorylation site function GelPrep assignments. UPDATED TO USE ONLY MISSENSE VARIANTS
multExpanded1_withDE_annotated$ClosestSNPtoSiteGelPrep <- mapply(DistToPhos, multExpanded1_withDE$ppProteinIDs, multExpanded1_withDE$ppPositionInProteins)


# calculate minimum of all the distances for each protein group
multExpanded1_withDE_annotated$ClosestSNPtoSiteMin <- sapply(multExpanded1_withDE_annotated$ClosestSNPtoSite, function(x){
  distances <- as.numeric(unlist(strsplit(x, ";")))
  if(!all(is.na(distances))){
    min(distances, na.rm = T)
  }else
  {NA}
})  

multExpanded1_withDE_annotated$ClosestSNPtoSiteMinGelPrep <- sapply(multExpanded1_withDE_annotated$ClosestSNPtoSiteGelPrep, function(x){
  distances <- as.numeric(unlist(strsplit(x, ";")))
  if(!all(is.na(distances))){
    min(distances, na.rm = T)
  }else
  {NA}
})  

# The minimal distance to the phosphorylation site significantly impacts variation in phosphorylation for that proximal site
GelPrep.distances <- multExpanded1_withDE_annotated[multExpanded1_withDE_annotated$GelPrepCovSubtoDE == "+",
                                                    c("GelPrepCovglobalFsig", "GelPrepCovFAdjPval", "ClosestSNPtoSiteMinGelPrep")]
y <- -log10(as.numeric(GelPrep.distances$GelPrepCovFAdjPval))
x <- log10(GelPrep.distances$ClosestSNPtoSiteMinGelPrep + 1)#

plot(x,y)
R <- cor(x, y, method = "pearson", use = "complete.obs")
R
cor.test(x,y)$p.value

#make and save plot
pdf("distance_pvalue_density.pdf", 7, 5)
smoothScatter(x,y, nbin = 150, bandwidth = 0.1,
              cex = .3,
              pch = 19, nrpoints = .15*length(x),
              colramp = colorRampPalette(c("white", "light gray", "dark gray", "red")),
              xlab = expression(log[10](AA~distance~between~SNP~and~phosphosite)),
              ylab = expression(-log[10](P~value)), lwd = 10,
              family = "serif"
)
reg.line <- lm(y~x, na.action = "na.omit")
abline(reg.line, lwd = 2, lty = 2)
text(3, 7.25, expression(R == -.12), col = "darkred", cex = 1, family = "serif") # rsquared and pvalue
text(3, 6.85, expression(p == 9.90e-06), col = "darkred", cex = 1, family = "serif")
dev.off()





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
multExpanded1_withDE_annotated$NsSnpPositive <- mapply(AnyMatchProtSNPeff, multExpanded1_withDE_annotated$Leading.proteins)
multExpanded1_withDE_annotated$GelPrepNsSnpPositive <- mapply(AnyMatchProtSNPeff, multExpanded1_withDE_annotated$ppMajorityProteinIDs)
multExpanded1_withDE_annotated$PhosPrepNsSnpPositive <- mapply(AnyMatchProtSNPeff, multExpanded1_withDE_annotated$PhosPrepMajorityProteinIDs)


#how many nonsyn snps per protein group assigned to a phosphopeptide?
NumMatchesProtSNPeff <- function(queryproteins){
  #This function counts the matches between queryproteins and the proteins containing a non-synonymous snp
  #queryproteins are proteins assigned to the peptide
  sum(unlist(strsplit(as.character(queryproteins), ";")) %in% SNPeffFinal$peptide)
}

#apply count function to all 'leading proteins' and 'majority protein ids'
multExpanded1_withDE_annotated$NsSnpCount <- mapply(NumMatchesProtSNPeff, multExpanded1_withDE_annotated$Leading.proteins)
multExpanded1_withDE_annotated$GelPrepNsSnpCount <- mapply(NumMatchesProtSNPeff, multExpanded1_withDE_annotated$ppMajorityProteinIDs)
multExpanded1_withDE_annotated$PhosPrepNsSnpCount <- mapply(NumMatchesProtSNPeff, multExpanded1_withDE_annotated$PhosPrepMajorityProteinIDs)


# Less than 1% of the phosphopeptides are mapped to a protein group that (collectively) contains > 1 snp #THIS NEEDS TO BE CORRECTED FOR UNIQUENESS
table(multExpanded1_withDE_annotated$NsSnpCount)

#For proteins mapped to peptides using protein prep data the picture is more complicated because more information is being used (more protein ids/peptide). That is the same snp is present in multiple isoforms, which cannot be disambiguated geven the shotgun level information. Proteotypic peptides are need for this.
table(multExpanded1_withDE_annotated$GelPrepNsSnpCount)
table(multExpanded1_withDE_annotated$PhosPrepNsSnpCount)



###pfam domain annotation ----
#use the pfam.boundary file created using interproscan program that calls HMMER3/pfam. Was searched on the proteome.

#read in and configure:
# Load output data from InterProScan. There should be 36 files!! For some reason cluster hangs on execution for a few of these...use "output2.tar.gz"
# lapply load into list and then do.call to make into a dataframe
if(!file.exists("./InterProScan/output")){dir.create("./InterProScan/output")}
untar("./InterProScan/output2.tar.gz", compressed = "g", exdir = "./InterProScan/output")
file.names <- list.files("./InterProScan/output", pattern = ".*.tsv")
pfam.boundary <- lapply(file.names, function(file){
  read.table(file.path(paste0(getwd(),"./InterProScan/output/", file)), sep = "\t", header = F, stringsAsFactors = F, quote = "", comment.char = "", skip = 3)
})
pfam.boundary <- do.call(rbind, pfam.boundary)
names(pfam.boundary) <- c("ENSPID", "MD5seq", "SeqLength", "Analysis", "PFamID", "Description", "Start", "Stop", "Evalue", "MatchStatus", "Date")

#subset snpeff dataframe s.t. only enspids matching the 'proteome' are returned
SNPeffFinal <- SNPeffFinal[SNPeffFinal$ProteomeIndex != "no match", ]



#binary assigment of snpeff AA mutation positions to in/out domain. For those within a domain a flag yes/no is given regarding phosphorylation relevant. A foreach solution: 
# col1: T/F in any domain 
# col2: T/F in any phospho relevant domain

#not easy to download so manually curated from pfam as of 8-24-15
pfam.phospho <- c("PF00498", "PF01846", "PF03166", "PF10401", "PF00244", "PF00533", "PF00400", "PF00659", "PF00397",#S/T
                  "PF00017", "PF08416", "PF00168",#Y
                  "PF00782", "PF00102", "PF13350", "PF06602", "PF04273", "PF03162", "PF14566", "PF14671", "PF04179", "PF05706", #phosphatase
                  "PF00069", "PF01636",  "PF07714", "PF03109", "PF03881", "PF06293", "PF01163", "PF01633", "PF10707", "PF06176", #kinase
                  "PF02958", "PF04655", "PF10009", "PF12260", "PF16474", "PF07914", "PF14531", "PF06734", "PF05445", "PF07387") #kinase



#trimmed to only pS/pT relevant. Removed Y binding and Y kinase. phosphatases include dual specificity members.
pfam.phospho.ST <- c("PF00498", "PF01846", "PF03166", "PF10401", "PF00244", "PF00533", "PF00400", "PF00659", "PF00397",#S/T
                     "PF00782", "PF06602", "PF04273", "PF14566", "PF14671", "PF04179", "PF05706", #phosphatase
                     "PF00069", "PF01636", "PF03109", "PF03881", "PF06293", "PF01163", "PF01633", "PF10707", "PF06176", #kinase
                     "PF02958", "PF04655", "PF10009", "PF12260", "PF16474", "PF07914", "PF14531", "PF06734", "PF05445", "PF07387") #kinase


cl <- makeCluster(cl)
registerDoParallel(cl)
snp.domain.boundary <- foreach(i = 1:length(SNPeffFinal$peptide), .combine = "rbind", .packages = "stringr") %dopar% {
  #retrieve protein and site
  protein <- SNPeffFinal$peptide[i]
  #alter to match 'pfam.boundary' DF
  protein <- gsub("[A-Z]", "", protein)
  protein <- gsub("(?<![0-9])0+", "", protein, perl = TRUE)
  #extract site information. Use the first number within the chacter string for insertions etc.
  site <- SNPeffFinal$aa[i]
  site <- as.numeric(str_extract(site, "[0-9]+"))
  #identify all domains matching the protein
  domains <- pfam.boundary[pfam.boundary$ENSPID == protein, c(5,7,8)]
  if(dim(domains)[1] != 0){
    #apply logical test of if protein/site combo is in ANY domain range (>=start and <=end)
    hits <- domains[site >= domains$Start & site <= domains$Stop, "PFamID"]
    in.domain <- length(hits) > 0
    phospho.relevant <- any(hits %in% pfam.phospho)
    phospho.relevant.ST <- any(hits %in% pfam.phospho.ST)
  }else{
    in.domain <- FALSE
    phospho.relevant <- FALSE
    phospho.relevant.ST <- FALSE
  }
  data.frame(variant.in.domain = in.domain, domain.phospho.relevant = phospho.relevant, domain.phospho.relevant.ST = phospho.relevant.ST)
}
stopCluster(cl)


SNPeffFinal <- cbind(SNPeffFinal, snp.domain.boundary)

#add these results to MEDF

cl <- makeCluster(5)
registerDoParallel(cl)
GelPrep.SNP.domain <- foreach(i = 1:length(multExpanded1_withDE_annotated$ppMajorityProteinIDs), .combine = "rbind") %dopar% {
  protein.group <- multExpanded1_withDE_annotated$ppMajorityProteinIDs[i]
  if(multExpanded1_withDE_annotated$GelPrepNsSnpPositive[i] == "+"){
    if(protein.group != ""){
      protein.group <- strsplit(protein.group, ";")
      protein.group <- as.character(unlist(protein.group))
      protein.group <- protein.group[!grepl("REV", protein.group)]#remove reverse entries
      #for each member of the group, does it contain a snp in a domain/phospho.domain?
      in.domain <- vector(mode = 'logical', length = length(protein.group))
      phospho.relevant <- vector(mode = 'logical', length = length(protein.group))
      phospho.relevant.ST <- vector(mode = 'logical', length = length(protein.group))
      for(protein in seq_along(protein.group)){
        in.domain[protein] <- ifelse(length(SNPeffFinal[SNPeffFinal$peptide == protein.group[protein], "variant.in.domain"]) > 0, 
                                     SNPeffFinal[SNPeffFinal$peptide == protein.group[protein], "variant.in.domain"], NA)
        
        phospho.relevant[protein] <- ifelse(length(SNPeffFinal[SNPeffFinal$peptide == protein.group[protein], "domain.phospho.relevant"]) > 0,
                                            SNPeffFinal[SNPeffFinal$peptide == protein.group[protein], "domain.phospho.relevant"], NA)
        
        phospho.relevant.ST[protein] <- ifelse(length(SNPeffFinal[SNPeffFinal$peptide == protein.group[protein],
                                                                  "domain.phospho.relevant.ST"]) > 0,
                                               SNPeffFinal[SNPeffFinal$peptide == protein.group[protein], "domain.phospho.relevant.ST"], NA)
      }
      data.frame(snp.in.domain = any(in.domain, na.rm = T), snp.domain.phospho.relevant = any(phospho.relevant, na.rm = T),
                 snp.domain.phospho.relevant.ST = any(phospho.relevant.ST, na.rm = T))
    }
  } else {
    data.frame(snp.in.domain = NA, snp.domain.phospho.relevant = NA, snp.domain.phospho.relevant.ST = NA)
  }
}
stopCluster(cl)

#quite a few with snps in domains
table(GelPrep.SNP.domain$snp.in.domain)
FALSE  TRUE 
3637   948 

#A large enough sample to perform enrichment analysis
table(GelPrep.SNP.domain$snp.domain.phospho.relevant)
FALSE  TRUE 
4522    63 

table(GelPrep.SNP.domain$snp.domain.phospho.relevant.ST)
FALSE  TRUE 
4544    41 

#add snp information back to mE_annotated table
multExpanded1_withDE_annotated <- cbind(multExpanded1_withDE_annotated, GelPrep.SNP.domain)


#and how many sites map to proteins with a snp
table(multExpanded1_withDE_annotated$GelPrepNsSnpPositive)
  -     + 
13189  4585 


#snps within a motif ----

#load the elm table with matching ENSPIDs (written out when adding annotations)
elm <- readRDS("./ELM/elmtable.rds")


cl <- makeCluster(5)
registerDoParallel(cl)
SNPeffFinal$snpInMotif <- foreach(i = 1:length(SNPeffFinal[[1]]), .combine = c, .packages = c("stringr", "seqinr") ) %dopar% {
  #retrieve protein and site
  protein <- SNPeffFinal$peptide[i]
  #extract site information. Use the first number within the chacter string for insertions etc.
  site <- SNPeffFinal$aa[i]
  site <- as.numeric(str_extract(site, "[0-9]+"))
  #protein match in the elm dataframe?
  index <- sapply(elm$ENSPID, function(x){
    elm.proteins <- unlist(strsplit(x, ";"))
    protein %in% elm.proteins    
  })
  match.elm <- elm[index, c("ENSPID", "motif")]
  if(dim(match.elm)[1] > 0){
    #retrieve the protein sequence and motif positions
    motif.hits <- match.elm$motif
    #where do these motifs reside within the matching enspid protein?
    proteome.index <- grep(protein, names(proteome))
    hit.sequence <- unlist(getSequence(proteome[[proteome.index]], as.string = T))
    #a matrix of start and stop positions for the motifs
    motif.indices <- do.call(rbind, str_locate_all(hit.sequence, motif.hits))
    # start end
    # [1,]   103 107
    # [2,]   106 110
    #Do these positions overlap with the SNP position?
    any(apply(motif.indices, 1, function(x) {
      site >= x[1] && site <= x[2]
    }))
  } else {
    FALSE
  }
}
stopCluster(cl)

#hardly any snps that disrupt a motif!
table(SNPeffFinal$snpInMotif)
FALSE  TRUE 
29412    32 

#add this info to the ME DF
cl <- makeCluster(5)
registerDoParallel(cl)
multExpanded1_withDE_annotated$GelPrepSNPInmotif <- foreach(i = 1:length(multExpanded1_withDE_annotated$ppMajorityProteinIDs), 
                                                            .combine = c) %dopar% {
  protein.group <- multExpanded1_withDE_annotated$ppMajorityProteinIDs[i]
  if(multExpanded1_withDE_annotated$GelPrepNsSnpPositive[i] == "+"){
    if(protein.group != ""){
      protein.group <- strsplit(protein.group, ";")
      protein.group <- as.character(unlist(protein.group))
      protein.group <- protein.group[!grepl("REV", protein.group)]#remove reverse entries
      #for each member of the group, does it contain a snp in a motif?
      snp.in.motif <- vector(mode = 'logical', length = length(protein.group))
      for(protein in seq_along(protein.group)){
        snp.in.motif[protein] <- any(SNPeffFinal[SNPeffFinal$peptide == protein.group[protein], "snpInMotif"])
      }
      snp.in.motif = any(snp.in.motif, na.rm = T)
    }
  } else {
    snp.in.motif = FALSE
  }
}
stopCluster(cl)

#only 21 phosphosites belong to a protein with a snp in an elm validated motif
table(multExpanded1_withDE_annotated$GelPrepSNPInmotif)

#from 5 unique proteins
length(unique(multExpanded1_withDE_annotated[multExpanded1_withDE_annotated$GelPrepSNPInmotif == T, "ppMajorityProteinIDs"]))


###SNP assignment relative to disordered regions -----

#Iupred amino acid level disorder 
source("./Disorder/iupredProcessing.R")#creates "Iupred" list of ENSPID dataframes with primary sequence disorder annotation. DisProb >.5 is considered disordered.


#two seconds faster to perform this with parallel processing
cl <- makeCluster(5)
registerDoParallel(cl)
system.time(SNPeffFinal$variant.in.disordered.region <- foreach(i = 1:length(SNPeffFinal[[1]]), .combine = c, .packages = "stringr") %dopar% {
  #retrieve protein and site
  protein <- SNPeffFinal$peptide[i]
  #extract site information. Use the first number within the chacter string for insertions etc.
  site <- SNPeffFinal$aa[i]
  site <- as.numeric(str_extract(site, "[0-9]+"))
  #find matching protein dataframe within the iupred list
  diso.pred <- Iupred[[which(names(Iupred)==protein)]]
  #find matching site and return ordered/disordered categorization. Note for stops etc the length will be 0.
  ifelse(length(diso.pred[diso.pred$Position == site, 3]) > 0, diso.pred[diso.pred$Position == site, 3] >= .5, NA)
})
stopCluster(cl)
  

#Add presence/absence of ANY SNP within a disordered region to the ME dataframe.
cl <- makeCluster(5)
registerDoParallel(cl)
multExpanded1_withDE_annotated$GelPrepAnySNPInDisorderedRegion  <- foreach(i = 1:length(multExpanded1_withDE_annotated$ppMajorityProteinIDs), .combine = c) %dopar% {
  protein.group <- multExpanded1_withDE_annotated$ppMajorityProteinIDs[i]
  if(multExpanded1_withDE_annotated$GelPrepNsSnpPositive[i] == "+"){
    if(protein.group != ""){
      protein.group <- strsplit(protein.group, ";")
      protein.group <- as.character(unlist(protein.group))
      protein.group <- protein.group[!grepl("REV", protein.group)]#remove reverse entries
      #for each member of the group, does it contain a snp in a disordered region?
      snp.in.disorder <- vector(mode = 'logical', length = length(protein.group))
      for(protein in seq_along(protein.group)){
        snp.in.disorder[protein] <- any(SNPeffFinal[SNPeffFinal$peptide == protein.group[protein], "variant.in.disordered.region"])
      }
      snp.in.disorder = any(snp.in.disorder, na.rm = T)
    }
  } else {
    snp.in.disorder = FALSE
  }
}
stopCluster(cl)



##############Enrichment tests ------
# 1)  Test for enrichment in diffphos. background is all sites subject to DiffPhos. Foreground is omnibus F significance. Category is 'with snp' or without snp at the phosphopeptide level. Contingency matrix is of the form:

#                            in category  not in category  
#                    DErow
#                    NotDErow
# 'in' category is any majority/leading protein(s) assigned to this phosphosite has a snp.  

SubtoDEGelProt <- multExpanded1_withDE_annotated[multExpanded1_withDE_annotated$GelPrepCovSubtoDE == "+",] #3257
SubtoDEConfounded <- multExpanded1_withDE_annotated[multExpanded1_withDE_annotated$ConfoundedSubtoDE == "+",] #4738
SubtoDEPhosProt <- multExpanded1_withDE_annotated[multExpanded1_withDE_annotated$PhosPrepCovSubtoDE == "+",] #1308

#GelPrep analysis using Zia's data
row1 <- c(nrow(SubtoDEGelProt[SubtoDEGelProt$GelPrepCovglobalFsig == "+" & SubtoDEGelProt$GelPrepNsSnpPositive == "+",]), 
          nrow(SubtoDEGelProt[SubtoDEGelProt$GelPrepCovglobalFsig == "+" & SubtoDEGelProt$GelPrepNsSnpPositive == "-",]))

row2 <- c(nrow(SubtoDEGelProt[SubtoDEGelProt$GelPrepCovglobalFsig == "-" & SubtoDEGelProt$GelPrepNsSnpPositive == "+",]), 
          nrow(SubtoDEGelProt[SubtoDEGelProt$GelPrepCovglobalFsig == "-" & SubtoDEGelProt$GelPrepNsSnpPositive == "-",]))

#FEtest
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g")
result$p.value
row1 
1.399411e-06 

# Threshold independent test of association using spearman rank cor coef. Here using nominal ps in the event I want to produce a qq plot
NSsnp.matrix <- SubtoDEGelProt[, c("GelPrepNsSnpPositive", "GelPrepCovFPval")]

#switch to 0/1 designation
NSsnp.matrix$GelPrepNsSnpPositive <- ifelse(NSsnp.matrix$GelPrepNsSnpPositive == "+", 1, 0)
NSsnp.matrix[] <- lapply(NSsnp.matrix, as.numeric)

plot(NSsnp.matrix[[1]], -log10(NSsnp.matrix[[2]]))
plot(-log10(NSsnp.matrix[[2]]), NSsnp.matrix[[1]])

# Working with negative transform where a positive association indicates enrichment.
cor(NSsnp.matrix[[1]], -log10(NSsnp.matrix[[2]]), method = "spearman")
cor(-log10(NSsnp.matrix[[2]]), NSsnp.matrix[[1]], method = "spearman")
0.09222988

#the correlation is significant.
cor.test(NSsnp.matrix[[1]], NSsnp.matrix[[2]], method = "spearman", exact = F)$p.value
[1] 1.343636e-07


#bootstrap based estimate (http://content.csbs.utah.edu/~rogers/datanal/labprj/bootstrap/index.html)
robs <- cor(-log10(NSsnp.matrix[[2]]), NSsnp.matrix[[1]], method = "spearman") # observed r
tail.prob <- 0
nreps <- 500
for(i in 1:nreps) {
  y <- sample(NSsnp.matrix[[1]])   # randomly reorder values of 2nd variable
  rsim <- cor(-log10(NSsnp.matrix[[2]]), y, method="spearman")     # simulated r
  if(abs(rsim) >= abs(robs)) {   # abs makes test 2-tailed
    tail.prob <- tail.prob + 1
  }
}
tail.prob <- tail.prob / nreps


# 2) snp in domain enrichment using snp positive as background


# Threshold independent test of association using spearman rank cor coef. Here using nominal ps in the event I want to produce a qq plot
NSsnp.domain.matrix <- SubtoDEGelProt[SubtoDEGelProt$GelPrepNsSnpPositive == "+", c("snp.in.domain", "GelPrepCovFPval")]

#switch to 0/1 designation. for now the NAs are a bug
NSsnp.domain.matrix$snp.in.domain <- ifelse(NSsnp.domain.matrix$snp.in.domain == T, 1, 0)
NSsnp.domain.matrix$snp.in.domain[is.na(NSsnp.domain.matrix$snp.in.domain)] <- 0

NSsnp.domain.matrix[] <- lapply(NSsnp.domain.matrix, as.numeric)

plot(NSsnp.domain.matrix[[1]], -log10(NSsnp.domain.matrix[[2]]))
plot(-log10(NSsnp.domain.matrix[[2]]), NSsnp.domain.matrix[[1]])

# Working with negative transform where a positive association indicates enrichment. Negative depletion. Here enrichment!
cor(NSsnp.domain.matrix[[1]], -log10(NSsnp.domain.matrix[[2]]), method = "spearman")
cor(-log10(NSsnp.domain.matrix[[2]]), NSsnp.domain.matrix[[1]], method = "spearman")
[1] 0.1158094 

#the correlation IS significant.
cor.test(NSsnp.domain.matrix[[1]], NSsnp.domain.matrix[[2]], method = "spearman", exact = F)$p.value
[1] 1.409406e-05

#Few proteins have at least one snp within a domain
table(NSsnp.domain.matrix$snp.in.domain)
FALSE  TRUE 
1111   288 


# Threshold independent test of association using spearman rank cor coef. Here using nominal ps in the event I want to produce a qq plot
NSsnp.phosphodomain.matrix <- SubtoDEGelProt[SubtoDEGelProt$GelPrepNsSnpPositive == "+" & SubtoDEGelProt$snp.in.domain == TRUE,
                                      c("snp.domain.phospho.relevant", "GelPrepCovFPval")]

#switch to 0/1 designation. for now the NAs are a bug
NSsnp.phosphodomain.matrix$snp.domain.phospho.relevant<- ifelse(NSsnp.phosphodomain.matrix$snp.domain.phospho.relevant == T, 1, 0)
NSsnp.phosphodomain.matrix$snp.domain.phospho.relevant[is.na(NSsnp.phosphodomain.matrix$snp.domain.phospho.relevant)] <- 0

NSsnp.phosphodomain.matrix[] <- lapply(NSsnp.phosphodomain.matrix, as.numeric)

plot(NSsnp.phosphodomain.matrix[[1]], -log10(NSsnp.phosphodomain.matrix[[2]]))
plot(-log10(NSsnp.phosphodomain.matrix[[2]]), NSsnp.phosphodomain.matrix[[1]])

# Working with negative transform where a positive association indicates enrichment. Negative depletion. Here depletion
cor(NSsnp.phosphodomain.matrix[[1]], -log10(NSsnp.phosphodomain.matrix[[2]]), method = "spearman")
cor(-log10(NSsnp.phosphodomain.matrix[[2]]), NSsnp.phosphodomain.matrix[[1]], method = "spearman")
-0.02818142 

#the negative correlation is significant at alpha  = .05; 
cor.test(NSsnp.phosphodomain.matrix[[1]], NSsnp.phosphodomain.matrix[[2]], method = "spearman", exact = F)$p.value
0.6338825





# Threshold independent test of association using spearman rank cor coef. Here using nominal ps in the event I want to produce a qq plot
NSsnp.phosphodomain.matrix <- SubtoDEGelProt[SubtoDEGelProt$GelPrepNsSnpPositive == "+" & SubtoDEGelProt$snp.in.domain == TRUE,
                                             c("snp.domain.phospho.relevant.ST", "GelPrepCovFPval")]

#switch to 0/1 designation. for now the NAs are a bug
NSsnp.phosphodomain.matrix$snp.domain.phospho.relevant.ST<- ifelse(NSsnp.phosphodomain.matrix$snp.domain.phospho.relevant.ST == T, 1, 0)
NSsnp.phosphodomain.matrix$snp.domain.phospho.relevant.ST[is.na(NSsnp.phosphodomain.matrix$snp.domain.phospho.relevant.ST)] <- 0

NSsnp.phosphodomain.matrix[] <- lapply(NSsnp.phosphodomain.matrix, as.numeric)

plot(NSsnp.phosphodomain.matrix[[1]], -log10(NSsnp.phosphodomain.matrix[[2]]))
plot(-log10(NSsnp.phosphodomain.matrix[[2]]), NSsnp.phosphodomain.matrix[[1]])

# Working with negative transform where a positive association indicates enrichment. Negative depletion. 
cor(NSsnp.phosphodomain.matrix[[1]], -log10(NSsnp.phosphodomain.matrix[[2]]), method = "spearman")
cor(-log10(NSsnp.phosphodomain.matrix[[2]]), NSsnp.phosphodomain.matrix[[1]], method = "spearman")
[1] 0.03691468


#the negative correlation is significant at alpha  = .05; 
cor.test(NSsnp.phosphodomain.matrix[[1]], NSsnp.phosphodomain.matrix[[2]], method = "spearman", exact = F)$p.value
0.5326597






#motifs.....Not a significant sampling
table(SubtoDEGelProt$GelPrepSNPInmotif)
FALSE  TRUE 
3252     5


#Test of ANY snp within a disordered region affecting variability
NSsnp.disorder.matrix <- SubtoDEGelProt[SubtoDEGelProt$GelPrepNsSnpPositive == "+" , c("GelPrepAnySNPInDisorderedRegion", "GelPrepCovFPval")]

#switch to 0/1 designation. for now the NAs are a bug
NSsnp.disorder.matrix$GelPrepAnySNPInDisorderedRegion <- ifelse(NSsnp.disorder.matrix$GelPrepAnySNPInDisorderedRegion == T, 1, 0)
NSsnp.disorder.matrix$GelPrepAnySNPInDisorderedRegion[is.na(NSsnp.disorder.matrix$GelPrepAnySNPInDisorderedRegion)] <- 0

NSsnp.disorder.matrix[] <- lapply(NSsnp.disorder.matrix, as.numeric)

plot(NSsnp.disorder.matrix[[1]], -log10(NSsnp.disorder.matrix[[2]]))
plot(-log10(NSsnp.disorder.matrix[[2]]), NSsnp.disorder.matrix[[1]])

# Working with negative transform where a positive association indicates enrichment. Negative depletion. Here enrichment!
cor(NSsnp.disorder.matrix[[1]], -log10(NSsnp.disorder.matrix[[2]]), method = "spearman")
cor(-log10(NSsnp.disorder.matrix[[2]]), NSsnp.disorder.matrix[[1]], method = "spearman")
[1] -0.002778683

#the correlation IS NOT significant.
cor.test(NSsnp.disorder.matrix[[1]], NSsnp.disorder.matrix[[2]], method = "spearman", exact = F)$p.value
[1] 0.9172971

#most genes have at least one snp within a disordered region
table(NSsnp.disorder.matrix$GelPrepAnySNPInDisorderedRegion)
0    1 
297 1102 




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

SubtoDEConfounded <- multExpanded1_withDE[multExpanded1_withDE$SubtoDEConfounded == "+",] #4738
SubtoDEGelProt <- multExpanded1_withDE[multExpanded1_withDE$SubtoDEGelProt == "+",] #3488
SubtoDEPhosProt <- multExpanded1_withDE[multExpanded1_withDE$SubtoDEPhosProt == "+",] #1308

#confounded analysis
row1 <- c(nrow(SubtoDEConfounded[SubtoDEConfounded$globalFsigConfounded == "+" & SubtoDEConfounded$NsSnpPositive == "+",]), 
          nrow(SubtoDEConfounded[SubtoDEConfounded$globalFsigConfounded == "+" & SubtoDEConfounded$NsSnpPositive == "-",]))

row2 <- c(nrow(SubtoDEConfounded[SubtoDEConfounded$globalFsigConfounded == "-" & SubtoDEConfounded$NsSnpPositive == "+",]), 
          nrow(SubtoDEConfounded[SubtoDEConfounded$globalFsigConfounded == "-" & SubtoDEConfounded$NsSnpPositive == "-",]))

#FEtest
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g")
result$p.value
6.107867e-06 


#GelPrep analysis using Zia's data
row1 <- c(nrow(SubtoDEGelProt[SubtoDEGelProt$globalFsigGelProt == "+" & SubtoDEGelProt$GelPrepNsSnpPositive == "+",]), 
          nrow(SubtoDEGelProt[SubtoDEGelProt$globalFsigGelProt == "+" & SubtoDEGelProt$GelPrepNsSnpPositive == "-",]))

row2 <- c(nrow(SubtoDEGelProt[SubtoDEGelProt$globalFsigGelProt == "-" & SubtoDEGelProt$GelPrepNsSnpPositive == "+",]), 
          nrow(SubtoDEGelProt[SubtoDEGelProt$globalFsigGelProt == "-" & SubtoDEGelProt$GelPrepNsSnpPositive == "-",]))

#FEtest
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g")
result$p.value
4.123789e-08 


#PhosPrep analysis
row1 <- c(nrow(SubtoDEPhosProt[SubtoDEPhosProt$globalFsigPhosProt == "+" & SubtoDEPhosProt$PhosPrepNsSnpPositive == "+",]), 
          nrow(SubtoDEPhosProt[SubtoDEPhosProt$globalFsigPhosProt == "+" & SubtoDEPhosProt$PhosPrepNsSnpPositive == "-",]))

row2 <- c(nrow(SubtoDEPhosProt[SubtoDEPhosProt$globalFsigPhosProt == "-" & SubtoDEPhosProt$PhosPrepNsSnpPositive == "+",]), 
          nrow(SubtoDEPhosProt[SubtoDEPhosProt$globalFsigPhosProt == "-" & SubtoDEPhosProt$PhosPrepNsSnpPositive == "-",]))

#FEtest
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g")
result$p.value
1.449032e-05 



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





