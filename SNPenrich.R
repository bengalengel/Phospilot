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
variables <- names(SNPeffFinal)[c(1:13, 17, 18)]
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
nrow(SNPeffFinal)##58078

# How many snps are only represented in the standard? 19238
sampleNames <- sampleNames[!grepl("NA19238",sampleNames)]
index <- apply(SNPeffFinal[,sampleNames], 1, function(x){
  any(hapTypes %in% x)})

#Because the standard line cannot contribute to the observed variation, variants unique to it are removed. (note if standard is homozygous positive for variant the peptide itself cannot be observed. See section on effect size estimates)
SNPeffFinal <- SNPeffFinal[index,]
nrow(SNPeffFinal)#51427


#subset to those variants where at least one non-standard line is different
index <- apply(SNPeffFinal[,sampleNames],1,function(x){
  x <- as.character(x)
  !all(sapply(x,identical, x[1]))
}
)
SNPeffFinal <- SNPeffFinal[index,]
nrow(SNPeffFinal) #46144



#19K coding variants in at least 1 line (not including standard)
length(unique(SNPeffFinal$snp))#19002
length(unique(SNPeffFinal$gene))#8656
length(unique(SNPeffFinal$peptide))#24036


## polyphen2 cleanup

# polyphen2 hdiv and hvar (humdiv model build from mendelian train set while humdiv was build from snp associations with disease training set).
# D - probably damaging
# P -  possibly damaging
# B - benign
# 
# The multiple scores correspond to multiple transcript forms (the annotation is based at the gene level). I take the the highest "score", which corresponds to
# D>P>B
# 
# A variant is assigned to the most damaging category



SNPeffFinal$Polyphen2_HDIV_pred <- sapply(SNPeffFinal$Polyphen2_HDIV_pred, function(x) {
  predictions <- unlist(strsplit(x, ","))
  if(any(predictions == "D")) {
    return("D") } else {
      if(any(predictions == "P")) {
        return("P") } else {
          if(any(predictions == "B")) {
            return("B") } else {
              return(NA)
            }
        }
    } 
} )


SNPeffFinal$Polyphen2_HVAR_pred <- sapply(SNPeffFinal$Polyphen2_HVAR_pred, function(x) {
  predictions <- unlist(strsplit(x, ","))
  if(any(predictions == "D")) {
    return("D") } else {
      if(any(predictions == "P")) {
        return("P") } else {
          if(any(predictions == "B")) {
            return("B") } else {
              return(NA)
            }
        }
    } 
} )



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
SNPeffFinalMSonly <- SNPeffFinal[SNPeffFinal$effect == "missense_variant" | SNPeffFinal$effect == "missense_variant&splice_region_variant",]#41,388 nonsynonymous snps
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
text(3, 7.25, expression(R == -.11), col = "darkred", cex = 1, family = "serif") # rsquared and pvalue
text(3, 6.85, expression(p == 4.37e-05), col = "darkred", cex = 1, family = "serif")
dev.off()


# subset to sites that map to disordered regions (GelPrep.Pos.Disorder = T)
GelPrep.distances.disorder <- multExpanded1_withDE_annotated[multExpanded1_withDE_annotated$GelPrepCovSubtoDE == "+" &
                                                               multExpanded1_withDE_annotated$GelPrep.Pos.Disorder == TRUE,
                                                    c("GelPrepCovglobalFsig", "GelPrepCovFAdjPval", "ClosestSNPtoSiteMinGelPrep")]
y.disorder <- -log10(as.numeric(GelPrep.distances.disorder$GelPrepCovFAdjPval))
x.disorder <- log10(GelPrep.distances.disorder$ClosestSNPtoSiteMinGelPrep + 1)#
disorder.line <- lm(y.disorder ~ x.disorder, na.action = "na.omit")
R.disorder <- cor(x.disorder, y.disorder, method = "pearson", use = "complete.obs")
R.disorder
cor.test(x.disorder,y.disorder)$p.value


GelPrep.distances.ordered <- multExpanded1_withDE_annotated[multExpanded1_withDE_annotated$GelPrepCovSubtoDE == "+" &
                                                               multExpanded1_withDE_annotated$GelPrep.Pos.Disorder == FALSE,
                                                             c("GelPrepCovglobalFsig", "GelPrepCovFAdjPval", "ClosestSNPtoSiteMinGelPrep")]
y.order <- -log10(as.numeric(GelPrep.distances.ordered$GelPrepCovFAdjPval))
x.order <- log10(GelPrep.distances.ordered$ClosestSNPtoSiteMinGelPrep + 1)#
order.line <- lm(y.order ~ x.order, na.action = "na.omit")
abline(order.line, lwd = 2, lty = 2)

R.order <- cor(x.order, y.order, method = "pearson", use = "complete.obs")
R.order
cor.test(x.order,y.order)$p.value

# plot with disordered phosphopeptides highlighted in red and ordered in black. Same with the regression lines.

pdf("distance_pvalue_cut_density.pdf", 7, 5)
plot(x, y,
     pch = 20,
     xlab = expression(log[10](AA~distance~between~SNP~and~phosphosite)),
     ylab = expression(-log[10](P~value)),
     family = "serif"
)

points(x.disorder, y.disorder,
       pch = 20,
       col = "red3"
)
abline(disorder.line, lwd = 2, lty = 2, col = "red3")
abline(order.line, lwd = 2, lty = 2, col = "black")
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




#### number of snps and phosphopeptide variation ----

GelPrep.snpcount <- multExpanded1_withDE_annotated[multExpanded1_withDE_annotated$GelPrepCovSubtoDE == "+",
                                                    c("GelPrepCovglobalFsig", "GelPrepCovFAdjPval", "GelPrepNsSnpCount")]
y <- -log10(as.numeric(GelPrep.snpcount$GelPrepCovFAdjPval))
x <- log10(GelPrep.snpcount$GelPrepNsSnpCount + 1)

# we are interested in a relationship specific to proteins that already have at least one snp.
GelPrep.snpcount2 <- GelPrep.snpcount[GelPrep.snpcount$GelPrepNsSnpCount > 0, ]
y <- -log10(as.numeric(GelPrep.snpcount2$GelPrepCovFAdjPval))
x <- log10(GelPrep.snpcount2$GelPrepNsSnpCount)

plot(x,y)
R <- cor(x, y, method = "pearson", use = "complete.obs")
R
cor.test(x,y)$p.value

#make and save plot if necessary
# pdf("nscount_pvalue_density.pdf", 7, 5)
# smoothScatter(x,y, nbin = 150, bandwidth = 0.1,
#               cex = .3,
#               pch = 19, nrpoints = .15*length(x),
#               colramp = colorRampPalette(c("white", "light gray", "dark gray", "red")),
#               xlab = expression(log[10](nsSNP~count)),
#               ylab = expression(-log[10](P~value)), lwd = 10,
#               family = "serif"
# )
# reg.line <- lm(y~x, na.action = "na.omit")
# abline(reg.line, lwd = 2, lty = 2)
# text(3, 7.25, expression(R == -.12), col = "darkred", cex = 1, family = "serif") # rsquared and pvalue
# text(3, 6.85, expression(p == 9.90e-06), col = "darkred", cex = 1, family = "serif")
# dev.off()







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



#Add presence/absence of ANY SNP (at least one) scored deleterious (polyphen HDIV) to the ME dataframe.
cl <- makeCluster(5)
registerDoParallel(cl)
multExpanded1_withDE_annotated$GelPrepAnySNPHDIVDeleterious  <- foreach(i = 1:length(multExpanded1_withDE_annotated$ppMajorityProteinIDs), .combine = c) %dopar% {
  protein.group <- multExpanded1_withDE_annotated$ppMajorityProteinIDs[i]
  if(multExpanded1_withDE_annotated$GelPrepNsSnpPositive[i] == "+"){
    if(protein.group != ""){
      protein.group <- strsplit(protein.group, ";")
      protein.group <- as.character(unlist(protein.group))
      protein.group <- protein.group[!grepl("REV", protein.group)]#remove reverse entries
      #for each member of the group, does it contain a snp scored HDIV deleterious?
      snp.hdiv.deleterious <- vector(mode = 'logical', length = length(protein.group))
      for(protein in seq_along(protein.group)){
        snp.hdiv.deleterious[protein] <- any(SNPeffFinal[SNPeffFinal$peptide == protein.group[protein], "Polyphen2_HDIV_pred"] == "D")
      }
      snp.hdiv.deleterious = any(snp.hdiv.deleterious, na.rm = T)
    }
  } else {
    snp.hdiv.deleterious = FALSE
  }
}
stopCluster(cl)


#Add presence/absence of ANY SNP (at least one) scored deleterious (polyphen HVAR) to the ME dataframe.
cl <- makeCluster(5)
registerDoParallel(cl)
multExpanded1_withDE_annotated$GelPrepAnySNPHVARDeleterious  <- foreach(i = 1:length(multExpanded1_withDE_annotated$ppMajorityProteinIDs), .combine = c) %dopar% {
  protein.group <- multExpanded1_withDE_annotated$ppMajorityProteinIDs[i]
  if(multExpanded1_withDE_annotated$GelPrepNsSnpPositive[i] == "+"){
    if(protein.group != ""){
      protein.group <- strsplit(protein.group, ";")
      protein.group <- as.character(unlist(protein.group))
      protein.group <- protein.group[!grepl("REV", protein.group)]#remove reverse entries
      #for each member of the group, does it contain a snp scored HVAR deleterious?
      snp.hvar.deleterious <- vector(mode = 'logical', length = length(protein.group))
      for(protein in seq_along(protein.group)){
        snp.hvar.deleterious[protein] <- any(SNPeffFinal[SNPeffFinal$peptide == protein.group[protein], "Polyphen2_HVAR_pred"] == "D")
      }
      snp.hvar.deleterious = any(snp.hvar.deleterious, na.rm = T)
    }
  } else {
    snp.hvar.deleterious = FALSE
  }
}
stopCluster(cl)


##############Enrichment tests ------
# 1)  Test for enrichment in diffphos. background is all sites subject to DiffPhos. Foreground is omnibus F significance. Category is 'with snp' or without snp at the phosphopeptide level. 

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
5.009443e-07 

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
[1] 0.0945307

#the correlation is significant.
cor.test(NSsnp.matrix[[1]], NSsnp.matrix[[2]], method = "spearman", exact = F)$p.value
[1] 6.483904e-08


# 2) snp in domain enrichment using snp positive as background

NSsnp.domain.matrix <- SubtoDEGelProt[SubtoDEGelProt$GelPrepNsSnpPositive == "+", c("snp.in.domain", "GelPrepCovFPval")]

#switch to 0/1 designation. for now the NAs are a bug
NSsnp.domain.matrix$snp.in.domain <- ifelse(NSsnp.domain.matrix$snp.in.domain == T, 1, 0)
NSsnp.domain.matrix$snp.in.domain[is.na(NSsnp.domain.matrix$snp.in.domain)] <- 0

NSsnp.domain.matrix[] <- lapply(NSsnp.domain.matrix, as.numeric)

plot(NSsnp.domain.matrix[[1]], -log10(NSsnp.domain.matrix[[2]]))
plot(-log10(NSsnp.domain.matrix[[2]]), NSsnp.domain.matrix[[1]])

# Working with negative transform where a positive association indicates enrichment. Negative depletion.
cor(NSsnp.domain.matrix[[1]], -log10(NSsnp.domain.matrix[[2]]), method = "spearman")
cor(-log10(NSsnp.domain.matrix[[2]]), NSsnp.domain.matrix[[1]], method = "spearman")
[1] 0.1208013

cor.test(NSsnp.domain.matrix[[1]], NSsnp.domain.matrix[[2]], method = "spearman", exact = F)$p.value
[1] 9.829453e-06


# 3) Is the phospho relevance of the domain driving the enrichment of domains within the ns.snp bg? (no)

# Threshold independent test of association using spearman rank cor coef. 
NSsnp.phosphodomain.matrix <- SubtoDEGelProt[SubtoDEGelProt$GelPrepNsSnpPositive == "+" & SubtoDEGelProt$snp.in.domain == TRUE,
                                      c("snp.domain.phospho.relevant", "GelPrepCovFPval")]

#switch to 0/1 designation.
NSsnp.phosphodomain.matrix$snp.domain.phospho.relevant<- ifelse(NSsnp.phosphodomain.matrix$snp.domain.phospho.relevant == T, 1, 0)
NSsnp.phosphodomain.matrix$snp.domain.phospho.relevant[is.na(NSsnp.phosphodomain.matrix$snp.domain.phospho.relevant)] <- 0

NSsnp.phosphodomain.matrix[] <- lapply(NSsnp.phosphodomain.matrix, as.numeric)

plot(NSsnp.phosphodomain.matrix[[1]], -log10(NSsnp.phosphodomain.matrix[[2]]))
plot(-log10(NSsnp.phosphodomain.matrix[[2]]), NSsnp.phosphodomain.matrix[[1]])

# Working with negative transform where a positive association indicates enrichment. Negative depletion.
cor(NSsnp.phosphodomain.matrix[[1]], -log10(NSsnp.phosphodomain.matrix[[2]]), method = "spearman")
cor(-log10(NSsnp.phosphodomain.matrix[[2]]), NSsnp.phosphodomain.matrix[[1]], method = "spearman")
-0.02818142 

#the negative correlation is not sig
cor.test(NSsnp.phosphodomain.matrix[[1]], NSsnp.phosphodomain.matrix[[2]], method = "spearman", exact = F)$p.value
0.6338825

##again,using only S/T specific domains

# Threshold independent test of association using spearman rank cor coef.
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
[1] 0.03404846

#the negative correlation is significant at alpha  = .05; 
cor.test(NSsnp.phosphodomain.matrix[[1]], NSsnp.phosphodomain.matrix[[2]], method = "spearman", exact = F)$p.value
[1] 0.5753744


# 4) Is the phospho relevancy of the protein an additional driver of variability?

# Threshold independent test of association using spearman rank cor coef. Here using nominal ps in the event I want to produce a qq plot
NSsnp.phosphoprotein.matrix <- SubtoDEGelProt[SubtoDEGelProt$GelPrepNsSnpPositive == "+",
                                             c("GelPrepPFamIDPhosphoST", "GelPrepCovFPval")]

#switch to 0/1 designation. for now the NAs are a bug
NSsnp.phosphoprotein.matrix$GelPrepPFamIDPhosphoST <- ifelse(NSsnp.phosphoprotein.matrix$GelPrepPFamIDPhosphoST == "yes", 1, 0)
NSsnp.phosphoprotein.matrix$GelPrepPFamIDPhosphoST[is.na(NSsnp.phosphoprotein.matrix$GelPrepPFamIDPhosphoST)] <- 0

NSsnp.phosphoprotein.matrix[] <- lapply(NSsnp.phosphoprotein.matrix, as.numeric)

plot(NSsnp.phosphoprotein.matrix[[1]], -log10(NSsnp.phosphoprotein.matrix[[2]]))
plot(-log10(NSsnp.phosphoprotein.matrix[[2]]), NSsnp.phosphoprotein.matrix[[1]])

# Working with negative transform where a positive association indicates enrichment. Negative depletion. 
cor(NSsnp.phosphoprotein.matrix[[1]], -log10(NSsnp.phosphoprotein.matrix[[2]]), method = "spearman")
cor(-log10(NSsnp.phosphoprotein.matrix[[2]]), NSsnp.phosphoprotein.matrix[[1]], method = "spearman")
[1] -0.006172978

#the negative correlation is significant at alpha  = .05; 
cor.test(NSsnp.phosphoprotein.matrix[[1]], NSsnp.phosphoprotein.matrix[[2]], method = "spearman", exact = F)$p.value
[1] 0.8219145

# 5) number of snps and phosphopeptide variability (see above)


#6) motifs.....Not a significant sampling
table(SubtoDEGelProt$GelPrepSNPInmotif)
FALSE  TRUE 
3252     5


#7) Does the presense of a SNP in a disordered region potentiate variability

#Test of ANY snp within a disordered region affecting variability
NSsnp.disorder.matrix <- SubtoDEGelProt[SubtoDEGelProt$GelPrepNsSnpPositive == "+" , c("GelPrepAnySNPInDisorderedRegion", "GelPrepCovFPval")]

NSsnp.disorder.matrix$GelPrepAnySNPInDisorderedRegion <- ifelse(NSsnp.disorder.matrix$GelPrepAnySNPInDisorderedRegion == T, 1, 0)
NSsnp.disorder.matrix$GelPrepAnySNPInDisorderedRegion[is.na(NSsnp.disorder.matrix$GelPrepAnySNPInDisorderedRegion)] <- 0

NSsnp.disorder.matrix[] <- lapply(NSsnp.disorder.matrix, as.numeric)

plot(NSsnp.disorder.matrix[[1]], -log10(NSsnp.disorder.matrix[[2]]))
plot(-log10(NSsnp.disorder.matrix[[2]]), NSsnp.disorder.matrix[[1]])

# Working with negative transform where a positive association indicates enrichment. Negative depletion.
cor(NSsnp.disorder.matrix[[1]], -log10(NSsnp.disorder.matrix[[2]]), method = "spearman")
cor(-log10(NSsnp.disorder.matrix[[2]]), NSsnp.disorder.matrix[[1]], method = "spearman")
[1] -0.008710016

#the correlation IS NOT significant.
cor.test(NSsnp.disorder.matrix[[1]], NSsnp.disorder.matrix[[2]], method = "spearman", exact = F)$p.value
[1] 0.7507934

#most genes have at least one snp within a disordered region
table(NSsnp.disorder.matrix$GelPrepAnySNPInDisorderedRegion)
0    1 
277 1055 

# 8) Do deleterious variants potentiate variability?
#Test of ANY snp HVAR deleterious affecting variability
NSsnp.hvar.matrix <- SubtoDEGelProt[SubtoDEGelProt$GelPrepNsSnpPositive == "+" , c("GelPrepAnySNPHVARDeleterious", "GelPrepCovFPval")]

#switch to 0/1 designation.
NSsnp.hvar.matrix$GelPrepAnySNPHVARDeleterious <- ifelse(NSsnp.hvar.matrix$GelPrepAnySNPHVARDeleterious == T, 1, 0)
NSsnp.hvar.matrix$GelPrepAnySNPHVARDeleterious[is.na(NSsnp.hvar.matrix$GelPrepAnySNPHVARDeleterious)] <- 0

NSsnp.hvar.matrix[] <- lapply(NSsnp.hvar.matrix, as.numeric)

plot(NSsnp.hvar.matrix[[1]], -log10(NSsnp.hvar.matrix[[2]]))
plot(-log10(NSsnp.hvar.matrix[[2]]), NSsnp.hvar.matrix[[1]])

# Working with negative transform where a positive association indicates enrichment. Negative depletion.
cor(NSsnp.hvar.matrix[[1]], -log10(NSsnp.hvar.matrix[[2]]), method = "spearman")
cor(-log10(NSsnp.hvar.matrix[[2]]), NSsnp.hvar.matrix[[1]], method = "spearman")
[1] 0.1235908

#the correlation is significant.
cor.test(NSsnp.hvar.matrix[[1]], NSsnp.hvar.matrix[[2]], method = "spearman", exact = F)$p.value
[1] 6.075066e-06

table(NSsnp.hvar.matrix$GelPrepAnySNPHVARDeleterious)
0    1 
1089  243

#Test of ANY snp HDIV deleterious affecting variability
NSsnp.hdiv.matrix <- SubtoDEGelProt[SubtoDEGelProt$GelPrepNsSnpPositive == "+" , c("GelPrepAnySNPHDIVDeleterious", "GelPrepCovFPval")]

#switch to 0/1 designation. for now the NAs are a bug
NSsnp.hdiv.matrix$GelPrepAnySNPHDIVDeleterious <- ifelse(NSsnp.hdiv.matrix$GelPrepAnySNPHDIVDeleterious == T, 1, 0)
NSsnp.hdiv.matrix$GelPrepAnySNPHDIVDeleterious[is.na(NSsnp.hdiv.matrix$GelPrepAnySNPHDIVDeleterious)] <- 0

NSsnp.hdiv.matrix[] <- lapply(NSsnp.hdiv.matrix, as.numeric)

plot(NSsnp.hdiv.matrix[[1]], -log10(NSsnp.hdiv.matrix[[2]]))
plot(-log10(NSsnp.hdiv.matrix[[2]]), NSsnp.hdiv.matrix[[1]])

# Working with negative transform where a positive association indicates enrichment. Negative depletion.
cor(NSsnp.hdiv.matrix[[1]], -log10(NSsnp.hdiv.matrix[[2]]), method = "spearman")
cor(-log10(NSsnp.hdiv.matrix[[2]]), NSsnp.hdiv.matrix[[1]], method = "spearman")
[1] 0.08658584

cor.test(NSsnp.hdiv.matrix[[1]], NSsnp.hdiv.matrix[[2]], method = "spearman", exact = F)$p.value
[1] 0.001561192

table(NSsnp.hdiv.matrix$GelPrepAnySNPHDIVDeleterious)



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

