#this script aims to assess if pQTLs are enriched in the differentially phosphorylated subset.

#load the pqtl file from Battle et al Science 2014

#requires the XLConnect package which also requires a JRE installation to work properly. This can be a slow solution apparently but will work on Unix machines
# http://www.oracle.com/technetwork/java/javase/downloads/jre8-downloads-2133155.html

install.packages("XLConnect")
library(XLConnect)
pqtl <- readWorksheetFromFile("E:/My Documents/Dropbox/Postdoc-Gilad/Zia/latest science/1260793_DatafileS1.xlsx", sheet = 4, header = T)


# retrieve enspid from pqtl ensgid ----------------------------------------

#how many of the ensgids are in the biomart ensembl75 GRCh37.p13/hg19 repository?
library("biomaRt")
library("seqinr")

#archive hsapiens mart
ensembl_75 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")


ensembl_75_pqtls <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_peptide_id','description','hgnc_id',
                                         'hgnc_symbol','gene_biotype', 'transcript_biotype', 'status', 'transcript_status', 'ccds'),
                            filters = "ensembl_gene_id", values = pqtl$ENSG, mart = ensembl_75)


#278 genes
length(unique(pqtl$ENSG)) - length(unique(ensembl_75_pqtls$ensembl_gene_id))#only one gene missing.

#from 278 pqtl genes, how many have at least one CCDs protein id and were therefore included in the database for the search?

ccds_index <- grep(pattern = "ccds", ensembl_75_pqtls$ccds, ignore.case = T)# this requires less writing than subsetting within the DF
ensembl_75_pqtlscds <- ensembl_75_pqtls[ccds_index,]
length(unique(ensembl_75_pqtlscds$ensembl_gene_id))
#275!!! Great ;)

##subset further based on gene/transcript biotype. start with 587 entries
#known or novel gene status index
kntranscript_index <- grep(pattern = "known|novel", ensembl_75_pqtlscds$transcript_status, ignore.case = T)#576 (removes putative transcripts)
ensembl_75_pqtlscds <- ensembl_75_pqtlscds[kntranscript_index,]

#subset to protein and transcript biotype - 'protein_coding' 
gene_biotype_index <- grep(pattern = "protein_coding", ensembl_75_pqtlscds$gene_biotype, ignore.case = T)#576 removes polymorphic pseudogenes
ensembl_75_pqtlscds <- ensembl_75_pqtlscds[gene_biotype_index,]

transcript_biotype_index <- grep(pattern = "protein_coding", ensembl_75_pqtlscds$transcript_biotype, ignore.case = T)#573 removes 3 NMD transcripts
ensembl_75_pqtlscds <- ensembl_75_pqtlscds[transcript_biotype_index,]

#still 275 possible hits
length(unique(ensembl_75_pqtlscds$ensembl_gene_id))



# adding pqtl presence/absence annotation to multExpanded1_withDE ---------

#match for any assigned protein from the protein group to a gene with a pqtl?
AnyMatchProtQTL <- function(queryproteins){
  #This function looks for a match between queryproteins and the proteins containing a non-synonymous snp
  #queryproteins are proteins assigned to the given phosphopeptide
  if(any(unlist(strsplit(as.character(queryproteins), ";")) %in% ensembl_75_pqtlscds$ensembl_peptide_id)){
    "+"
  }else{
    "-"
  }
}
#apply anymatch function to all 'leading proteins' and 'majority protein ids'
multExpanded1_withDE$pQTLPositive <- mapply(AnyMatchProtQTL, multExpanded1_withDE$Leading.proteins)
multExpanded1_withDE$GelPreppQTLPositive <- mapply(AnyMatchProtQTL, multExpanded1_withDE$ppMajorityProteinIDs)

#how many ENSPids from a pQTL gene per protein group are assigned to a phosphopeptide?
NumMatchesProtQTL <- function(queryproteins){
  #This function counts the matches between queryproteins and the proteins containing a non-synonymous snp
  #queryproteins are proteins assigned to the peptide
  sum(unlist(strsplit(as.character(queryproteins), ";")) %in% ensembl_75_pqtlscds$ensembl_peptide_id)
}

#apply count function to all 'leading proteins' and 'majority protein ids'
multExpanded1_withDE$pQTLCount <- mapply(NumMatchesProtQTL, multExpanded1_withDE$Leading.proteins)
multExpanded1_withDE$GelPreppQTLCount <- mapply(NumMatchesProtQTL, multExpanded1_withDE$ppMajorityProteinIDs)

# This is likely due to multiple protein isoforms corresponding to a single gene. How many pQTLs correspond are mapped to multiple genes (paralogs)?
table(multExpanded1_withDE$pQTLCount)
table(multExpanded1_withDE$GelPreppQTLCount)


# Diff Phos enrichment test -----------------------------------------------

# 1)  Test for enrichment in diffphos. background is all sites subject to DiffPhos. Foreground is omnibus F significance. Category is 'with pQTL' or without pQTL at the phosphopeptide level. That is the same protein is counted as many times as it has an identified phosphopeptide in the subtodiffphos datafram. Contingency matrix is of the form:

#                            in category  not in category  
#                    DErow
#                    NotDErow
# 'in' category is any majority/leading protein(s) assigned to this phosphosite has a pQTL.  

subtoDE <- multExpanded1_withDE[multExpanded1_withDE$SubtoDE == "+",] #4738
subtoDEpn <- multExpanded1_withDE[multExpanded1_withDE$SubtoDEpn == "+",] #3488

#confounded analysis
row1 <- c(nrow(subtoDE[subtoDE$globalFsig == "+" & subtoDE$pQTLPositive == "+",]), 
          nrow(subtoDE[subtoDE$globalFsig == "+" & subtoDE$pQTLPositive == "-",]))

row2 <- c(nrow(subtoDE[subtoDE$globalFsig == "-" & subtoDE$pQTLPositive == "+",]), 
          nrow(subtoDE[subtoDE$globalFsig == "-" & subtoDE$pQTLPositive == "-",]))

#FEtest
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g")
result$p.value
row1 
0.1776086  


#protnormalized analysis using Zia's data (much less significant and speaks to penetrance at the post translational level)
row1 <- c(nrow(subtoDEpn[subtoDEpn$globalFsigpn == "+" & subtoDEpn$GelPreppQTLPositive == "+",]), 
          nrow(subtoDEpn[subtoDEpn$globalFsigpn == "+" & subtoDEpn$GelPreppQTLPositive == "-",]))

row2 <- c(nrow(subtoDEpn[subtoDEpn$globalFsigpn == "-" & subtoDEpn$GelPreppQTLPositive == "+",]), 
          nrow(subtoDEpn[subtoDEpn$globalFsigpn == "-" & subtoDEpn$GelPreppQTLPositive == "-",]))

#FEtest
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g")
result$p.value
0.0006978427

#here is the relative proportion difference
apply(contmatrix,1,function(x) x[1]/sum(x))

row1       row2 
0.06484425 0.04020888  


# pQTL enrichment where alleles differ ------------------------------------
#subset the pQTL analysis such that there is a difference in alleles across the samples. One wouldn't expect a difference to emerge if the variants are the same for the three individuals analyzed here.

###########First I need to convert the hg18 coordinates used for pqtls to hg19 coordinates
#I will use the 'liftover' tool from ucsc browser track infrastructure via the 'liftOver' function within the package Rtracklayer
library(rtracklayer)
library(GenomicRanges)

#A 'chain file' is required to convert genomic coordinates from one reference assembly to another
#create 'ChainFiles' directory if it doesn't already exist
dir.create(file.path(getwd(), "ChainFiles"))

#download file and save the date (first pass on 5/1/15)
url <- "http://hgdownload-test.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz"
destfile <- "ChainFiles/hg18ToHg19.over.chain.gz"
download.file(url, destfile)
date_downld <- date()

#file unzipped locally with bash (gzfile doesn't work)
ch <- import.chain("./ChainFiles/hg18ToHg19.over.chain")

#create Granges object for use with liftover function
pqtl18 <- data.frame(chr = pqtl$chr, start = pqtl$hg18.pos, end = pqtl$hg18.pos)
pqtl18 <- makeGRangesFromDataFrame(pqtl18)
pqtl19 <- liftOver(pqtl18, ch)#GRangeslist object of 278 Granges objects. Spot check confirms concordance with online webapp

#convert to data frame
pqtl19 <- as.data.frame(pqtl19, row.names = NULL, optional = FALSE)

#append hg19 coordinates to pqtl file
pqtl$hg19.pos <- pqtl19$start


#######Identify pqtls that have differing alleles across any of the three individuals studied here
head(SNPeffFinal)
dim(SNPeffFinal)
length(unique(SNPeffFinal$snp))#17963 unique snps in the dataset where at least one of the three individuals has a minor alelle copy

#how many pqtls are found in at least one of the three individuals with the minor allele? 
sum(pqtl$hg19.pos %in% SNPeffFinal$pos)
# Fuck only 2!! 

# Do they make sense? That is do the chromosomes match? YES
pqtl[which(pqtl$hg19.pos %in% SNPeffFinal$pos),]
#first match chromosome verified (16)
SNPeffFinal[SNPeffFinal$pos == 89167094, c("chr", "pos")]
#second match chromosome verified (11)
SNPeffFinal[SNPeffFinal$pos == 499120, c("chr", "pos")]

#The low number of variants present here may be due in part to the 'snpefffinal' df being nonsynonymous variants only.

# Is there value in relying the top-n most significant snps for a gene as opposed to the very top snp?

#how many pqtls (of any type) are present in the three yoruba with varying genotypes?




