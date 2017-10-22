#this script aims to assess if pQTLs are enriched in the differentially phosphorylated subset.

#John's new file. Modified with excel (EM). Needs to have '#' removed in one of the column headers! 
pqtl <- read.table("D:/EMpqtl-genos.txt", sep = "\t", header = T, stringsAsFactors = F, quote = "")#278
pqtl <- read.table("E:/My Documents/Pilot/EMpqtl-genos.txt", sep = "\t", header = T, stringsAsFactors = F, quote = "")#278

#subset to those carried over to hg19
pqtl <- pqtl[!is.na(pqtl$CHROM),]

#subset to those present in any of our samples

#subset to 3 samples of interest: 18486, 18862, 19160, (19238 is the standard so is not included)
sampleNames <- grep("18486|18862|19160", names(pqtl), value = T)
variables <- names(pqtl)[1:15]
pqtl <- pqtl[,c(variables,sampleNames)]

#subset to those variants present in at least one of the three lines
hapTypes <- c("0|1", "1|0", "1|1")
#get logical index for subsetting
index <- apply(pqtl[,sampleNames], 1, function(x){
  any(hapTypes %in% x)})
pqtl <- pqtl[index,]
nrow(pqtl)#199


#subset to those variants where at least one line is different
index <- apply(pqtl[,sampleNames],1,function(x){
  x <- as.character(x)
  !all(sapply(x,identical, x[1]))
    }
  )
pqtl <- pqtl[index,]
nrow(pqtl) #181


### Assign pqtls to phosphosites ----

# I will use biomart to asssign ENSGids to ENSPids and from there to the phosphosites
library("biomaRt")
library("seqinr")

#archive hsapiens mart
ensembl_75 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

# Biomart entries for the pQTL genes
ensembl_75_pqtls <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_peptide_id','description','hgnc_id',
                                         'hgnc_symbol','gene_biotype', 'transcript_biotype', 'status', 'transcript_status', 'ccds'),
                          filters = "ensembl_gene_id", values = pqtl$ENSG, mart = ensembl_75)


length(unique(pqtl$ENSG)) == length(unique(ensembl_75_pqtls$ensembl_gene_id))#no genes missing. (181 of 181)

#from 181 pqtl genes, how many have at least one CCDs protein id and were therefore included in the database for the search?
ccds_index <- grep(pattern = "ccds", ensembl_75_pqtls$ccds, ignore.case = T)
ensembl_75_pqtlscds <- ensembl_75_pqtls[ccds_index,]
length(unique(ensembl_75_pqtlscds$ensembl_gene_id))
#179

##subset further based on gene/transcript biotype. start with 587 entries
#known or novel gene status index
kntranscript_index <- grep(pattern = "known|novel", ensembl_75_pqtlscds$transcript_status, ignore.case = T)#576 (removes putative transcripts)
ensembl_75_pqtlscds <- ensembl_75_pqtlscds[kntranscript_index,]

#subset to protein and transcript biotype - 'protein_coding' 
gene_biotype_index <- grep(pattern = "protein_coding", ensembl_75_pqtlscds$gene_biotype, ignore.case = T)#576 removes polymorphic pseudogenes
ensembl_75_pqtlscds <- ensembl_75_pqtlscds[gene_biotype_index,]

transcript_biotype_index <- grep(pattern = "protein_coding", ensembl_75_pqtlscds$transcript_biotype, ignore.case = T)#573 removes 3 NMD transcripts
ensembl_75_pqtlscds <- ensembl_75_pqtlscds[transcript_biotype_index,]

#still 179 possible hits
length(unique(ensembl_75_pqtlscds$ensembl_gene_id))




#apply a match test for any assigned protein from the protein group to a gene with a pqtl.
AnyMatchProtQTL <- function(queryproteins){
  #This function looks for a match between queryproteins and the proteins containing a non-synonymous snp
  #queryproteins are proteins assigned to the given phosphopeptide
  if(any(unlist(strsplit(as.character(queryproteins), ";")) %in% ensembl_75_pqtlscds$ensembl_peptide_id)){
    "+"
  }else{
    "-"
  }
}

#apply anymatch function to all 'leading proteins' and 'majority protein ids'. around 340 sites match to a protein with a pQTL
multExpanded1_withDE_annotated$pQTLPositive <- mapply(AnyMatchProtQTL, multExpanded1_withDE_annotated$Leading.proteins)
multExpanded1_withDE_annotated$GelPreppQTLPositive <- mapply(AnyMatchProtQTL, multExpanded1_withDE_annotated$ppMajorityProteinIDs)
multExpanded1_withDE_annotated$PhosPreppQTLPositive <- mapply(AnyMatchProtQTL, multExpanded1_withDE_annotated$PhosPrepMajorityProteinIDs)

#339 sites map to a pQTL protein
table(multExpanded1_withDE_annotated$GelPreppQTLPositive)

#80 proteins
length(unique(multExpanded1_withDE_annotated[multExpanded1_withDE_annotated$GelPreppQTLPositive == "+", "ppMajorityProteinIDs"]))


#how many sites/proteins subjected to diffphos?

#38 proteins
length(unique(multExpanded1_withDE_annotated[multExpanded1_withDE_annotated$GelPreppQTLPositive == "+" &
                                               multExpanded1_withDE_annotated$GelPrepCovSubtoDE == "+", 
                                             "ppMajorityProteinIDs"]))
#100 sites
nrow(multExpanded1_withDE_annotated[multExpanded1_withDE_annotated$GelPreppQTLPositive == "+" &
                                               multExpanded1_withDE_annotated$GelPrepCovSubtoDE == "+", ])



# Diff Phos enrichment tests -----------------------------------------------

# 1)  Test for enrichment in diffphos. background is all sites subject to DiffPhos. Foreground is omnibus F significance. Category is 'with pQTL' or without pQTL at the phosphopeptide level.

SubtoDEGelProt <- multExpanded1_withDE_annotated[multExpanded1_withDE_annotated$GelPrepCovSubtoDE == "+",] #3257
SubtoDEConfounded <- multExpanded1_withDE_annotated[multExpanded1_withDE_annotated$ConfoundedSubtoDE == "+",] #4738
SubtoDEPhosProt <- multExpanded1_withDE_annotated[multExpanded1_withDE_annotated$PhosPrepCovSubtoDE == "+",] #1308

#GelPrep analysis using Zia's data
row1 <- c(nrow(SubtoDEGelProt[SubtoDEGelProt$GelPrepCovglobalFsig == "+" & SubtoDEGelProt$GelPreppQTLPositive == "+",]), 
          nrow(SubtoDEGelProt[SubtoDEGelProt$GelPrepCovglobalFsig == "+" & SubtoDEGelProt$GelPreppQTLPositive == "-",]))

row2 <- c(nrow(SubtoDEGelProt[SubtoDEGelProt$GelPrepCovglobalFsig == "-" & SubtoDEGelProt$GelPreppQTLPositive == "+",]), 
          nrow(SubtoDEGelProt[SubtoDEGelProt$GelPrepCovglobalFsig == "-" & SubtoDEGelProt$GelPreppQTLPositive == "-",]))

#FEtest
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g")
result$p.value
# .00384 (10/18/2015)


# Threshold independent test of association ----

# Regress p-values against binary assignment of sites. Assess if significant spearman rank based correlation.
pqtl.matrix <- SubtoDEGelProt[, c("GelPreppQTLPositive", "GelPrepCovFPval")]

#switch to 0/1 designation
pqtl.matrix$GelPreppQTLPositive <- ifelse(pqtl.matrix$GelPreppQTLPositive == "+", 1, 0)
pqtl.matrix[] <- lapply(pqtl.matrix, as.numeric)

cor(pqtl.matrix, method = "spearman")
cor.test(pqtl.matrix[[1]], pqtl.matrix[[2]], method = "spearman", exact = F)$p.value
# [1] 0.02667636 (10/18/2015)


