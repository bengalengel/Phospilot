BuildProt <- function(protein1){
#build the protein table required for pQTL analysis

#required packages
  require(plyr)
  require(limma)
  require(biomaRt)
  require(iterators)
  require(foreach)
  require(doParallel)
  require(parallel)
  library(compiler); setCompilerOptions(optimize=3); enableJIT(3)


#### add gene (ENSG) column for pc correction. ----------
#Protein groups with majority protein ids mapping to separate genes are removed. This approach leads to a reasonable confidence that the quantifications I am using derive from a single gene.

#the approach is to make a list of ENGIDS. The name of each list object are the enspids from the majority proteins column name.


#grab the ensembl75 database
#archive hsapiens mart
ensembl_75 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

#get the gene, transcript and protein ids
ensembl_75_CCDS <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_peptide_id'), 
                         filters = 'with_ccds', values = T, mart = ensembl_75)

#grab the majority proteins
MajorityProteins <- unique(as.character(protein1$Majority.protein.IDs))
MajorityProteins <- as.list(MajorityProteins)
ENSGids <- vector(mode = 'list', length = length(MajorityProteins))#vector function flexible for preallocation

cl <- makeCluster(4)#4 cores works well
registerDoParallel(cl)
#for protein groups with majority ids consisting wholly of isoforms (same gene) that ENSGID is returned. IF not, then 'NA' is returned.
ENSGids <- foreach(i=1:length(MajorityProteins)) %dopar% {
  ENSPIDs <- strsplit(MajorityProteins[[i]], ";")
  ENSPIDs <- as.character(unlist(ENSPIDs))
  ENSPIDs <- paste(ENSPIDs, collapse = "|")
  index <- grep(ENSPIDs, ensembl_75_CCDS$ensembl_peptide_id)
  ids <- ensembl_75_CCDS$ensembl_gene_id[index]
  ids <- unique(ids)
  if(length(ids) == 1 & ids[1] != "") {
    ids[1]
  }else {
    NA
  }
}
stopCluster(cl)

#How many 'NA'?
sum(unlist(lapply(ENSGids, is.na)))#96 lost

#assign gene id to ENSG column of protein1
protein1$ENSG <- unlist(ENSGids)

#subset to 'proteingene'
ProteinGene <- protein1[!is.na(protein1$ENSG),]#6117

# 72 genes were mapped to multiple (one mapped to three) protein groups with maj ids from the same gene. 
# These are unique quants of separate isoforms from the same genes 
# These are removed for now
length(which(duplicated(ProteinGene$ENSG)))
dup.genes <- ProteinGene[which(duplicated(ProteinGene$ENSG)), "ENSG"]

#for each gene in dup.genes keep the isoform with the most unique plus razor ids
#use rowname list to subset away
discard.list <- c()
for(gene in dup.genes){
  tmp <- ProteinGene[ProteinGene$ENSG == gene, c("Razor...unique.peptides", "ENSG")]
  #good luck future self :)
  tmp2 <- row.names(tmp[is.na(match(seq(1:nrow(tmp)), which.max(tmp$Razor...unique.peptides))),])
  discard.list <- c(discard.list,tmp2)
}
ProteinGene <- ProteinGene[!row.names(ProteinGene) %in% discard.list,]#73 duplicates are removed

####Use gene models to asssign gene stop and start regions for downstream pQTL based PC regression optimization--------
#assign transcription start and stop along with chromosome number using gene models
#gene models from http://eqtl.uchicago.edu/RNA_Seq_data/ensembl_4_30.gz 

## Load file with gene models
## The gene models used are here: http://eqtl.uchicago.edu/RNA_Seq_data/ensembl_4_30.gz
gene.models <- read.delim(gzfile("./PCregression/ensembl_4_30.gz"), stringsAsFactors=F, header=F)
## Set relevant column names.
names(gene.models) <- c("V1", "ENST", "chr", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "V12", "ENSG", "V14", "V15", "V16")

# Build a table with gene name, chromosome, transcript start, end, and strand.
# cl <- makeCluster(4)#4 cores works well
# registerDoParallel(cl)
gene.start.end <- function(cur.ensg) {  c(txStart=min(cur.ensg$txStart,na.rm=T), txEnd=max(cur.ensg$txEnd,na.rm=T) )}
gene.table <- ddply(gene.models, c("ENSG"), gene.start.end, .parallel=F)
gene.chr <- function(cur.ensg) {  c(chr=cur.ensg$chr[1]) }
gene.chr.table <- ddply(gene.models, c("ENSG"), gene.chr, .parallel=F)
# stopCluster(cl)

#Here I lose 126 protein groups. But this is a single merge double merge, much like a double rainbow. So thats cool.
ProteinGene <- merge(ProteinGene, merge(gene.table, gene.chr.table, by="ENSG"), by="ENSG")#5991

#remove entries that map to genes that map to X chrom or 
## Remove peptides that map to chrX or unassembled regions.
# ProteinGene <- ProteinGene[-grep("_", ProteinGene$chr, fixed=T),]#all mapped to assembled region given these are CCDSs
ProteinGene <- ProteinGene[-grep("chrX", ProteinGene$chr, fixed=T),]

#extract pheno.pro dataframe to run through Zia's pipeline
stuff <- ProteinGene[,c("txStart", "txEnd", "chr")]
goodstuff <- ProteinGene[,1:61]
pheno.pro <- cbind(goodstuff,stuff)

#invert and log transform pheno columns. Heavy is standard in this dataset
pheno.cols <- grep("HL", colnames(ProteinGene))
logratios <- -log2(pheno.pro[,pheno.cols])
pheno.pro[,pheno.cols] <- logratios

#replace HL with GM
names(pheno.pro) <- gsub("HL", "GM", names(pheno.pro))

## Normalize to mean/median (here I am using mean like my boy zia does)
pheno.cols <- grep("GM", names(pheno.pro), fixed=T)
for(pcol in pheno.cols) {
  mu <- mean(pheno.pro[,pcol], na.rm=T)
  pheno.pro[,pcol] <-  pheno.pro[,pcol] - mu
}

# #note the differenet sample
# boxplot(pheno.pro[,pheno.cols])
# 
# #what is up with that? which one is it? This is the same with pview processing. Moving on
# which.min(apply(pheno.pro[,pheno.cols], 2, function(x) median(x, na.rm = T)))
# # GM18871 
# 24 

#rename and fix duplicated sample. NAs from duplicated sample are simply inserted. One is chosen arbitrarily to represent the individual? The mean is not taken across culture replicates:

## This sample might have been mislabeled.
names(pheno.pro)[which(names(pheno.pro) == "GM19012")] <- "GM19102"

## GM18355 is actually GM18855.
if(sum(names(pheno.pro) == "GM18355") > 0) {
  run1 <- pheno.pro$GM18855
  run2 <- pheno.pro$GM18355
  for(i in 1:length(run1)) {
    if(is.na(run1[i]) & (is.na(run2[i]) == F)) {
      run1[i] <- run2[i]
    }
  }
  pheno.pro$GM18855 <- run1 - mean(run1, na.rm=T)
  pheno.pro <- pheno.pro[,-grep("GM18355", names(pheno.pro), fixed=T)]
  pheno.cols <- grep("GM", names(pheno.pro), fixed=T)
}

#remove if not in at least half the samples
pheno.cols <- grep("GM", names(pheno.pro), fixed=T)
cnt.na <- function(i) {
  sum(is.na(pheno.pro[i,pheno.cols]) == F)
}
N.cnt <- unlist(lapply(1:dim(pheno.pro)[1], cnt.na))
pheno.pro <- pheno.pro[N.cnt >= length(pheno.cols)/2,]

#A bit of sample name changes were performed in Zia's script....!? Mine are labelled correctly.

## Sort table columns by name.
do.table.sort <- function(pheno.table) {
  ## Sort phenotype columns so they match genotype columns
  pheno.cols <- grep("GM", names(pheno.table), fixed=T)
  pheno.sort <- data.frame(ids=names(pheno.table)[pheno.cols], idx=pheno.cols, stringsAsFactors=F)
  pheno.cols.sort <- pheno.sort[sort(pheno.sort$ids, index.return=T)$ix,]$idx
  pheno.table[, pheno.cols] <- pheno.table[,pheno.cols.sort]
  names(pheno.table)[pheno.cols] <- names(pheno.table[pheno.cols.sort])
  pheno.cols <- grep("GM", names(pheno.table), fixed=T)
  return(pheno.table)
}
pheno.pro <- pheno.pro[sort(pheno.pro$ENSG,index.return=T)$ix,]
pheno.pro <- do.table.sort(pheno.pro)

save(pheno.pro, file="./PCregression/PRO_raw.RData")
}

