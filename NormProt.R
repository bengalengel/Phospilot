NormProt <- function(directory){
  #this function loads MQoutput from SILAC quantitative proteomic analysis of 60 LCLsm reduces it to the three samples of interest, median and quantile normalizes it. Returns a list of three data frames. MQoutput, median normalized, and quantile normalized. 
  source("loadMQZ.R")
  require(plyr)
  require(limma)
  require(biomaRt)
  require(iterators)
  require(foreach)
  require(doParallel)
  require(parallel)
  library(compiler); setCompilerOptions(optimize=3); enableJIT(3)
  
  # load protein files with particular variables populated using "loadMQ"
  protein <- load.MQZ(directory)#7120 protein groups
  
  # remove contaminants and reverse database hits
  protein <- protein[(protein$Potential.contaminant != "+" & protein$Reverse != "+"),]#6768
  
  # "only identified by site" hits CAN BE removed because they tend to have lower PEPs (wouldn't pass the FDR TH anyway) and can't be quantified since they are not idd by non-modified peptides. 
  protein1 <- protein[(protein$Only.identified.by.site != "+"),]#6465
  
  colnames(protein1)<- gsub(colnames(protein1), pattern = "Ratio.H.L.normalized.", replacement = "HL") ##remove redundant information
  
  #some strangeness sometimes there are two extra rows!
  #protein1 <- protein1[2:length(protein1)]
  
  #remove proteins if not quantified in at least one sample
  expCol <- grep("HL(.*)", colnames(protein1))
  
  ##removes rows containing all NA's using the sums of the logical per row 
  protein1 <- protein1[rowSums(is.na(protein1[,expCol]))!=length(expCol),]#6213
  
  #### add gene (ENSG) column for pc correction. ----------
#   Protein groups with majority protein ids mapping to separate genes are removed. This approach leads to a reasonable confidence that the quantifications I am using derive from a single gene.
  
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

####Use gene models to asssign gene stop and start regions for downstream pQTL based PC regression optimization--------
  #assign transcription start and stop along with chromosome number using gene models
  #gene models from http://eqtl.uchicago.edu/RNA_Seq_data/ensembl_4_30.gz 

## Load file with gene models
## The gene models used are here: http://eqtl.uchicago.edu/RNA_Seq_data/ensembl_4_30.gz
gene.models <- read.delim(gzfile("ensembl_4_30.gz"), stringsAsFactors=F, header=F)
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

#Here I loose 126 protein groups. But this is a single merge double merge, much like a double rainbow. So thats cool.
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

save(pheno.pro, file="PRO_raw.RData")

# PC regression based pQTL optimization ---------------------
#I will source the necessary Zia scripts to build the genotype files, regress PCs, and calculate pQTLs. 

#Files needed: 
## Get iso data files.
mean.genotype.files <- Sys.glob("E:/My Documents/R/ZiaScriptsEdited/data/genotypes_imputed/*.mean.genotype.txt")
snpdata.files <- Sys.glob("E:/My Documents/R/ZiaScriptsEdited/data/genotypes_imputed/*.mean.genotype.txt")

pull.chr <- function(x) {
  res <- gsub(".+chr", "chr", x,  perl=T)
  res <- gsub("\\..+", "", res, perl=T)
  res
}
snpdata.files <- data.frame(chr=unlist(lapply(snpdata.files, pull.chr)), snpdata.file=snpdata.files, stringsAsFactors=F)
mean.genotype.files <- data.frame(chr=unlist(lapply(mean.genotype.files, pull.chr)), mean.genotype.file=mean.genotype.files, stringsAsFactors=F)

genotype.files <- merge(snpdata.files, mean.genotype.files)

## Load individuals by column
con  <- file("E:/My Documents/R/ZiaScriptsEdited/data/genotypes_imputed/individual.order.txt", open = "r")
ind.line <- readLines(con, n=1)
ind.line <- gsub("IND ", "", ind.line, fixed=T)
close(con)

load.genotypes <- function(i)  {
  snpdata.file <- genotype.files[i,]$snpdata.file
  mean.genotype.file <- genotype.files[i,]$mean.genotype.file
  chr <- genotype.files[i,]$chr
  
  mean.genotype <- read.delim(mean.genotype.file, header=F, stringsAsFactors=F)
  
  ## Get rid of over hanging column.
  last.column <- dim(mean.genotype)[2]
  mean.genotype <- mean.genotype[,1:(last.column-1)]
  gc()
  
  ## Use individual names for columns.
  sample.names <- gsub("NA", "GM", unlist(strsplit(ind.line, split=",", fixed=T)), fixed=T)
  names(mean.genotype)[4:(last.column-1)] <- sample.names
  
  ## Rename initial columns
  names(mean.genotype)[1:3] <- c("snp.name", "A.allele", "B.allele")
  
  ## Load SNP position data and rename columns set chromosome number
  snp.data <- read.delim(snpdata.file, stringsAsFactors=F, skip=1)
  names(snp.data) <- c("snp.name", "A.allele.snpdata", "B.allele.snpdata", "minor.allele.freq", "chromosome", "hg18.pos")
  snp.data$chromosome <- chr
  
  ## 0-0.5] = BB allele
  ## (0.5, 1.5] = heterozygous, AB alleles
  ## (1.5, 2.5] = AA allele
  
  print(dim(mean.genotype))
  print(names(mean.genotype))
  ## Merge data sets on SNP identifier and sort the rows/SNPs based on position in the genome.
  genotype.data <- merge(snp.data, mean.genotype, by="snp.name")
  print(names(genotype.data))  
  rm(mean.genotype)
  rm(snp.data)
  gc()
  genotype.data <- genotype.data[sort(genotype.data$hg18.pos, index.return=T)$ix,]
  print(dim(genotype.data))
  if(!file.exists("data/geno_data")){
    dir.create(file.path(getwd(), "data/geno_data"))
  }
  save(genotype.data, file=paste("data/geno_data/",chr, ".RData", sep=""), compress="gzip", compression_level=1)
  NA
}

blank <- unlist(lapply(1:dim(genotype.files)[1], load.genotypes))# will not work with 6GB of memory.
# A foreach analog
# cl <- makeCluster(2)#This will depend on available memory (these files are rather large...). I can do this on laptop perhaps but not at home.
# registerDoParallel(cl)
# blank <- foreach(i = 1:dim(genotype.files)[1], .combine = 'c') %dopar% load.genotypes(i)
# stopCluster(cl)














  data <- protein1[,expCol]#60 cell lines
  
  row.names(data) <- protein1$id
  data <- 1/data   #inverse becausebecause the Heavy sample is the standard for the proteomics work
  data <- log2(data)
  
  #subset data to only the samples of interest
  data <- data[,c("HL18862","HL18486","HL19160")]
  
  ##change column names to match inversion
  colnames(data)<- gsub(colnames(data), pattern = "HL", replacement = "LH")
  
  #median normalize
  names <- colnames(data)
  median.subtract <- function(x){ x - median(x, na.rm = TRUE)}##create a function for median subtraction
  MedianNorm <- colwise(median.subtract, names)(data) #create median subtracted data but loose intensity and the row names
  
  #add back protien ids
  row.names(MedianNorm) <- protein1$id
  
  #summaries
  summary(MedianNorm)
  boxplot(MedianNorm)
  
  #remove if protein group not found in all samples
  MedianNorm <- na.omit(MedianNorm)#4270
  boxplot(MedianNorm)#differences in distribution shape for sure with HL18486 and HL19160
  par(mfrow = c(1, 1))
  for (i in 1:(ncol(MedianNorm))){
    if(i==1) plot(density(MedianNorm[, i], na.rm=T), col = i, ylim = c(0,2))
    else lines(density(MedianNorm[, i], na.rm=T), col = i)
  }
  #quantile normalize data being compared
  quantiled <- normalizeQuantiles(MedianNorm,ties = T)#ties are all assigned the same value for the common quantile
  summary(quantiled)
  boxplot(data)
  boxplot(quantiled)
  # density plots all look the same now of course
  plot.new()
  par(mfrow = c(1, 1))
  for (i in 1:(ncol(quantiled))){
    if(i==1) plot(density(quantiled[, i], na.rm=T), col = i, ylim = c(0,1.9))
    else lines(density(quantiled[, i], na.rm=T), col = i)
  }
  DFs <- list(protein1,data,MedianNorm,quantiled)
  return(DFs)
}
  
  