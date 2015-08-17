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
# gene.models <- read.delim(gzfile("ensembl_4_30.gz"), stringsAsFactors=F, header=F)#home
gene.models <- read.delim(gzfile("D:/for_brett/ensembl_4_30.gz"), stringsAsFactors=F, header=F)#laptop
#cluster
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

save(pheno.pro, file="PRO_raw.RData")

# PC regression based pQTL optimization ---------------------
#I will source the necessary Zia scripts to build the genotype files, regress PCs, and calculate pQTLs. 

#Files needed: 
## Get iso data files.
# mean.genotype.files <- Sys.glob("E:/My Documents/R/ZiaScriptsEdited/data/genotypes_imputed/*.mean.genotype.txt")
# snpdata.files <- Sys.glob("E:/My Documents/R/ZiaScriptsEdited/data/genotypes_imputed/*.snpdata.txt")
mean.genotype.files <- Sys.glob("D:/for_brett/data/genotypes_imputed/*.mean.genotype.txt")
snpdata.files <- Sys.glob("D:/for_brett/data/genotypes_imputed/*.snpdata.txt")

pull.chr <- function(x) {
  res <- gsub(".+chr", "chr", x,  perl=T)
  res <- gsub("\\..+", "", res, perl=T)
  res
}
snpdata.files <- data.frame(chr=unlist(lapply(snpdata.files, pull.chr)), snpdata.file=snpdata.files, stringsAsFactors=F)
mean.genotype.files <- data.frame(chr=unlist(lapply(mean.genotype.files, pull.chr)), mean.genotype.file=mean.genotype.files, stringsAsFactors=F)

genotype.files <- merge(snpdata.files, mean.genotype.files)

## Load individuals by column
# con  <- file("E:/My Documents/R/ZiaScriptsEdited/data/genotypes_imputed/individual.order.txt", open = "r")
con  <- file("D:/for_brett/data/genotypes_imputed/individual.order.txt", open = "r")

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
  if(!file.exists("data")){
    dir.create(file.path(getwd(), "data"))
  }  
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

#alter genotype cache st MAF > X and samples the same as those identified via proteomics
pheno.data <- pheno.pro; suffix <- "pro"
## Get table with genotype files.
genotype.files <- Sys.glob("data/geno_data/chr*.RData")
pull.chr <- function(x) {
  res <- gsub(".+chr", "chr", x,  perl=T)
  res <- gsub("\\..+", "", res, perl=T)
  res
}
genotype.files <- data.frame(chr=unlist(lapply(genotype.files, pull.chr)), genotype.file=genotype.files, stringsAsFactors=F)

## Loads data for analyzed cell lines and filters based on MAF.
build.geno.cache <- function(g) {
  ## Load genotype data for given chromosome.
  chr <- genotype.files$chr[g]
  ##cache.file <- paste("data/geno_data_cache_", suffix, "_maf5/", chr, ".RData", sep="")
  cache.file <- paste("data/geno_data_cache_", suffix, "/", chr, ".RData", sep="")
  data.file <- genotype.files$genotype.file[g]
  load(data.file)
  geno.cols <- grep("GM", names(genotype.data), fixed=T)
  pheno.cols <- grep("GM", names(pheno.data), fixed=T)
  
  keep.geno <- merge(data.frame(ids=names(pheno.data)[pheno.cols]),
                     data.frame(ids=names(genotype.data)[geno.cols], idx=geno.cols))
  
  ## Remove columns that don't intersect phenotype data.
  remove.cols <- setdiff(geno.cols, keep.geno$idx)
  if(length(remove.cols) > 0) {
    genotype.data <- genotype.data[,-remove.cols]
  }
  ## Sort the columns based on name so they match pheno data.
  geno.cols <- grep("GM", names(genotype.data), fixed=T)
  geno.sort <- data.frame(ids=names(genotype.data)[geno.cols], idx=geno.cols, stringsAsFactors=F)
  geno.cols.sort <- geno.sort[sort(geno.sort$ids, index.return=T)$ix,]$idx
  genotype.data[, geno.cols] <- genotype.data[,geno.cols.sort]
  ## Don't forget to update column names after sort.
  names(genotype.data)[geno.cols] <- names(genotype.data[geno.cols.sort])
  geno.cols <- grep("GM", names(genotype.data), fixed=T)
  gc()
  ## Get final genotype matrix.
  geno.matrix <- as.matrix(genotype.data[,geno.cols])
  
  ## Compute minor allele frequency of each SNP.
  compute.maf <- function(r) {
    genotypes <- geno.matrix[r,]
    ## Genotype codes.
    ## 0-0.5] = BB allele
    ## (0.5, 1.5] = heterozygous, AB alleles
    ## (1.5, 2.5] = AA allele
    major.cnt <- sum(0 <= genotypes & genotypes <= 0.5)
    het.cnt <- sum(0.5 < genotypes & genotypes <= 1.5)
    minor.cnt <- sum(1.5 < genotypes & genotypes <= 2.5)
    N <- dim(geno.matrix)[2] 
    major.freq <- (2*major.cnt + het.cnt) / (2*N)
    minor.freq <- (2*minor.cnt + het.cnt) / (2*N)
    min(major.freq, minor.freq)
  }
  maf <- unlist(lapply(1:dim(geno.matrix)[1], compute.maf))
  
  ## Test SNPs only with cutoff MAF
  before.MAF.cut <- dim(genotype.data)[1]
  ##genotype.data <- genotype.data[maf > 0.05,]
  genotype.data <- genotype.data[maf > 0.1,]
  if(!file.exists("data/geno_data_cache_pro")){
    dir.create(file.path(getwd(), "data/geno_data_cache_pro"))
  }  
  save(genotype.data, file=cache.file)
  NA
}
unlist(lapply(1:dim(genotype.files)[1], build.geno.cache))

########


#make the PC regression and imputation function
library(impute)

load.pro <- function(numPCs, keep.NAs=TRUE, standardize=TRUE) {
  require(impute)
  load("PRO_raw.RData")
  if(numPCs == -1) {  return(pheno.pro)  }
  
  pheno.cols <- grep("GM", names(pheno.pro), fixed=T)
  pheno.matrix <- as.matrix(pheno.pro[,pheno.cols])
  
  ## Gene-level standardization
  ## Make each row Normal(0,1).
  if(standardize) {
    pheno.matrix <- t(apply(pheno.matrix,1,function(row){(row-mean(row,na.rm=T))/sd(row,na.rm=T)})) ##sd uses (n-1))
  }
  ## Column normalization by qqnorm, 
  pheno.matrix <- apply(pheno.matrix,2,function(kk){qqnorm(kk, plot.it=F)$x})
  pheno.matrix <- as.matrix(pheno.matrix)
  pheno.matrix.with.NA <- pheno.matrix
  
  ## Impute missing values by KNN.
  if(sum(is.na(pheno.matrix)) > 0 ) {
    pheno.matrix <- impute.knn(pheno.matrix)$data
  }
  if(numPCs > 0) {
    PCAobject <- prcomp(t(pheno.matrix))
    expected <- t(PCAobject$x[,1:numPCs]  %*% t(PCAobject$rotation[,1:numPCs]))
    pheno.matrix <- pheno.matrix - expected
  }
  
  if(keep.NAs) {
    for(i in 1:dim(pheno.matrix)[1]) {
      NA.idx <- which(is.na(pheno.matrix.with.NA[i,]))
      pheno.matrix[i,NA.idx] <- NA
    }
  }
  
  pheno.pro[,pheno.cols] <- pheno.matrix
  return(pheno.pro)
}


########run the regressions-----------
reg.name <- "pro"

## Sample code for computing correlation p-value.
## n <- 25
## x <- rnorm(n)
## y <- -x + rnorm(n)
## R <- cor(x,y)
## p.value <- 2*pt(-abs(R * sqrt(n-2) / sqrt(1-R*R)),df=n-2)
## cor.test(x,y,)$p.value

## Get table with genotype files.
genotype.files <- Sys.glob(paste("data/geno_data_cache_", reg.name, "/chr*.RData", sep=""))
pull.chr <- function(x) {
  res <- gsub(".+chr", "chr", x,  perl=T)
  res <- gsub("\\..+", "", res, perl=T)
  res
}
genotype.files <- data.frame(chr=unlist(lapply(genotype.files, pull.chr)), genotype.file=genotype.files, stringsAsFactors=F)


#error in test.matrix[not.na.idx,] subscript out of bounds

##for(numPCs in 0:30) {
##for(numPCs in seq(0,30,2)) {
##for(numPCs in seq(1,30,2)) {
##for(numPCs in seq(0,25,1)) {

# for(numPCs in seq(0,20,1)) {
#a foreach solution?
cl <- makeCluster(2)#4 cores works well
registerDoParallel(cl)
#for protein groups with majority ids consisting wholly of isoforms (same gene) that ENSGID is returned. IF not, then 'NA' is returned.
Blank <- foreach(numPCs = 0:1) %dopar% {
## load protien data
  pheno.data <- load.pro(numPCs)#PC corrected/regressed

  process.chr <- function(g) {        
    ## Load genotype data for given chromosome.
    chr <- genotype.files$chr[g]
    cache.file <- paste("data/geno_data_cache_", reg.name, "/", chr, ".RData", sep="")
    print(cache.file)
    
    pheno.cols <- grep("GM", names(pheno.data), fixed=T)
    
    if(file.exists(cache.file)) {
      load(cache.file)
    }else {
      stop("no cache file")
    }
    print(chr)
    geno.cols <- grep("GM", names(genotype.data), fixed=T)
    
    ## Convert data to matricies for fast access.
    geno.matrix <- as.matrix(genotype.data[,geno.cols])
    row.names(geno.matrix) <- genotype.data$snp.name
    snp.hg18.pos <- genotype.data$hg18.pos
    
    ## Get protein measurements for chromosome.
    pheno.data.chr <- pheno.data[pheno.data$chr == chr,]
    pheno.data.chr <- pheno.data.chr[sort(pheno.data.chr$ENSG, index.return=T)$ix,]  ## Sort so matched to mRNA
    pheno.matrix.data <- as.matrix(pheno.data.chr[,pheno.cols]) #60 entries
    
    run.permutation <- function(i) {
      pheno.cur <- pheno.matrix.data[i,]
      n.dim <- length(pheno.cur)
      n <- sum(is.na(pheno.cur) == F); nm2 <- n-2; sqrt.nm2 <- sqrt(nm2)
      n.iter <- 10000
      null.P.data <- numeric(n.iter)
      ## Get genotypes in a given window.
      txStart <- pheno.data.chr$txStart[i]; txEnd <- pheno.data.chr$txEnd[i]  
      snp.rows <- which(txStart - 20000 <= snp.hg18.pos & snp.hg18.pos <= txEnd + 20000)
      ##snp.rows <- which(txStart - 100000 <= snp.hg18.pos & snp.hg18.pos <= txStart + 100000)
      if(length(snp.rows) > 0) {
        if(length(snp.rows) == 1) {
          test.matrix <- as.matrix(geno.matrix[snp.rows,]) ## transpose for faster correlation computation
        }else {
          test.matrix <- t(geno.matrix[snp.rows,]) ## transpose for faster correlation computation
        }
        ## Get corresponding SNPs.
        test.snps <- genotype.data$snp.name[snp.rows]
        
        ## NA genotypes for NA phenotypes.
        test.matrix[which(is.na(pheno.cur)),] <- NA
        
        ## Recalculate MAF for each SNP. 
        compute.maf <- function(r) {
          genotypes <- test.matrix[,r]
          ## Genotype codes.
          ## 0-0.5] = BB allele
          ## (0.5, 1.5] = heterozygous, AB alleles
          ## (1.5, 2.5] = AA allele
          major.cnt <- sum(0 <= genotypes & genotypes <= 0.5, na.rm=T)
          het.cnt <- sum(0.5 < genotypes & genotypes <= 1.5, na.rm=T)
          minor.cnt <- sum(1.5 < genotypes & genotypes <= 2.5, na.rm=T)
          ##N <- dim(geno.matrix)[2]
          N <- sum(is.na(genotypes) == F)
          major.freq <- (2*major.cnt + het.cnt) / (2*N)
          minor.freq <- (2*minor.cnt + het.cnt) / (2*N)
          min(major.freq, minor.freq)
        }
        maf <- unlist(lapply(1:dim(test.matrix)[2], compute.maf))
        maf.cut <- 0.1
        maf.OK <- which(maf > maf.cut)
        test.matrix <- as.matrix(test.matrix[, maf.OK])
        test.snps <- test.snps[maf.OK]
        
        if(length(test.snps) > 0) { 
          ## Get protein p-values.
          not.na.idx <- which(is.na(pheno.cur) == F)
          test.matrix <- test.matrix[not.na.idx,]
          pheno.cur <- pheno.cur[not.na.idx]
          
          R.data <- cor(pheno.cur, test.matrix)
          
          absT.data <- abs(R.data * sqrt.nm2 / sqrt(1-R.data*R.data))
          p.values.data <- 2*pt(-absT.data,df=nm2)
          maxT.data.idx <- which.max(absT.data)
          maxT.data <- absT.data[maxT.data.idx]
          snpname <- test.snps[maxT.data.idx]
          snp.pvalue <- p.values.data[maxT.data.idx]
          snp.R.value <- R.data[maxT.data.idx]
          
          for(iter in 1:n.iter) {    ## Use pearson correlation p-values.
            R.data <- cor(pheno.cur[sample(n)], test.matrix)
            null.P.data[iter] <-  max(abs(R.data * sqrt.nm2 / sqrt(1-R.data*R.data)), na.rm=T)
          }
          adj.gene.pvalue <- mean(null.P.data > maxT.data, na.rm=T)
          list(p.value=adj.gene.pvalue,snp.name=snpname, snp.pvalue=snp.pvalue, snp.R.value=snp.R.value)
        }
        else {
          list(p.value=NA, snp.name=NA, snp.pvalue=NA, snp.R.value=NA)
        }
      }
      else {
        list(p.value=NA, snp.name=NA, snp.pvalue=NA, snp.R.value=NA)
      }
    }
    perm.results <- lapply(1:dim(pheno.data.chr)[1], run.permutation)
    #perm.results <- mclapply(1:dim(pheno.data.chr)[1], run.permutation, mc.cores=6, mc.preschedule=F)
    ##perm.results <- mclapply(1:dim(pheno.data.chr)[1], run.permutation, mc.cores=6, mc.preschedule=T)
    
    p.values <- unlist(lapply(perm.results, function(x) { x$p.value } ))
    snp.name <- unlist(lapply(perm.results, function(x) { x$snp.name } ))
    snp.pvalue <- unlist(lapply(perm.results, function(x) { x$snp.pvalue } ))
    snp.R.value <- unlist(lapply(perm.results, function(x) { x$snp.R.value } ))
    
    res <- data.frame(ENSG=pheno.data.chr$ENSG,chr=chr,
                      p.values=p.values,                      
                      snp.name=snp.name,
                      snp.pvalue=snp.pvalue,
                      snp.R.value=snp.R.value,
                      stringsAsFactors=F)
    return(res)
  }
  
  #     df.list <- mclapply(1:dim(genotype.files)[1], process.chr, mc.cores=12, mc.preschedule=F)
  df.list <- lapply(1:dim(genotype.files)[1], process.chr)
  
  combined <- df.list[[1]]
  if(length(df.list) > 1) {
    for(i in 2:length(df.list)) {
      combined <- rbind(combined, df.list[[i]])
    }
  }
  
  file.name <- paste("combined_", reg.name, "_numPCs_", numPCs ,".RData", sep="")
  save(combined, file=file.name)
}
stopCluster(cl)

#read in all PC regressed SNP pQTL association files and compute number of pQTLs (at given pvalue threshold).
##threshold is empirical (permutation derived) gene level p value based on 10K permutations at .01
pQTL.files <- Sys.glob("./combined_pro_numPCs*.RData")
pQTL <- c()
for (file in pQTL.files){
  load(file)
  #sort and truncate combined by "p.values" to 1% level.
  combined <- combined[order(combined$p.values),]
  combined <- combined[combined$p.values <= .50,]
  pQTL <- c(pQTL, nrow(combined))
}
#plot results.
plot(0:20, pQTL, xlab = "# PCs regressed from protein data", ylab = "pQTLs at .50 threshold")


#trying using the snp.pvalue column with p.adjust function
pQTL.files <- Sys.glob("./combined_pro_numPCs*.RData")
pQTL <- c()
for (file in pQTL.files){
  load(file)
  #sort and truncate combined by "p.values" to 1% level.
  combined <- combined[order(combined$snp.pvalue),]
  combined$padjust <- p.adjust(combined$snp.pvalue, method = "BH")
  combined <- combined[combined$padjust <= .1,]
  pQTL <- c(pQTL, nrow(combined))
}
#plot results.
plot(0:20, pQTL, xlab = "# PCs regressed from protein data", ylab = "pQTLs at .1 'snp.pvalue' + BH threshold")









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
  
  