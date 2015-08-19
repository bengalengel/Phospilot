#alter genotype cache st MAF > .10 and samples the same as those identified via proteomics

require(parallel); require(plyr)
library(compiler)
enableJIT(3)

#load the protein file
load("./PCregression/PRO_raw.RData")
pheno.data <- pheno.pro; suffix <- "pro"
## Get table with genotype files.
genotype.files <- Sys.glob("./PCregression/data/geno_data/chr*.RData")
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
  cache.file <- paste("./PCregression/data/geno_data_cache_", suffix, "/", chr, ".RData", sep="")
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
  if(!file.exists("./PCregression/data/geno_data_cache_pro")){
    dir.create(file.path(getwd(), "./PCregression/data/geno_data_cache_pro"))
  }  
  save(genotype.data, file=cache.file)
  NA
}
unlist(lapply(1:dim(genotype.files)[1], build.geno.cache))
