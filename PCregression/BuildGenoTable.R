
## Reads in genotype data. Renames the columns and saves as .RData files.
require(parallel)
library(compiler); setCompilerOptions(optimize=3); enableJIT(3)


# A foreach analog may work here
# cl <- makeCluster(2)#This will depend on available memory (these files are rather large...). I can do this on laptop perhaps but not at home.
# registerDoParallel(cl)
# blank <- foreach(i = 1:dim(genotype.files)[1], .combine = 'c') %dopar% load.genotypes(i)
# stopCluster(cl)

#Files needed: 
mean.genotype.files <- Sys.glob("./PCregression/data/genotypes_imputed/*.mean.genotype.txt")
snpdata.files <- Sys.glob("./PCregression/data/genotypes_imputed/*.snpdata.txt")

pull.chr <- function(x) {
  res <- gsub(".+chr", "chr", x,  perl=T)
  res <- gsub("\\..+", "", res, perl=T)
  res
}
snpdata.files <- data.frame(chr=unlist(lapply(snpdata.files, pull.chr)), snpdata.file=snpdata.files, stringsAsFactors=F)
mean.genotype.files <- data.frame(chr=unlist(lapply(mean.genotype.files, pull.chr)), mean.genotype.file=mean.genotype.files, stringsAsFactors=F)

genotype.files <- merge(snpdata.files, mean.genotype.files)

## Load individuals by column
con  <- file("./PCregression/data/genotypes_imputed/individual.order.txt", open = "r")

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
  if(!file.exists("./PCregression/data")){
    dir.create(file.path(getwd(), "./PCregression/data"))
  }  
  if(!file.exists("./PCregression/data/geno_data")){
    dir.create(file.path(getwd(), "PCregression/data/geno_data"))
  }
  save(genotype.data, file=paste("./PCregression/data/geno_data/",chr, ".RData", sep=""), compress="gzip", compression_level=1)
  NA
}

blank <- unlist(lapply(1:dim(genotype.files)[1], load.genotypes))