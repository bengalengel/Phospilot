##run the regressions on Zia's data...
rm(list=ls(all=TRUE)) #start with empty workspace

#required packages and libraries
require(plyr)
require(limma)
require(biomaRt)
require(iterators)
require(foreach)
require(doParallel)
require(parallel)
library(compiler); setCompilerOptions(optimize=3); enableJIT(3)


#change the working directory
setwd("D:/for_brett/")


#estimate and regress away PCs from protein matrix
#make the PC regression and imputation function


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
  combined <- combined[combined$p.values <= .05,]
  pQTL <- c(pQTL, nrow(combined))
}
#plot results.
plot(c(0,1,8), pQTL, xlab = "# PCs regressed from protein data", ylab = "pQTLs at .05 threshold")


