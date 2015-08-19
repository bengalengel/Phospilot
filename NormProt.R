NormProt <- function(directory){
  #this function loads MQoutput from SILAC quantitative proteomic analysis of 60 LCLsm reduces it to the three samples of interest, median and quantile normalizes it. Returns a list of three data frames. MQoutput, median normalized, and quantile normalized. 
  source("loadMQZ.R")
  require(plyr)
  require(limma)
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
  
  
  # PC regression based pQTL optimization ---------------------
  #note that this requires "genotypes_imputed" files and ensembl_4_30.gz from Jack. See readme file in genotypes imputed folder. 
  #I will source alterations of Zia's pQTL mapping scripts to build the genotype files, regress PCs, and calculate pQTLs. 
  
  #check for output from PC regression regression in local directory. run required scripts or prompt if files are not present
  
  #1 - check for protein table with single gene assignment and chromosome and transcription start/stop. Build if necessary
  #note uses ensembl_4_30_gz file from eQTL website that Jack produced
  if(!file.exists("./PCregression/PRO_raw.RData")) {
    source("./PCregression/BuildProteinTable.R")
    BuildProt(protein1)
  }
  
  #2 - check for genotype data and cahe files. Build if necessary. Note using Jack's imputed data. Note it takes awhile to make cache files. Here I am using MAF of .10
  if(length(Sys.glob("./PCregression/data/geno_data/*.RData")) == 0) {
    source("./PCregression/BuildGenoTable.R")
  }
  if(length(Sys.glob("./PCregression/data/geno_data_cache_pro/*.RData")) == 0) {
    source("./PCregression/BuildCache.R")
  }
  
  #check for regression output. Stop function if not present and prompt user to run pc regression on cluster.
  if(length(Sys.glob("./PCregression/Output/*.RData")) == 0){
    stop("Run PC_regress on cluster. See the readme file")
  }
 
  #Find max pQTL count  
  
  #Output from batcharray job should be placed in output directory. Run to 30 PCs. Ensure all 31 output files
#   present.
  
  #read in all PC regressed SNP pQTL association files and compute number of pQTLs (at given pvalue threshold).
  ##threshold is empirical (permutation derived) gene level p value based on 10K permutations
  pQTL.files <- Sys.glob("./PCregression/Output/combined_pro_numPCs*.RData")
  pQTL.files <- pQTL.files[order(nchar(pQTL.files),pQTL.files)]
  pQTL <- c()
  for (file in pQTL.files){
    require(qvalue)
    load(file)
    #sort and truncate combined by "p.values" to 1% level.
    #   combined <- combined[order(combined$p.values),]
    #   combined <- combined[combined$p.values <= .50,]
    #   not.NA <- which(is.na(combined$p.values) == F)
    #   qobj <- qvalue(combined$p.values[not.NA])
    #   combined$q.values <- NA
    #   combined$q.values[not.NA] <- qobj$qvalues
    qobj <- qvalue(combined$p.values)
    combined$q.values <- qobj$qvalues
    combined <- combined[combined$q.values <= .1,]
    print(nrow(combined))
    pQTL <- c(pQTL, nrow(combined))
  }
  #plot results.
  plot(0:30, pQTL, xlab = "# PCs regressed from protein data", ylab = "pQTLs at .10 threshold")
  #which is the max?
  numPCs <- which.max(pQTL) - 1
  
  
  # define load.pro function
  load.pro <- function(numPCs, keep.NAs=TRUE, standardize=TRUE) {
    require(impute)
    load("./PCregression/PRO_raw.RData")
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
  
  # Return and save the matrix with the maximum number of pQTL calls
  GelProtPCRegressed <- load.pro(numPCs)#PC corrected/regressed
  row.names(GelProtPCRegressed) <- GelProtPCRegressed$ENSG
  
  file.name <- paste("GelProt_", numPCs ,"_PCsRegressed", ".RData", sep="")
  save(GelProtPCRegressed, file = file.name)
  
  #subsetted to present in all samples of interest
  #subset data to only the samples of interest
  RegressedCommon <- GelProtPCRegressed[,c("GM18862","GM18486","GM19160")]
  RegressedCommon <- na.omit(RegressedCommon)#3896! Only lose 374 proteins compared to 'data' below!
  file.name <- paste("GelProt_", numPCs ,"_PCsRegressedCommon", ".RData", sep="")
  save(RegressedCommon, file = file.name)
  
  ####now for the nonregressed dataframes ------
  
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
  DFs <- list(protein1, data, MedianNorm, quantiled, GelProtPCRegressed, RegressedCommon)
  return(DFs)
}
  
  