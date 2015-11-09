BatchNorm <- function(multExpanded1){
  #this function accepts the expanded phosphosite file and performs normalization, batch correction
  # EDA and returns the 'pilot' dataframe. The pilot dataframe is the per biological replicate average
  # across all phosphopeptides. See 'combat tests.R' for a more complete version of this.
  require(sva)
  require(plyr)
  require(swamp)
  require(limma)
  library(gplots)
  require(RColorBrewer)
  
  #make the sparse matrix from real data (all quants of class 1 here) 
  expCol <- grep("HL(.*)", colnames(multExpanded1))
  data <- multExpanded1[,expCol]
  
  # add row names with site id and multiplicity designation 
  # row.names(data) <- multExpanded$id#note this won't work because of the multiplicity issue
  idmult <- paste(multExpanded1$id, multExpanded1$multiplicity, sep="_")
  row.names(data) <- idmult
  multExpanded1 <- cbind(multExpanded1,idmult)
  data <- log2(data)
  
  #summaries
  # summary(data)
  # boxplot(data[,1:12])#extreme outlier in HL18486_1_2
  # boxplot(data[,1:12], ylim=c(-4,4))
  
  #remove outlier which is 3X greater than any other datapoint. Not present in ensembl data.
#   data[,2][which.max(data[,2])] <- NaN
  boxplot(data[,1:12])#extreme outlier in HL18486_1_2 removed
  
  RawRatios <- data
  
  #median normalize
  names <- colnames(data)[1:12]
  median.subtract <- function(x){ x - median(x, na.rm = TRUE)}##create a wrapper for median subtraction
  data <- colwise(median.subtract, names)(data)
  row.names(data) <- idmult##add back the row names
  
  MedianNorm <- data
 
  
  # quantile normalization. from normalize.quantiles {preprocessCore}  
  # "This functions will handle missing data (ie NA values), based on the assumption that the data is missing at random."
  
  quantiled <- normalizeQuantiles(data,ties = T)#ties are all assigned the same value for the common quantile
  # summary(quantiled)
  # boxplot(data)
  boxplot(quantiled)
  # density plots all look the same
  plot.new()
  par(mfrow = c(1, 1))
  for (i in 1:(ncol(quantiled))){
    if(i==1) plot(density(quantiled[, i], na.rm=T), col = i, ylim = c(0,.55))
    else lines(density(quantiled[, i], na.rm=T), col = i)
  }
  
  # remove exp obs if not observed in each sample quantiled
  quantiled2 <- quantiled[rowSums(is.na(quantiled[ , 1:4])) < 4 & rowSums(is.na(quantiled[ , 5:8])) < 4 & rowSums(is.na(quantiled[ , 9:12])) < 4,]    
  
  # remove exp obs if not observed in each batch
  quantiled3 <- quantiled[rowSums(is.na(quantiled[ , c("HL18486_1_1", "HL18486_1_2", "HL18862_1_1", "HL18862_1_2", "HL19160_1_1", "HL19160_1_2")])) < 6 
                          & rowSums(is.na(quantiled[, c("HL18486_2_1", "HL18486_2_2", "HL18862_2_1", "HL18862_2_2", "HL19160_2_1", "HL19160_2_2")])) < 6,]    
  
  # remove exp obs if not observed two or more times in each batch to ensure a variance measurement
  quantiled4 <- quantiled[rowSums(is.na(quantiled[ , c("HL18486_1_1", "HL18486_1_2", "HL18862_1_1", "HL18862_1_2", "HL19160_1_1", "HL19160_1_2")])) < 5 
                          & rowSums(is.na(quantiled[, c("HL18486_2_1", "HL18486_2_2", "HL18862_2_1", "HL18862_2_2", "HL19160_2_1", "HL19160_2_2")])) < 5,] 
  quantiled5 <- na.omit(quantiled)##common across all
  
  
  #batch effect identification and adjustment using swamp/combat*************************
  swamp <- as.matrix(quantiled5)
  swamp <- swamp[,1:12]
  ##### sample annotations (data.frame)
  set.seed(50)
  o1<-data.frame(Batch = factor(rep(c("A","A","B","B"),3)),
                 Random = rnorm(12), row.names=colnames(swamp))
  
  # PCA analysis
  res1<-prince(swamp,o1,top=10,permute=T)
  str(res1)
  a <- res1$linp#plot p values
  b <- res1$linpperm#plot p values for permuted data
  pdf("PC_batch_correlation.pdf", 11.5, 8)
  prince.plot(prince=res1, key = T, note = T, breaks = 16, 
              col =  rev(colorRampPalette(c("white", "light gray", "dark gray", "red"))(15)),
#               col =  heat.colors(15)
              )
  dev.off()
  #There is a batch effect associated with the process date.
  # I must combat this
  ##batch adjustment using quantiled4
  swamp <- as.matrix(quantiled4)
  ##### sample annotations (data.frame)
  set.seed(50)
  o1<-data.frame(Batch=factor(rep(c("A","A","B","B"),3)),
                 Random=rnorm(12),row.names=colnames(swamp))
  
  
  com2<-combat(swamp,o1$Batch,batchcolumn=1) #WORKS AFTER ENSURING AT LEAST TWO IN A BATCH.
  # Found 2 batches
  # Found 0 covariate(s)
  # Found 25032 Missing Data Values
  # Standardizing Data across genes
  # Fitting L/S model and finding priors
  # Finding parametric adjustments
  # Adjusting the Data
  
  ##batch effect correction using sva combat and 'covariate' matrix
  
  # now for the full dataset n=8560
  cdata <- na.omit(com2)#also used below
  prince.plot(prince(cdata,o1,top=10)) #huzzah!
  
  # PCA analysis
  res1<-prince(cdata,o1,top=10,permute=T)
  #str(res1)
  c <- res1$linp#plot p values
  d <- res1$linpperm#plot p values for permuted data
  out <- rbind(a,b,c,d)
  write.table(out, "PC_ba_batch.csv", sep=',', col.names=T, row.names=F) #PCs before and after batch correction.
  
  
  ##batch corrected EDA********************************************************************************
  par(mfrow = c(1, 1))
  boxplot(com2, cex.axis = 1, cex.names = .5, cex.lab = .5, las=2)#fix the margins later
  summary(com2)
  # density plots
  plot.new()
  for (i in 1:(ncol(com2))){
    if(i==1) plot(density(com2[, i], na.rm=T), col = i, ylim = c(0,.75))
    else lines(density(com2[, i], na.rm=T), col = i)
  }
  # I am not going to normalize again after batch correction.
  
  
  
  # now with missing data removed perform the clustering and heatmaps*******************************************
  dataZ <- scale(cdata)##Z-scored column wise
  
  # now all data excepting complete cases (note that the sample dendograms look the same)
  #hist(dataZ[,6], breaks = 100)
  
  # dendogram using euclidian distance (default) and ward or complete agglomeration
  dend.ward<- as.dendrogram(hclust(dist(t(dataZ)),method="ward"))
  dend.complete<- as.dendrogram(hclust(dist(t(dataZ))))
  
  ward.o<- order.dendrogram(dend.ward)
  complete.o<- order.dendrogram(dend.complete)
  
  plot(dend.complete,ylab="height", main = "Euclidian/Complete")
  plot(dend.ward, leaflab = "perpendicular", ylab = "height", main = "Euclidian/Ward")
  
  plot.new()##produces a blank canvas
  
  # Cluster using euclidian distance and ward linkage for both sites(rows) and samples (columns)
  # Note that both dendograms are created independently and row Z scores are presented in the heatmap
  
  # row scaled
  r <- t(scale(t(cdata)))#transpose to zscale the rows then transpose back to original format
  
  # sample scaled
  c <- scale(cdata)
 
  # Create dendrogram using the data without NAs
  feature.dend<- as.dendrogram(hclust(dist(r),method="ward"))
  sample.dend<- as.dendrogram(hclust(dist(t(c)),method="ward"))##note that dist caclculates distance between rows by default
  
  
  ##produce the heatmap. Note that the help page has a nice section on identifying subregions by color. Although I will likely have to cut the dendogram to id clusters of interest
  
  heatmap.2(
    r,#row Z scores
    Colv=sample.dend,
    Rowv=feature.dend,
    col=bluered(25),
    scale="none",
    trace="none",
    density.info="none",
    key.xlab = "Row Z scores", key.ylab=NULL, key.title = "",
    srtCol=45,  ,adjCol = c(1,1),
    margins = c(6,5),
    cexCol=1,
    labRow = NA#remove row labels
  )
  
  
  #PCA analysis  - cdata = combat corrected confounded phospho data. 
  x <- t(cdata)#samples are the rows of the column matrix
  pc <- prcomp(x, scale = T, center = T) #I am now scaling and centering x

  names(pc)

  cols <- as.factor(substr(colnames(cdata), 3, 7))##check me out. use 5 digit exp name.
  pdf("PCA_BEcorrect_confounded.pdf")
  plot(pc$x[, 1], pc$x[, 2], col = as.numeric(cols), main = "Batch Effect Corrected PCA", xlab = "PC1", ylab = "PC2", pch = 1, 
     cex = 1.5, family = "serif")
  legend("bottomleft", levels(cols), col = seq(along=levels(cols)), pch = 1, cex = 1.25)
  dev.off()

  summary(pc)
  
  #SVD for calculating variance explained;
  cx <- sweep(x, 2, colMeans(x), "-")
  sv <- svd(cx)
  names(sv)
  plot(sv$u[, 1], sv$u[, 2], col = as.numeric(cols), main = "SVD", xlab = "U1", ylab = "U2")  
  plot(sv$d^2/sum(sv$d^2), xlim = c(1, 12), type = "b", pch = 16, xlab = "principal components", 
       ylab = "variance explained")
  
  
  
  
  # AND NOW FOR ALL THE DATA around 5K*****************************************
  adata <- com2[rowSums(is.na(com2[ , 1:2])) < 2 & rowSums(is.na(com2[ , 3:4])) < 2 & rowSums(is.na(com2[ , 5:6])) < 2 
                & rowSums(is.na(com2[ , 7:8])) < 2 & rowSums(is.na(com2[ , 9:10])) < 2 & rowSums(is.na(com2[ , 11:12])) < 2,] #there needs to be at least one measurement per biological replicate                   
  
  # Produce dataframe from sample means ignoring missing data
  
  HL18486_1 <- rowMeans(adata[,1:2], na.rm = T)
  HL18486_2 <- rowMeans(adata[,3:4], na.rm = T)
  HL18862_1 <- rowMeans(adata[,5:6], na.rm = T)
  HL18862_2 <- rowMeans(adata[,7:8], na.rm = T)
  HL19160_1 <- rowMeans(adata[,9:10], na.rm = T)
  HL19160_2 <- rowMeans(adata[,11:12], na.rm = T)
  
  
  pilot <- cbind(HL18486_1, HL18486_2, HL18862_1, HL18862_2, HL19160_1, HL19160_2)#4,996 class 1 measurements with at least one quant in each biological replicate
  
  boxplot(pilot)
  
  #write the normalized and batch corrected "pilot" matrix
  write.csv(pilot,"pilot_dataframe.csv", row.names = T)
  
  #data frames to be returned; quantiled1-5, com2, adata, pilot
  DFs <- list(RawRatios, MedianNorm, quantiled, quantiled2, quantiled3, quantiled4, quantiled5, com2, adata, pilot)
  saveRDS(DFs, file = "./CorrectedData.rds")
  return(DFs)
}