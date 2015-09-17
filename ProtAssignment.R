ProtAssignment <- function(protein, proteome, multExpanded1){
  ## This function performs a loop to assign protein ids, H/L values and ibaq values to phosphosites for normalization using "protein groups" file produced from the phospho workup. The 'proteome' FASTA file is used to assign modification positions to the ratios used for quantification. Intellegent naming of the columns is used here as well. 
  
  #inputs needed are ME/ME1 and 'protein' from SCX-IMAC protein groups file WITHOUT SUBSETTING BY ONLY IDENTIFIED BY SITE!! (that is not protein1)
  
  ## I require protein copynumber and/or iBAQ estimates to test for enrichment of expression level in differential phosphorylation subsets or NRE model quandrants
  
  #When there are multiple protein groups that could match the phosphosite I choose the one with the most unique/razor ids
  require(qdapRegex)
  require(sva)
  require(plyr)
  require(swamp)
  require(limma)
  library(gplots)
  
  #########################
  #collect ids from multexpanded1 DF. Results will be appended to this dataframe with appropriate headers. 
  ids <- multExpanded1$Protein.group.IDs
  ids <- as.character(ids)
  
  #for those ids that contain multiple groups. pick the one with the most unique plus razor 
  index <- grep(";", ids)#location within vector of multiple group hits
  MultiMatch <- ids[index]
  
  #for each of the phosphopeptide ids with multiple protein ids, choose as a match the proteingroup id with the most razor + unique peptides
  proteinids <- protein$id#for the loops
  FinalIDMultMatch <- sapply(MultiMatch,function(x){
    y <- as.numeric(unlist(strsplit(x,split = ";")))
    #below sapply works! but ddply doesn't?
    hits <- sapply(proteinids, function(z) any(y %in% z))
    matches <- protein[hits, c("id","Razor...unique.peptides")] 
    #assign id with greatest amount of razor + unique peptides
    finalid <- matches$id[which.max(matches$Razor...unique.peptides)]
  })
  
  #assign these final ids to the main id vector
  ids[index] <- as.numeric(FinalIDMultMatch) #changes 'ids' to class list
  # ids2 <- as.numeric(as.character(ids2))#but why?...length 1 loss with unlist
  #one na introduced by coercion
  
  #add this list of ids to me1df
  multExpanded1$PhosPrepMatchProteinGroupID <- ids
  
  #locate the positions of the phosphosites within the matched protein groups. Because only one protein group is chosen, the original position assignments may include non-matched group assignments and therefore should be corrected.
  PositionInProteinsFinal <- c()
  for(i in seq_along(ids)){
    #     for(i in 1:50){
    ##retrieve the peptide
    #get the peptide to search
    peptide <- as.character(multExpanded1$Phospho..STY..Probabilities[i])
    peptide <- gsub(pattern = " *\\(.*?\\) *", replacement = "", peptide)
    
    ##return the ENSPids for each group assigned to a peptide.
    matches <- as.character(protein[protein$id == ids[i], "Protein.IDs"])
    matches <- unlist(strsplit(as.character(matches), ";"))

    #Retrieve position(s) of phosphosite within protein match(s). Where a protein group member doesn't match an 'NA' is returned. Oftentimes these are contaminant proteins.
    PositionInProteins <- c()
    for(j in seq_along(matches)){
      index <- grep(matches[j], names(proteome))
      seq <- unlist(seqinr::getSequence(object = proteome[index], as.string = T))
      #position of peptide within protein
      protpos <- regexpr(peptide,seq)
      #now I need the position of the modified site within the peptide
      windows <- as.character(multExpanded1$Sequence.window[i])
      windows <- unlist(strsplit(as.character(windows), ";"))
      window <- windows[1]
      sitepos <- regexpr(peptide,window)
      modsiteinpeptide <- 16 - sitepos[1]#modified site is always the 16th position within the 'sequence window'.
      ##add this number to protein position!!!! Hurray
      tmp <- protpos[1] + modsiteinpeptide    
      PositionInProteins <- c(PositionInProteins, tmp)
    }
    
    #add the position of the modified site information based on matching protein ids
    PositionInProteins <- paste(PositionInProteins, collapse = ";")
    PositionInProteinsFinal <- c(PositionInProteinsFinal, PositionInProteins)
  }
    
   #Check the length of the PositionInProteinFinal vector
    length(PositionInProteinsFinal)
    
  #add ids and phosphosite positions to multexpanded df
  multExpanded1$PhosPrepPositionInProteins <- PositionInProteinsFinal
  
  #retrieve the H/L, ibaq and majority protein id information from protein1
  ibaqNames <- grep("iBAQ.H.*",names(protein), value = T)
  ratios <- grep("Ratio*",names(protein), value = T)
  info <- protein[,c("id", "Protein.IDs", "Majority.protein.IDs", "Only.identified.by.site", ratios, ibaqNames)]
  names(info)[1] <- "PhosPrepProteinGroupID"
  names(info)[2] <- "PhosPrepProteinID"
  names(info)[3] <- "PhosPrepMajorityProteinIDs"
  names(info)[4] <- "PhosPrepProteinOnlyIDdbySite"
  
  #make protein names informative
  ibaqindex <- grep("iBAQ.H.*",names(info))
  ratioindex <- grep("Ratio*",names(info))
  names(info)[ratioindex] <- gsub("Ratio.H.L.normalized.", "PhosPrepProteinHL", names(info)[ratioindex])
  names(info)[ibaqindex] <- paste0("PhosPrepProtein", names(info)[ibaqindex])
  
  #merge by 'matchids' from ME1 and protein group 'id' from protein groups file. 
  
  # merge into ME DF
  multExpanded1 <- merge(multExpanded1, info, by.x = "PhosPrepMatchProteinGroupID", by.y = "PhosPrepProteinGroupID")
  
  #remove the 'NAs' from the 'PositionsInProteins' and the 'ProteinID' columns. The ProteinID columns may contain reverse or contaminant entries. Note that there are no mapping issues in 'Majority protein" columns because these identifications are from the same raw data. Note that this mapply call was the neatest solution. sapply follwed by unlist returned one huge logical vector. It would be nice to know how to unlist each vector element of a list! I could paste/collapse each element first to make it a character vector of length one I suppose, but the mapply call did this in a neater fashion.
  TruncateCommaDelim <- function(Positions,TruncationTarget){
    index <- unlist(strsplit(Positions,";")) != "NA"
    Target <- unlist(strsplit(TruncationTarget,";"))
    Truncated <- Target[index]
    paste(Truncated, collapse = ";")
  }
  
  multExpanded1$PhosPrepProteinIDTruncated <- mapply(TruncateCommaDelim, multExpanded1$PhosPrepPositionInProteins, multExpanded1$PhosPrepProteinID)
  multExpanded1$PhosPrepPositionInProteinsTruncated <- mapply(TruncateCommaDelim, multExpanded1$PhosPrepPositionInProteins, multExpanded1$PhosPrepPositionInProteins)
  
  # Extraction and normalization of PhosPrep protein ratios -----------------
  require(limma)
  #Retrieve and transform data
  expCol <- grep("PhosPrepProteinHL(.*)", colnames(multExpanded1))
  RawRatios <- multExpanded1[,expCol]
  
  # add row names with site id and multiplicity designation 
  row.names(RawRatios) <- multExpanded1$idmult
  RawRatios <- log2(RawRatios)
  
  #median normalize
  names <- colnames(RawRatios)
  median.subtract <- function(x){ x - median(x, na.rm = TRUE)}##create a wrapper for median subtraction
  MedianNorm <- colwise(median.subtract, names)(RawRatios)
  row.names(MedianNorm) <- multExpanded1$idmult##add back the row names
  
  # quantile normalization. from normalize.quantiles {preprocessCore}  
  # "This functions will handle missing data (ie NA values), based on the assumption that the data is missing at random."
  
  quantiled <- normalizeQuantiles(MedianNorm,ties = T)#ties are all assigned the same value for the common quantile
  # summary(quantiled)
  # boxplot(data)
  boxplot(quantiled)
  # density plots all look the same of course
  par(mfrow = c(1, 1))
  for (i in 1:(ncol(quantiled))){
    if(i==1) plot(density(quantiled[, i], na.rm=T), col = i, ylim = c(0,.85))
    else lines(density(quantiled[, i], na.rm=T), col = i)
  }
  #subsets of quantiled
  #clean up names from quantiled
  names(quantiled) <- gsub("PhosPrepProtein", replacement = "", names(quantiled))
  # remove exp obs if not observed in each sample quantiled
  quantiled2 <- quantiled[rowSums(is.na(quantiled[ , 1:4])) < 4 & rowSums(is.na(quantiled[ , 5:8])) < 4 & rowSums(is.na(quantiled[ , 9:12])) < 4,]    
  
  # remove exp obs if not observed in each batch
  quantiled3 <- quantiled[rowSums(is.na(quantiled[ , c("HL18486_1_1", "HL18486_1_2", "HL18862_1_1", "HL18862_1_2", "HL19160_1_1", "HL19160_1_2")])) < 6 
                          & rowSums(is.na(quantiled[, c("HL18486_2_1", "HL18486_2_2", "HL18862_2_1", "HL18862_2_2", "HL19160_2_1", "HL19160_2_2")])) < 6,]    
  
  # remove exp obs if not observed two or more times in each batch to ensure a variance measurement
  quantiled4 <- quantiled[rowSums(is.na(quantiled[ , c("HL18486_1_1", "HL18486_1_2", "HL18862_1_1", "HL18862_1_2", "HL19160_1_1", "HL19160_1_2")])) < 5 
                          & rowSums(is.na(quantiled[, c("HL18486_2_1", "HL18486_2_2", "HL18862_2_1", "HL18862_2_2", "HL19160_2_1", "HL19160_2_2")])) < 5,] 
  quantiled5 <- na.omit(quantiled)##common across all
  
  #once in each biological replicate
  quantiledBio <- quantiled[rowSums(is.na(quantiled[ , 1:2])) < 2 & rowSums(is.na(quantiled[ , 3:4])) < 2 & rowSums(is.na(quantiled[ , 5:6])) < 2 
                & rowSums(is.na(quantiled[ , 7:8])) < 2 & rowSums(is.na(quantiled[ , 9:10])) < 2 & rowSums(is.na(quantiled[ , 11:12])) < 2,] #there needs to be at
  
  #combat batch correction and EDA validation
  
  #batch effect identification and adjustment using swamp/combat*************************
  swamp <- as.matrix(quantiled5)
  swamp <- swamp[,1:12]
  ##### sample annotations (data.frame)
  set.seed(50)
  o1<-data.frame(Factor1=factor(rep(c("A","A","B","B"),3)),
                 Numeric1=rnorm(12),row.names=colnames(swamp))
  
  # PCA analysis
  res1<-prince(swamp,o1,top=10,permute=T)
  str(res1)
  a <- res1$linp#plot p values
  b <- res1$linpperm#plot p values for permuted data
  prince.plot(prince=res1)
  
  #There is a batch effect associated with the process date.
  # I must combat this
  ##batch adjustment using quantiled4
  swamp <- as.matrix(quantiled4)
  ##### sample annotations (data.frame)
  set.seed(50)
  o1<-data.frame(Factor1=factor(rep(c("A","A","B","B"),3)),
                 Numeric1=rnorm(12),row.names=colnames(swamp))
  
  
  com2<-combat(swamp,o1$Factor1,batchcolumn=1) #plots could look better...
#   Found 2 batches
#   Found 0 covariate(s)
#   Found 13002 Missing Data Values
#   Standardizing Data across genes
#   Fitting L/S model and finding priors
#   Finding parametric adjustments
# #   Adjusting the Data

par(mfrow = c(1, 1))
for (i in 1:(ncol(quantiled4))){
  if(i==1) plot(density(quantiled4[, i], na.rm=T), col = i, ylim = c(0,.85))
  else lines(density(quantiled4[, i], na.rm=T), col = i)
}
  
  ##batch effect correction using sva combat and 'covariate' matrix
  
  # now for the full dataset n=8560
  cdata <- na.omit(com2)#also used below
  prince.plot(prince(cdata,o1,top=10)) #huzzah!\
  
  
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
  
  plot(dend.complete,ylab="height", main = "PhosProt Combat Euclidian/Complete")
  plot(dend.ward, leaflab = "perpendicular", ylab = "height", main = "PhosProt Combat Euclidian/Ward")
  
  plot.new()##produces a blank canvas
  
  # Cluster using euclidian distance and ward linkage for both sites(rows) and samples (columns)
  # Note that both dendograms are created independently and row Z scores are presented in the heatmap
  
  # row scaled
  r <- t(scale(t(cdata)))#transpose to zscale the rows then transpose back to original format
  
  # sample scaled
  c <- scale(cdata)
  
  
  # install heatmap.2 package
  # install.packages("gplots")
  
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
  
  
  # plot.new()
  
  #PCA analysis 
  x <- t(cdata)#samples are the rows of the column matrix
  pc <- prcomp(x)#scale = T, center = T) as of now I am not scaling
  
  names(pc)
  
  cols <- as.factor(substr(colnames(cdata), 3, 7))##check me out. use 5 digit exp name.
  plot(pc$x[, 1], pc$x[, 2], col=as.numeric(cols), main = "PhosProt PCA", xlab = "PC1", ylab = "PC2")
  legend("bottomleft", levels(cols), col = seq(along=levels(cols)), pch = 1)
  
  
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
    
  
  #At least two calls in each batch (for Combat)
  DFs <- list(multExpanded1, RawRatios, MedianNorm, quantiled, quantiled2, quantiled3, quantiled4, quantiled5, quantiledBio,
              com2, adata, pilot)
  saveRDS(DFs, file = "./PhosPrepMatrices.rds")
  return(DFs)
}



  





  