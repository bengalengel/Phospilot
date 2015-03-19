NormProt <- function(directory){
  #this function loads MQoutput from SILAC quantitative proteomic analysis of 60 LCLsm reduces it to the three samples of interest, median and quantile normalizes it. Returns a list of three data frames. MQoutput, median normalized, and quantile normalized. 
  source("loadMQZ.R")
  require(plyr)
  require(limma)
  
  # load protein files with particular variables populated using "loadMQ"
  protein <- load.MQZ(directory)#7710 protein groups/81 variables
  
  # remove contaminants and reverse database hits
  protein <- protein[(protein$Potential.contaminant != "+" & protein$Reverse != "+"),]#7354/78
  
  # "only identified by site" hits CAN BE removed because they tend to have lower PEPs (wouldn't pass the FDR TH anyway) and can't be quantified since they are not idd by non-modified peptides. 
  protein1 <- protein[(protein$Only.identified.by.site != "+"),]#6705
  
  colnames(protein1)<- gsub(colnames(protein1), pattern = "Ratio.H.L.normalized.", replacement = "HL") ##remove redundant information
  
  #some strangeness sometimes there are two extra rows!
  #protein1 <- protein1[2:length(protein1)]
  
  #remove proteins if not quantified in at least one sample
  expCol <- grep("HL(.*)", colnames(protein1))
  
  protein1 <- protein1[rowSums(is.na(protein1[,expCol]))!=length(expCol),]##removes rows containing all 
  #'NA's using the sums of the logical per row 
  #6421
  
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
  
  