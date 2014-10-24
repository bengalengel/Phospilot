phoscount <- function(phospho,phospho1,multExpanded,multExpanded1){

  ## total sites, quantified sites, observations of the quantified site, and protein groups associated with the site for all and class 1
  #if passed
  
  ##number of unique sites
  sites <- nrow(phospho)
  
  ##number of quantified sites
  quantsites <- sum(!is.na(phospho$Ratio.H.L))
  
  # number of proteins mapped from quant sites
  tmp <- phospho[!is.na(phospho$Ratio.H.L),]
  proteinquantsites <- nrow(table(tmp$Protein))
  
  # number of leading proteins and protein groups
  leadproteinquantsites <- nrow(table(tmp$Leading.proteins))
  pgroupsquantsites <- nrow(table(tmp$Proteins))
  
  
  ##number of unique class 1 sites
  class1 <- nrow(phospho1)
  #number of quantified class 1 sites
  quantClass1 <- sum(!is.na(phospho1$Ratio.H.L))
  
  # number of proteins mapped from quant sites
  tmp <- phospho1[!is.na(phospho$Ratio.H.L),]
  c1proteinquantsites <- nrow(table(tmp$Protein))
  
  # number of leading proteins and protein groups
  c1leadproteinquantsites <- nrow(table(tmp$Leading.proteins))
  c1pgroupsquantsites <- nrow(table(tmp$Proteins))
  
  
  # output for table
  output <- data.frame(sites,quantsites, proteinquantsites, leadproteinquantsites, pgroupsquantsites,class1,quantClass1,
                       c1proteinquantsites,c1leadproteinquantsites,c1pgroupsquantsites)
  
  write.table(output,"phos_counts.csv",sep=',',col.names=T,row.names=F)
  
  
  
  
  # Overlap per sample at the phosphosite level
  # number of unique ids per experiment (CERTAINLY A MUCH BETTER WAY TO DO THIS!)
  
  nmeasure <- function(x) sum(!(is.na(x)))#function to count valid values
  
  #nmeasure(multExpanded$HL16770_1)##works
  #nmeasure(multExpanded[,expCol]) ##sums over all the columns! must use colwise (plyr) see below
  
  #t <- colwise(nmeasure,newnames)(multExpanded)#same as unique events above
  
  newnames <- colnames(multExpanded)[grep("_", colnames(multExpanded))]
  
  idBreakdown <- ddply(multExpanded,.(id), colwise(nmeasure,newnames))##breaks down by id number of mults observed per sample
  
  data <- idBreakdown
  
  # remove exp obs if not observed in each sample 
  data2 <- data[rowSums(data[ , 2:5]) >= 1 & rowSums(data[ , 6:9]) >= 1 & rowSums(data[ , 10:13]) >= 1,]    
  
  # remove exp obs if not observed in each bio replicate (not sure how to automate this for larger datasets)
  data3 <- data[rowSums(data[ , 2:3]) >= 1 & rowSums(data[ , 4:5]) >= 1 & rowSums(data[ , 6:7]) >= 1 & rowSums(data[ , 8:9]) >= 1 & 
                  rowSums(data[ , 10:11]) >= 1 & rowSums(data[ , 12:13]) >= 1,]                    
  
  # remove exp obs if not observed in every experiment
  data4 <- data[rowSums(data[, 2:13]) >= 12,]
  
  quantified <- nrow(data)
  quantifiedsample <- nrow(data2)
  quantifiedbio <- nrow(data3)
  quantifiedall <- nrow(data4)
  #overlap summary stats for the  
  output <- data.frame(quantified, quantifiedsample, quantifiedbio, quantifiedall)
  write.table(output, "phosquantoverlap.csv", sep = ",", col.names = T, row.names = F)
  
  
  
  
  # Overlap per sample at the phosphosite level C1
  # number of unique ids per experiment (CERTAINLY A MUCH BETTER WAY TO DO THIS!)
  
  newnames <- colnames(multExpanded1)[grep("_", colnames(multExpanded1))]
  
  idBreakdown <- ddply(multExpanded1,.(id), colwise(nmeasure,newnames))##breaks down by id number of mults observed per sample
  
  data <- idBreakdown
  
  # remove exp obs if not observed in each sample 
  data2 <- data[rowSums(data[ , 2:5]) >= 1 & rowSums(data[ , 6:9]) >= 1 & rowSums(data[ , 10:13]) >= 1,]    
  
  # remove exp obs if not observed in each bio replicate (not sure how to automate this for larger datasets)
  data3 <- data[rowSums(data[ , 2:3]) >= 1 & rowSums(data[ , 4:5]) >= 1 & rowSums(data[ , 6:7]) >= 1 & rowSums(data[ , 8:9]) >= 1 & 
                  rowSums(data[ , 10:11]) >= 1 & rowSums(data[ , 12:13]) >= 1,]                    
  
  # remove exp obs if not observed in every experiment
  data4 <- data[rowSums(data[, 2:13]) >= 12,]
  
  quantified <- nrow(data)
  quantifiedsample <- nrow(data2)
  quantifiedbio <- nrow(data3)
  quantifiedall <- nrow(data4)
  #overlap summary stats for the  
  output <- data.frame(quantified, quantifiedsample, quantifiedbio, quantifiedall)
  write.table(output, "c1phosquantoverlap.csv", sep = ",", col.names = T, row.names = F)
  
  
  
  #sample overlap at the observation level
  ##gives string of experiment and replicate
  experiments <- colnames(multExpanded)[grep("_", colnames(multExpanded))]
  
  data <- multExpanded[,experiments]
  
  # remove exp obs if not observed in each sample 
  data2 <- data[rowSums(is.na(data[ , 1:4])) < 4 & rowSums(is.na(data[ , 5:8])) < 4 & rowSums(is.na(data[ , 9:12])) < 4,]    
  
  # remove exp obs if not observed in each bio replicate (not sure how to automate this for larger datasets)
  data3 <- data[rowSums(is.na(data[ , 1:2])) < 2 & rowSums(is.na(data[ , 3:4])) < 2 & rowSums(is.na(data[ , 5:6])) < 2 & 
                  rowSums(is.na(data[ , 7:8])) < 2 & rowSums(is.na(data[ , 9:10])) < 2 & rowSums(is.na(data[ , 11:12])) < 2,]                    
  
  # remove exp obs if not observed in every experiment
  data4 <- na.omit(data)
  
  quantified <- nrow(data)
  quantifiedsample <- nrow(data2)
  quantifiedbio <- nrow(data3)
  quantifiedall <- nrow(data4)
  #overlap summary stats for the  
  output <- data.frame(quantified, quantifiedsample, quantifiedbio, quantifiedall)
  write.table(output, "phosquantoverlap_obs.csv", sep = ",", col.names = T, row.names = F)
  
  #class1 sample overlap at the observation level for class 1
  ##gives string of experiment and replicate
  experiments <- colnames(multExpanded1)[grep("_", colnames(multExpanded))]
  
  data <- multExpanded1[,experiments]
  
  # remove exp obs if not observed in each sample 
  data2 <- data[rowSums(is.na(data[ , 1:4])) < 4 & rowSums(is.na(data[ , 5:8])) < 4 & rowSums(is.na(data[ , 9:12])) < 4,]    
  
  # remove exp obs if not observed in each bio replicate (not sure how to automate this for larger datasets)
  data3 <- data[rowSums(is.na(data[ , 1:2])) < 2 & rowSums(is.na(data[ , 3:4])) < 2 & rowSums(is.na(data[ , 5:6])) < 2 & 
                  rowSums(is.na(data[ , 7:8])) < 2 & rowSums(is.na(data[ , 9:10])) < 2 & rowSums(is.na(data[ , 11:12])) < 2,]                    
  
  # remove exp obs if not observed in every experiment
  data4 <- na.omit(data)
  
  quantified <- nrow(data)
  quantifiedsample <- nrow(data2)
  quantifiedbio <- nrow(data3)
  quantifiedall <- nrow(data4)
  #overlap summary stats for the  
  output <- data.frame(quantified, quantifiedsample, quantifiedbio, quantifiedall)
  write.table(output, "c1phosquantoverlap_obs.csv", sep = ",", col.names = T, row.names = F)
  
  
  
}