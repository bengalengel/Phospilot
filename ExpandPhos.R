ExpandPhos <- function(phospho){
  ##now to ensure each multiplicity has its own row and variables are condensed. I am converting from wide to long format to ensure that each
  # observation is uniquely represented in the data. Then I cast back and remove empty rows
  # This function accepts a dataframe 'phospho' loaded from MQ output
  require(plyr)
  require(reshape2)
  
  expression <- phospho[,grep("Ratio.H.L.normalized(.*)_[12]___", colnames(phospho))]
  
  melted <- melt(phospho, measure.vars = names(expression))
  
  # Here I split (that is add a variable identifier) the melted 'variable' column so that the sample, replicate, 
  # and multiplicity are now explicit
  melted <- cbind(melted, colsplit(melted$variable, "_", c("sample", "bio_tech_mult"))) ##first split
  melted <- cbind(melted, colsplit(melted$bio_tech_mult, "_", c("bio", "tech_mult"))) ##second split
  melted <- cbind(melted, colsplit(melted$tech_mult, "___", c("tech","multiplicity"))) ##third split
  melted$sample <- gsub(melted$sample, pattern = "Ratio.H.L.normalized.", replacement = "") ##remove redundant information next 3 lines
  drop <- c("bio_tech_mult","tech_mult")
  melted <- melted[,!(names(melted) %in% drop)]
  
  ##cast data so that each unique 'sample/replicate' combination has its own column populated by the measurement 'currently the value column'.  
  casted <- dcast(melted, ... ~ sample + bio + tech, value.var="value") ##close but creates extra rows
  
  
  ##produce the multiplicity explicit table **********************************
  
  ##gives index of experiment and replicate
  data <- grep("^[0-9]+_[12]_[12]", colnames(casted))
  
  ##gives string of experiment and replicate
  data2 <- colnames(casted)[grep("^[0-9]+_[12]_[12]", colnames(casted))]
  
  #produces a new string with proper alpha leading R variable names
  newnames <- paste0("HL",data2)
  
  # rename the experiment variables within the dataframe
  colnames(casted)[data] <- newnames
  
  ## columnwise application of mean to condense the dataframe.
  out <- ddply(casted, .(id, multiplicity), colwise(mean,newnames,na.rm=T))
  
  #merge with identifying information by id to produce the multiplicity expanded table (each obs has a row)
  other_data <- phospho[, -which(names(phospho) %in% names(expression))]         
    
  multExpanded <- merge(other_data, out, by="id")
  
  ## remove rows with only NAs in expression columns
  expCol <- grep("HL(.*)", colnames(multExpanded))
  
  multExpanded <- multExpanded[rowSums(is.na(multExpanded[,expCol]))!=length(expCol),]##removes rows containing all 
  #'NA's using the sums of the logical per row
  
  return(multExpanded)
  }


