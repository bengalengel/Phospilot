ProtAssignment <- function(protein, proteome, multExpanded1){
  ## This function performs a loop to assign protein ids, H/L values and ibaq values to phosphosites for normalization using "protein groups" file produced from the phospho workup. The 'proteome' FASTA file is used to assign modification positions to the ratios used for quantification. Intellegent naming of the columns is used here as well. 
  
  #inputs needed are ME/ME1 and 'protein' from SCX-IMAC protein groups file WITHOUT SUBSETTING BY ONLY IDENTIFIED BY SITE!! (that is not protein1)
  
  ## I require protein copynumber and/or iBAQ estimates to test for enrichment of expression level in differential phosphorylation subsets or NRE model quandrants
  
  #When there are multiple protein groups that could match the phosphosite I choose the one with the most unique/razor ids
  require(qdapRegex)
  
  
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
    for(i in seq_along(matches)){
      index <- grep(matches[i], names(proteome))
      seq <- unlist(getSequence(object = proteome[index], as.string = T))
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
  
  #At least two calls in each batch (for Combat)
  DFs <- list(multExpanded1, RawRatios, MedianNorm, quantiled, quantiled2, quantiled3, quantiled4, quantiled5, quantiledBio)
  
  return(DFs)
}



  





  