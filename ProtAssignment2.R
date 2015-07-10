ProtAssignment2 <- function(proteinfull, proteinnorm, multExpanded1_withDE, phosphonorm, proteome){
  require(qdapRegex)
  ##This program assigns protein groups (from previous proteomic analysis) to phosphosites from the SCX-TiO2 workflow for normalization. It also brings along the ibaq estimates from those groups in the pertinant samples. 
  
  #   It also does the same thing for the protein groups identified/quantified via SCX-TiO2.
  
  ##proteinfull is the 60 estimates from zias work. proteinnorm is the quantile normalized protein data.  phosphonorm is the normalized/batch corrected /confounded phospho data.multExpanded1_withDE is the parent dataframe for the class one sites with DiffPhos annotation from using the confounded data. Proteome is the passed fasta file used for the database search. 
  
  #subset proteinfull by protein group id found in the three normalized samples of interest. 
  proteinfull <- proteinfull[proteinfull$id %in% row.names(proteinnorm),]
  
  #revert names of proteinnorm back to 'HL' for continuity with the structure below. (Changed back at the end)
  names(proteinnorm) <- gsub(names(proteinnorm), pattern = "LH", replacement = "HL")
  
  ##combine the normalized protein information and the annotation data into a common dataframe for matching to the phospho data
  #   ibaq values to append
  ibaq <- grep("ibaq(.*)18862|ibaq(.*)18486|ibaq(.*)19160", names(proteinfull), ignore.case = T, value = T)
  
  
  datacomp <- cbind(proteinnorm,proteinfull[c("Protein.IDs","Majority.protein.IDs","Sequence.coverage....",
                                              "Number.of.proteins", "Sequence.length", "Sequence.lengths", "Peptides", 
                                              "Razor...unique.peptides", "Unique.peptides",
                                              "Razor...unique.peptides.18862", "Razor...unique.peptides.18486", 
                                              "Razor...unique.peptides.19160", ibaq)])
  
  #Descriptive
  #how many unique proteins in the phospho data are identified and subjected to DE analysis? Below are counts with and without isoform designation.
  ######################
  SubtoDE <- multExpanded1_withDE[as.character(multExpanded1_withDE$SubtoDE) == "+",]
  
  #remove the reverse (ME should be passed this way) and add 'AllPhos', a data frame for all phosphopeptides.
  SubtoDE <- SubtoDE[!grepl(SubtoDE$Protein, pattern = "REV"),]#the logical subsets the DF by row. 4991 observations.
  AllPhos <- multExpanded1_withDE[!grepl(multExpanded1_withDE$Protein, pattern = "REV"),] #18238
  
  #How many 'protiens' subject to DE? These are majority protein ids from the protein group assigned to the phosphopeptide.
  SubtoDEtable <- as.matrix(table(SubtoDE$Protein))
  SubtoDEtable <- SubtoDEtable[SubtoDEtable!=0,,drop=F]#1991 proteins
  
  #How many 'proteins' in the entire class1 dataset?
  AllPhostable <- as.matrix(table(AllPhos$Protein))
  AllPhostable <- AllPhostable[AllPhostable!=0,,drop=F]#4184 proteins
  
  #Remove isoform designation (this should be noted in the discussion of the paper. Unless there is a specific sequence aligning the peptide to a particular isoform the peptide could reasonably belong to any of the isoforms. This workflow (like those that use only the annotated uniprot databases) assumes this is not the case. Only a minor issue using this approach anyway 1991 to 1968 and 4184 to 4127.
  
  row.names(SubtoDEtable) <- substr(row.names(SubtoDEtable),1,6)#eventually 
  SubtoDEproteins <- row.names(SubtoDEtable)#1991
  SubtoDEproteins <- unique(SubtoDEproteins)#1968 unique proteins (excluding isoforms) subjected to DE
  
  row.names(AllPhostable) <- substr(row.names(AllPhostable),1,6)#eventually 
  AllPhosproteins <- row.names(AllPhostable)#4184
  AllPhosproteins <- unique(AllPhosproteins)#4127 unique proteins (excluding isoforms) subjected to DE
  ###############
  
  #Descriptive
  #How many proteins subjected to DE analysis are also IDd and quantified by proteomic analysis (Zia)? 
  ###############################################
  #Here I will use the majority protein IDs. I will use the majority protein IDs for each protein group quantification (these proteins have at least half the peptides of the leading protein within the group)
  Ziaproteins <- datacomp$Majority.protein.IDs
  
  #unparse and turn into a long string for comparison
  #now I need to design a loop to dig into the phospho file 
  pcount <- strsplit(as.character(Ziaproteins), ";")
  pcount <- as.character(unlist(pcount))
  pcount <- unique(pcount)
  any(duplicated(pcount))#12642 proteins identified with at least 1/2 of the peptides of the majority protein(s) in Zias work within 4270 protein groups
  
  pcount <- substr(pcount,1,6)#now I have some duplicates due to isoform designation (3500 duplicates)
  pcount <- unique(pcount) #Now I have 8885 unique protein with at least 1/2 of the peptides of the majority protein(s) within 4270 groups.
  
  table(SubtoDEproteins%in%pcount)#isoform free
  #1206 of 1968  sub to diffphos (61.2% are quantified in all three samples in Zia's work)
  
  table(AllPhosproteins%in%pcount)
  #2117 of 4127 (51% of proteins)
  #####################
  
  #Assigning protein group quantifications to phosphopeptide quantifications
  #####################################################################
  
  ##find the peptide sequence in the fasta file, assign it to multiple proteins. Find out how many protein groups (from proteome workup) contain the matching proteins. Assign peptide to protein group with the most unique plus razor ids. 
  
  #after I do this I will want to spot check the assignments to see if the protein assigned by MQ in the phospho sample corresponds to the protein assigned using only the proteome sample
  
  
  #############
  #how many fasta sequences are in this file?
  #     length(proteome)#
  #     [1] 88993
  #     
  
  # Each element is a sequence object of the class SeqFastadna or SeqFastaAA.
  # is.SeqFastaAA(proteome[[1]])
  # [1] TRUE
  
  #from interesting list i = 2623
  ################
  
  protein_norm <- data.frame()
  
  #   ibaq values to append
  ibaq <- grep("ibaq(.*)18862|ibaq(.*)18486|ibaq(.*)19160", names(proteinfull), ignore.case = T, value = T)
  
  for(i in seq_along(multExpanded1_withDE[,1])){
#     for(i in 1:50){
      
    
    #get the peptide to search
    peptide <- as.character(multExpanded1_withDE$Phospho..STY..Probabilities[i])
    peptide <- gsub(pattern = " *\\(.*?\\) *", replacement = "", peptide)
    
    #search the fasta database and retrieve uniprot/ensp ids corresponding to matches.
    matches <- grep(peptide, proteome)
    
    #Retrieve position(s) of phosphosite within protein match(s).
    PositionInProteins <- c()
    for(i in seq_along(matches)){
    seq <- unlist(getSequence(object = proteome[[matches[1]]], as.string = T))
    #position of peptide within protein
    protpos <- regexpr(peptide,seq)
    #now I need the position of the modified site within the peptide
    windows <- as.character(multExpanded1_withDE$Sequence.window[i])
    windows <- unlist(strsplit(as.character(windows), ";"))
    window <- windows[1]
    sitepos <- regexpr(peptide,window)
    modsiteinpeptide <- 16 - sitepos[1]#modified site is always the 16th position within the 'sequence window'.
    ##add this number to protein position!!!! Hurray
    tmp <- protpos[1] + modsiteinpeptide    
    PositionInProteins <- c(PositionInProteins, tmp)
    }
    #produce 'matches' character vector of utility
    matches <- names(proteome)[matches]
    matches <- sub(pattern = "\\|", replacement = "(", matches)
    matches <- sub(pattern = "\\|", replacement = ")", matches)
    matches <- sapply(rm_round(matches, extract = T), paste, collapse = "")
    
    #search the protein file from Zia's dataset to find protein group(s) containing at least one of these identifiers containing the peptide
    #Note that there will always be something toMatch because the peptides were identified using this database.     
    toMatch <- paste(matches, collapse = "|")
    matchingGroups <- grep(toMatch, as.character(datacomp$Protein.IDs))
    matchingGroups <- datacomp[matchingGroups,]
    
    #assign the matching group's ids and values for samples of interest to protein_norm data frame
    if(nrow(matchingGroups) > 1){#which has the most counts?
      counts <- matchingGroups$Razor...unique.peptides
      assignedProtein <- matchingGroups[which.max(counts),] #need to fix this to return desired numbers
      tmp <- assignedProtein[c("Protein.IDs","Majority.protein.IDs","Sequence.coverage....", "Sequence.length", "Sequence.lengths", 
                               "HL18862", "HL18486", "HL19160", ibaq)]
      #protein_norm <- rbind(protein_norm,tmp)
    }
    if(nrow(matchingGroups) == 1){
      tmp <- matchingGroups[c("Protein.IDs","Majority.protein.IDs","Sequence.coverage....", "Sequence.length", "Sequence.lengths", 
                              "HL18862", "HL18486", "HL19160", ibaq)]
      #protein_norm <- rbind(protein_norm,tmp)
    }
    if(nrow(matchingGroups)==0){#set names of dataframe in the event the first loop doesn't produce a match
      tmp <- as.data.frame(t(rep(NA,17)))
      names(tmp) <- c("Protein.IDs","Majority.protein.IDs","Sequence.coverage....", "Sequence.length", "Sequence.lengths", 
                      "HL18862", "HL18486", "HL19160", ibaq)
      #protein_norm <- rbind(protein_norm,tmp)
    }
    
    # prune the protein_norm dataframe s.t. proteins that can't contain the phosphopeptide are removed. This is for annotation enrichment.
    # I am first assuming that I can work with tmp outside of the if statements. I really need to learn about environments.
    
    Protein.IDs <- strsplit(as.character(tmp$Protein.IDs), ";")
    Protein.IDs <- unlist(Protein.IDs)
    Protein.IDsF <- Protein.IDs[Protein.IDs %in% matches]
    Protein.IDs <- paste(Protein.IDsF, collapse = ";")
    tmp$Protein.IDs <- Protein.IDs
    
    #add the position of the modified site information based on matching protein ids
    PositionInProteins <- PositionInProteins[matches %in% Protein.IDsF]
    tmp$PositionInProteins <- paste(PositionInProteins, collapse = ";")
    
    Majority.protein.IDs <- strsplit(as.character(tmp$Majority.protein.IDs), ";")
    Majority.protein.IDs <- unlist(Majority.protein.IDs)
    Majority.protein.IDs <- Majority.protein.IDs[Majority.protein.IDs %in% matches]
    Majority.protein.IDs <- paste(Majority.protein.IDs, collapse = ";")
    tmp$Majority.protein.IDs <- Majority.protein.IDs
    
    #bind the dataframe
    protein_norm <- rbind(protein_norm,tmp)
  }
  
  #make the protein_norm names specific to protein workup and be explicit about the ratios
  #affix string to the beginning of each element in the character vector
  
  ibaqnames <- paste("pp",ibaq, sep = "")
  
  names(protein_norm) <- c("ppProteinIDs", "ppMajorityProteinIDs", "ppSequenceCoverage", "ppSequence.length", "ppSequence.lengths",
                           "LH18862", "LH18486", "LH19160", "ppPositionInProteins", ibaqnames)
  
  #link the protein quants to the phospho ids to make a dataframe with normalized protein quants appended. Note "REV_" entries are removed again within this function in case they were passed accidentally.
  AllPhos <- cbind(AllPhos, protein_norm)
  
  #############################
  
  #Normalize the passed phospho dataframe for return and produce EDA plots on normalized dataframe.
  ################################################################
  
  
  #subset phosphonorm s.t. all phosphopeptides were mapped to a protein in Zia's data
  
  #first generate the ids for subsetting
  normphos <- AllPhos[,c("idmult", "Protein", "Leading.proteins")]
  normphos <- cbind(normphos,protein_norm)
  MappedPhos <- na.omit(normphos)#12033 phosphopeptides are mapped to a protein in zia's dataset
  
  #subset phosphonorm
  phosphonorm <- phosphonorm[row.names(phosphonorm) %in% MappedPhos$idmult, ]
  
  #now subset the MappedPhospho
  MappedPhos <- MappedPhos[MappedPhos$idmult %in% row.names(phosphonorm), ]
  
  #combine the DFs
  phosphonorm <- cbind(phosphonorm,MappedPhos)
  
  #normalize and subject to EDA and varcomp. protein will be run as a covariate for DiffPhos.
  
  ##now subset and normalize
  expCol <- grep("HL(.*)|LH(.*)", colnames(phosphonorm))
  data <- phosphonorm[,expCol]
  row.names(data) <- phosphonorm$idmult
  
  #perform the normalization
  
  #perform the normalization
  HL18486_1_1norm <- data$HL18486_1_1-data$LH18486
  HL18486_1_2norm <- data$HL18486_1_2-data$LH18486
  HL18486_2_1norm <- data$HL18486_2_1-data$LH18486
  HL18486_2_2norm <- data$HL18486_2_2-data$LH18486
  HL18862_1_1norm <- data$HL18862_1_1-data$LH18862
  HL18862_1_2norm <- data$HL18862_1_2-data$LH18862
  HL18862_2_1norm <- data$HL18862_2_1-data$LH18862
  HL18862_2_2norm <- data$HL18862_2_2-data$LH18862
  HL19160_1_1norm <- data$HL19160_1_1-data$LH19160
  HL19160_1_2norm <- data$HL19160_1_2-data$LH19160
  HL19160_2_1norm <- data$HL19160_2_1-data$LH19160
  HL19160_2_2norm <- data$HL19160_2_2-data$LH19160
  
  ProtNormalized <- cbind(HL18486_1_1norm, HL18486_1_2norm, HL18486_2_1norm, HL18486_2_2norm, HL18862_1_1norm, HL18862_1_2norm, HL18862_2_1norm,
                          HL18862_2_2norm, HL19160_1_1norm, HL19160_1_2norm, HL19160_2_1norm, HL19160_2_2norm)
  row.names(ProtNormalized) <- row.names(data)
  ##############
  
  #summary plots
  ###############
  boxplot(ProtNormalized)
  par(mfrow = c(1, 1))
  for (i in 1:(ncol(ProtNormalized))){
    if(i==1) plot(density(ProtNormalized[, i], na.rm=T), col = i, ylim = c(0,1.5))
    else lines(density(ProtNormalized[, i], na.rm=T), col = i)
  }
  
  #clustering and EDA using complete cases
  ProtNormalized2 <- na.omit(ProtNormalized)
  
  dataZ <- scale(ProtNormalized2)##Z-scored column wise the complete data matrix
  
  # dendogram using euclidian distance (default) and ward or complete agglomeration
  dend.ward<- as.dendrogram(hclust(dist(t(dataZ)),method="ward"))
  dend.complete<- as.dendrogram(hclust(dist(t(dataZ))))
  
  ward.o<- order.dendrogram(dend.ward)
  complete.o<- order.dendrogram(dend.complete)
  
  plot(dend.complete,ylab="height", main = "Euclidian/Complete")
  plot(dend.ward, leaflab = "perpendicular", ylab = "height", main = "Euclidian/Ward")
  
  # row scaled
  r <- t(scale(t(ProtNormalized2)))#transpose to zscale the rows then transpose back to original format
  
  # sample scaled
  c <- scale(ProtNormalized2)
  
  require(gplots)
  
  # Create dendrogram using the data without NAs
  feature.dend<- as.dendrogram(hclust(dist(r),method="ward"))
  sample.dend<- as.dendrogram(hclust(dist(t(c)),method="ward"))##note that dist caclculates distance between rows by default
  
  
  ##produce the heatmap. 
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
    margins = c(7,6),
    cexCol=1,
    labRow = NA#remove row labels
  )
  # plot.new()
  
  #PCA analysis 
  x <- t(ProtNormalized2)#samples are the rows of the column matrix
  pc <- prcomp(x)#scale = T, center = T) as of now I am not scaling
  
  cols <- as.factor(substr(colnames(ProtNormalized2), 3, 7))##use 5 digit exp name.
  plot(pc$x[, 1], pc$x[, 2], col=as.numeric(cols), main = "PCA", xlab = "PC1", ylab = "PC2")
  legend("bottomleft", levels(cols), col = seq(along=levels(cols)), pch = 1)
  
  
  summary(pc)
  
  #SVD for calculating variance explained; see Rafa's notes for an explaination
  cx <- sweep(x, 2, colMeans(x), "-")
  sv <- svd(cx)
  names(sv)
  #plot(sv$u[, 1], sv$u[, 2], col = as.numeric(cols), main = "SVD", xlab = "U1", ylab = "U2")
  plot(sv$d^2/sum(sv$d^2), xlim = c(1, 12), type = "b", pch = 16, xlab = "principal components", 
       ylab = "variance explained")
  
  
  DFs <- list(AllPhos,ProtNormalized)
  
  return(DFs)
}