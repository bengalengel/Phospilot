

require(foreach)
require(doParallel)
require(qdapRegex)
require(iterators)

cl <- makeCluster(5)
registerDoParallel(cl)

protein_norm <- data.frame()
#   for(i in seq_along(multExpanded1_withDE[,1])){
foreach(i=1:50, .packages = "qdapRegex") %do% {
  
  #get the peptide to search
  peptide <- as.character(multExpanded1_withDE$Phospho..STY..Probabilities[i])
  peptide <- gsub(pattern = " *\\(.*?\\) *", replacement = "", peptide)
  
  #search the fasta database and retrieve uniprot ids corresponding to matches. Soetimes regex is a mfer
  matches <- grep(peptide, proteome)
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
  %:% when(nrow(matchingGroups) > 1) %do% {#which has the most counts?
    counts <- matchingGroups$Razor...unique.peptides
    assignedProtein <- matchingGroups[which.max(counts),] #need to fix this to return desired numbers
    tmp <- assignedProtein[c("Protein.IDs","Majority.protein.IDs","Sequence.coverage....", "HL18862", "HL18486", "HL19160")]
    #protein_norm <- rbind(protein_norm,tmp)
  }
  %:% when(nrow(matchingGroups) == 1) %do% {
    tmp <- matchingGroups[c("Protein.IDs","Majority.protein.IDs","Sequence.coverage....", "HL18862", "HL18486", "HL19160")]
    #protein_norm <- rbind(protein_norm,tmp)
  }
  %:% when(nrow(matchingGroups)==0) %do% {#set names of dataframe in the event the first loop doesn't produce a match
    tmp <- as.data.frame(t(rep(NA,6)))
    names(tmp) <- c("Protein.IDs","Majority.protein.IDs","Sequence.coverage....", "HL18862", "HL18486", "HL19160")
    #protein_norm <- rbind(protein_norm,tmp)
  }
  
  # prune the protein_norm dataframe s.t. proteins that can't contain the phosphopeptide are removed. This is for annotation enrichment.
  # I am first assuming that I can work with tmp outside of the if statements. I really need to learn about environments.
  
  Protein.IDs <- strsplit(as.character(tmp$Protein.IDs), ";")
  Protein.IDs <- unlist(Protein.IDs)
  Protein.IDs <- Protein.IDs[Protein.IDs %in% matches]
  Protein.IDs <- paste(Protein.IDs, collapse = ";")
  tmp$Protein.IDs <- Protein.IDs
  
  Majority.protein.IDs <- strsplit(as.character(tmp$Majority.protein.IDs), ";")
  Majority.protein.IDs <- unlist(Majority.protein.IDs)
  Majority.protein.IDs <- Majority.protein.IDs[Majority.protein.IDs %in% matches]
  Majority.protein.IDs <- paste(Majority.protein.IDs, collapse = ";")
  tmp$Majority.protein.IDs <- Majority.protein.IDs
  
  #bind the dataframe
  protein_norm <- rbind(protein_norm,tmp)
}

stopCluster(cl)



