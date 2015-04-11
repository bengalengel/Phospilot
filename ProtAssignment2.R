##find the peptide sequence in the fasta file, assign it to multiple proteins. Find out how many protein groups contain the matching proteins. Assign peptide to protein group with the most ids. 

#after I do this I will want to spot check the assignments to see if the protein assigned by MQ in the phospho sample corresponds to the protein assigned using only the proteome sample

require(seqinr)
##read in the proteome fasta file
proteome <- read.fasta( file = "E:/My Documents/MQ SS/fasta/HUMAN.fasta", seqtype = "AA", as.string = TRUE)#will have to update

proteome2 <- read.fasta( file = "E:/My Documents/MQ SS/fasta/HUMAN.fasta", seqtype = "AA", as.string = TRUE, set.attributes = FALSE)#no attributes

#how many fasta sequences are in this file?
#     length(proteome)#86725
#     [1] 86725
#     

# Each element is a sequence object of the class SeqFastadna or SeqFastaAA.
# is.SeqFastaAA(proteome[[1]])
# [1] TRUE

require(qdapRegex)

for(i in seq_len(rows(multExpanded1_withDE)))

#from interesting list i = 2623
  
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
toMatch <- paste(matches, collapse = "|")
matchingGroups <- grep(toMatch, as.character(datacomp$Protein.IDs))
matchingGroups <- datacomp[matchingGroups,]

if(nrow(matchingGroups) != 0){
  if(nrow(matchingGroups) > 1){#which has the most counts?
    counts <- matchingGroups$Razor...unique.peptides
    matchingGroup <- matchingGroups[which.max(counts),] #need to fix this to return desired numbers
    tmp <- datacomp[c("Majority.protein.IDs","HL18862", "HL18486", "HL19160")]
    
    
    
    else{
      
    }
  }
}

if more than one match pick the one with the most ids




#retrieve the sequence for a given uniprot ID. This will be within the loop
toMatch <- unlist(proteins)
toMatch <- paste(toMatch, collapse = "|")
matches <- (grep(toMatch,names(proteome), value=TRUE))
tt <- getSequence(proteome[matches], as.string = T)
tt <- sapply(tt,as.character)
