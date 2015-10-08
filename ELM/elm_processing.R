# This script will annotate the phosphosite table with positional overlap to the ~ 1500 annotted elm instances.

#elm instances. around 1500 elm instances in the database. limit to 'true positives'. These instances need to be aligned to the proteins in the fasta file. I need the name of the elm instance and the map 

require(seqinr)
require(foreach)
require(iterators)
require(biomaRt)
require(stringr)

#download all annotated motif instances from human with true positive logic instance
if(!file.exists("./ELM")){
  dir.create("./ELM")
}
if(!file.exists("./ELM/ELM_instances_hs.tsv")){
  url <- "http://elm.eu.org/instances.tsv?q=*&instance_logic=true+positive&taxon=Homo+sapiens"  
  download.file(url, "./ELM/ELM_instances_hs.tsv")
  date.download <- date()
  write.csv(date.download, "./ELM/download.date.csv")
}

#do the same for the fasta sequences in order to retrieve actual motif string
if(!file.exists("./ELM/ELM_instances.fasta")){
  url <- "http://elm.eu.org/instances.fasta?q=*&instance_logic=true+positive&taxon=Homo+sapiens"  
  download.file(url, "./ELM/ELM_instances.fasta")
  date.download <- date()
  write.csv(date.download, "./ELM/fasta_download.date.csv")
}

##read in table
elm <- read.table("./ELM/ELM_instances_hs.tsv", sep = "\t", header = T, stringsAsFactors = F, skip = 5)

#read in fasta file. These correspond to the number of unique primary uniprot accessions (850).
elmfasta <- read.fasta("./ELM/ELM_instances.fasta", seqtype = "AA", as.string = T)


#retrieve motif substring for each entry


elm$motif <- foreach(i = 1:length(elm[[1]]), .combine = c) %do% {
  index <- grep(elm$Primary_Acc[i], names(elmfasta))
  # Work from fasta file without hyphens. i.e. use primary uniprot accession because hyphenated accessions will match ONLY to that particular
#   value in the fasta file
  if(length(index) > 1){
    # Note apply function returns the wrong answer
    subindex <- vector(mode = "logical", length = length(index))
    for(j in seq_along(index)){
      header <- names(elmfasta[index[j]])
      header <- unlist(strsplit(header,split = "\\|"))[2]
      subindex[j] <- !grepl("-", header)
    }
    index <- index[subindex]
  }
  motif <- getSequence.character(elmfasta[[index]])[elm$Start[i]:elm$End[i]]
  paste0(motif, collapse = "")
}

## a spot check confirms the accuracy of the specific motif assignments

##Match the ENSP ids to uniprot ids and add these to the elm table. Note that there are no hyphenated entries in this table
ensembl_75_CCDS_uniprot <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_peptide_id',"uniprot_swissprot_accession"), 
                              filters = 'with_ccds', values = T, mart = ensembl_75)


#add matched ENSPID to elm table
elm$ENSPID <- foreach(i = 1:length(elm[[1]]), .combine = c) %do% {
#   collect the uniprot accessions
  uni.accessions <- elm$Accessions[i]
  uni.accessions <- unlist(str_split(uni.accessions, " "))
  #find matching uniprot ids and return linked ENSPIDs
  ENSPID.match <- ensembl_75_CCDS_uniprot[ensembl_75_CCDS_uniprot$uniprot_swissprot_accession %in% uni.accessions, "ensembl_peptide_id"]
  paste0(ENSPID.match, collapse = ";")
}

#add indicator of phosphosite overlap with ELM annotated motif instance

all.proteins <- unlist(strsplit(multExpanded1_withDE_annotated$ppMajorityProteinIDs, ";"))
all.elm <- unlist(strsplit(elm$ENSPID, ";"))
length(intersect(all.proteins, all.elm))
#426 overlapping unique ENSPIDs. This is going to be close and there likely will not be many sites within elm instances

#now I want to align the uniprot sequences over the ensembl peptide sequences so that I can map the positions of the motif within the uniprot sequence onto the positions within the ensembl peptide id sequence

test <- foreach(i = 1:100, .combine = rbind) %do% {
  proteins <- multExpanded1_withDE_annotated$ppMajorityProteinIDs[i]
  proteins <- unlist(strsplit(proteins, ";"))
  proteins <- proteins[!grepl("REV", proteins)]#these are omitted when finding the phosphoposition
  sites <- multExpanded1_withDE_annotated$ppPositionInProteins[i]
  sites <- strsplit(as.character(sites), ";")
  sites <- as.character(unlist(sites))
  #find matches in the elm data frame
  #for each protein
#   elm.dataframes <- vector(mode = "list", length = length(proteins))
  position.overlap <- vector(mode = "logical", length = length(proteins))
for(j in seq_along(proteins)){
  #indentify an index to subset the elm data frame when ANY enspid matches 
  index <- sapply(elm$ENSPID, function(x){
    elm.proteins <- unlist(strsplit(x, ";"))
    proteins[j] %in% elm.proteins    
  })
  #subset to get elm.dataframe[j] of motifs
#   elm.dataframes[[j]] <- elm[index,]
  match.elm <- elm[index, c("ENSPID", "motif")]
  motif.hits <- match.elm$motif
  #where do these motifs reside within the matching enspid protein?
  proteome.index <- grep(proteins[j], names(proteome))
  hit.sequence <- unlist(getSequence(proteome[[proteome.index]], as.string = T))
  #a matrix of start and stop positions
  motif.indices <- do.call(rbind, str_locate_all(hit.sequence, motif.hits))
# start end
# [1,]   103 107
# [2,]   106 110
  #do these positions overlap with the paired phosphosite position?
  position.overlap[j] <- any(apply(motif.indices, 1, function(x) {
    sites[j] >= x[1] && sites[j] <= x[2]
  }))
}
return(any(position.overlap))
}
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  }
  motif.match <- sum(sapply(elm.dataframes, dim)[1,]) > 0
  
    
  elm.proteins <- elm[any(unlist(strsplit(elm$ENSPID, ";")) %in% proteins),]
  
  elm.proteins <- unlist(strsplit(elm.proteins, ";"))
  if(any(proteins %in% elm.proteins)


for each peptide
do any **gelprot** proteins match elm table (match key will be ENSP id)
for the matches find the region within the protein that matches the motif (start/stop)
Do any of the ENSPID linked phosphosite positions fall within the motif regions of the protein?
return answer to "any matches" as binary TRUE or FALSE
return matching elm accession number

test.match  <- "ENSP00000260731"
doublematch "two matches"
i <- 17 for multiple proteins



