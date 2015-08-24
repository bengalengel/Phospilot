AddAnnotation <- function(multExpanded1_withDE){
#This function adds the annotations from GO, reactome, entrez and HGNC. Also may add phosphositeplus, and corum? ME DF passed with DiffPhos analysis already performed on confounded and non-confounded datasets. 

#curated kinase annotation and known modifications from phosphositeplus downloaded manually
# #files downloaded on 3/31/15
# PSPdwnload <- date()

#assign to multExpanded1 using leading proteins and majority protein ids for protein normalized data.
#GO terms and reactome added using R annotation packages


#source("http://bioconductor.org/biocLite.R")
#biocLite("OrganismDbi")
#biocLite("Homo.sapiens")
require(Homo.sapiens)
require(GO.db)
require(reactome.db)#only entrezids to reactome ids!
require(iterators)
require(foreach)
require(doParallel)
require(biomaRt)

#Examples using db objects
##############
# To list the kinds of things that can be retrieved, use the columns method.
# columns(Homo.sapiens)
# # To list the kinds of things that can be used as keys we can use the keytypes method
# keytypes(Homo.sapiens)
# # And to extract viable keys of a particular kind (keytype), we can use the keys method.
# idtest <- head(keys(Homo.sapiens, keytype="GOID"),25)
# # Since the keys method can tell us specific things that can be used as keys, here we will use it to extract a few ids to use for demonstrating the fourth method type.
# ids = head(keys(Homo.sapiens, keytype="UNIPROT"),25)
# # Once you have some ids that you want to look up data for, the select method allows you to map these ids as long as you use the columns argument to indicate what you need to know and the keytype argument to specify what kind of keys they are.
# select(Homo.sapiens, keys=ids, columns="GOID", keytype="UNIPROT")
# head(select(Homo.sapiens, keys=ids, columns="GOALL", keytype="UNIPROT"))#above is more direct
######################


# Adding annotations
######################
# use the ensembl 75 data wherever possible because depreciated identifiers lead to issues when using updated ensembl annotations. 
# Can't get reactome from ensembl 75 so will use latest reactome.db to retireve ids using HGNC identifiers.

#majority protein ids used for protein preps. In both cases this is the 'middle' ground between conservative and liberal.

#grab the ensembl75 database
#archive hsapiens mart
ensembl_75 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

#lets see how many we get and how many have an HGNC identifier etc etc.
ensembl_75_CCDS <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_peptide_id','description','hgnc_id','hgnc_symbol','gene_biotype',
                                        'transcript_biotype', 'status', 'transcript_status', 'ccds', "go_id"), 
                         filters = 'with_ccds', values = T, mart = ensembl_75)

ensembl_75_CCDS_EG <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_peptide_id',"entrezgene"), 
                        filters = 'with_ccds', values = T, mart = ensembl_75)#four values without entrez ids

ensembl_75_CCDS_pfam <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_peptide_id',"pfam"), 
                            filters = 'with_ccds', values = T, mart = ensembl_75)#four values without entrez ids


#GOIDs, (hugo gene nomenclature committee) hgncid, hgncnam and associated description

### Define function for adding GOIDs, description, hgncid, hgncsymbol, entrezgene ----------

annotate <- function(proteins, ensembl_75_CCDS, ensembl_75_CCDS_EG, ensembl_75_CCDS_pfam){
# This function accepts a LIST of unique protein groups and assigns annotation using passed annotation DF
  #not easy to download so manually curated from pfam as of 8-24-15
  pfam.phospho <- c("PF00498", "PF01846", "PF03166", "PF10401", "PF00244", "PF00533", "PF00400", "PF00659", "PF00397",#S/T
                    "PF00017", "PF08416", "PF00168",#Y
                    "PF00782", "PF00102", "PF13350", "PF06602", "PF04273", "PF03162", "PF14566", "PF14671", "PF04179", "PF05706", #phosphatase
                    "PF00069", "PF01636",  "PF07714", "PF03109", "PF03881", "PF06293", "PF01163", "PF01633", "PF10707", "PF06176", #kinase
                    "PF02958", "PF04655", "PF10009", "PF12260", "PF16474", "PF07914", "PF14531", "PF06734", "PF05445", "PF07387") #kinase
  
  
  
  #vector function flexible for preallocation
  GOIDs <- vector(mode = 'list', length = length(proteins))
  description <- vector(mode = 'list', length = length(proteins))
  hgncid <- vector(mode = 'list', length = length(proteins))
  hgncsymbol <- vector(mode = 'list', length = length(proteins))
  entrezgene <- vector(mode = 'list', length = length(proteins))
  ReactIDs <- vector(mode = 'list', length = length(proteins))
  pfamIDs <- vector(mode = 'list', length = length(proteins))
  pfamIDphospho <- vector(mode = 'list', length = length(proteins))
  
  cl <- makeCluster(5)#I have 8 cores but had a crash when using all 8
  registerDoParallel(cl)
  GOIDs <- foreach(i=1:length(proteins)) %dopar% {
    ENSPIDs <- strsplit(proteins[[i]], ";")
    ENSPIDs <- as.character(unlist(ENSPIDs))
    ENSPIDs <- paste(ENSPIDs, collapse = "|")
    index <- grep(ENSPIDs, ensembl_75_CCDS$ensembl_peptide_id)
    GOIDs[[i]] <- ensembl_75_CCDS$go_id[index]
  }
  stopCluster(cl)
  
  cl <- makeCluster(5)#I have 8 cores but had a crash when using all 8
  registerDoParallel(cl)
  description <- foreach(i=1:length(proteins)) %dopar% {
    ENSPIDs <- strsplit(proteins[[i]], ";")
    ENSPIDs <- as.character(unlist(ENSPIDs))
    ENSPIDs <- paste(ENSPIDs, collapse = "|")
    index <- grep(ENSPIDs, ensembl_75_CCDS$ensembl_peptide_id)
    description[[i]] <- unique(ensembl_75_CCDS$description[index])
  }
  stopCluster(cl)
  
  cl <- makeCluster(5)#I have 8 cores but had a crash when using all 8
  registerDoParallel(cl)
  hgncid <- foreach(i=1:length(proteins)) %dopar% {
    ENSPIDs <- strsplit(proteins[[i]], ";")
    ENSPIDs <- as.character(unlist(ENSPIDs))
    ENSPIDs <- paste(ENSPIDs, collapse = "|")
    index <- grep(ENSPIDs, ensembl_75_CCDS$ensembl_peptide_id)
    hgncid[[i]] <- unique(ensembl_75_CCDS$hgnc_id[index])
  }
  stopCluster(cl)
  
  cl <- makeCluster(5)#I have 8 cores but had a crash when using all 8
  registerDoParallel(cl)
  hgncsymbol <- foreach(i=1:length(proteins)) %dopar% {
    ENSPIDs <- strsplit(proteins[[i]], ";")
    ENSPIDs <- as.character(unlist(ENSPIDs))
    ENSPIDs <- paste(ENSPIDs, collapse = "|")
    index <- grep(ENSPIDs, ensembl_75_CCDS$ensembl_peptide_id)
    hgncsymbol[[i]] <- unique(ensembl_75_CCDS$hgnc_symbol[index])
  }
  stopCluster(cl)
  
  cl <- makeCluster(5)#I have 8 cores but had a crash when using all 8
  registerDoParallel(cl)
  entrezgene <- foreach(i=1:length(proteins)) %dopar% {
    ENSPIDs <- strsplit(proteins[[i]], ";")
    ENSPIDs <- as.character(unlist(ENSPIDs))
    ENSPIDs <- paste(ENSPIDs, collapse = "|")
    index <- grep(ENSPIDs, ensembl_75_CCDS_EG$ensembl_peptide_id)
    entrezgene[[i]] <- unique(ensembl_75_CCDS_EG$entrezgene[index])
  }
  stopCluster(cl)
  
  #now for the reactomeIds. Entrezgene id is the only valid keytype
  cl <- makeCluster(5)#I had a crash when using all 8
  registerDoParallel(cl)
  ReactIDs <- foreach(i=seq_along(entrezgene), .packages = "reactome.db") %dopar% {
    #define if hgnc is valid keytype for reactomedb
    pos_err <- tryCatch(select(reactome.db, keys=as.character(entrezgene[[i]]), columns="REACTOMEID", keytype="ENTREZID"),error=function(e) e)
    if(!inherits(pos_err, "error")){
      tmp <- select(reactome.db, keys=as.character(entrezgene[[i]]), columns="REACTOMEID", keytype="ENTREZID")
      ReactIDs[[i]] <- as.character(tmp$REACTOMEID)
    }
  }
  stopCluster(cl)
  
  
  #pfam ids
  cl <- makeCluster(5)
  registerDoParallel(cl)
  pfamIDs <- foreach(i=1:length(proteins)) %dopar% {
    ENSPIDs <- strsplit(proteins[[i]], ";")
    ENSPIDs <- as.character(unlist(ENSPIDs))
    ENSPIDs <- paste(ENSPIDs, collapse = "|")
    index <- grep(ENSPIDs, ensembl_75_CCDS_pfam$ensembl_peptide_id)
    tmp <- unique(ensembl_75_CCDS_pfam$pfam[index])
    if(tmp != ""){
      pfamIDs[[i]] <- unique(ensembl_75_CCDS_pfam$pfam[index])
    }else {
      pfamIDs[[i]] <- NA
    }
  }
  stopCluster(cl)
  
  #binary assignment of pfam id as phospho 'relevant'
  cl <- makeCluster(5)
  registerDoParallel(cl)
  pfamIDphospho <- foreach(i=1:length(pfamIDs)) %dopar% {
    if(any(pfamIDs[[i]] %in% pfam.phospho)){
      pfamIDphospho[[i]] <- "yes"
    } else {
      pfamIDphospho[[i]] <- "no"
    }
  }
  stopCluster(cl)
  
  
  proteins <- as.character(proteins)
  names(GOIDs) <- proteins
  names(description) <- proteins
  names(hgncid) <- proteins
  names(hgncsymbol) <- proteins
  names(entrezgene) <- proteins
  names(ReactIDs) <- proteins
  names(pfamIDs) <- proteins
  names(pfamIDphospho) <- proteins
  
  annotation.results <- list(GOIDs, description, hgncid, hgncsymbol, entrezgene, ReactIDs, pfamIDs, pfamIDphospho)
  names(annotation.results) <-  c("GOID", "Description", "HGNCID", "HGNCSymbol", "EntrezGene", "ReactIDs", "PFamIDs", "PFamIDPhospho")
  
  return(annotation.results)
}

###Bulk annotation add ----------

#confounded analysis
#'leading proteins' are used for the phospho alone annotated data. Leading proteins are the FIRST ENTRY from each matching protein group assigned to the peptide. Within a protein group ties are ignored. 'proteins' is all identified (potentially only by a modified site) proteins the peptide is associated with. 

proteins <- unique(multExpanded1_withDE$Leading.proteins)#note that all members of the protein group are not used for annotation.
proteins <- as.character(proteins)
proteins <- proteins[proteins != ""]#some phosphosites are not assigned to a quantified protein
proteins <- as.list(proteins)

confounded.annotation <- annotate(proteins, ensembl_75_CCDS, ensembl_75_CCDS_EG, ensembl_75_CCDS_pfam)

###phosprep annotation
#Using protein ids containing the phosphopeptide from the protein group with the most razor and unique peptides
#note that this caN be used for confounded enrichment analyses as well
proteins <- unique(multExpanded1_withDE$PhosPrepMajorityProteinIDs)
proteins <- as.character(proteins)
proteins <- proteins[proteins != ""]#some phosphosites are not assigned to a quantified protein
proteins <- as.list(proteins)

PhosPrep.annotation <- annotate(proteins, ensembl_75_CCDS, ensembl_75_CCDS_EG, ensembl_75_CCDS_pfam)

###gelprep annotation
#Using protein ids containing the phosphopeptide from the protein group with the most razor and unique peptides
proteins <- unique(multExpanded1_withDE$ppMajorityProteinIDs)
proteins <- as.character(proteins)
proteins <- proteins[proteins != ""]#some phosphosites are not assigned to a quantified protein
proteins <- as.list(proteins)

GelPrep.annotation <- annotate(proteins, ensembl_75_CCDS, ensembl_75_CCDS_EG, ensembl_75_CCDS_pfam)


#Add these annotation IDs to the original ME dataframe. 
############ 

process.anno <- function(Anno, proteins){
  #Identify protein matches from annotation list
  Match <- sapply(proteins, FUN = function(x) {
    if(is.na(x)){"NA"}
    else{
    hits = which(names(Anno) %in% x)
    Anno[hits]
    }
  })
  Match <- unlist(lapply(Match, function(x) paste(x, collapse = ";")))
}

#Note that each element of each annotation list is itself a list
#process confounded. 
proteins <- as.character(multExpanded1_withDE$Leading.proteins)
confounded <- lapply(confounded.annotation, function(x) process.anno(Anno = x, proteins))
confounded <- do.call(cbind.data.frame, confounded)
names(confounded) <- paste("confounded", names(confounded), sep = "")

#process phosprep
proteins <- as.character(multExpanded1_withDE$PhosPrepMajorityProteinIDs)
PhosPrep <- lapply(PhosPrep.annotation, function(x) process.anno(Anno = x, proteins))
PhosPrep <- do.call(cbind.data.frame, PhosPrep)
names(PhosPrep) <- paste("PhosPrep", names(PhosPrep), sep = "")

#process gelprep
proteins <- as.character(multExpanded1_withDE$ppMajorityProteinIDs)
index <- which(proteins == "")
proteins[index] <- NA #some very strange behavior in recognizing ' "" '. had to do this first
GelPrep <- lapply(GelPrep.annotation, function(x) process.anno(Anno = x, proteins))
GelPrep <- do.call(cbind.data.frame, GelPrep)
names(GelPrep) <- paste("GelPrep", names(GelPrep), sep = "")


multExpanded1_withDE <- cbind(multExpanded1_withDE, confounded, PhosPrep, GelPrep)

save(multExpanded1_withDE, file = "./multExpanded_withDE_annotated.RData")
return(multExpanded1_withDE)
}

