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

#normal for loop takes way too long
##################
#for each unique set of proteins assigned to a site retrieve GOIDs.
# uniIDs <- c()
# tmp <- c()
# proteins <- unique(multExpanded1_withDE$Proteins)
# proteins <- as.character(proteins)
# GOIDs <- vector(mode = 'list', length = length(unique(multExpanded1_withDE$Proteins)))#vector function flexible for preallocation
# for(i in 1:length(proteins)){
#   #unique protein groups require annotation assignment
#   names(GOIDs)[i] <- proteins[i]
#   uniIDs <- strsplit(proteins[i], ";")
#   uniIDs <- as.character(unlist(uniIDs))
#   uniIDs <- substr(uniIDs,1,6)#de-isoform becuase if not it throws an error
#   pos_err <- tryCatch(select(Homo.sapiens, keys=uniIDs, columns="GOID", keytype="UNIPROT"),error=function(e) e)
#   if(!inherits(pos_err, "error")){
#     tmp <- select(Homo.sapiens, keys=uniIDs, columns="GOID", keytype="UNIPROT")#retrieve GOIDs I should use the list form!
#     GOIDs[[i]] <- as.character(tmp$GOID)
#   }
# }
#############################

# Adding annotations
######################
# use the ensembl 75 data wherever possible because depreciated identifiers lead to issues when using updated ensembl annotations. 
# Can't get reactome from ensembl 75 so will use latest reactome.db to retireve ids using HGNC identifiers.

#I will use the 'leading proteins' for the phospho alone annotated data. Leading proteins are the FIRST ENTRY from each matching protein group assigned to the peptide. Within a protein group ties are ignored. 'proteins' is all identified (potentially only by a modified site) proteins the peptide is associated with. 

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


#GOIDs, (hugo gene nomenclature committee) hgncid, hgncnam and associated description
####################
# proteins <- unique(multExpanded1_withDE$Leading.proteins)
# proteins <- as.character(proteins)
# proteins <- as.list(proteins)
# GOIDs <- vector(mode = 'list', length = length(unique(multExpanded1_withDE$Leading.proteins)))#vector function flexible for preallocation
# description <- vector(mode = 'list', length = length(unique(multExpanded1_withDE$Leading.proteins)))#vector function flexible for preallocation
# hgncid <- vector(mode = 'list', length = length(unique(multExpanded1_withDE$Leading.proteins)))#vector function flexible for preallocation
# hgncsymbol <- vector(mode = 'list', length = length(unique(multExpanded1_withDE$Leading.proteins)))#vector function flexible for preallocation
# #using 7 cores
# cl <- makeCluster(5)#I have 8 cores but had a crash when using all 8
# registerDoParallel(cl)
# system.time(foreach(i=1:10) %dopar% {
#   ENSPIDs <- strsplit(proteins[i], ";")
#   ENSPIDs <- as.character(unlist(ENSPIDs))
#   ENSPIDs <- paste(ENSPIDs, collapse = "|")
#   index <- grep(ENSPIDs, ensembl_75_CCDS$ensembl_peptide_id)
#   GOIDs[[i]] <- ensembl_75_CCDS$go_id[index]
# #   description[[i]] <- unique(ensembl_75_CCDS$description[index])
# #   hgncid[[i]] <- unique(ensembl_75_CCDS$hgnc_id[index])
# #   hgncsymbol[[i]] <- unique(ensembl_75_CCDS$hgnc_symbol[index])
#   }
# )
# # user  system elapsed 
# # 4.61    1.89  313.11 #huge (20X) improvement using the dataframe as opposed to calling 'homo.sapiens' and 3X improvement over lapply (below commented out)
# stopCluster(cl)
# getDoParName()
# names(GOIDs) <- proteins
# # user  system elapsed 
# 7.30    0.67 5016.08 

# #make proteins into a list. The search using lapply. multiple functions with lapply?...later :(
# proteins <- unique(multExpanded1_withDE$Proteins)
# proteins <- as.character(proteins)
# proteins <- as.list(proteins)
# system.time(GOIDs <- lapply(proteins, function(x) {
#   ENSPIDs <- strsplit(x, ";")
#   ENSPIDs <- as.character(unlist(ENSPIDs))
#   ENSPIDs <- paste(ENSPIDs, collapse = "|")
#   index <- grep(ENSPIDs, ensembl_75_CCDS$ensembl_peptide_id)
#   x <- ensembl_75_CCDS$go_id[index]
# }
# ))
# # user  system elapsed 
# # 1180.42    0.01 1181.32 
################
proteins <- unique(multExpanded1_withDE$Leading.proteins)
proteins <- as.character(proteins)
proteins <- as.list(proteins)
GOIDs <- vector(mode = 'list', length = length(unique(multExpanded1_withDE$Leading.proteins)))#vector function flexible for preallocation
description <- vector(mode = 'list', length = length(unique(multExpanded1_withDE$Leading.proteins)))#vector function flexible for preallocation
hgncid <- vector(mode = 'list', length = length(unique(multExpanded1_withDE$Leading.proteins)))#vector function flexible for preallocation
hgncsymbol <- vector(mode = 'list', length = length(unique(multExpanded1_withDE$Leading.proteins)))#vector function flexible for preallocation
entrezgene <- vector(mode = 'list', length = length(unique(multExpanded1_withDE$Leading.proteins)))#vector function flexible for preallocation


cl <- makeCluster(5)#I have 8 cores but had a crash when using all 8
registerDoParallel(cl)
GOIDs <- foreach(i=1:length(proteins)) %dopar% {
  ENSPIDs <- strsplit(proteins[[i]], ";")
  ENSPIDs <- as.character(unlist(ENSPIDs))
  ENSPIDs <- paste(ENSPIDs, collapse = "|")
  index <- grep(ENSPIDs, ensembl_75_CCDS$ensembl_peptide_id)
  GOIDs[[i]] <- ensembl_75_CCDS$go_id[index]
  #   description[[i]] <- unique(ensembl_75_CCDS$description[index])
  #   hgncid[[i]] <- unique(ensembl_75_CCDS$hgnc_id[index])
  #   hgncsymbol[[i]] <- unique(ensembl_75_CCDS$hgnc_symbol[index])
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
  #   hgncid[[i]] <- unique(ensembl_75_CCDS$hgnc_id[index])
  #   hgncsymbol[[i]] <- unique(ensembl_75_CCDS$hgnc_symbol[index])
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
  #   hgncsymbol[[i]] <- unique(ensembl_75_CCDS$hgnc_symbol[index])
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

proteins <- as.character(proteins)
names(GOIDs) <- proteins
names(description) <- proteins
names(hgncid) <- proteins
names(hgncsymbol) <- proteins
names(entrezgene) <- proteins


#repeat for the protein prep IDs assigned to phosphosites using Zia's data. Majority protein IDs used here.

# replace the name of the ProtPrep column
PPproteins <- unique(multExpanded1_withDE$ppMajorityProteinIDs)
PPproteins <- as.character(PPproteins)
PPproteins <- PPproteins[PPproteins != ""]#some phosphosites are not assigned to a quantified protein
PPproteins <- as.list(PPproteins)
ppGOIDs <- vector(mode = 'list', length = length(PPproteins))#vector function flexible for preallocation
ppdescription <- vector(mode = 'list', length = length(PPproteins))#vector function flexible for preallocation
pphgncid <- vector(mode = 'list', length = length(PPproteins))#vector function flexible for preallocation
pphgncsymbol <- vector(mode = 'list', length = length(PPproteins))#vector function flexible for preallocation
ppentrezgene <- vector(mode = 'list', length = length(PPproteins))#vector function flexible for preallocation

cl <- makeCluster(5)#I have 8 cores but had a crash when using all 8
registerDoParallel(cl)
ppGOIDs <- foreach(i=1:length(PPproteins)) %dopar% {
  ppENSPIDs <- strsplit(PPproteins[[i]], ";")
  ppENSPIDs <- as.character(unlist(ppENSPIDs))
  ppENSPIDs <- paste(ppENSPIDs, collapse = "|")
  index <- grep(ppENSPIDs, ensembl_75_CCDS$ensembl_peptide_id)
  ppGOIDs[[i]] <- ensembl_75_CCDS$go_id[index]
  #   description[[i]] <- unique(ensembl_75_CCDS$description[index])
  #   hgncid[[i]] <- unique(ensembl_75_CCDS$hgnc_id[index])
  #   hgncsymbol[[i]] <- unique(ensembl_75_CCDS$hgnc_symbol[index])
}
stopCluster(cl)

cl <- makeCluster(5)#I have 8 cores but had a crash when using all 8
registerDoParallel(cl)
ppdescription <- foreach(i=1:length(PPproteins)) %dopar% {
  ppENSPIDs <- strsplit(PPproteins[[i]], ";")
  ppENSPIDs <- as.character(unlist(ppENSPIDs))
  ppENSPIDs <- paste(ppENSPIDs, collapse = "|")
  index <- grep(ppENSPIDs, ensembl_75_CCDS$ensembl_peptide_id)
  ppdescription[[i]] <- unique(ensembl_75_CCDS$description[index])
  #   hgncid[[i]] <- unique(ensembl_75_CCDS$hgnc_id[index])
  #   hgncsymbol[[i]] <- unique(ensembl_75_CCDS$hgnc_symbol[index])
}
stopCluster(cl)

cl <- makeCluster(5)#I have 8 cores but had a crash when using all 8
registerDoParallel(cl)
pphgncid <- foreach(i=1:length(PPproteins)) %dopar% {
  ENSPIDs <- strsplit(PPproteins[[i]], ";")
  ENSPIDs <- as.character(unlist(ENSPIDs))
  ENSPIDs <- paste(ENSPIDs, collapse = "|")
  index <- grep(ENSPIDs, ensembl_75_CCDS$ensembl_peptide_id)
  pphgncid[[i]] <- unique(ensembl_75_CCDS$hgnc_id[index])
  #   hgncsymbol[[i]] <- unique(ensembl_75_CCDS$hgnc_symbol[index])
}
stopCluster(cl)

cl <- makeCluster(5)#I have 8 cores but had a crash when using all 8
registerDoParallel(cl)
pphgncsymbol <- foreach(i=1:length(PPproteins)) %dopar% {
  ENSPIDs <- strsplit(PPproteins[[i]], ";")
  ENSPIDs <- as.character(unlist(ENSPIDs))
  ENSPIDs <- paste(ENSPIDs, collapse = "|")
  index <- grep(ENSPIDs, ensembl_75_CCDS$ensembl_peptide_id)
  pphgncsymbol[[i]] <- unique(ensembl_75_CCDS$hgnc_symbol[index])
}
stopCluster(cl)

cl <- makeCluster(5)#I have 8 cores but had a crash when using all 8
registerDoParallel(cl)
ppentrezgene <- foreach(i=1:length(PPproteins)) %dopar% {
  ENSPIDs <- strsplit(proteins[[i]], ";")
  ENSPIDs <- as.character(unlist(ENSPIDs))
  ENSPIDs <- paste(ENSPIDs, collapse = "|")
  index <- grep(ENSPIDs, ensembl_75_CCDS_EG$ensembl_peptide_id)
  ppentrezgene[[i]] <- unique(ensembl_75_CCDS_EG$entrezgene[index])
}
stopCluster(cl)


#assign names 
PPproteins <- as.character(PPproteins)
names(ppGOIDs) <- PPproteins
names(ppdescription) <- PPproteins
names(pphgncid) <- PPproteins
names(pphgncsymbol) <- PPproteins
names(ppentrezgene) <- PPproteins


#now for the reactomeIds. Entrezgene id is the only valid keytype
ReactIDs <- vector(mode = 'list', length = length(unique(multExpanded1_withDE$Proteins)))#vector function flexible for preallocation
#using 7 cores
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
names(ReactIDs) <- proteins


ppReactIDs <- vector(mode = 'list', length = length(PPproteins))#vector function flexible for preallocation
cl <- makeCluster(5)
registerDoParallel(cl)
ppReactIDs <- foreach(i=seq_along(ppentrezgene), .packages = "reactome.db") %dopar% {
  pos_err <- tryCatch(select(reactome.db, keys=as.character(ppentrezgene[[i]]), columns="REACTOMEID", keytype="ENTREZID"),error=function(e) e)
  if(!inherits(pos_err, "error")){
    tmp <- select(reactome.db, keys=as.character(ppentrezgene[[i]]), columns="REACTOMEID", keytype="ENTREZID")
    ppReactIDs[[i]] <- as.character(tmp$REACTOMEID)
  }
}
stopCluster(cl)
names(ppReactIDs) <- PPproteins





############

# #first is GOIDs
# proteins <- unique(multExpanded1_withDE$Proteins)
# proteins <- as.character(proteins)
# uniIDs <- c()
# tmp <- c()
# GOIDs <- vector(mode = 'list', length = length(unique(multExpanded1_withDE$Proteins)))#vector function flexible for preallocation
# names(GOIDs) <- proteins
# #using 7 cores
# cl <- makeCluster(5)#I have 8 cores but had a crash when using all 8
# registerDoParallel(cl)
# system.time(GOIDs <- foreach(i=1:length(proteins), .packages = "Homo.sapiens") %dopar% {
#   ENSPIDs <- strsplit(proteins[i], ";")
#   ENSPIDs <- as.character(unlist(ENSPIDs))
#   pos_err <- tryCatch(select(ensembl_75, keys=ENSPIDs, columns="go_id", keytype="ensembl_peptide_id"),error=function(e) e)
#   if(!inherits(pos_err, "error")){
#     tmp <- select(ensembl_75, keys=ENSPIDs, columns="go_id", keytype="ensembl_peptide_id")#retrieve GOIDs I should use the list form!
#     GOIDs[[i]] <- as.character(tmp$GOID)
#   }
# }
# )
# stopCluster(cl)
# getDoParName()
# names(GOIDs) <- proteins
# # user  system elapsed 
# # 7.30    0.67 5016.08 
# 
# 
# #for the protein normalized data
# 
# # replace the name of the ProtPrep column
# names(multExpanded1_withDE)[grep(names(multExpanded1_withDE), pattern = "ProtPrep")] <- "ProtPrepIds"
# 
# proteins <- unique(multExpanded1_withDE$ProtPrepIds)
# proteins <- as.character(proteins)
# proteins <- na.omit(proteins)
# uniIDs <- c()
# tmp <- c()
# ppGOIDs <- vector(mode = 'list', length = length(unique(multExpanded1_withDE$ProtPrepIds)))#vector function flexible for preallocation
# names(ppGOIDs) <- proteins
# #using 7 cores
# cl <- makeCluster(5)#I have 8 cores but had a crash when using all 8
# registerDoParallel(cl)
# system.time(ppGOIDs <- foreach(i=1:length(proteins), .packages = "Homo.sapiens") %dopar% {
#   uniIDs <- strsplit(proteins[i], ";")
#   uniIDs <- as.character(unlist(uniIDs))
#   uniIDs <- substr(uniIDs,1,6)#de-isoform becuase if not it throws an error
#   pos_err <- tryCatch(select(Homo.sapiens, keys=uniIDs, columns="GOID", keytype="UNIPROT"),error=function(e) e)
#   if(!inherits(pos_err, "error")){
#     tmp <- select(Homo.sapiens, keys=uniIDs, columns="GOID", keytype="UNIPROT")#retrieve ppGOIDs I should use the list form!
#     ppGOIDs[[i]] <- as.character(tmp$GOID)
#   }
# }
# )
# stopCluster(cl)
# getDoParName()
# 
# names(ppGOIDs) <- proteins



#Produce a reactome list in a similar fashion. first EntrezIDs
# uniIDs <- c()
# tmp <- c()
# EntrezIDs <- vector(mode = 'list', length = length(unique(multExpanded1_withDE$Proteins)))#vector function flexible for
# #using 7 cores
# cl <- makeCluster(5)#I had a crash when using all 8
# registerDoParallel(cl)
# system.time(EntrezIDs <- foreach(i=1:length(proteins), .packages = "Homo.sapiens") %dopar% {
#   uniIDs <- strsplit(proteins[i], ";")
#   uniIDs <- as.character(unlist(uniIDs))
#   uniIDs <- substr(uniIDs,1,6)#de-isoform becuase if not it throws an error
#   pos_err <- tryCatch(select(Homo.sapiens, keys=uniIDs, columns="ENTREZID", keytype="UNIPROT"),error=function(e) e)
#   if(!inherits(pos_err, "error")){
#     tmp <- select(Homo.sapiens, keys=uniIDs, columns="ENTREZID", keytype="UNIPROT")#retrieve Entrezids I should use the list form!
#     EntrezIDs[[i]] <- as.character(tmp$ENTREZID)
#   }
# }
# )
# stopCluster(cl)
# names(EntrezIDs) <- proteins
# 
# 
# # ppEntrezIDs
# uniIDs <- c()
# tmp <- c()
# ppEntrezIDs <- vector(mode = 'list', length = length(unique(multExpanded1_withDE$ProtPrepIds)))#vector function flexible for
# #using 7 cores
# cl <- makeCluster(5)#I had a crash when using all 8
# registerDoParallel(cl)
# system.time(ppEntrezIDs <- foreach(i=1:length(proteins), .packages = "Homo.sapiens") %dopar% {
#   uniIDs <- strsplit(proteins[i], ";")
#   uniIDs <- as.character(unlist(uniIDs))
#   uniIDs <- substr(uniIDs,1,6)#de-isoform becuase if not it throws an error
#   pos_err <- tryCatch(select(Homo.sapiens, keys=uniIDs, columns="ENTREZID", keytype="UNIPROT"),error=function(e) e)
#   if(!inherits(pos_err, "error")){
#     tmp <- select(Homo.sapiens, keys=uniIDs, columns="ENTREZID", keytype="UNIPROT")#retrieve ppEntrezIDs I should use the list form!
#     ppEntrezIDs[[i]] <- as.character(tmp$ENTREZID)
#   }
# }
# )
# stopCluster(cl)
# names(ppEntrezIDs) <- proteins
# 
# 
# 
# #now for the reactomeIds
# ReactIDs <- vector(mode = 'list', length = length(unique(multExpanded1_withDE$Proteins)))#vector function flexible for preallocation
# #using 7 cores
# cl <- makeCluster(5)#I had a crash when using all 8
# registerDoParallel(cl)
# system.time(ReactIDs <- foreach(i=seq_along(EntrezIDs), .packages = "reactome.db") %dopar% {
#   pos_err2 <- tryCatch(select(reactome.db, keys=EntrezIDs[[i]], columns="REACTOMEID", keytype="ENTREZID"),error=function(e) e)
#   if(!inherits(pos_err2, "error")){
#     tmp <- select(reactome.db, keys=EntrezIDs[[i]], columns="REACTOMEID", keytype="ENTREZID")
#     ReactIDs[[i]] <- as.character(tmp$REACTOMEID)
#   }
# }
# )
# stopCluster(cl)
# names(ReactIDs) <- proteins
# 
# 
# ppReactIDs <- vector(mode = 'list', length = length(unique(multExpanded1_withDE$ProtPrepIds)))#vector function flexible for preallocation
# #using 7 cores
# cl <- makeCluster(5)#I had a crash when using all 8
# registerDoParallel(cl)
# system.time(ppReactIDs <- foreach(i=seq_along(EntrezIDs), .packages = "reactome.db") %dopar% {
#   pos_err2 <- tryCatch(select(reactome.db, keys=EntrezIDs[[i]], columns="REACTOMEID", keytype="ENTREZID"),error=function(e) e)
#   if(!inherits(pos_err2, "error")){
#     tmp <- select(reactome.db, keys=EntrezIDs[[i]], columns="REACTOMEID", keytype="ENTREZID")
#     ppReactIDs[[i]] <- as.character(tmp$REACTOMEID)
#   }
# }
# )
# stopCluster(cl)
# names(ppReactIDs) <- proteins

#################################################




#Add these annotation IDs to the original ME dataframe.

#GOIDS 
############ 

#Identify protein matches from GOIDs list to add to DF
GOmatch <- sapply(as.character(multExpanded1_withDE$Leading.proteins), FUN = function(x) {
  match = which(names(GOIDs) %in% x)
  GOIDs[[match]]
})
#I need GOmatch as a character vector of length one so that I may create a data.frame/matrix
GOmatch <- lapply(GOmatch, function(x) paste(x,collapse = ";"))
GOmatch <- as.matrix(GOmatch)#as data frame doesn't work becuase can't have duplicate row names.
#add to the dataframe
multExpanded1_withDE$GOIDs <- as.character(GOmatch[,1])


#Identify protein matches from ppGOIDs list

ppMajorityProteinIDs <- as.character(multExpanded1_withDE$ppMajorityProteinIDs)
index <- which(ppMajorityProteinIDs=="")
ppMajorityProteinIDs[index] <- NA #some very strange behavior in recognizing ' "" '. had to do this first
GOmatch <- sapply(ppMajorityProteinIDs, FUN = function(x) {
  if(is.na(x)){"NA"}
  else{
    match = which(names(ppGOIDs) %in% x)
    ppGOIDs[[match]]
  }
})
#I need GOmatch as a character vector of length one so that I may create a data.frame/matrix
GOmatch <- lapply(GOmatch, function(x) paste(x,collapse = ";"))
GOmatch <- as.matrix(GOmatch)#as data frame doesn't work becuase can't have duplicate row names.
#add to the dataframe
multExpanded1_withDE$ppGOIDs <- as.character(GOmatch[,1])

######################

#HGNCids
############ 

#Identify leading protein matches from hgncids list to add to DF
hgncidmatch <- sapply(as.character(multExpanded1_withDE$Leading.proteins), FUN = function(x) {
  match = which(names(hgncid) %in% x)
  hgncid[[match]]
})
#I need hgncidmatch as a character vector of length one so that I may create a data.frame/matrix
hgncidmatch <- lapply(hgncidmatch, function(x) paste(x,collapse = ";"))
hgncidmatch <- as.matrix(hgncidmatch)#as data frame doesn't work becuase can't have duplicate row names.
#add to the dataframe
multExpanded1_withDE$hgncid <- as.character(hgncidmatch[,1])


#Identify protein matches from ppGOIDs list

ppMajorityProteinIDs <- as.character(multExpanded1_withDE$ppMajorityProteinIDs)
index <- which(ppMajorityProteinIDs=="")
ppMajorityProteinIDs[index] <- NA #some very strange behavior in recognizing ' "" '. had to do this first
hgncidmatch <- sapply(ppMajorityProteinIDs, FUN = function(x) {
  if(is.na(x)){"NA"}
  else{
    match = which(names(pphgncid) %in% x)
    pphgncid[[match]]
  }
})
#I need hgncidmatch as a character vector of length one so that I may create a data.frame/matrix
hgncidmatch <- lapply(hgncidmatch, function(x) paste(x,collapse = ";"))
hgncidmatch <- as.matrix(hgncidmatch)#as data frame doesn't work becuase can't have duplicate row names.
#add to the dataframe
multExpanded1_withDE$pphgncid <- as.character(hgncidmatch[,1])

######################

#HGNCsymbol
############ 

#Identify leading protein matches from hgncids list to add to DF
hgncSM <- sapply(as.character(multExpanded1_withDE$Leading.proteins), FUN = function(x) {
  match = which(names(hgncsymbol) %in% x)
  hgncsymbol[[match]]
})
#I need hgncidmatch as a character vector of length one so that I may create a data.frame/matrix
hgncSM <- lapply(hgncSM, function(x) paste(x,collapse = ";"))
hgncSM <- as.matrix(hgncSM)#as data frame doesn't work becuase can't have duplicate row names.
#add to the dataframe
multExpanded1_withDE$hgncsymbol <- as.character(hgncSM[,1])


#Identify protein matches from ppGOIDs list

ppMajorityProteinIDs <- as.character(multExpanded1_withDE$ppMajorityProteinIDs)
index <- which(ppMajorityProteinIDs=="")
ppMajorityProteinIDs[index] <- NA #some very strange behavior in recognizing ' "" '. had to do this first
hgncSM <- sapply(ppMajorityProteinIDs, FUN = function(x) {
  if(is.na(x)){"NA"}
  else{
    match = which(names(pphgncsymbol) %in% x)
    pphgncsymbol[[match]]
  }
})
#I need hgncidmatch as a character vector of length one so that I may create a data.frame/matrix
hgncSM <- lapply(hgncSM, function(x) paste(x,collapse = ";"))
hgncSM <- as.matrix(hgncSM)#as data frame doesn't work becuase can't have duplicate row names.
#add to the dataframe
multExpanded1_withDE$pphgncsymbol <- as.character(hgncSM[,1])

######################


#description
############ 

#Identify leading protein matches from hgncids list to add to DF
dMatch <- sapply(as.character(multExpanded1_withDE$Leading.proteins), FUN = function(x) {
  match = which(names(description) %in% x)
  description[[match]]
})
#I need hgncidmatch as a character vector of length one so that I may create a data.frame/matrix
dMatch <- lapply(dMatch, function(x) paste(x,collapse = ";"))
dMatch <- as.matrix(dMatch)#as data frame doesn't work becuase can't have duplicate row names.
#add to the dataframe
multExpanded1_withDE$description <- as.character(dMatch[,1])


#Identify protein matches from ppGOIDs list
ppMajorityProteinIDs <- as.character(multExpanded1_withDE$ppMajorityProteinIDs)
index <- which(ppMajorityProteinIDs=="")
ppMajorityProteinIDs[index] <- NA #some very strange behavior in recognizing ' "" '. had to do this first
dMatch <- sapply(ppMajorityProteinIDs, FUN = function(x) {
  if(is.na(x)){"NA"}
  else{
    match = which(names(ppdescription) %in% x)
    ppdescription[[match]]
  }
})
#I need hgncidmatch as a character vector of length one so that I may create a data.frame/matrix
dMatch <- lapply(dMatch, function(x) paste(x,collapse = ";"))
dMatch <- as.matrix(dMatch)#as data frame doesn't work becuase can't have duplicate row names.
#add to the dataframe
multExpanded1_withDE$ppdescription <- as.character(dMatch[,1])

######################


#entrezid
############ 

#Identify leading protein matches from hgncids list to add to DF
entmatch <- sapply(as.character(multExpanded1_withDE$Leading.proteins), FUN = function(x) {
  match = which(names(entrezgene) %in% x)
  entrezgene[[match]]
})
#I need hgncidmatch as a character vector of length one so that I may create a data.frame/matrix
entmatch <- lapply(entmatch, function(x) paste(x,collapse = ";"))
entmatch <- as.matrix(entmatch)#as data frame doesn't work becuase can't have duplicate row names.
#add to the dataframe
multExpanded1_withDE$entrezid <- as.character(entmatch[,1])


#Identify protein matches from ppGOIDs list
ppMajorityProteinIDs <- as.character(multExpanded1_withDE$ppMajorityProteinIDs)
index <- which(ppMajorityProteinIDs=="")
ppMajorityProteinIDs[index] <- NA #some very strange behavior in recognizing ' "" '. had to do this first
entmatch <- sapply(ppMajorityProteinIDs, FUN = function(x) {
  if(is.na(x)){"NA"}
  else{
    match = which(names(ppentrezgene) %in% x)
    ppentrezgene[[match]]
  }
})
#I need hgncidmatch as a character vector of length one so that I may create a data.frame/matrix
entmatch <- lapply(entmatch, function(x) paste(x,collapse = ";"))
entmatch <- as.matrix(entmatch)#as data frame doesn't work becuase can't have duplicate row names.
#add to the dataframe
multExpanded1_withDE$ppentrezid <- as.character(entmatch[,1])

######################



#reactomeid
############ 

#Identify leading protein matches from hgncids list to add to DF
reactmatch <- sapply(as.character(multExpanded1_withDE$Leading.proteins), FUN = function(x) {
  match = which(names(ReactIDs) %in% x)
  ReactIDs[[match]]
})
#I need hgncidmatch as a character vector of length one so that I may create a data.frame/matrix
reactmatch <- lapply(reactmatch, function(x) paste(x,collapse = ";"))
reactmatch <- as.matrix(reactmatch)#as data frame doesn't work becuase can't have duplicate row names.
#add to the dataframe
multExpanded1_withDE$ReactIDs <- as.character(reactmatch[,1])


#Identify protein matches from ppGOIDs list
ppMajorityProteinIDs <- as.character(multExpanded1_withDE$ppMajorityProteinIDs)
index <- which(ppMajorityProteinIDs=="")
ppMajorityProteinIDs[index] <- NA #some very strange behavior in recognizing ' "" '. had to do this first
reactmatch <- sapply(ppMajorityProteinIDs, FUN = function(x) {
  if(is.na(x)){"NA"}
  else{
    match = which(names(ppReactIDs) %in% x)
    ppReactIDs[[match]]
  }
})
#I need hgncidmatch as a character vector of length one so that I may create a data.frame/matrix
reactmatch <- lapply(reactmatch, function(x) paste(x,collapse = ";"))
reactmatch <- as.matrix(reactmatch)#as data frame doesn't work becuase can't have duplicate row names.
#add to the dataframe
multExpanded1_withDE$ppReactIDs <- as.character(reactmatch[,1])






#######################
# #Identify protein matches from EntrezIDs list
# Entmatch <- sapply(as.character(multExpanded1_withDE$Proteins), FUN = function(x) {
#   match = which(names(EntrezIDs) %in% x)
#   EntrezIDs[[match]]
# })
# #I need Entmatch as a character vector of length one so that I may create a data.frame/matrix
# Entmatch <- lapply(Entmatch, function(x) paste(x,collapse = ";"))
# Entmatch <- as.matrix(Entmatch)#as data frame doesn't work becuase can't have duplicate row names.
# #add to the dataframe
# multExpanded1_withDE$EntrezIDs <- as.character(Entmatch[,1])
# 
# 
# 
# #Identify protein matches from ppEntrezIDs list
# Entmatch <- sapply(as.character(multExpanded1_withDE$ProtPrepIds), FUN = function(x) {
#   if(is.na(x)){"NA"}
#   else{
#     match = which(names(ppEntrezIDs) %in% x)
#     ppEntrezIDs[[match]]
#   }
# })
# #I need Entmatch as a character vector of length one so that I may create a data.frame/matrix
# Entmatch <- lapply(Entmatch, function(x) paste(x,collapse = ";"))
# Entmatch <- as.matrix(Entmatch)#as data frame doesn't work becuase can't have duplicate row names.
# #add to the dataframe
# multExpanded1_withDE$ppEntrezIDs <- as.character(Entmatch[,1])
# 
# 
# 
# 
# #Identify protein matches from ReactIDs list
# Rmatch <- sapply(as.character(multExpanded1_withDE$Proteins), FUN = function(x) {
#   match = which(names(ReactIDs) %in% x)
#   ReactIDs[[match]]
# })
# #I need Rmatch as a character vector of length one so that I may create a data.frame/matrix
# Rmatch <- lapply(Rmatch, function(x) paste(x,collapse = ";"))
# Rmatch <- as.matrix(Rmatch)#as data frame doesn't work becuase can't have duplicate row names.
# #add to the dataframe
# multExpanded1_withDE$ReactIDs <- as.character(Rmatch[,1])
# 
# 
# 
# #Identify protein matches from ReactIDs list
# Rmatch <- sapply(as.character(multExpanded1_withDE$ProtPrepIds), FUN = function(x) {
#   if(is.na(x)){"NA"}
#   else{
#     match = which(names(ppReactIDs) %in% x)
#     ppReactIDs[[match]]
#   }
# })
# #I need Rmatch as a character vector of length one so that I may create a data.frame/matrix
# Rmatch <- lapply(Rmatch, function(x) paste(x,collapse = ";"))
# Rmatch <- as.matrix(Rmatch)#as data frame doesn't work becuase can't have duplicate row names.
# #add to the dataframe
# multExpanded1_withDE$ppReactIDs <- as.character(Rmatch[,1])

return(multExpanded1_withDE)
}







