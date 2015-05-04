AddAnnotation <- function(multExpanded1_withDE){
#This function adds the annotations from GO, reactome, phosphositeplus, and corum? ME DF passed with DiffPhos analysis already performed on confounded and non-confounded datasets. 

#curated kinase annotation and known modifications from phosphositeplus downloaded manually
# #files downloaded on 3/31/15
# PSPdwnload <- date()

#assign to multExpanded1 (proteins - that is the majority proteins from each of the protein groups assigned to a phosphopeptide. Perhaps the number of isoforms has to do with the lack of tie-breaker approach used at the protein group level. That is, simply the highest ranked protein is returned.) 
  
#GO terms and reactome added using R annotation packages


#source("http://bioconductor.org/biocLite.R")
#biocLite("OrganismDbi")
#biocLite("Homo.sapiens")
require(Homo.sapiens)
require(GO.db)
require(reactome.db)
require(iterators)
require(foreach)
require(doParallel)


#Examples
##############
# To list the kinds of things that can be retrieved, use the columns method.
columns(Homo.sapiens)
# To list the kinds of things that can be used as keys we can use the keytypes method
keytypes(Homo.sapiens)
# And to extract viable keys of a particular kind (keytype), we can use the keys method.
idtest <- head(keys(Homo.sapiens, keytype="GOID"),25)
# Since the keys method can tell us specific things that can be used as keys, here we will use it to extract a few ids to use for demonstrating the fourth method type.
ids = head(keys(Homo.sapiens, keytype="UNIPROT"),25)
# Once you have some ids that you want to look up data for, the select method allows you to map these ids as long as you use the columns argument to indicate what you need to know and the keytype argument to specify what kind of keys they are.
select(Homo.sapiens, keys=ids, columns="GOID", keytype="UNIPROT")
head(select(Homo.sapiens, keys=ids, columns="GOALL", keytype="UNIPROT"))#above is more direct
######################


#retrieve the GOIDs from the 'proteins' column of the ME DF for the confounded data and the "ProtPrep Majority Protein IDs" from the protein normalized dataset. The annotation should be performed for each 'ontology' separately unfortunately.

#preliminary question 1 - Do isoforms retrieve different GOIDs than the parent uniprot identifiers? Isoforms may be uniquely assigned to different biological processes.
##################
MEtest <- multExpanded1_withDE
MEtest$Protein <- as.character(MEtest$Protein)
uniIDS <- MEtest[nchar(MEtest$Protein) > 6, ]
uniIDs <- head(uniIDS$Protein, 25)
uniIDs <- c(uniIDs, substr(uniIDs,1,6))

isotest <- select(Homo.sapiens, keys=uniIDs, columns="GOID", keytype="UNIPROT")#hmm, no isoforms produce good maps

uniIDs2 <- uniIDs[1:25]#just isoforms
isotest2 <- select(Homo.sapiens, keys=uniIDs2, columns="GOID", keytype="UNIPROT")#isoforms in fact throw an error in uniprot mapping
uniIDs3 <- uniIDs[26:length(uniIDs)]
isotest3 <- select(Homo.sapiens, keys=uniIDs3, columns="GOID", keytype="UNIPROT")#warning message says that duplicate query results in 1:many mapping issues. however I don't see this.
#are duplicate uniprot entries mapped to duplicate GO entries?
uniIDs4 <- unique(uniIDs3)
isotest4 <- select(Homo.sapiens, keys=uniIDs4, columns="GOID", keytype="UNIPROT")

#No, this is not the case
dim(isotest3)==dim(isotest4)
#[1] TRUE TRUE
##########################
#ANSWER - isoforms do not retrieve GO IDs at all and duplicate entries return the same number of GOIDs as unique entries despite warning message.

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

# Trying parallelization with 'foreach' 'iterators' and 'doParallel'

# Adding annotations
############

#first is GOIDs
proteins <- unique(multExpanded1_withDE$Proteins)
proteins <- as.character(proteins)
uniIDs <- c()
tmp <- c()
GOIDs <- vector(mode = 'list', length = length(unique(multExpanded1_withDE$Proteins)))#vector function flexible for preallocation
names(GOIDs) <- proteins
#using 7 cores
cl <- makeCluster(5)#I have 8 cores but had a crash when using all 8
registerDoParallel(cl)
system.time(GOIDs <- foreach(i=1:length(proteins), .packages = "Homo.sapiens") %dopar% {
  uniIDs <- strsplit(proteins[i], ";")
  uniIDs <- as.character(unlist(uniIDs))
  uniIDs <- substr(uniIDs,1,6)#de-isoform becuase if not it throws an error
  pos_err <- tryCatch(select(Homo.sapiens, keys=uniIDs, columns="GOID", keytype="UNIPROT"),error=function(e) e)
  if(!inherits(pos_err, "error")){
    tmp <- select(Homo.sapiens, keys=uniIDs, columns="GOID", keytype="UNIPROT")#retrieve GOIDs I should use the list form!
    GOIDs[[i]] <- as.character(tmp$GOID)
  }
}
)
stopCluster(cl)
getDoParName()
names(GOIDs) <- proteins
# user  system elapsed 
# 7.30    0.67 5016.08 


#for the protein normalized data

# replace the name of the ProtPrep column
names(multExpanded1_withDE)[grep(names(multExpanded1_withDE), pattern = "ProtPrep")] <- "ProtPrepIds"

proteins <- unique(multExpanded1_withDE$ProtPrepIds)
proteins <- as.character(proteins)
proteins <- na.omit(proteins)
uniIDs <- c()
tmp <- c()
ppGOIDs <- vector(mode = 'list', length = length(unique(multExpanded1_withDE$ProtPrepIds)))#vector function flexible for preallocation
names(ppGOIDs) <- proteins
#using 7 cores
cl <- makeCluster(5)#I have 8 cores but had a crash when using all 8
registerDoParallel(cl)
system.time(ppGOIDs <- foreach(i=1:length(proteins), .packages = "Homo.sapiens") %dopar% {
  uniIDs <- strsplit(proteins[i], ";")
  uniIDs <- as.character(unlist(uniIDs))
  uniIDs <- substr(uniIDs,1,6)#de-isoform becuase if not it throws an error
  pos_err <- tryCatch(select(Homo.sapiens, keys=uniIDs, columns="GOID", keytype="UNIPROT"),error=function(e) e)
  if(!inherits(pos_err, "error")){
    tmp <- select(Homo.sapiens, keys=uniIDs, columns="GOID", keytype="UNIPROT")#retrieve ppGOIDs I should use the list form!
    ppGOIDs[[i]] <- as.character(tmp$GOID)
  }
}
)
stopCluster(cl)
getDoParName()

names(ppGOIDs) <- proteins



#Produce a reactome list in a similar fashion. first EntrezIDs
uniIDs <- c()
tmp <- c()
EntrezIDs <- vector(mode = 'list', length = length(unique(multExpanded1_withDE$Proteins)))#vector function flexible for
#using 7 cores
cl <- makeCluster(5)#I had a crash when using all 8
registerDoParallel(cl)
system.time(EntrezIDs <- foreach(i=1:length(proteins), .packages = "Homo.sapiens") %dopar% {
  uniIDs <- strsplit(proteins[i], ";")
  uniIDs <- as.character(unlist(uniIDs))
  uniIDs <- substr(uniIDs,1,6)#de-isoform becuase if not it throws an error
  pos_err <- tryCatch(select(Homo.sapiens, keys=uniIDs, columns="ENTREZID", keytype="UNIPROT"),error=function(e) e)
  if(!inherits(pos_err, "error")){
    tmp <- select(Homo.sapiens, keys=uniIDs, columns="ENTREZID", keytype="UNIPROT")#retrieve Entrezids I should use the list form!
    EntrezIDs[[i]] <- as.character(tmp$ENTREZID)
  }
}
)
stopCluster(cl)
names(EntrezIDs) <- proteins


# ppEntrezIDs
uniIDs <- c()
tmp <- c()
ppEntrezIDs <- vector(mode = 'list', length = length(unique(multExpanded1_withDE$ProtPrepIds)))#vector function flexible for
#using 7 cores
cl <- makeCluster(5)#I had a crash when using all 8
registerDoParallel(cl)
system.time(ppEntrezIDs <- foreach(i=1:length(proteins), .packages = "Homo.sapiens") %dopar% {
  uniIDs <- strsplit(proteins[i], ";")
  uniIDs <- as.character(unlist(uniIDs))
  uniIDs <- substr(uniIDs,1,6)#de-isoform becuase if not it throws an error
  pos_err <- tryCatch(select(Homo.sapiens, keys=uniIDs, columns="ENTREZID", keytype="UNIPROT"),error=function(e) e)
  if(!inherits(pos_err, "error")){
    tmp <- select(Homo.sapiens, keys=uniIDs, columns="ENTREZID", keytype="UNIPROT")#retrieve ppEntrezIDs I should use the list form!
    ppEntrezIDs[[i]] <- as.character(tmp$ENTREZID)
  }
}
)
stopCluster(cl)
names(ppEntrezIDs) <- proteins



#now for the reactomeIds
ReactIDs <- vector(mode = 'list', length = length(unique(multExpanded1_withDE$Proteins)))#vector function flexible for preallocation
#using 7 cores
cl <- makeCluster(5)#I had a crash when using all 8
registerDoParallel(cl)
system.time(ReactIDs <- foreach(i=seq_along(EntrezIDs), .packages = "reactome.db") %dopar% {
  pos_err2 <- tryCatch(select(reactome.db, keys=EntrezIDs[[i]], columns="REACTOMEID", keytype="ENTREZID"),error=function(e) e)
  if(!inherits(pos_err2, "error")){
    tmp <- select(reactome.db, keys=EntrezIDs[[i]], columns="REACTOMEID", keytype="ENTREZID")
    ReactIDs[[i]] <- as.character(tmp$REACTOMEID)
  }
}
)
stopCluster(cl)
names(ReactIDs) <- proteins


ppReactIDs <- vector(mode = 'list', length = length(unique(multExpanded1_withDE$ProtPrepIds)))#vector function flexible for preallocation
#using 7 cores
cl <- makeCluster(5)#I had a crash when using all 8
registerDoParallel(cl)
system.time(ppReactIDs <- foreach(i=seq_along(EntrezIDs), .packages = "reactome.db") %dopar% {
  pos_err2 <- tryCatch(select(reactome.db, keys=EntrezIDs[[i]], columns="REACTOMEID", keytype="ENTREZID"),error=function(e) e)
  if(!inherits(pos_err2, "error")){
    tmp <- select(reactome.db, keys=EntrezIDs[[i]], columns="REACTOMEID", keytype="ENTREZID")
    ppReactIDs[[i]] <- as.character(tmp$REACTOMEID)
  }
}
)
stopCluster(cl)
names(ppReactIDs) <- proteins

#################################################




#Add these annotation IDs to the original ME dataframe.

#Identify protein matches from GOIDs list
GOmatch <- sapply(as.character(multExpanded1_withDE$Proteins), FUN = function(x) {
  match = which(names(GOIDs) %in% x)
  GOIDs[[match]]
})
#I need GOmatch as a character vector of length one so that I may create a data.frame/matrix
GOmatch <- lapply(GOmatch, function(x) paste(x,collapse = ";"))
GOmatch <- as.matrix(GOmatch)#as data frame doesn't work becuase can't have duplicate row names.
#add to the dataframe
multExpanded1_withDE$GOIDs <- as.character(GOmatch[,1])


#Identify protein matches from ppGOIDs list
GOmatch <- sapply(as.character(multExpanded1_withDE$ProtPrepIds), FUN = function(x) {
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


#Identify protein matches from EntrezIDs list
Entmatch <- sapply(as.character(multExpanded1_withDE$Proteins), FUN = function(x) {
  match = which(names(EntrezIDs) %in% x)
  EntrezIDs[[match]]
})
#I need Entmatch as a character vector of length one so that I may create a data.frame/matrix
Entmatch <- lapply(Entmatch, function(x) paste(x,collapse = ";"))
Entmatch <- as.matrix(Entmatch)#as data frame doesn't work becuase can't have duplicate row names.
#add to the dataframe
multExpanded1_withDE$EntrezIDs <- as.character(Entmatch[,1])



#Identify protein matches from ppEntrezIDs list
Entmatch <- sapply(as.character(multExpanded1_withDE$ProtPrepIds), FUN = function(x) {
  if(is.na(x)){"NA"}
  else{
    match = which(names(ppEntrezIDs) %in% x)
    ppEntrezIDs[[match]]
  }
})
#I need Entmatch as a character vector of length one so that I may create a data.frame/matrix
Entmatch <- lapply(Entmatch, function(x) paste(x,collapse = ";"))
Entmatch <- as.matrix(Entmatch)#as data frame doesn't work becuase can't have duplicate row names.
#add to the dataframe
multExpanded1_withDE$ppEntrezIDs <- as.character(Entmatch[,1])




#Identify protein matches from ReactIDs list
Rmatch <- sapply(as.character(multExpanded1_withDE$Proteins), FUN = function(x) {
  match = which(names(ReactIDs) %in% x)
  ReactIDs[[match]]
})
#I need Rmatch as a character vector of length one so that I may create a data.frame/matrix
Rmatch <- lapply(Rmatch, function(x) paste(x,collapse = ";"))
Rmatch <- as.matrix(Rmatch)#as data frame doesn't work becuase can't have duplicate row names.
#add to the dataframe
multExpanded1_withDE$ReactIDs <- as.character(Rmatch[,1])



#Identify protein matches from ReactIDs list
Rmatch <- sapply(as.character(multExpanded1_withDE$ProtPrepIds), FUN = function(x) {
  if(is.na(x)){"NA"}
  else{
    match = which(names(ppReactIDs) %in% x)
    ppReactIDs[[match]]
  }
})
#I need Rmatch as a character vector of length one so that I may create a data.frame/matrix
Rmatch <- lapply(Rmatch, function(x) paste(x,collapse = ";"))
Rmatch <- as.matrix(Rmatch)#as data frame doesn't work becuase can't have duplicate row names.
#add to the dataframe
multExpanded1_withDE$ppReactIDs <- as.character(Rmatch[,1])

return(multExpanded1_withDE)
}




# Enrichment below. Will be a separate function
##########################

#First is GO. 

#get background and correct for multiple observations of the same peptides (multiplicities)
backgroundGO <- multExpanded1_withDE[multExpanded1_withDE$SubtoDE == "+",]
#remove duplicate entries of the same peptide (multiplicities). Subset by unique id. id refers to a phosphosite
backgroundGO <- backgroundGO[!duplicated(backgroundGO$id),]#the duplicated returns a logical identifying the first instance of duplication from an ordererd query. Therefore it is safe to subset using 'duplicated'
backgroundGO <- backgroundGO$GOIDs
BGGO <- sapply(backgroundGO, function(x) strsplit(x, ";"))
BGGO <- unlist(BGGO)
BGGO <- BGGO[BGGO!='NA']#to be passed

enrichedGO <- multExpanded1_withDE[multExpanded1_withDE$globalFsig == "+",]
enrichedGO <- enrichedGO[!duplicated(enrichedGO$id),]
enrichedGO <- enrichedGO$GOIDs
DEGO <- sapply(enrichedGO, function(y) strsplit(y, ";"))
DEGO <- unlist(DEGO)
DEGO <- DEGO[DEGO!='NA']#to be passed

#now for the protein normalized DE
backgroundGO <- multExpanded1_withDE[multExpanded1_withDE$SubtoDEpn =="+",]
#remove duplicate entries of the same peptide (multiplicities). Subset by unique id. id refers to a phosphosite
backgroundGO <- backgroundGO[!duplicated(backgroundGO$id),]#the duplicated returns a logical identifying the first instance of duplication from an ordererd query. Therefore it is safe to subset using 'duplicated'
backgroundGO <- backgroundGO$ppGOIDs
BGGOpn <- sapply(backgroundGO, function(x) strsplit(x, ";"))
BGGOpn <- unlist(BGGOpn)
BGGOpn <- BGGOpn[BGGOpn!='NA']#to be passed

enrichedGO <- multExpanded1_withDE[multExpanded1_withDE$globalFsigpn == "+",]
enrichedGO <- enrichedGO[!duplicated(enrichedGO$id),]
enrichedGO <- enrichedGO$GOIDs
DEGOpn <- sapply(enrichedGO, function(y) strsplit(y, ";"))
DEGOpn <- unlist(DEGOpn)
DEGOpn <- DEGOpn[DEGOpn!='NA']#to be passed



Enrich <- function(x,y){
  #This function accepts lists of kinases found in a DE subset and background subset identified by the networkin program and return a list of adjusted pvalues for categorical enrichment using a fisher's exact test.x=background and y=enriched
  require(plyr)
  BGtable <- as.matrix(table(x))#the parenthetical should be a factor vector of ids passed to this function
  #remove entries with 0
  BGtable <- BGtable[BGtable!=0,,drop=F]
  DEtable <- as.matrix(table(y))#the parenthetical should be a passed DE factor vector of networKin output
  DEtable <- as.matrix(DEtable[row.names(DEtable) %in% row.names(BGtable),])#removing zeros and all factors not present in BG data
  #note the use of the %in% statement to control for rownames in DeTable but not in background
  #subset the background table in a similar way to ensure we are making the proper comparisons
  BGtable <- as.matrix(BGtable[row.names(BGtable) %in% row.names(DEtable),])
  NotDE <- BGtable-DEtable
  facttemp <- as.factor(row.names(DEtable))
  pvals <- c()
  frequency <- c()
  ids <- c()
  #for each unique kinase in the DE1 kinases table I need to make a contingency table
  for(i in levels(facttemp)){#this would be the DE dataframe
    #make the first row of the contingency table
    if(DEtable[as.character(i),] & NotDE[as.character(i),] >= 0){#this condition should always be true
      DErow <- c(DEtable[as.character(i),],(sum(DEtable)-DEtable[as.character(i),]))
      NotDErow <- c(NotDE[as.character(i),],(sum(NotDE)-NotDE[as.character(i),]))
      contmatrix <- rbind(DErow,NotDErow)
      tmp <- fisher.test(contmatrix, alternative = "g")
      pvals <- c(pvals,tmp$p.value)
      ids <- c(ids,i)
      tmp2 <- paste(as.character(contmatrix[1,1]), as.character(colSums(contmatrix)[1]), sep = "/")
      frequency <- c(frequency, tmp2)
    }
  }
  ##multiple testing correction for pvals
  adjps <- p.adjust(pvals,method="BH")#none pass significance in this instance
  #put into a dataframe object, sort by adjusted p values
  DF <- data.frame('ID' = ids, 'Frequency' = frequency, 'pvalue' = pvals, 'adjpvalue' = adjps)
  DF <- arrange(DF, adjpvalue)#arrange
  DF <- DF[DF$adjpvalue<=.05,]
  return(DF)
}

#run the function
GODF <- Enrich(BGGO,DEGO)#from above
#assign the annotation
require(GO.db)
GOannotation <- select(GO.db, keys = as.character(GODF$ID), keytype = "GOID", columns = c("TERM", "ONTOLOGY", "DEFINITION"))
GOenrichment <- cbind(GODF,GOannotation)

#run the function
GODFpn <- Enrich(BGGOpn,DEGOpn)#from above
#assign the annotation
require(GO.db)
GOannotationpn <- select(GO.db, keys = as.character(GODFpn$ID), keytype = "GOID", columns = c("TERM", "ONTOLOGY", "DEFINITION"))
GOenrichmentpn <- cbind(GODFpn,GOannotationpn)


#the two are quite different. Note that for proteins I am using the majority protein IDs from zia's protein prep. This may be incorrect given that a protein within the maj. protein ids. group may not have the same sequence as the phosphorylated peptide...

#the above is true iff the phosphopeptide is a unique sequence relative all other sequences defined by the protein group in Zia's dataset.

#two steps for this
# 1) match using 'proteins' for phoshphopeptides to majority proteins in zia's prep. multiple matches are parsed by occam's razor.
# 2) Matching protein group from Zia's prep is pruned for proteins that don't contain the matching phosphopeptide (nonphosphorylated form)

x <- GOenrichment$GOID#114
y <- GOannotationpn$GOID#82

length(intersect(y,x))#22
length(setdiff(y,x))#60


################################################

#same steps for reactome
###################
#get background and correct for multiple observations of the same peptides (multiplicities)
backgroundReact <- multExpanded1_withDE[multExpanded1_withDE$SubtoDE=="+",]
#remove duplicate entries of the same peptide (multiplicities). Subset by unique id. id refers to a phosphosite
backgroundReact <- backgroundReact[!duplicated(backgroundReact$id),]#the duplicated returns a logical identifying the first instance of duplication from an ordererd query. Therefore it is safe to subset using 'duplicated'
backgroundReact <- backgroundReact$ReactIDs
BGRO <- sapply(backgroundReact, function(x) strsplit(x, ";"))
BGRO <- unlist(BGRO)
BGRO <- BGRO[BGRO!='NA']#to be passed

enrichedReact <- multExpanded1_withDE[multExpanded1_withDE$globalFsig == "+",]
enrichedReact <- enrichedReact[!duplicated(enrichedReact$id),]
enrichedReact <- enrichedReact$ReactIDs
DERO <- sapply(enrichedReact, function(y) strsplit(y, ";"))
DERO <- unlist(DERO)
DERO <- DERO[DERO!='NA']#to be passed

#now for the protein normalized DE
#get background and correct for multiple observations of the same peptides (multiplicities)
backgroundReact <- multExpanded1_withDE[multExpanded1_withDE$SubtoDE=="+",]
#remove duplicate entries of the same peptide (multiplicities). Subset by unique id. id refers to a phosphosite
backgroundReact <- backgroundReact[!duplicated(backgroundReact$id),]#the duplicated returns a logical identifying the first instance of duplication from an ordererd query. Therefore it is safe to subset using 'duplicated'
backgroundReact <- backgroundReact$ppReactIDs
pnBGRO <- sapply(backgroundReact, function(x) strsplit(x, ";"))
pnBGRO <- unlist(pnBGRO)
pnBGRO <- pnBGRO[pnBGRO!='NA']#to be passed

enrichedReact <- multExpanded1_withDE[multExpanded1_withDE$globalFsig == "+",]
enrichedReact <- enrichedReact[!duplicated(enrichedReact$id),]
enrichedReact <- enrichedReact$ppReactIDs
pnDERO <- sapply(enrichedReact, function(y) strsplit(y, ";"))
pnDERO <- unlist(pnDERO)
pnDERO <- pnDERO[pnDERO!='NA']#to be passed


#run the function
RODF <- Enrich(BGRO,DERO)
#assign the annotation
require(reactome.db)

#assign the annotation
ROannotation <- select(reactome.db, keys = as.character(RODF$ID), keytype = "REACTOMEID", columns = "PATHNAME")
ROenrichment <- cbind(RODF,ROannotation)

#protein normalized data
RODFpn <- Enrich(pnBGRO,pnDERO)

#assign the annotation
ROannotationpn <- select(reactome.db, keys = as.character(RODFpn$ID), keytype = "REACTOMEID", columns = "PATHNAME")
ROenrichmentpn <- cbind(RODFpn,ROannotationpn)


##############################






