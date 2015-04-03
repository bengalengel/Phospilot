AddAnnotation <- function(multExpanded1_withDE){
#This funciotn adds the annotations from GO, reactome, phosphositeplus, and corum? ME DF passed with DiffPhos analysis performed on confounded and non-confounded datasets. 

#curated kinase annotation and known modifications from phosphositeplus? can't get access so downloaded manually
# #files downloaded on 3/31/15
# PSPdwnload <- date()

#assign to multExpanded1 (proteins - that is the majority proteins from each of the protein groups assigned to a phosphopeptide. Perhaps the number of isoforms has to do with the lack of tie-breaker approach used at the protein group level. That is, simply the highest ranked protein is returned.) GO terms using R annotation packages


# Here I attempt to create my own organismDb object that combines the OrgDb and OntDb objects into one object 

#source("http://bioconductor.org/biocLite.R")
#biocLite("OrganismDbi")
#biocLite("Homo.sapiens")
library(Homo.sapiens)
library(GO.db)



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
#idtest <- head(keys(Homo.sapiens, keytype="GOID"),25)
#write.csv(file = "idtest.csv", idtest)

#load the GO.db and reactome.db packages
library(GO.db)
columns(GO.db)#requires GOID to map to ontology term

#retrieve the GOIDs from the 'proteins' column of the ME DF for the confounded data and the "ProtPrep Majority Protein IDs" from the protein normalized dataset. Should the annotation be performed for each 'ontology' separately?

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

#for each unique set of proteins assigned to a site retrieve GOIDs.
uniIDs <- c()
tmp <- c()
proteins <- unique(multExpanded1_withDE$Proteins)
proteins <- as.character(proteins)
GOIDs <- vector(mode = 'list', length = length(unique(multExpanded1_withDE$Proteins)))#vector function flexible for preallocation
for(i in 1:length(proteins)){
  #unique protein groups require annotation assignment
  names(GOIDs)[i] <- proteins[i]
  uniIDs <- strsplit(proteins[i], ";")
  uniIDs <- as.character(unlist(uniIDs))
  uniIDs <- substr(uniIDs,1,6)#de-isoform becuase if not it throws an error
  pos_err <- tryCatch(select(Homo.sapiens, keys=uniIDs, columns="GOID", keytype="UNIPROT"),error=function(e) e)
  if(!inherits(pos_err, "error")){
    tmp <- select(Homo.sapiens, keys=uniIDs, columns="GOID", keytype="UNIPROT")#retrieve GOIDs I should use the list form!
    GOIDs[[i]] <- as.character(tmp$GOID)
  }
}

#takes way too long. Trying parallelization with 'foreach' 'iterators' and 'doParallel'
library(iterators)
library(foreach)
# install.packages("doParallel")
library(doParallel)

proteins <- unique(multExpanded1_withDE$Proteins)
proteins <- as.character(proteins)
uniIDs <- c()
tmp <- c()
GOIDs <- vector(mode = 'list', length = length(unique(multExpanded1_withDE$Proteins)))#vector function flexible for preallocation
names(GOIDs) <- proteins
#using 7 cores
cl <- makeCluster(7)#I have 8 cores but had a crash when using all 8
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
# user  system elapsed 
# 7.30    0.67 5016.08 

stopCluster(cl)
getDoParName()




#load the reactome package
library(reactome.db)
columns(reactome.db)#requires entrez id to map to reactomeid,pathname, pathid (reactomeid vs pathid? pathid shorter)



#get the Enrezids from a set of uniprot ids
Entids <- select(Homo.sapiens, keys = ids, keytype = "UNIPROT", columns = "ENTREZID")

##note that sometimes the same entrezid is returned for different uniprot ids. I need to correct for this somehow. Can use annotation status on uniprot website.


#retrieve the reactome ids
reactids <- select(reactome.db, keys = Entids$ENTREZID, keytype = "ENTREZID", columns = "REACTOMEID")
reactids <- select(reactome.db, keys = Entids$ENTREZID, keytype = "ENTREZID", columns = c("PATHNAME","PATHID"))

#must perform enrichment using only the unique sites per protein and not the multiplicities. Also need to subset background to ensure that sites mapped to uniprot ids without annotation in a given database are not counted.





## ----echo=FALSE----------------------------------------------------------
suppressPackageStartupMessages(library(org.Hs.eg.db))

#load the human gene centric OrgDb package.
## ------------------------------------------------------------------------
library(org.Hs.eg.db)

#I can create my own OrganismDb package or use an existing one with diverse annotations

#annotationhub packages allow me to get access to a LARGE variety of external databases.

library(AnnotationHub)

#I have to create my own ah object
ah = AnnotationHub()#must be online. This is pretty powerful


#biomaRt exposes a huge family of online annotation resources called marts. Load the backage and decide which 'mart' I want to use. useMart() selects the annoation mart.

library("biomaRt")
listMarts()#need the interwebs for this. Data is current and should work from this locally once it is downloaded so that things are reproducible.

#anyway after I choose a mart I need to decide on the dataset. use the listDatasets(martobject) function
#use the useMart() function again to decide on a dataset.

#using a combination of filters and attributes I can retrieve stuff using 'getBM'

#I NEED TO PAY ATTENTION TO THE GENOME BUILD WHEN RETRIEVING ANNOTATIONS!!

#BSgenome objects exist for extracting sequence information
library(AnnotationDbi)
help(package="AnnotationDbi")

!!!!!!!!!!!!!!!!!!!!!!
#variantannotation package to assess if something is within a coding region or what. I should produce a figure using this package. 

# OrganismDb packages are named for the species they represent (such as the Homo.sapiens package). These packages contain references to other key annotations packages and can thus represent all the underlying data as if it were coming from one place. So for example, the Homo.sapiens package can allow you to retrieve data about the ranges of a genes transcripts at the same time that you extract it's gene name because it represents both a the transcriptome and the relevant org package for Homo sapiens. These can be generated using functions in the OrganismDbi package if you have specific packages that you want to link together.












