AddAnnotation <- function(multExpanded1_withDE){
  #This function adds the annotations from GO, reactome, entrez and HGNC, PFAM, iupred, rapid, biogrid, and phosphosite (at least).
  #Get cozy I am not a computer scientist


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
require(plyr)

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
################## Adding annotations (biogrid, rapid, iupred, phosphosite)



##BIOGRID----
# BIOGRID REST service for individual queries - http://wiki.thebiogrid.org/doku.php/biogridrest 
# access key - 87836af89316ba4f59bcb1e8d81df273

#get human interactions from the biogrid server.
if(!file.exists("./BIOGRID")){
  dir.create("./BIOGRID")
}

if(!file.exists("./BIOGRID/BIOGRID-ORGANISM-3.4.127.tab2.zip")){
url <- 'http://thebiogrid.org/downloads/archives/Release%20Archive/BIOGRID-3.4.127/BIOGRID-ORGANISM-3.4.127.tab2.zip'
destfile <- "./BIOGRID/BIOGRID-ORGANISM-3.4.127.tab2.zip"
download.file(url, destfile)
date_downld <- date()
write.csv(date_downld, file = "./BIOGRID/BIOGRID_download_date.csv")
}

#unzip and load the human interactions
if(!file.exists("./BIOGRID/BIOGRID-ORGANISM-Homo_sapiens-3.4.127.tab2.txt")){
human <- unzip("./BIOGRID/BIOGRID-ORGANISM-3.4.127.tab2.zip", list = TRUE)[23,1]#this is the human file
unzip("./BIOGRID/BIOGRID-ORGANISM-3.4.127.tab2.zip", files = human, exdir = "./BIOGRID")
}
Biogrid <- read.table("./BIOGRID/BIOGRID-ORGANISM-Homo_sapiens-3.4.127.tab2.txt", sep = "\t", header = T, stringsAsFactors = F, quote = "", comment.char = "")

#subset to physical interactions within human. keep only the entrez gene ids for both interactants
Biogrid <- Biogrid[Biogrid$Organism.Interactor.A == 9606 & Biogrid$Organism.Interactor.B == 9606,]
Biogrid <- Biogrid[Biogrid$Organism.Interactor.A == 9606 & Biogrid$Organism.Interactor.B == 9606,]
Biogrid <- Biogrid[Biogrid$Experimental.System.Type == "physical",]
Biogrid <- Biogrid[,c("Entrez.Gene.Interactor.A", "Entrez.Gene.Interactor.B")]

###RAPID protein level disorder data frame ----
Rapid <- read.table("./Disorder/Rapid/results_fullnames.csv", sep = ",", header = T, stringsAsFactors = F, quote = "")

#Iupred amino acid level disorder 
source("./Disorder/iupredProcessing.R")#creates list of dataframes

###PHOSPHOSITE positions of modifications.
#curated kinase annotation and known modifications from phosphositeplus downloaded manually
# #files downloaded anew on 9/8/15. These are the august files. unarchived using git/bash/bin/tar. each file is read as a .gz file
# write.csv(downld.date, file = "./PSP/PSP_download_date.csv")

#read all 'dataset' files into memory as a list. Includes all mods and curated K-S data
file.names <- grep( "dataset", list.files("./PSP/"), ignore.case = T, value = T)
file.names <- file.names[c(1,3:8)]#remove kinase substrate for now
phosphosite <- lapply(file.names, function(file){
  read.table(gzfile(file.path(paste0(getwd(),"./PSP/", file))), sep = "\t", header = T, stringsAsFactors = F, quote = "", comment.char = "", skip = 3)
})
names(phosphosite) <- gsub(".gz", "", file.names)

#remove empty assignments and subset to human
human.subset <- function(x){
  x <- x[which(x$ORGANISM == "human"),]
  x <- x[which(x$GENE != ""),]
}

phosphosite <- lapply(phosphosite, human.subset)

#first thing is to annotate novelty of the identified residues but this requires sequence conversion to uniprot identifiers.
#gene(s) level assignment of the number of modifications for now.

####INTERPROSCAN FOR PFAM BOUNDARIES -----
# This data will be used to cross ref with annotated pfam domain to assign boundaries. These boundaries will then be used to assign binary indicators of sites to 'within domain' and 'within phospho associated domain'

#Load output data from InterProScan. There should be 36 files!! For some reason cluster hangs on execution for a few of these...use "output2.tar.gz"
# lapply load into list and then do.call to make into a dataframe
if(!file.exists("./InterProScan/output")){dir.create("./InterProScan/output")}
untar("./InterProScan/output2.tar.gz", compressed = "g", exdir = "./InterProScan/output")
file.names <- list.files("./InterProScan/output", pattern = ".*.tsv")
pfam.boundary <- lapply(file.names, function(file){
  read.table(file.path(paste0(getwd(),"./InterProScan/output/", file)), sep = "\t", header = F, stringsAsFactors = F, quote = "", comment.char = "", skip = 3)
})
pfam.boundary <- do.call(rbind, pfam.boundary)
names(pfam.boundary) <- c("ENSPID", "MD5seq", "SeqLength", "Analysis", "PFamID", "Description", "Start", "Stop", "Evalue", "MatchStatus", "Date")


###################### Ensembl Derived annotations -----

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
                        filters = 'with_ccds', values = T, mart = ensembl_75)#only four values without entrez ids

ensembl_75_CCDS_pfam <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_peptide_id',"pfam"), 
                            filters = 'with_ccds', values = T, mart = ensembl_75)

#uniprot information is required when after SEQUENCE SPECIFIC information from phosphosite plus data! As of now this is on hold 
# ensembl_75_CCDS_uniprot <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_peptide_id',"uniprot_swissprot_accession"), 
#                               filters = 'with_ccds', values = T, mart = ensembl_75)


#GOIDs, (hugo gene nomenclature committee) hgncid, hgncnam and associated description

### Define function for assigning annotation to protein group (protein level) ----------

annotate <- function(proteins, ensembl_75_CCDS, ensembl_75_CCDS_EG, ensembl_75_CCDS_pfam, Biogrid, Rapid, pfam.phospho){
# This function accepts a LIST of unique protein groups and assigns annotation using passed annotation DF  
  
  
  #vector function flexible for preallocation
  GOIDs <- vector(mode = 'list', length = length(proteins))
  description <- vector(mode = 'list', length = length(proteins))
  hgncid <- vector(mode = 'list', length = length(proteins))
  hgncsymbol <- vector(mode = 'list', length = length(proteins))
  entrezgene <- vector(mode = 'list', length = length(proteins))
  ReactIDs <- vector(mode = 'list', length = length(proteins))
  pfamIDs <- vector(mode = 'list', length = length(proteins))
  pfamIDphospho <- vector(mode = 'list', length = length(proteins))
  IntCount <- vector(mode = 'list', length = length(proteins))
  Disorder <- vector(mode = 'list', length = length(proteins))
  
  cl <- makeCluster(5)#I have 8 cores but had a crash when using all 8
  registerDoParallel(cl)
  GOIDs <- foreach(i=1:length(proteins)) %dopar% {
    ENSPIDs <- strsplit(proteins[[i]], ";")
    ENSPIDs <- as.character(unlist(ENSPIDs))
    ENSPIDs <- paste(ENSPIDs, collapse = "|")
    index <- grep(ENSPIDs, ensembl_75_CCDS$ensembl_peptide_id)
    GOIDs[[i]] <- unique(ensembl_75_CCDS$go_id[index])
  }
  stopCluster(cl)
  
  cl <- makeCluster(5)
  registerDoParallel(cl)
  description <- foreach(i=1:length(proteins)) %dopar% {
    ENSPIDs <- strsplit(proteins[[i]], ";")
    ENSPIDs <- as.character(unlist(ENSPIDs))
    ENSPIDs <- paste(ENSPIDs, collapse = "|")
    index <- grep(ENSPIDs, ensembl_75_CCDS$ensembl_peptide_id)
    description[[i]] <- unique(ensembl_75_CCDS$description[index])
  }
  stopCluster(cl)
  
  cl <- makeCluster(5)
  registerDoParallel(cl)
  hgncid <- foreach(i=1:length(proteins)) %dopar% {
    ENSPIDs <- strsplit(proteins[[i]], ";")
    ENSPIDs <- as.character(unlist(ENSPIDs))
    ENSPIDs <- paste(ENSPIDs, collapse = "|")
    index <- grep(ENSPIDs, ensembl_75_CCDS$ensembl_peptide_id)
    hgncid[[i]] <- unique(ensembl_75_CCDS$hgnc_id[index])
  }
  stopCluster(cl)
  
  cl <- makeCluster(5)
  registerDoParallel(cl)
  hgncsymbol <- foreach(i=1:length(proteins)) %dopar% {
    ENSPIDs <- strsplit(proteins[[i]], ";")
    ENSPIDs <- as.character(unlist(ENSPIDs))
    ENSPIDs <- paste(ENSPIDs, collapse = "|")
    index <- grep(ENSPIDs, ensembl_75_CCDS$ensembl_peptide_id)
    hgncsymbol[[i]] <- unique(ensembl_75_CCDS$hgnc_symbol[index])
  }
  stopCluster(cl)
  
  cl <- makeCluster(5)
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
  cl <- makeCluster(5)
  registerDoParallel(cl)
  ReactIDs <- foreach(i=seq_along(entrezgene), .packages = "reactome.db") %dopar% {
    #define if hgnc is valid keytype for reactomedb
    pos_err <- tryCatch(select(reactome.db, keys=as.character(entrezgene[[i]]), columns="REACTOMEID", keytype="ENTREZID"),error=function(e) e)
    if(!inherits(pos_err, "error")){
      tmp <- select(reactome.db, keys=as.character(entrezgene[[i]]), columns="REACTOMEID", keytype="ENTREZID")
      ReactIDs[[i]] <- unique(as.character(tmp$REACTOMEID))
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
  
  #assign max interaction count among all genes within a protein group using the entrez gene id
  cl <- makeCluster(5)
  registerDoParallel(cl)
  IntCount <- foreach(i=seq_along(entrezgene)) %dopar% {
    #for each entrezgene entry (note multiple entries for single gene) find it in either of the two biogrid columns
    hits <- Biogrid[Biogrid[,1] %in% entrezgene[[i]] | Biogrid[,2] %in% entrezgene[[i]],]
    #return the number of unique entries not equal to the search term
    hits <- c(hits, recursive = T)
    hits <- unique(hits)
    IntCount[[i]] <- sum(!hits %in% entrezgene[[i]])
  }
  stopCluster(cl)
  
  #assign % disorder to each protein using the 'RAPID' dataframe
  cl <- makeCluster(5)
  registerDoParallel(cl)
  Disorder <- foreach(i=seq_along(proteins)) %dopar% {
    ENSPIDs <- strsplit(proteins[[i]], ";")
    ENSPIDs <- as.character(unlist(ENSPIDs))
    ENSPIDs <- gsub("[A-Z]", "", ENSPIDs)
    ENSPIDs <- gsub("(?<![0-9])0+", "", ENSPIDs, perl = TRUE)#rapid removes leading zeros. This is annoying.
    ENSPIDs <- paste(ENSPIDs, collapse = "|")
    index <- grep(ENSPIDs, Rapid$Prot..ID)
    Disorder[[i]] <- Rapid$Disorder.Content..[index]
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
  names(IntCount) <- proteins
  names(Disorder) <- proteins
  
  annotation.results <- list(GOIDs, description, hgncid, hgncsymbol, entrezgene, ReactIDs, pfamIDs, pfamIDphospho, IntCount, Disorder)
  names(annotation.results) <-  c("GOID", "Description", "HGNCID", "HGNCSymbol", "EntrezGene", "ReactIDs", "PFamIDs", 
                                  "PFamIDPhospho", "InteractCount", "PercentDisorder")
  
  return(annotation.results)
}

###Bulk annotation add ----------

#not easy to download so manually curated from pfam as of 8-24-15
pfam.phospho <- c("PF00498", "PF01846", "PF03166", "PF10401", "PF00244", "PF00533", "PF00400", "PF00659", "PF00397",#S/T
                  "PF00017", "PF08416", "PF00168",#Y
                  "PF00782", "PF00102", "PF13350", "PF06602", "PF04273", "PF03162", "PF14566", "PF14671", "PF04179", "PF05706", #phosphatase
                  "PF00069", "PF01636",  "PF07714", "PF03109", "PF03881", "PF06293", "PF01163", "PF01633", "PF10707", "PF06176", #kinase
                  "PF02958", "PF04655", "PF10009", "PF12260", "PF16474", "PF07914", "PF14531", "PF06734", "PF05445", "PF07387") #kinase


#confounded analysis
#'leading proteins' are used for the phospho alone annotated data. Leading proteins are the FIRST ENTRY from each matching protein group assigned to the peptide. Within a protein group ties are ignored. 'proteins' is all identified (potentially only by a modified site) proteins the peptide is associated with. 
proteins <- unique(multExpanded1_withDE$Leading.proteins)#note that all members of the protein group are not used for annotation.
proteins <- as.character(proteins)
proteins <- proteins[proteins != ""]
proteins <- as.list(proteins)

confounded.annotation <- annotate(proteins, ensembl_75_CCDS, ensembl_75_CCDS_EG, ensembl_75_CCDS_pfam, Biogrid, Rapid, pfam.phospho)


###phosprep annotation
#Using protein ids containing the phosphopeptide from the protein group with the most razor and unique peptides
#note that this caN be used for confounded enrichment analyses as well
proteins <- unique(multExpanded1_withDE$PhosPrepMajorityProteinIDs)
proteins <- as.character(proteins)
proteins <- proteins[proteins != ""]#some phosphosites are not assigned to a quantified protein
proteins <- as.list(proteins)

PhosPrep.annotation <- annotate(proteins, ensembl_75_CCDS, ensembl_75_CCDS_EG, ensembl_75_CCDS_pfam, Biogrid, Rapid, pfam.phospho)

###gelprep annotation
#Using protein ids containing the phosphopeptide from the protein group with the most razor and unique peptides
proteins <- unique(multExpanded1_withDE$ppMajorityProteinIDs)
proteins <- as.character(proteins)
proteins <- proteins[proteins != ""]#some phosphosites are not assigned to a quantified protein
proteins <- as.list(proteins)

GelPrep.annotation <- annotate(proteins, ensembl_75_CCDS, ensembl_75_CCDS_EG, ensembl_75_CCDS_pfam, Biogrid, Rapid, pfam.phospho)


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


###### Phosphosite protein level modification counts using hgnc symbol to match ----
#max on gene level is returned. For paralogous site matches the idea is that the dominant paralog will have more annotations. For example laminin1/2 or Abl1/2. Certainly a caveat but better than throwing them away.

#a function to collapse each data frame and return a data frame with two columns, the hgncname and the number of modifications/gene
count.mod <- function(x) {c(n.sites=dim(x)[1])}
phosphosite.counts <- lapply(phosphosite, function(x) ddply(x, c("GENE"), count.mod))

#merge them all and count the total number of mods
# test <- join_all(phosphosite.counts, by = "GENE", type = 'full', match = 'all'). I don't know why this doesn't work...

merged.counts <- phosphosite.counts[[1]]
for (i in 2:length(phosphosite.counts)){
  merged.counts <- merge(merged.counts, phosphosite.counts[[i]], by = "GENE", all = TRUE)
  names(merged.counts)[c(i, i + 1)] <- names(phosphosite.counts)[c(i-1, i)]
}
names(merged.counts) <- gsub("site_dataset", "site.count", names(merged.counts))
#rowsums
merged.counts$total.mod.count <- rowSums(merged.counts[,2:8], na.rm = T)
names(merged.counts)[1] <- "HGNCSymbol"

### Add this data to the MEDF
# test <- merge(multExpanded1_withDE, merged.counts, by.x = "confoundedHGNCSymbol", by.y = "HGNCSymbol", all.x = TRUE)#this will not work with multiple hgnc ids..

# a foreach solution works nicely and can be run in parallel. -Inf entries signify nomatch in database

#confounded
HGNCsymbol <- as.character(multExpanded1_withDE$confoundedHGNCSymbol)
cl <- makeCluster(5)
registerDoParallel(cl)
confounded.mod.counts <- foreach(i = 1:length(HGNCsymbol), .combine = "rbind") %dopar% {
   hgnc.symbol <- strsplit(HGNCsymbol[i], ";")
   hgnc.symbol <- as.character(unlist(hgnc.symbol))
   matches <- merged.counts[merged.counts$HGNCSymbol %in% hgnc.symbol,]
   matches <- matches[,2:9]
   if(length(matches) > 1){
     matches <- apply(matches, 2, max)
     matches <- as.data.frame(t(matches))
   }
}
#PhosPrep
HGNCsymbol <- as.character(multExpanded1_withDE$PhosPrepHGNCSymbol)
phosprep.mod.counts <- foreach(i = 1:length(HGNCsymbol), .combine = "rbind") %dopar% {
  hgnc.symbol <- strsplit(HGNCsymbol[i], ";")
  hgnc.symbol <- as.character(unlist(hgnc.symbol))
  matches <- merged.counts[merged.counts$HGNCSymbol %in% hgnc.symbol,]
  matches <- matches[,2:9]
  if(length(matches) > 1){
    matches <- apply(matches, 2, max)
    matches <- as.data.frame(t(matches))
  }
}
#GelPrep
HGNCsymbol <- as.character(multExpanded1_withDE$GelPrepHGNCSymbol)
gelprep.mod.counts <- foreach(i = 1:length(HGNCsymbol), .combine = "rbind") %dopar% {
  hgnc.symbol <- strsplit(HGNCsymbol[i], ";")
  hgnc.symbol <- as.character(unlist(hgnc.symbol))
  matches <- merged.counts[merged.counts$HGNCSymbol %in% hgnc.symbol,]
  matches <- matches[,2:9]
  if(length(matches) > 1){
    matches <- apply(matches, 2, max)
    matches <- as.data.frame(t(matches))
  }
}
stopCluster(cl)

#fix names and append to MEDF
names(confounded.mod.counts) <- paste(names(confounded.mod.counts), ".confounded", sep = "")
names(phosprep.mod.counts) <- paste(names(phosprep.mod.counts), ".PhosPrep", sep = "")
names(gelprep.mod.counts) <- paste(names(gelprep.mod.counts), ".GelPrep", sep = "")


multExpanded1_withDE <- cbind(multExpanded1_withDE, confounded.mod.counts, phosprep.mod.counts, gelprep.mod.counts)


### Amino acid position specific annotations ----
# local disorder, within pfam domains, and phosphosite level annotations

#first identify positions for 'leading proteins' in the confounded dataset
#function retrieve matching peptide location
calc.lead.prot.pos <- function(leading.prot, proteins, positions){
  leading.prot <- strsplit(as.character(leading.prot), ";")
  leading.prot <- as.character(unlist(leading.prot))
  proteins <- strsplit(as.character(proteins), ";")
  proteins <- as.character(unlist(proteins))
  index <- proteins %in% leading.prot
  #subset and return the positions vector. ensure single element passed back
  positions <- strsplit(as.character(positions), ";")
  positions <- as.character(unlist(positions))
  positions <- positions[index]
  positions <- paste(positions, collapse = ";")
}

multExpanded1_withDE$Leading.proteins.position <- mapply(calc.lead.prot.pos, multExpanded1_withDE$Leading.proteins, multExpanded1_withDE$Proteins, multExpanded1_withDE$Positions.within.proteins)



###pfam domain annotation ----
#use the pfam.boundary file created using interproscan program that calls HMMER3/pfam. Was searched on the proteome.

#binary assigment of these positions to in/out domain. If any of the sites in a protein group are in a domain it is considered 'within domain'. For those within a domain a flag yes/no is given regarding phosphorylation relevant. A foreach solution: nas when no protein group assigned to site or no domains in protein

# Confounded data
proteins <- as.character(multExpanded1_withDE$Leading.proteins)
positions <- as.character(multExpanded1_withDE$Leading.proteins.position)

cl <- makeCluster(5)
registerDoParallel(cl)
confounded.domain.boundary <- foreach(i = 1:length(proteins), .combine = "rbind") %dopar% {
  protein.group <- strsplit(proteins[i], ";")
  protein.group <- as.character(unlist(protein.group))
  protein.group <- protein.group[!grepl("REV", protein.group)]#these are omitted when finding the phosphoposition
  #Entries in pfam.boundary with no alpha or leading zeros. Altered for matching
  protein.group <- gsub("[A-Z]", "", protein.group)
  protein.group <- gsub("(?<![0-9])0+", "", protein.group, perl = TRUE)
  sites <- strsplit(as.character(positions[i]), ";")
  sites <- as.character(unlist(sites))
  #if the protein contains domains, ascertain if site is within and, if so if that domain is associated with phoshporylation
  if(length(protein.group) > 0){#not all sites have been matched to a protein group
    #preallocate in.domain and phospho.relevant
    in.domain <- vector(mode = 'logical', length = length(protein.group))
    phospho.relevant <- vector(mode = 'logical', length = length(protein.group))
    for(j in seq_along(protein.group)){
      domains <- pfam.boundary[pfam.boundary$ENSPID == protein.group[j], c(5,7,8)]
      if(dim(domains)[1] != 0){
        #apply logical test of if protein/site combo is in ANY domain range (>=start and <=end)
        hits <- domains[sites[j] >= domains$Start & sites[j] <= domains$Stop, "PFamID"]
        in.domain[j] <- length(hits) > 0
        phospho.relevant[j] <- any(hits %in% pfam.phospho)
      }else{
        in.domain[j] <- NA
        phospho.relevant[j] <- NA
      }
    }
    data.frame(site.in.domain = any(in.domain), phospho.relevant = any(phospho.relevant))
  }else{
    data.frame(site.in.domain = NA, phospho.relevant = NA)
  }
}


# PhosPrep data
proteins <- as.character(multExpanded1_withDE$PhosPrepMajorityProteinIDs)
positions <- as.character(multExpanded1_withDE$PhosPrepPositionInProteinsTruncated)

PhosPrep.domain.boundary <- foreach(i = 1:length(proteins), .combine = "rbind") %dopar% {
  protein.group <- strsplit(proteins[i], ";")
  protein.group <- as.character(unlist(protein.group))
  protein.group <- protein.group[!grepl("REV", protein.group)]#these are omitted when finding the phosphoposition
  #Entries in pfam.boundary with no alpha or leading zeros. Altered for matching
  protein.group <- gsub("[A-Z]", "", protein.group)
  protein.group <- gsub("(?<![0-9])0+", "", protein.group, perl = TRUE)
  sites <- strsplit(as.character(positions[i]), ";")
  sites <- as.character(unlist(sites))
  #if the protein contains domains, ascertain if site is within and, if so if that domain is associated with phoshporylation
  if(length(protein.group) > 0){#not all sites have been matched to a protein group
    #preallocate in.domain and phospho.relevant
    in.domain <- vector(mode = 'logical', length = length(protein.group))
    phospho.relevant <- vector(mode = 'logical', length = length(protein.group))
    for(j in seq_along(protein.group)){
      domains <- pfam.boundary[pfam.boundary$ENSPID == protein.group[j], c(5,7,8)]
      if(dim(domains)[1] != 0){
        #apply logical test of if protein/site combo is in ANY domain range (>=start and <=end)
        hits <- domains[sites[j] >= domains$Start & sites[j] <= domains$Stop, "PFamID"]
        in.domain[j] <- length(hits) > 0
        phospho.relevant[j] <- any(hits %in% pfam.phospho)
      }else{
        in.domain[j] <- NA
        phospho.relevant[j] <- NA
      }
    }
    data.frame(site.in.domain = any(in.domain), phospho.relevant = any(phospho.relevant))
  }else{
    data.frame(site.in.domain = NA, phospho.relevant = NA)
  }
}

        
# GelPrep data
proteins <- as.character(multExpanded1_withDE$ppMajorityProteinIDs)
positions <- as.character(multExpanded1_withDE$ppPositionInProteins)

GelPrep.domain.boundary <- foreach(i = 1:length(proteins), .combine = "rbind") %dopar% {
  protein.group <- strsplit(proteins[i], ";")
  protein.group <- as.character(unlist(protein.group))
  protein.group <- protein.group[!grepl("REV", protein.group)]#these are omitted when finding the phosphoposition
  #Entries in pfam.boundary with no alpha or leading zeros. Altered for matching
  protein.group <- gsub("[A-Z]", "", protein.group)
  protein.group <- gsub("(?<![0-9])0+", "", protein.group, perl = TRUE)
  sites <- strsplit(as.character(positions[i]), ";")
  sites <- as.character(unlist(sites))
  #if the protein contains domains, ascertain if site is within and, if so if that domain is associated with phoshporylation
  if(length(protein.group) > 0){#not all sites have been matched to a protein group
    #preallocate in.domain and phospho.relevant
    in.domain <- vector(mode = 'logical', length = length(protein.group))
    phospho.relevant <- vector(mode = 'logical', length = length(protein.group))
    for(j in seq_along(protein.group)){
      domains <- pfam.boundary[pfam.boundary$ENSPID == protein.group[j], c(5,7,8)]
      if(dim(domains)[1] != 0){
        #apply logical test of if protein/site combo is in ANY domain range (>=start and <=end)
        hits <- domains[sites[j] >= domains$Start & sites[j] <= domains$Stop, "PFamID"]
        in.domain[j] <- length(hits) > 0
        phospho.relevant[j] <- any(hits %in% pfam.phospho)
      }else{
        in.domain[j] <- NA
        phospho.relevant[j] <- NA
      }
    }
    data.frame(site.in.domain = any(in.domain), phospho.relevant = any(phospho.relevant))
  }else{
    data.frame(site.in.domain = NA, phospho.relevant = NA)
  }
}
stopCluster(cl)

#make names distinctive then append to dataframe
names(confounded.domain.boundary) <- paste0(names(confounded.domain.boundary), ".confounded")
names(PhosPrep.domain.boundary) <- paste0(names(PhosPrep.domain.boundary), ".PhosPrep")
names(GelPrep.domain.boundary) <- paste0(names(GelPrep.domain.boundary), ".GelPrep")


#append annotation to MEDF
multExpanded1_withDE <- cbind(multExpanded1_withDE, confounded.domain.boundary, PhosPrep.domain.boundary, GelPrep.domain.boundary)



### Local disorder using iupred -----


#binary assigment of these positions to disorder/order using iupred list of dataframes. If all of the phosphorylation positions within a protein group are in a disordered region the site is considered disordered. #Note negligible difference when calling disorder when any protein group members are assigned >.5! Most sites are in disordered regions.

disorder.assignment <- function(proteins, positions){
#   for each protein position pair, assign to order/disorder
  proteins <- strsplit(as.character(proteins), ";")
  proteins <- as.character(unlist(proteins))
  proteins <- proteins[!grepl("REV", proteins)]#these are omitted when finding the phosphoposition
  positions <- strsplit(as.character(positions), ";")
  positions <- as.character(unlist(positions))
  if(length(proteins) > 0){
  tmp <- c()
  for(i in seq_along(proteins)){
    diso.pred <- Iupred[[which(names(Iupred)==proteins[i])]]
    tmp <- c(tmp, diso.pred[diso.pred$Position == positions[i],3] >= .5)
  }
  all(tmp==TRUE)
  } else{
    NA
  }
}

multExpanded1_withDE$Confounded.Pos.Disorder <- mapply(disorder.assignment, as.character(multExpanded1_withDE$Leading.proteins), as.character(multExpanded1_withDE$Leading.proteins.position))


#phosprep. No NAs becasuse every site is assigned to a protein, although perhaps not quantified
multExpanded1_withDE$PhosPrep.Pos.Disorder <- mapply(disorder.assignment, as.character(multExpanded1_withDE$PhosPrepMajorityProteinIDs), as.character(multExpanded1_withDE$PhosPrepPositionInProteinsTruncated))#truncated?


#gelprep. NAs because some sites do not have matching protein quant estimates and therefore were not matched
multExpanded1_withDE$GelPrep.Pos.Disorder <- mapply(disorder.assignment, as.character(multExpanded1_withDE$ppMajorityProteinIDs), as.character(multExpanded1_withDE$ppPositionInProteins))

## Phosphosite annotations. on hold ---- 

# # Is residue annotated on phosphosite? Using Hugo nomenclature ids. this is really slow. need to import phosphosite with human already subsetted away.
# is.present <- function(hgnc.symbol, positions){
#   hgnc.symbol <- strsplit(as.character(hgnc.symbol), ";")
#   hgnc.symbol <- as.character(unlist(hgnc.symbol))
#   annotated.sites <- phosphosite$Phosphorylation_site_dataset[phosphosite$Phosphorylation_site_dataset$ORGANISM == "human" &
#                                                           phosphosite$Phosphorylation_site_dataset$GENE %in% hgnc.symbol,  ]
#   annotated.sites <- annotated.sites$MOD_RSD
#   annotated.sites <- gsub("[^0-9]", "", annotated.sites)
#   positions <- strsplit(as.character(positions), ";")
#   positions <- as.character(unlist(positions))
#   any(!positions %in% annotated.sites)
#   #well if only the major isoform is annotated this will give many hits! I likely need to sync this annotation with the site.group.id...
#   #uniprot ids must be used. 
# #   Positions must match the uniprot ids. Note that for a tryptic phosphopeptide matching multiple proteins will have the same site group id.
# }
# 
# #confounded
# test <- mapply(is.present, as.character(multExpanded1_withDE$confoundedHGNCSymbol), as.character(multExpanded1_withDE$Leading.proteins.position))
# table(test)


save(multExpanded1_withDE, file = "./multExpanded_withDE_annotated.RData")
return(multExpanded1_withDE)
}

