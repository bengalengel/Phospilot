#this script produces the FASTA file used for proteome database search and downstream processing. It is ugly. Sorry future self.

#biomart practice. I want to see if I can pull the peptide fasta sequences for transcripts that meet a certain level of annotation confidence. I also want the fastas for the ENSP ids.Some ENSP transcripts have ambiguous start and stops. (5' and 3' truncations). I don't want to search a protein database with truncated proteins.
#options are CCDS and GENCODE basic. The latter may include truncated transcripts. Therefore the CCDs SET IS USED subject to a few filters.

# http://www.ensembl.org/Help/Glossary?id=500
library("biomaRt")
library("seqinr")

#archive hsapiens mart
ensembl_75 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

#what is this thing?
ensembl_75
# Object of class 'Mart':
#   Using the ENSEMBL_MART_ENSEMBL BioMart database
# Using the hsapiens_gene_ensembl dataset

#filters are the things with values that restrict what returns and attributes are what is returned. Can I return with 5' CDS annotation flags as a 'filter'? NO!
keytypes(ensembl_75)
columns(ensembl_75)
filters <- listFilters(ensembl_75)
attributes <- listAttributes(ensembl_75)


#return all ENSG/T/Ps that are part of the GENCODE basic set and the CCDS set. Then look for overlap. 

#CCDS filters
grep(pattern = "ccds", filters$name, ignore.case = T, value = T)
# [1] "with_ccds" "ccds"


#lets see how many we get and how many have an HGNC identifier etc etc.
ensembl_75_CCDS <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_peptide_id','description','hgnc_id','hgnc_symbol','gene_biotype',
                                        'transcript_biotype', 'status', 'transcript_status', 'ccds'), filters = 'with_ccds', values = T, mart = ensembl_75)

#spot check confirms all transcripts contain a cCDS identifier. 



#subset to gene and transcript status: known or novel. https://gencodegenes.wordpress.com/toolbox/
# Specifically, a Known CDS is 100% identical to a RefSeq NP or Swissprot / Uniprot entry along its length, a Novel CDS shares >60% length with a Known CDS in the same gene and uses the same initiation and termination codons, or has a corresponding known paralog or ortholog, and a Putative CDS shares <60% length with a Known CDS or has an alternative first or last coding exon.

#known or novel gene status index
kngenes_index <- grep(pattern = "known|novel", ensembl_75_CCDS$status, ignore.case = T)# this requires less writing than subsetting within the DF
ensembl_75_CCDS <- ensembl_75_CCDS[kngenes_index,]

#known or novel gene status index
kntranscript_index <- grep(pattern = "known|novel", ensembl_75_CCDS$transcript_status, ignore.case = T)# removes 1373 putative transcripts
ensembl_75_CCDS <- ensembl_75_CCDS[kntranscript_index,]

#subset to protein and transcript biotype - 'protein_coding' 
gene_biotype_index <- grep(pattern = "protein_coding", ensembl_75_CCDS$gene_biotype, ignore.case = T)#removes polymorphic pseudogenes (8 of them)
ensembl_75_CCDS <- ensembl_75_CCDS[gene_biotype_index,]

transcript_biotype_index <- grep(pattern = "protein_coding", ensembl_75_CCDS$transcript_biotype, ignore.case = T)#removes 343 NMD transcripts
ensembl_75_CCDS <- ensembl_75_CCDS[transcript_biotype_index,]


#checks
table(ensembl_75_CCDS$status)
table(ensembl_75_CCDS$transcript_status)
table(ensembl_75_CCDS$gene_biotype)
table(ensembl_75_CCDS$transcript_biotype)

#total of 35587 transcripts


####subset the ENSP fasta s.t. only ENSP ids within ensembl_75_CCDS set are retained -----

#upload FASTA from archive site;
#create FASTA directory if it doesn't already exist
if(!file.exists(file.path(getwd(), "FASTA"))){
  dir.create(file.path(file.path(getwd(), "FASTA"))
}  

#download file and save the date (first pass on 5/1/15)
if(!file.exists(file.path(getwd(), "FASTA", "Homo_sapiens.GRCh37.75.pep.all.fa.gz"))){
url <- "ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/pep/Homo_sapiens.GRCh37.75.pep.all.fa.gz"
destfile <- "FASTA/Homo_sapiens.GRCh37.75.pep.all.fa.gz"
download.file(url, destfile)
date_downld <- date()
save(date_downld, )
write.csv(date_downld, file = "./FASTA/FASTA_download_date.csv")
}

#read in the FASTA file using 'read.fasta' and 'gzfile' to uncompress the file
EnsProt <- read.fasta(file = gzfile("./FASTA/Homo_sapiens.GRCh37.75.pep.all.fa.gz"), seqtype = "AA", as.string = T)

# length(EnsProt)
# #100778
# is.SeqFastaAA(EnsProt[[1]])
# 
# names(EnsProt[1])#gives ENSP identifier
# 
# attributes(EnsProt[[1]])# each list object has an 'Annot' attribute associated with it. This has all the fasta header information
# attr(EnsProt[[1]], 'Annot')


# parse FASTA; 
ENSPids <- ensembl_75_CCDS$ensembl_peptide_id
ENSPidsFASTA <- names(EnsProt)
ENSPids_index <- ENSPidsFASTA%in%ENSPids
EnsProt <- EnsProt[ENSPids_index]#35585 close enough

#Change headers for compatibility with Andromeda output tables.  

#really need to learn how to vectorize things. pick either plyr or s/lapply
test <- getAnnot(EnsProt[[1]])
names(EnsProt[1])

#below is what I want for a FASTA header. Analagous to uniprot header scheme and will allow MQ to fill in all relevant columns
">ensembl_75|ENSPID|HGNCid description(from HGNC) OS=HomoSapiens GN=hgnc_symbol"

#create a list of character vectors to serve as fasta headers for new outputted fasta file
FastaHeaders2 = lapply(names(EnsProt), function(x){
  #extract the ensp id
  #   EnspId <- names(x)
  
  #extract annotation from the ensembl table
  annot <- ensembl_75_CCDS[ensembl_75_CCDS$ensembl_peptide_id == x,]
  HGNCid <- annot$hgnc_id
  description <- annot$description
  hgnc_symbol <- annot$hgnc_symbol
  hgnc_symbol <- paste0("GN=",hgnc_symbol)
  
  #make the header file
  header <- paste("ensembl_75", x, sep = "|")
  header <- paste(header, HGNCid, sep = "|")
  header <- paste(header, description, sep = " ")
  header <- paste(header, "OS=HomoSapiens", sep = " ")
  header <- paste(header, hgnc_symbol, sep = " ")
})

#write parsed fasta with new headers
filepath <- file.path(getwd(),"FASTA/Homo_sapiens.GRCh37.75.pep.all.parsedCCDS.fa")
write.fasta(EnsProt, names = FastaHeaders2, file = filepath, as.string=T)

#spotcheck for id matching positive.
EnsProt[1:10]





