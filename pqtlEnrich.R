#this script aims to assess if pQTLs are enriched in the differentially phosphorylated subset.

#load the pqtl file from Battle et al Science 2014

#requires the XLConnect package which also requires a JRE installation to work properly. This can be a slow solution apparently but will work on Unix machines
# http://www.oracle.com/technetwork/java/javase/downloads/jre8-downloads-2133155.html

install.packages("XLConnect")
library(XLConnect)
pqtl <- readWorksheetFromFile("E:/My Documents/Dropbox/Postdoc-Gilad/Zia/latest science/1260793_DatafileS1.xlsx", sheet = 4, header = T)


# retrieve enspid from pqtl ensgid ----------------------------------------

#how many of the ensgids are in the biomart ensembl75 GRCh37.p13/hg19 repository?
library("biomaRt")
library("seqinr")

#archive hsapiens mart
ensembl_75 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")


ensembl_75_pqtls <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_peptide_id','description','hgnc_id',
                                         'hgnc_symbol','gene_biotype', 'transcript_biotype', 'status', 'transcript_status', 'ccds'),
                            filters = "ensembl_gene_id", values = pqtl$ENSG, mart = ensembl_75)


#278 genes
length(unique(pqtl$ENSG)) - length(unique(ensembl_75_pqtls$ensembl_gene_id))#only one gene missing.

#from 278 pqtl genes, how many have at least one CCDs protein id and were therefore included in the database for the search?

ccds_index <- grep(pattern = "ccds", ensembl_75_pqtls$ccds, ignore.case = T)# this requires less writing than subsetting within the DF
ensembl_75_pqtlscds <- ensembl_75_pqtls[ccds_index,]
length(unique(ensembl_75_pqtlscds$ensembl_gene_id))
#275!!! Great ;)

##subset further based on gene/transcript biotype. start with 587 entries
#known or novel gene status index
kntranscript_index <- grep(pattern = "known|novel", ensembl_75_pqtlscds$transcript_status, ignore.case = T)#576 (removes putative transcripts)
ensembl_75_pqtlscds <- ensembl_75_pqtlscds[kntranscript_index,]

#subset to protein and transcript biotype - 'protein_coding' 
gene_biotype_index <- grep(pattern = "protein_coding", ensembl_75_pqtlscds$gene_biotype, ignore.case = T)#576 removes polymorphic pseudogenes
ensembl_75_pqtlscds <- ensembl_75_pqtlscds[gene_biotype_index,]

transcript_biotype_index <- grep(pattern = "protein_coding", ensembl_75_pqtlscds$transcript_biotype, ignore.case = T)#573 removes 3 NMD transcripts
ensembl_75_pqtlscds <- ensembl_75_pqtlscds[transcript_biotype_index,]

#still 275 possible hits
length(unique(ensembl_75_pqtlscds$ensembl_gene_id))



# adding pqtl presence/absence annotation to multExpanded1_withDE ---------

#match for any assigned protein from the protein group to a gene with a pqtl?
AnyMatchProtQTL <- function(queryproteins){
  #This function looks for a match between queryproteins and the proteins containing a non-synonymous snp
  #queryproteins are proteins assigned to the given phosphopeptide
  if(any(unlist(strsplit(as.character(queryproteins), ";")) %in% ensembl_75_pqtlscds$ensembl_peptide_id)){
    "+"
  }else{
    "-"
  }
}
#apply anymatch function to all 'leading proteins' and 'majority protein ids'
multExpanded1_withDE$pQTLPositive <- mapply(AnyMatchProtQTL, multExpanded1_withDE$Leading.proteins)
multExpanded1_withDE$GelPreppQTLPositive <- mapply(AnyMatchProtQTL, multExpanded1_withDE$ppMajorityProteinIDs)

#how many ENSPids from a pQTL gene per protein group are assigned to a phosphopeptide?
NumMatchesProtQTL <- function(queryproteins){
  #This function counts the matches between queryproteins and the proteins containing a non-synonymous snp
  #queryproteins are proteins assigned to the peptide
  sum(unlist(strsplit(as.character(queryproteins), ";")) %in% ensembl_75_pqtlscds$ensembl_peptide_id)
}

#apply count function to all 'leading proteins' and 'majority protein ids'
multExpanded1_withDE$pQTLCount <- mapply(NumMatchesProtQTL, multExpanded1_withDE$Leading.proteins)
multExpanded1_withDE$GelPreppQTLCount <- mapply(NumMatchesProtQTL, multExpanded1_withDE$ppMajorityProteinIDs)

# This is likely due to multiple protein isoforms corresponding to a single gene. How many pQTLs correspond are mapped to multiple genes (paralogs)?
table(multExpanded1_withDE$pQTLCount)
table(multExpanded1_withDE$GelPreppQTLCount)



