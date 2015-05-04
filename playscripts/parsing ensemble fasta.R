#parsing FASTA files from Ensembl for protein database search. Downloaded on 4/29/15. Ensembl version 75 Grch37/hg19
library(seqinr)

#load the ENSP FASTA database
EnsProt <- read.fasta(file = "C:/Users/Brett/Dropbox/Postdoc-Gilad/Ensembl Database/Homo_sapiens.GRCh38.pep.all.fa", seqtype = "AA", as.string = T)

length(EnsProt)
#100778
is.SeqFastaAA(EnsProt[[1]])

names(EnsProt[1])#gives ENSP identifier

attributes(EnsProt[[1]])# each list object has an 'Annot' attribute associated with it. This has all the fasta header information
attr(EnsProt[[1]], 'Annot')

#subset list st genes are only 'known' and 'novel'. Only those list objects with 'pep:novel' or 'pep:known'. This is incorrect. Some genes have entries that are listed both pep:known and pep:putative etc. Some 5' CDS incomplete are annotated as 'known'
genetype <- "pep:novel|pep:known"

#fist get the annotation
annotation <- getAnnot(EnsProt)

#count number of known or novel genes
genetype_index <- grep(genetype, annotation)
length(genetype_index)
#823337
#10665 novel
#71672 known

#subset fasta
EnsProt <- EnsProt[genetype_index]


#subset FASTA st gene biotypes are 'protein_coding' or 'polymorphic_pseudogene'
genebiotype <- "gene_biotype:protein_coding|gene_biotype:polymorphic_pseudogene"

#retrieve updated annotation
annotation <- getAnnot(EnsProt)

genebiotype_index <- grep(genebiotype, annotation)
length(genebiotype_index)
#81878

#subset
EnsProt <- EnsProt[genebiotype_index]


#subset FASTA st transcript biotypes are "protein_coding"
tranbiotype <- "transcript_biotype:protein_coding"
#retrieve updated annotation
annotation <- getAnnot(EnsProt)

tranbiotype_index <- grep(tranbiotype, annotation)
length(tranbiotype_index)
#67993
#surprising drop but ah well

#subset
EnsProt <- EnsProt[tranbiotype_index]


#subset out seqences with an 'X' 5' UTR CDS truncation and a '.' 3' UTR truncation
start_trunc <- grep("X", EnsProt)
stop_trunc <- grep("\.", EnsProt)#nothing for the 3' UTR

#subset
EnsProt <- EnsProt[!start_trunc]

#before I write out the FASTA file I would like to add the HGNC gene symbol and approved name. The 'homo.sapiens' package can do this. What about biomaRt?

#extract for each 




#write the new fasta file
write.fasta(EnsProt, names = getAnnot(EnsProt), "C:/Users/Brett/Dropbox/Postdoc-Gilad/Ensembl Database/Homo_sapiens.GRCh38.pep.all.parsed.fa", as.string=T)









grepl(genetype, attr(EnsProt[[1]], 'Annot'))

Store the fasta header
annotation <- getAnnot(ncrna)
#count the number of piRNAs
pirna_index <- grep("piRNA",annotation,ignore.case=T)
length(pirna_index)
[1] 174724
pirna <- ncrna[pirna_index]



lapply(test,function(x) attr(EnsProt[[x]], 'Annot'))


mylist <- list(1:5, 6:10, 11:15)
sapply(mylist, "[", c(2,3))


test <- EnsProt[1:10]

tout <- lapply(test, function(x) grep(genetype, attr(EnsProt[[x]], 'Annot')))


tout <- lapply(test, function(x) grepl(genetype, attr(EnsProt[[x]], 'Annot')))

list.condition <- sapply(input.list, function(x) class(x)=="desired.class")
output.list  <- input.list[list.condition]



#############
#how many fasta sequences are in this file?
#     length(proteome)#
#     [1] 88993
#     

# Each element is a sequence object of the class SeqFastadna or SeqFastaAA.
# is.SeqFastaAA(proteome[[1]])
# [1] TRUE

#from interesting list i = 2623
################