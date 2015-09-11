###output CCDS fasta for RAPID protein % intrinsic disorder prediction ------------
require("seqinr")
require("foreach")
require("iterators")
require("doParallel")

#load proteome. working directory is phospilot.
proteome <- read.fasta( file = "./FASTA/Homo_sapiens.GRCh37.75.pep.all.parsedCCDS.fa", seqtype = "AA", as.string = TRUE)

#make new headers to the fasta file using only the numeric regions of the ENSPIDs (<12 chars)
headers <- lapply(names(proteome), function(x){
  id <- unlist(strsplit(x, "[|]"))[2]
  id <- gsub("[A-z]","", id)#keep last 11 characters
})

#write parsed fasta with new headers
filepath <- file.path(getwd(),"Disorder/Rapid/Homo_sapiens.GRCh37.75.pep.all.CCDS.ENSPIDs.fa")
write.fasta(proteome, names = headers, file = filepath, as.string=T)


###output CCDS individual protein fastas for Iupred predicitons -------

#Iupred parsed FASTA files
#create Indfiles directory and run this script to create parsed fasta files. This could use a foreach solution
if(!file.exists("./Disorder/IUPred/Indfiles")) dir.create("./Disorder/IUPred/Indfiles")
cl <- makeCluster(5)
registerDoParallel(cl)
foreach(i=seq_along(proteome), .packages = "seqinr") %dopar% {
  filepath <- file.path(getwd(), paste("Disorder/IUPred/Indfiles/",headers[i], ".fa", sep = ""))
  write.fasta(proteome[i], names = headers[i], file = filepath, as.string=T)
}
stopCluster(cl)

# create tar.gz the individual fasta file directory for sharing using shell

###Interpro FASTAs ------

#break proteome into regular intervals and write to separate files using truncated headers
# Files of the type 'Interpro.$SGE_TASK_ID.fa'

#SPLIT 'proteome' and 'headers' into groups of 1000. NOTE THE NUMBER OF ELEMENTS FOR PASSING TO SGE!!
max <- 1000
x <- seq_along(proteome)
proteome.parsed <- split(proteome, ceiling(x/max))
headers.parsed <- split(headers, ceiling(x/max))

#create directory and truncated FASTAs
if(!file.exists("./InterProScan/FASTAs/")) dir.create("./InterProScan/FASTAs/")
cl <- makeCluster(5)
registerDoParallel(cl)
foreach(i=seq_along(proteome.parsed), .packages = "seqinr") %dopar% {
  filepath <- file.path(getwd(), paste("InterProScan/FASTAs/","Interpro.", i, ".fa", sep = ""))
  write.fasta(proteome.parsed[[i]], names = headers.parsed[[i]], file = filepath, as.string=T)
}
stopCluster(cl)

# create tar.gz the individual fasta file directory (using the shell) to share
tar('./InterProScan/FASTAs.tar.gz', './InterProScan/FASTAs/', compression = "gzip")














