#read iupred results into a list. unpack the tar.gz file (tarball) first
require("foreach")
require("iterators")
require("doParallel")


#get a vector with all the filenames
iupred.files <- list.files("./Disorder/Iupred/Results/Indfiles", full.names = T)#can sysglob here as well

#read these into a list of dataframes using foreach parallel processing.
cl <- makeCluster(5)
registerDoParallel(cl)
Iupred <- foreach(i=seq_along(iupred.files)) %dopar% {
  read.table(iupred.files[i], sep = "", col.names = c("Position", "AA", "DisProb"), stringsAsFactors = F)
}
stopCluster(cl)


#add enspid as names to this file
Iupred.names <- list.files("./Disorder/Iupred/Results/Indfiles", full.names = F)#can sysglob here as well
Iupred.names <- gsub("[^0-9]","", Iupred.names)
Iupred.names <- paste0("ENSP", Iupred.names)
names(Iupred) <- Iupred.names

