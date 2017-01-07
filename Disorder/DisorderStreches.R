#This script calculates the longest contiguous stretch of disordered residues for each protein within the Iupred list
Iupred.disorder.contiguous.max <- sapply(Iupred, function(x) {
  prediction <- x$DisProb >= .5
  run.lengths <-  rle(prediction)
  with(run.lengths, {
    max(lengths[values == T])
  })
})


# For each proteingroup assigned to a phosphosite, assign the max disorder stretch. 

disorder.length.assignment <- function(proteins){
proteins <- strsplit(as.character(proteins), ";")
proteins <- as.character(unlist(proteins))
proteins <- proteins[!grepl("REV", proteins)]
tmp <- c()
for(i in seq_along(proteins)){
  diso.length <- Iupred.disorder.contiguous.max[[which(names(Iupred.disorder.contiguous.max)==proteins[i])]]
  tmp <- c(tmp, diso.length)
}
max(tmp)
}


#gelprep assignment (can be included with annannotation script)
multExpanded1_withDE_annotated$GelPrep.Disorder.MaxLength <- mapply(disorder.length.assignment, as.character(multExpanded1_withDE_annotated$ppMajorityProteinIDs))