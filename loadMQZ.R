load.MQZ <- function(directory) {
  ##This function returns the maxquant file requested with the relevant columns.
  #Note that it also returns the NORMALIZED column values!!!
  
  ## 'directory' is a character vector of length 1 indicating
  ## the location of the maxquant files. 
  ## Note this must be passed to the function using ""
  
  ## all column headers are for v1.5.0.3 

    ##construct the file path
    filepath <- file.path(directory,"proteinGroups.txt")
    
    ## open the file
    data <- read.table(file=filepath, sep = "\t", header=T, fill = T, quote = "")
    
    #select the columns of interest 
    vars <- c("id", "Protein.IDs", "Majority.protein.IDs",  "Protein.names", "Gene.names", "Number.of.proteins", "Peptides", 
              "Razor...unique.peptides", "Unique.peptides", "Sequence.coverage....", "Mol..weight..kDa.", "Sequence.length",
              "PEP", "Peptide.IDs", "Mod..peptide.IDs", "Only.identified.by.site", "Potential.contaminant", "Reverse",
              "Razor...unique.peptides.18862", "Razor...unique.peptides.18486", "Razor...unique.peptides.19160")
    
    other_data <- data[,vars]
    
    ##dataframe that collects only the relevent expression columns. NOTE THE NEED TO USE REP!!!!!
    ##The sample number precedes 'Rep' (technical replicate) and the triple underscore denotes the multiplicity 
    expression <- data[,grep("Ratio.H.L.normalized.(.*)", colnames(data))]
    
    ##combine the two
    data <- cbind(expression,other_data)
  ## return the opened MQ file and have a nice day
  
  return(data)
}