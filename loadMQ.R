load.MQ <- function(type, directory) {
  ##This function returns the maxquant file requested with the relevant columns.
  #Note that it also returns the NORMALIZED column values!!!   
  
  ## 'type' is a character vector of length 1 indicating
  ## the type of file to be opened. The choices are 'phospho', 'protein', and 'evidence'
  ## Note this must be passed to the function using ""
  
  ## 'directory' is a character vector of length 1 indicating
  ## the location of the maxquant files. 
  ## Note this must be passed to the function using ""
  
  ## all column headers are for v1.5.0.3 
  
  if(type=="phospho"){
    ##construct the file path
    filepath <- file.path(directory,"Phospho (STY)sites.txt")
    
    ## open the file
    data <- read.table(file=filepath, sep = "\t", header=T, fill = T, quote = "")
    
    #select the columns of interest     
    vars <- c("id","Amino.acid","Charge","Reverse","Potential.contaminant","Proteins","Positions.within.proteins","Leading.proteins",
              "Sequence.window","Phospho..STY..Probabilities","Localization.prob","PEP", "Score", "Delta.score", "Score.for.localization", 
              "Mass.error..ppm.", "Intensity", "Intensity.L", "Intensity.H", "Position", "Number.of.Phospho..STY.", 
              "Protein.group.IDs", "Protein")
    
    other_data <- data[,vars]
    
    ##dataframe that collects only the relevent expression columns. NOTE THE NEED TO USE REP!!!!!
    ##The sample number precedes 'Rep' (technical replicate) and the triple underscore denotes the multiplicity 
    expression <- data[,grep("Ratio.H.L.normalized(.*)_[12]___", colnames(data))]
    
    ##combine the two
    data <- cbind(expression,other_data)
  }
  
  if(type=="protein"){
    ##construct the file path
    filepath <- file.path(directory,"ProteinGroups.txt")
    
    ## open the file
    data <- read.table(file=filepath, sep = "\t", header=T, fill = T, quote = "")
    
    #select the columns of interest 
    vars <- c("id", "Protein.IDs", "Majority.protein.IDs",  "Protein.names", "Gene.names", "Number.of.proteins", "Peptides", 
              "Razor...unique.peptides", "Unique.peptides", "Sequence.coverage....", "Mol..weight..kDa.", "Sequence.length", "PEP", 
              "Peptide.IDs", "Mod..peptide.IDs", "Phospho..STY..site.IDs")
    
    other_data <- data[,vars]
    
    ##dataframe that collects only the relevent expression columns. NOTE THE NEED TO USE REP!!!!!
    ##The sample number precedes 'Rep' (technical replicate) and the triple underscore denotes the multiplicity 
    expression <- data[,grep("Ratio.H.L.normalized(.*)_[12]___", colnames(data))]
    
    ##combine the two
    data <- cbind(expression,other_data)
  }
  ## return the opened MQ file and have a nice day
  
  return(data)
}
  
    
    
    
    
    
    
    
