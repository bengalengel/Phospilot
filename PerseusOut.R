PerseusOut <- function(directory){
#   This function produces the data matrix output required for the 'proteomic ruler' based absolute quantification of proteins within Perseus software
#   Input is directory to zia's protien groups 
#   Total, Heavy, and Light intensities will be used for the proteins estimated by Zia's approach to start
#   WILL ADD THE PROTEIN ESTIMATES FROM THE PHOSPHO DATA LATER? 
  
  #read in the required data from zia proteome and remove contaminants and other things from the dataset!
  
  ## 'directory' is a character vector of length 1 indicating
  ## the location of the maxquant files. 
  ## Note this must be passed to the function using ""
  
  ## all column headers are for v1.5.0.3 
  
  ##construct the file path
  filepath <- file.path(directory,"proteinGroups.txt")
  
  ## open the file
  data <- read.table(file=filepath, sep = "\t", header=T, fill = T, quote = "")
  
  #select the columns of interest 
  vars <- c("id", "Protein.IDs", "Majority.protein.IDs", "Number.of.proteins", "Peptides", 
            "Razor...unique.peptides", "Unique.peptides", "Sequence.coverage....", "Sequence.length", "Sequence.lengths",
            "Mol..weight..kDa.", "PEP", "Peptide.IDs", "Mod..peptide.IDs", "Only.identified.by.site", 
            "Potential.contaminant", "Reverse", "Razor...unique.peptides.18862", "Razor...unique.peptides.18486", 
            "Razor...unique.peptides.19160")
  
  other_data <- data[,vars]
  
  ##dataframe that collects only the relevent expression columns. NOTE THE NEED TO USE REP!!!!!
  ##The sample number precedes 'Rep' (technical replicate) and the triple underscore denotes the multiplicity 
  intensity <- data[,grep("Intensity.*18486|Intensity.*18862|Intensity.*19160", colnames(data))]
  
#   expression <- data[,grep("Ratio.H.L.normalized.(.*)", colnames(data))]
  
  ##combine the two
  ZiaIntensities <- cbind(intensity,other_data)

  ##remove contaminants and reverse entries
#   remove contaminants and reverse database hits
  ZiaIntensities <- ZiaIntensities[(ZiaIntensities$Potential.contaminant != "+" & ZiaIntensities$Reverse != "+"),]#7354/78

  # "only identified by site" hits CAN BE removed because they tend to have lower PEPs (wouldn't pass the FDR TH anyway) and can't be quantified   since they are not idd by non-modified peptides. 
  ZiaIntensities <- ZiaIntensities[(ZiaIntensities$Only.identified.by.site != "+"),]#6705

  #keep complete cases of intensity in both channels
  #remove proteins if not quantified in at least one sample
  expCol <- grep("Intensity(.*)", colnames(ZiaIntensities))

  ZiaIntensities <- ZiaIntensities[rowSums(ZiaIntensities[,expCol] > 0) == length(expCol),]##removes rows containing any 0s (4928 complete cases)
  
  #upload FASTA from archive site; 
  #create FASTA directory if it doesn't already exist
  dir.create(file.path(getwd(), "Perseus"))

  write.csv(ZiaIntensities, "./Perseus/ZiaIntensities.csv", row.names=F)
}





  
  