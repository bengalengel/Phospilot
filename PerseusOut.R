PerseusOut <- function(directory, dataset = c("Zia","Brett")){
#   This function produces the data matrix output required for the 'proteomic ruler' based absolute quantification of proteins within Perseus software
  #   Input is directory to Brett's or zia's protien groups 
  
  #required data is read and contaminants and reverse hits are removed from the dataset!
  
  ## 'directory' is a character vector of length 1 indicating
  ## the location of the maxquant files. 
  ## Note this must be passed to the function using ""
  
  ## all column headers are for v1.5.0.3 
  if(dataset == "Zia"){
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
    
    #create Perseus directory if it doesn't already exist and write the file
    dir.create(file.path(getwd(), "Perseus"))
    write.csv(ZiaIntensities, "./Perseus/ZiaIntensities.csv", row.names=F)
  }
  
  if(dataset == "Brett"){
    ##construct the file path
    filepath <- file.path(directory,"proteinGroups.txt")
    
    ## open the file
    data <- read.table(file=filepath, sep = "\t", header=T, fill = T, quote = "")
    
    #select the columns of interest 
    vars <- c("id", "Protein.IDs", "Majority.protein.IDs", "Number.of.proteins", "Peptides", 
              "Razor...unique.peptides", "Unique.peptides", "Sequence.coverage....", "Sequence.length", "Sequence.lengths",
              "Mol..weight..kDa.", "PEP", "Peptide.IDs", "Mod..peptide.IDs", "Only.identified.by.site", 
              "Potential.contaminant", "Reverse")
    
    other_data <- data[,vars]
    
    ##dataframe that collects only the relevent expression columns. NOTE THE NEED TO USE REP!!!!!
    ##The sample number precedes 'Rep' (technical replicate) and the triple underscore denotes the multiplicity 
    intensity <- data[,grep("Intensity.*18486|Intensity.*18862|Intensity.*19160", colnames(data))]
    
    ##combine the two
    BrettIntensities <- cbind(intensity,other_data)
    
    ##remove contaminants and reverse entries
    #   remove contaminants and reverse database hits
    BrettIntensities <- BrettIntensities[(BrettIntensities$Potential.contaminant != "+" & BrettIntensities$Reverse != "+"),]#5961
    
    # "only identified by site" hits CAN BE removed because they tend to have lower PEPs (wouldn't pass the FDR TH anyway) and can't be quantified   since they are not idd by non-modified peptides. 
    BrettIntensities <- BrettIntensities[(BrettIntensities$Only.identified.by.site != "+"),]#5127
    
    #keep complete cases of intensity in both channels
    #remove proteins if not quantified in at least one sample. For this run I will use
    expCol <- grep("Intensity.[0-9]+_[12]_[12]", colnames(BrettIntensities))
    
    BrettIntensities <- BrettIntensities[rowSums(BrettIntensities[,expCol][1:4] > 0) >= 1 & rowSums(BrettIntensities[,expCol][5:8] > 0) >= 1 & rowSums(BrettIntensities[,expCol][9:12] > 0) >= 1,]##removes rows containing all 0s in any of the samples(3036 complete cases)
    
    #caluculate the average intensity for each line for H,L, and total intensity
    
    #replace 0s with NAs
    BrettIntensities[BrettIntensities == 0] <- NA
    
    #expCol gives total intensity
    BrettIntensities$IntensityMean18486 <- rowMeans(BrettIntensities[,expCol[1:4]], na.rm = T)
    BrettIntensities$IntensityMean18862 <- rowMeans(BrettIntensities[,expCol[5:8]], na.rm = T)
    BrettIntensities$IntensityMean19160 <- rowMeans(BrettIntensities[,expCol[9:12]], na.rm = T)
    
    #expColL and expColH
    expColL <- grep("Intensity.L.[0-9]+_[12]_[12]", colnames(BrettIntensities))
    expColH <- grep("Intensity.H.[0-9]+_[12]_[12]", colnames(BrettIntensities))
    
    BrettIntensities$IntensityLightMean18486 <- rowMeans(BrettIntensities[,expColL[1:4]], na.rm = T)
    BrettIntensities$IntensityLightMean18862 <- rowMeans(BrettIntensities[,expColL[5:8]], na.rm = T)
    BrettIntensities$IntensityLightMean19160 <- rowMeans(BrettIntensities[,expColL[9:12]], na.rm = T)
    
    BrettIntensities$IntensityHeavyMean18486 <- rowMeans(BrettIntensities[,expColH[1:4]], na.rm = T)
    BrettIntensities$IntensityHeavyMean18862 <- rowMeans(BrettIntensities[,expColH[5:8]], na.rm = T)
    BrettIntensities$IntensityHeavyMean19160 <- rowMeans(BrettIntensities[,expColH[9:12]], na.rm = T)
    
    
    #create Perseus directory if it doesn't already exist and write the file
    dir.create(file.path(getwd(), "Perseus"))
    write.csv(BrettIntensities, "./Perseus/BrettIntensities.csv", row.names=F)
  }
}  
  
  
  
  
  
  