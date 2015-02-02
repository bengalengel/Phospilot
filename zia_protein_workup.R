##zia protein workup. goal is to get protein measurements for all 60 human samples. Normalize them against each other (median then quantile normalized) and then subtract (log scale) the relative protein concentration measurements from the relative phospho measurements. Issues will be in isoform quantification and batch correction across the protein measurements.

#these measurements will then be subjected to the same workflow as the confounded phosphomeasurements. 

rm(list=ls(all=TRUE)) #start with empty workspace


# First perform all processing steps using plyr and related tools.
# load required libraries
library(reshape2)
library(stringr)
library(plyr)
require(limma)
require(sva)
require(swamp)
require(statmod)
source("loadMQZ.R")
#source("ExpandPhos.R")
#source("counts.R")
#source("breakdown.R")

# load protein files with particular variables populated using "loadMQ"
protein <- load.MQZ(directory = "D:/November Zia MBR MQ analysis/txt/")#7710 protein groups/81 variables

# load protein files with particular variables populated using "loadMQ" at home
#protein <- load.MQZ(directory = "E:/My Documents/Pilot/November Zia MBR MQ analysis/txt/")#7,710 protein groups


# remove contaminants and reverse database hits
protein <- protein[(protein$Potential.contaminant != "+" & protein$Reverse != "+"),]#7354/78


# "only identified by site" hits CAN BE removed because they tend to have lower PEPs (wouldn't pass the FDR TH anyway) and can't be quantified since they are not idd by non-modified peptides. 
# Note there are some high probability proteins here given some proteins are idd by 20+ phosphopeptides.
# eg is A6NKT7 (PEP = 2.23E-70)
protein1 <- protein[(protein$Only.identified.by.site != "+"),]#6705 (really getting down there...)


colnames(protein1)<- gsub(colnames(protein1), pattern = "Ratio.H.L.normalized.", replacement = "HL") ##remove redundant information 

#some strangeness sometimes there are two extra rows!
#protein1 <- protein1[2:length(protein1)]

#remove proteins if not quantified in at least one sample
expCol <- grep("HL(.*)", colnames(protein1))

protein1 <- protein1[rowSums(is.na(protein1[,expCol]))!=length(expCol),]##removes rows containing all 
#'NA's using the sums of the logical per row 
#6413 now I have 6421 (previous workflow uploaded incorrect column)

data <- protein1[,expCol]#60 cell lines

row.names(data) <- protein1$id
data <- log2(data)

#18871 is all over the place. actually its 18871?
boxplot(data)
which.max(sapply(data,median,na.rm=T))

#I am going to remove it for now
drop <- names(which.max(sapply(data,median,na.rm=T)))
data <- data[,!(names(data) %in% drop)]

#median normalize
names <- colnames(data)
median.subtract <- function(x){ x - median(x, na.rm = TRUE)}##create a wrapper for median subtraction
data <- colwise(median.subtract, names)(data) #create median subtracted data but loose intensity and the row names here...

#add back protien ids
row.names(data) <- protein1$id

#summaries
summary(data)
boxplot(data)#

#subset data to only the samples of interest
data <- data[,c("HL18862","HL18486","HL19160")]

#remove if protein group not found in all samples
data <- na.omit(data)#3,925 now this is 4270? Confirm from home
boxplot(data)#differences in distribution shape for sure with HL18486 and HL19160
par(mfrow = c(1, 1))
for (i in 1:(ncol(data))){
  if(i==1) plot(density(data[, i], na.rm=T), col = i, ylim = c(0,1.5))
  else lines(density(data[, i], na.rm=T), col = i)
}


#quantile normalize using all data for now
quantiled <- normalizeQuantiles(data,ties = T)#ties are all assigned the same value for the common quantile
summary(quantiled)
boxplot(data)
boxplot(quantiled)
# density plots all look the same
plot.new()
par(mfrow = c(1, 1))
for (i in 1:(ncol(quantiled))){
  if(i==1) plot(density(quantiled[, i], na.rm=T), col = i, ylim = c(0,1.9))
  else lines(density(quantiled[, i], na.rm=T), col = i)
}


##normalized protein ratios from the three samples of interest
#quantiled <- quantiled[,c("HL18862","HL18486","HL19160")]

#and with quants in all three
#datacomp <- na.omit(quantiled)#3,925 now this is 4270? Confirm from home
#boxplot(datacomp)

#add an indicator to protein1 and then subset the gene names or protein names or whatever to get a count of of how many proteins I can use to normalize the phosphosites with...

#subset protein1 by id
tmp <- protein1[protein1$id %in% row.names(quantiled),]

##add the gene names etc of interest

datacomp <- cbind(quantiled,tmp[c("Protein.IDs","Majority.protein.IDs","Protein.names","Gene.names","Sequence.coverage....",
                                 "Number.of.proteins", "Peptides", "Razor...unique.peptides", "Unique.peptides",
                                 "Razor...unique.peptides.18862", "Razor...unique.peptides.18486", 
                                 "Razor...unique.peptides.19160")])


#Protein quants need to be assigned to phosphosites.
#*******************************************************************************************************************

#now how many unique proteins in the phospho data are subjected to DE analysis? 

#subset of multExpanded that are subjected to DE (load a local file with DE information)
multExpanded1_withDE <- read.csv("multExpanded1.csv", header=T)

#names(multExpanded1_withDE) <- sub("X.","",names(multExpanded1_withDE))#replace Xs
#names(multExpanded1_withDE) <- sub(".$","",names(multExpanded1_withDE))#replace dot at the end


SubtoDE <- multExpanded1_withDE[as.character(multExpanded1_withDE$SubtoDE) == "+",]

#DEany <- multExpanded1_withDE[multExpanded1_withDE$DEcont1=="+"|multExpanded1_withDE$DEcont2=="+"|multExpanded1_withDE$DEcont3=="+",]

#remove the reverse for now
SubtoDE <- SubtoDE[SubtoDE$Protein[!grepl("REV",SubtoDE$Protein)],]#4991 observations

#How many protiens subject to DE?
SubtoDEtable <- as.matrix(table(SubtoDE$Protein))
SubtoDEtable <- SubtoDEtable[SubtoDEtable!=0,,drop=F]#1122 proteins

#Remove isoform designation (this should be noted in the discussion of the paper. Unless there is a specific sequence aligning the peptide to a particular isoform the peptide could reasonably belong to any of the isoforms. Is there a designator in the table indicating if this is so?)

row.names(SubtoDEtable) <- substr(row.names(SubtoDEtable),1,6)#eventually 
SubtoDEproteins <- row.names(SubtoDEtable)
SubtoDEproteins <- unique(SubtoDEproteins)#1115 unique proteins (excluding isoforms) subjected to DE

#*******************************************************************************

#table for comparison with protein. 

# I need protein ID without isoform and ID_multiplicity for phosphopeptide. For matches with Ziaproteins, a new column will be added to that datatable, "norm_id_mult" that maps to the ids this protein group will normalize. For proteins that map to multiple ID_mults they will be added as a semicolon separated list. 

# This Ziaprotein table will be crossed with the original multexpanded table to identify any sites with matching ids. if there is a match, the three measurments of that protein group will be appended to the data table. Three new columns will then be populated containing the "normalized" measurements. 

#How many protiens subject to DE?
SubtoDEtable <- SubtoDE[c("Protein","idmult")]
SubtoDEtable$Protein <- substr(SubtoDEtable$Protein,1,6)#remove isoform designation






#****************************************************************************************************************
#How many proteins that are subjected to DE analysis are also IDd and quantified by proteomic analysis (Zia)? Here I will use the majority protein IDs. I will use the majority protein IDs for each protein group quantification (these proteins have at least half the peptides of the leading protein within the group)
Ziaproteins <- datacomp$Majority.protein.IDs

#unparse and turn into a long string for comparison
#now I need to design a loop to dig into the phospho file 
test <- strsplit(as.character(Ziaproteins), ";")
test <- as.character(unlist(test))
test <- unique(test)
any(duplicated(test))#no duplicates with 11500 proteins (now 12642 proteins with proper workflow)

test <- substr(test,1,6)#now I have some duplicates (3500 duplicates)
test <- unique(test)#around 8116 unique proteins with at least 1/2 of the peptides of the majority protein within 3925 protien groups. Now I have 8885 unique protein with at least 1/2 of the peptides of the majority proteins within 4270 groups.

table(SubtoDEproteins%in%test)
#171 of 270 (63.3% are quantified through Zia's work) (of the DE in any contrast phosphosites using the flawed workflow).
#713 of 1115 (63.9% are quantified in all three samples in Zia's work)




#how many of the overlapping proteins map to unique protein groups in the Zia dataset? (the same run gets around 40% of the protein groups ) I think this should be all of them.

#from Mann silac ratios 2010 paper:
# "For the phosphopeptides shared between multiple protein identifiers, the
# identifier with the maximum occupancy stoichiometry was used in the
# analysis."


#look for a uniuque identifier within a list of colon separated identifiers (grep?)
#row.names(DEtable)%in%as.character(Ziaproteins)#fix isoforms

# Ziaproteins2 <- substr(Ziaproteins,1,6)
# row.names(DEtable)%in%as.character(Ziaproteins2)#better but removes the other majority proteins
# table(row.names(DEtable)%in%as.character(Ziaproteins2))#

# Ziaproteins <- datacomp[c("Majority.protein.IDs"v)]
# Ziaproteins3 <- gsub("-.", "", Ziaproteins)#removes the isoform indicator
# Ziaproteins$id <- row.names(datacomp)
# Ziaproteins <- cbind(Ziaproteins,datacomp[c("Razor...unique.peptides", "Unique.peptides"             "Razor...unique.peptides.18862", "Razor...unique.peptides.18486", "Razor...unique.peptides.19160")]

Ziaproteins <- datacomp[c("Majority.protein.IDs","Razor...unique.peptides", "Unique.peptides", "Razor...unique.peptides.18862", "Razor...unique.peptides.18486", "Razor...unique.peptides.19160")]
Ziaproteins$id <- row.names(datacomp)
Ziaproteins$Majority.protein.IDs <- gsub("-.", "", Ziaproteins$Majority.protein.IDs)#removes the isoform indicator

                     

                     # table(row.names(DEtable)%in%as.character(Ziaproteins3))#hmm need it to match ANY of the semicolon separated values within an element



# a quick for loop for each level of rownames that returns the number of hits and the protein ids that they match
facttemp <- as.factor(row.names(SubtoDEtable))
proteinindex <- c()
morethan1 <- c()

## I need protein ID without isoform and ID_multiplicity from the phosphopeptide table (multexpanded). For matches with Ziaproteins, a new column will be added to that datatable, "norm_id_mult" that maps to the ids this protein group will normalize. For Ziaproteins that map to multiple phosphoobservations they will be added as a semicolon separated list. 

#for every protein linked to an id_mult from the phosphotable, a paired protein group from the ziaproteins table is found (if present) using any of the majority protein ids within that group. If the phospho id maps to multiple protein groups, the one with the most peptides is used.

#are there any problems with this?

for(i in seq_along(SubtoDEtable[,1])){
  tmp <- grep(SubtoDEtable$Protein[i], Ziaproteins$Majority.protein.IDs)
  #more than one value? This can happen with isoforms
  if(length(tmp)>1){
    #compare the razor plus unique count across the two matches and choose the one with the most   matches
    counts <- Ziaproteins$Razor...unique.peptides[tmp]
    proteinindex[i] <- tmp[which.max(counts)]
    morethan1 <- c(morethan1,tmp)
  }
  #if length of tmp >0 add to index
  if(length(tmp) == 1){
    proteinindex[i] <- tmp
  }
}

SubtoDE$ziaindex <- proteinindex

##for each index value add the three values and the majority ids. next for loop performs the normalization.

protein_norm <- data.frame()
for(i in seq_along(proteinindex)){
  if(!is.na(proteinindex[i])){
    tmp <- datacomp[c("Majority.protein.IDs","HL18862", "HL18486", "HL19160")][proteinindex[i],]
    protein_norm <- rbind(protein_norm,tmp)
  }
  
  if(is.na(proteinindex[i])){
    tmp <- rep("NA",4)
    protein_norm <- rbind(protein_norm,tmp)
  }
}    

#make values numeric
protein_norm$HL18862 <- as.numeric(protein_norm$HL18862)
protein_norm$HL19160 <- as.numeric(protein_norm$HL19160)
protein_norm$HL18486 <- as.numeric(protein_norm$HL18486)

#add values to subtoDE dataframe
SubtoDE$Majority_proteins_Zia <- protein_norm$Majority.protein.IDs
SubtoDE$protein_ratio_HL18862 <- protein_norm$HL18862
SubtoDE$protein_ratio_HL18486 <- protein_norm$HL18486
SubtoDE$protein_ratio_HL19160 <- protein_norm$HL19160

#since the normalized and batch treated samples are also not included they need to be added.




# When there are two matches this tends to mean that here is unique evidence for an isoform which is quantified separately. I will use the measurements from the (major) isoform with the most razor + unique peptide for normalization. (for the DE peptides only two fulfill this criteria)









mt1test <- 
proteinindex <- c(proteinindex,tmp)
    DErow <- c(DEtable[as.character(i),],(sum(DEtable)-DEtable[as.character(i),]))
    NotDErow <- c(NotDE[as.character(i),],(sum(NotDE)-NotDE[as.character(i),]))
    contmatrix <- rbind(DErow,NotDErow)
    tmp <- fisher.test(contmatrix, alternative = "g")
    pvals <- c(pvals,tmp$p.value)
  
}








pmatch(row.names(DEtable),as.character(Ziaproteins3))

lapply(grep(row.names(DEtable), as.character(Ziaproteins3)))
grep(row.names(DEtable), as.character(Ziaproteins3))

apply(row.names()








