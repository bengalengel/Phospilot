##zia protein workup. goal is to get protein measurements for all 60 human samples. Normalize them against each other (median then quantile normalized) and then subtract (log scale) the relative protein concentration measurements from the relative phospho measurements. Issues will be in isoform quantification and batch correction across the protein measurements.

#these measurements will then be subjected to the same workflow as the confounded phosphomeasurements. 

rm(list=ls(all=TRUE)) #start with empty workspace


# First perform all processing steps using plyr and related tools.
# load required libraries
library(reshape2)
library(stringr)
library(plyr)
source("loadMQZ.R")
source("ExpandPhos.R")
source("counts.R")
source("breakdown.R")

# load protein files with particular variables populated using "loadMQ"
protein <- load.MQZ(directory = "C:/Users/Brett/Documents/Pilot/10_9_14/txt/", type = "protein")

# load protein files with particular variables populated using "loadMQ" at home
protein <- load.MQZ(directory = "E:/My Documents/Pilot/November Zia MBR MQ analysis/txt/")#7,710 protein groups


# remove contaminants and reverse database hits
protein <- protein[(protein$Potential.contaminant != "+" & protein$Reverse != "+"),]#7351


# "only identified by site" hits CAN BE removed because they tend to have lower PEPs (wouldn't pass the FDR TH anyway) and can't be quantified since they are not idd by non-modified peptides. 
# Note there are some high probability proteins here given some proteins are idd by 20+ phosphopeptides.
# eg is A6NKT7 (PEP = 2.23E-70)
protein1 <- protein[(protein$Only.identified.by.site != "+"),]#6705 (really getting down there...)


colnames(protein1)<- gsub(colnames(protein1), pattern = "Ratio.H.L.normalized.", replacement = "HL") ##remove redundant information 

#some strangeness
protein1 <- protein1[3:length(protein1)]

#remove proteins if not quantified in at least one sample
expCol <- grep("HL(.*)", colnames(protein1))

protein1 <- protein1[rowSums(is.na(protein1[,expCol]))!=length(expCol),]##removes rows containing all 
#'NA's using the sums of the logical per row 
#6413

data <- protein1[,expCol]

row.names(data) <- protein1$id
data <- log2(data)

#18870 is all over the place 
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

require(limma)
require(sva)
require(plyr)
require(swamp)
require(statmod)


#quantile normalize
quantiled <- normalizeQuantiles(data,ties = T)#ties are all assigned the same value for the common quantile
summary(quantiled)
boxplot(data)
boxplot(quantiled)
# density plots all look the same
plot.new()
par(mfrow = c(1, 1))
for (i in 1:(ncol(quantiled))){
  if(i==1) plot(density(quantiled[, i], na.rm=T), col = i, ylim = c(0,.9))
  else lines(density(quantiled[, i], na.rm=T), col = i)
}


##lets gut us some protein from the three samples of interest
data <- data[,c("HL18862","HL18486","HL19160")]

#and with quants in all three
datacomp <- na.omit(data)#3,925

#add an indicator to protein1 and then subset the gene names or protein names or whatever to get a count of of how many proteins I can use to normalize the phosphosites with...

#subset protein1 by id
tmp <- protein1[protein1$id %in% row.names(datacomp),]

##add the gene names etc of interest

datacomp <- cbind(datacomp,tmp[c("Protein.IDs","Majority.protein.IDs","Protein.names","Gene.names","Sequence.coverage....")])

#now how many of these are in the phospho data subject to DE?

#subset of multExpanded that are subjected to DE
DEany <- multExpanded1_withDE[multExpanded1_withDE$DEcont1=="+"|multExpanded1_withDE$DEcont2=="+"|multExpanded1_withDE$DEcont3=="+",]

#remove the reverse for now
DEany <- DEany[DEany$Protein[!grepl("REV",DEany$Protein)],]#1547 observations

#How many protiens subject to DE?
DEtable <- as.matrix(table(DEany$Protein))#the parenthetical should be a passed DE factor vector of networKin output
DEtable <- DEtable[DEtable!=0,,drop=F]#270 proteins

#Remove isoform designation (this is simply a guess anyway and should be noted in the discussion of the paper)
row.names(DEtable) <- substr(row.names(DEtable),1,6)

#I will use the majority protein IDs for each protein group quantification (these proteins have at least half the peptides of the leading protein within the group)
Ziaproteins <- datacomp$Majority.protein.IDs

#unparse and turn into a long string for comparison
#now I need to design a loop to dig into the phospho file 
test <- strsplit(as.character(Ziaproteins), ";")
test <- as.character(unlist(test))
test <- unique(test)
any(duplicated(test))#no duplicates with 11500 proteins

test <- substr(test,1,6)#now I have some duplicates (3500 duplicates)
test <- unique(test)#around 8116 unique proteins with at least 1/2 of the peptides of the majority protein within 3925 protien groups

table(row.names(DEtable)%in%test)
#171 of 270 (63.3% are quantified through Zia's work)

#how many of the overlapping proteins map to unique protein groups in the Zia dataset? (the same run gets around 40% of the protein groups ) I think this should be all of them.

#from Mann silac ratios 2010 paper:
# "For the phosphopeptides shared between multiple protein identifiers, the
# identifier with the maximum occupancy stoichiometry was used in the
# analysis."


#look for a uniuque identifier within a list of colon separated identifiers (grep?)
row.names(DEtable)%in%as.character(Ziaproteins)#fix isoforms

Ziaproteins2 <- substr(Ziaproteins,1,6)
row.names(DEtable)%in%as.character(Ziaproteins2)#better but removes the other majority proteins
table(row.names(DEtable)%in%as.character(Ziaproteins2))#


Ziaproteins3 <- gsub("-.", "", Ziaproteins)#removes the isoform indicator

table(row.names(DEtable)%in%as.character(Ziaproteins3))#hmm need it to match ANY of the semicolon separated values within an element



# a quick for loop for each level of rownames that returns the number of hits and the protein ids that they match
facttemp <- as.factor(row.names(DEtable))
proteinindex <- c()
morethan1 <- c()
#for each unique kinase in the DE1 kinases table I need to make a contingency table
for(i in seq_along(DEtable)){
  tmp <- grep(row.names(DEtable)[i], Ziaproteins3)
  #more than one value? this should not happen
  if(length(tmp)>1){
    morethan1 <- c(morethan1,tmp)
  }
  #if length of tmp >0 add to index
  if(length(tmp)>0){
    proteinindex <- c(proteinindex,tmp)
  }
}

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








