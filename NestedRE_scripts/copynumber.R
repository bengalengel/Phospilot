####################perseus
#output data matrix for perseus absolute quantification
PerseusOut(directory = "D:/EnsemblDBProteome/txt/", dataset = "Zia")
PerseusOut(directory = "D:/EnsemblDBPhospho/PilotPhosphoensemblDB/combined/txt/", dataset = "Brett")

#Protein expression level estimation approach within perseus:

#'raw' light intensity values from MQ for each of the three Zia samples.
#'raw' heavy intensity values from MQ for each of the three Zia samples.

#plugin assumes linearity between the Intensity values and the cumulative mass of each protein. In other words, the molecular mass of each protein serves as the correction factor if one wants to calculate copy numbers. HERE I USED THE MEDIAN MW OF EACH OF THE MAJORITY PROTEINS. 'AVERAGE MOLECULR MASS' column in Perseus output. An additional file is produced with a 'detectability correction', or the theoretical number of tryptic peptides per protein. For each file a thresholding of quantification accuracy is performed using the number of peptides observed, the percentage of those peptides that are razor or unique, and the number of theoretical peptides/100 amino acids for each protien group.

#(NOTE THIS IS DIFFERENT THAN HOW MQ DOES IBAQ. IT USES MASS FROM LEADING PROTEIN divided by the NUMBER OF THEORETICAL PEPTIDES FROM LEADING PROTEIN.)

# The total protein amount per cell was used as a scaling factor. In this case the total intensity will be considered proportional to the total protein amount per cell. 200 pg of protein/cell and total cellular protein concentration of 200 g/l is used here.


#input and explore Zia's copynumber data. Include 'corrected' data as well.
Zia <- read.table(file="./Perseus/ZiaOut/TPZia.txt", sep = "\t", header=T, fill = T, quote = "")#column independent (no normalization)
Ziacorrected <- read.table(file="./Perseus/ZiaOut/TPCorrectedZia.txt", sep = "\t", header=T, fill = T, quote = "")#column independent (no normalization)

keepers <- grep("copy.", names(Zia), value = T, ignore.case = T)
keepers <- keepers[c(1,4,7)]# only copynumbers kept
copynumbers <- Zia[,names(Zia) %in% keepers]
names(copynumbers) <- c("cn18486","cn18862","cn19160")
copynumbers <- log10(copynumbers)
boxplot(copynumbers, ylab = "log10(copynumber)")#looks better than cns1

#corrected data
keepers <- grep("copy.", names(Ziacorrected), value = T, ignore.case = T)
keepers <- keepers[c(1,4,7)]# only copynumbers kept
copynumbers2 <- Ziacorrected[,names(Ziacorrected) %in% keepers]
names(copynumbers2) <- c("cn19160corrected","cn18862corrected","cn18486corrected")
copynumbers2 <- log10(copynumbers2)
boxplot(copynumbers2, ylab = "log10(copynumber)")

#combine the two.. these look the same
copynumbersall <- cbind(copynumbers,copynumbers2)
boxplot(copynumbersall)

# The spearman correlation across samples is good. Scatterplots reveal (as expected) most of variation at low intensity end of distribution.

require(gplots)
require(Hmisc)
pairs(copynumbersall)
cn.corr <- rcorr(as.matrix(copynumbersall), type = "spearman")
heatmap.2(
  cn.corr$r,
  key = FALSE,
  dendrogram = "none",
  trace = "none",
  cellnote = round(cn.corr$r,3),
  notecex = 1.25,
  notecol = "black",
  col=bluered(10),
  cexCol=1.254,
  cexRow=1.254,
  margins = c(10,10)
)
#because it doesn't seem to matter I will be using the 'non corrected' values from here on out.

#note this SO article for orienting heatmap using 'layout' http://stackoverflow.com/questions/16683026/r-gplots-removing-white-space-in-heatmap-2-when-key-false


##make a 2D density plot of the median copynumber and look for 1D enrichment of the clusters? Otherwise I will have to categorize/threshold the copynumbers, and this is somewhat arbitrary. I would then have to test the enrichment values by bootstrapping. Easiest way to do this is with the S curve.

#make an S curve of rank vs median(logcn). subset this to those values subject to varcomp with protein normalization. Overlay with color the four quadrants. Perform the same analysis with omnibus enrichment values. 

#get median copynumber and median rank, calculate relative rank.
cnMedian <- apply(as.matrix(copynumbers), 1, median)
Rank <- rank(1/cnMedian)
RelRank <- sapply(Rank,function(x) x/max(Rank))
#combine to create DF
MedCNDF <- data.frame(copynumber = cnMedian, Rank = Rank, RelRank = RelRank)
row.names(MedCNDF) <- Zianonorm$Majority.protein.IDs
par(mfrow = c(1,1))
plot(MedCNDF$RelRank, MedCNDF$copynumber)
dim(MedCNDF)#4928

#now subset and look for bias in cns for the varcomp quadrants

#protein prep sub to varcomp
proteinvarcompNames <- multExpanded1_withDE[multExpanded1_withDE$ppSubtoVarcomp == "+", "ppMajorityProteinIDs"]

##subset copynumber DF st only those subjected to NRE model are considered
MedCNDFvc <- MedCNDF[row.names(MedCNDF) %in% proteinvarcompNames, ]#1152

##VarcompProt formatting fix
VarcompProt$individual <- as.numeric(as.vector(VarcompProt$individual))
VarcompProt$biorep <- as.numeric(as.vector(VarcompProt$biorep))
VarcompProt$residual <- as.numeric(as.vector(VarcompProt$residual))

#add ppMajority IDs and copynumber information to VarcompProt
VarcompProt$ppMajorityIDs <- multExpanded1_withDE[multExpanded1_withDE$idmult %in% row.names(VarcompProt), "ppMajorityProteinIDs"]

##merge by medCNDFvc rownames
#rowNames <- row.names(VarcompProt)
VarcompProtMerged <- merge(VarcompProt,MedCNDFvc, by.x = "ppMajorityIDs", by.y = "row.names")
# row.names(VarcompProt) <- rowNames#lost rownames somehow. A few missing here.

#using red gray blue color gradient to represent copynumber relative rank
rgbpal <- colorRampPalette(c("Red","Gray","Blue"))

#plot final choice. RGB using relcnrank using 1 row and two columns
par(mfcol = c(1,2))
index <- length(unique(VarcompProtMerged$ppMajorityIDs))
plot(VarcompProtMerged$RelRank, VarcompProtMerged$copynumber, 
     col = rgbpal(index)[findInterval(VarcompProtMerged$RelRank,
                                      seq(range(VarcompProtMerged$RelRank)[1], 
                                          range(VarcompProtMerged$RelRank)[2], 
                                          length.out = index), rightmost.closed = T)], 
     pch = 19, ylab = "log10(copynumber)", xlab = "Relative Rank")

#quadrant plot with copynumber RelRank heatmap overlay
plot(log10(VarcompProtMerged$individual),log10(VarcompProtMerged$biorep),
     col = rgbpal(index)[findInterval(VarcompProtMerged$RelRank,
                                  seq(range(VarcompProtMerged$RelRank)[1], 
                                      range(VarcompProtMerged$RelRank)[2], 
                                      length.out = index), rightmost.closed = T)], 
     pch = 19, ylab = "log10 BioRep VarComp", xlab = "log10 Ind VarComp")



##repeat for confounded data. I SHOULD HAVE LET PERSEUS CALCULATE INFORMATION THEN TAKEN THE MEDIAN HERE! NOTE THAT 18862 HAS AN 'NA'!
#######################

#input and explore Brett's copynumber data. Include 'corrected' data as well.
Brett <- read.table(file="./Perseus/BrettOut/TPBrett.txt", sep = "\t", header=T, fill = T, quote = "")#column independent (no normalization)
Brettcorrected <- read.table(file="./Perseus/BrettOut/TPTheoreticalPepCorrectBrett.txt", sep = "\t", header=T, fill = T, quote = "")#column independent (no normalization)

#remove strange NA values from 18862
Brett <- Brett[!is.na(Brett$HeavyMean18862),]
Brettcorrected <- Brettcorrected[!is.na(Brettcorrected$HeavyMean18862),]

keepers <- grep("copy.", names(Brett), value = T, ignore.case = T)
keepers <- keepers[c(1,4,7)]# only copynumbers kept
copynumbers <- Brett[,names(Brett) %in% keepers]
names(copynumbers) <- c("cn18486","cn18862","cn19160")
copynumbers <- log10(copynumbers)
boxplot(copynumbers, ylab = "log10(copynumber)")#looks better than cns1

#corrected data
keepers <- grep("copy.", names(Brettcorrected), value = T, ignore.case = T)
keepers <- keepers[c(1,4,7)]# only copynumbers kept
copynumbers2 <- Brettcorrected[,names(Brettcorrected) %in% keepers]
names(copynumbers2) <- c("cn19160corrected","cn18862corrected","cn18486corrected")
copynumbers2 <- log10(copynumbers2)
boxplot(copynumbers2, ylab = "log10(copynumber)")

#combine the two.. these look the same
copynumbersall <- cbind(copynumbers,copynumbers2)
boxplot(copynumbersall)

# The spearman correlation across samples is good but less than the proteome samples as expected. Scatterplots reveal (as expected) most of variation at low intensity end of distribution.

require(gplots)
require(Hmisc)
pairs(copynumbersall)
cn.corr <- rcorr(as.matrix(copynumbersall), type = "spearman")
heatmap.2(
  cn.corr$r,
  key = FALSE,
  dendrogram = "none",
  trace = "none",
  cellnote = round(cn.corr$r,3),
  notecex = 1.25,
  notecol = "black",
  col=bluered(10),
  cexCol=1.254,
  cexRow=1.254,
  margins = c(10,10)
)
#because it doesn't seem to matter I will be using the 'non corrected' values from here on out.

#note this SO article for orienting heatmap using 'layout' http://stackoverflow.com/questions/16683026/r-gplots-removing-white-space-in-heatmap-2-when-key-false


##make a 2D density plot of the median copynumber and look for 1D enrichment of the clusters? Otherwise I will have to categorize/threshold the copynumbers, and this is somewhat arbitrary. I would then have to test the enrichment values by bootstrapping. Easiest way to do this is with the S curve.

#make an S curve of rank vs median(logcn). subset this to those values subject to varcomp with protein normalization. Overlay with color the four quadrants. Perform the same analysis with omnibus enrichment values. 

#get median copynumber and median rank, calculate relative rank.
cnMedian <- apply(as.matrix(copynumbers), 1, median)
Rank <- rank(1/cnMedian)
RelRank <- sapply(Rank,function(x) x/max(Rank))
#combine to create DF
MedCNDF <- data.frame(copynumber = cnMedian, Rank = Rank, RelRank = RelRank)
row.names(MedCNDF) <- Brett$Majority.protein.IDs
par(mfrow = c(1,1))
plot(MedCNDF$RelRank, MedCNDF$copynumber)
dim(MedCNDF)#3035

#now subset and look for bias in cns for the varcomp quadrants

#there is no protein identifier from the confounded protein groups data frame added to multExpanded1! I am going to remedy this tomorrow!

##subset copynumber DF st only those subjected to NRE model are considered
proteinvarcompNames <- multExpanded1_withDE[multExpanded1_withDE$SubtoVarcomp == "+", "MajorityProteinIDs"]
MedCNDFvc <- MedCNDF[row.names(MedCNDF) %in% proteinvarcompNames, ]#1152

##VarcompProt formatting fix
VarcompProt$individual <- as.numeric(as.vector(VarcompProt$individual))
VarcompProt$biorep <- as.numeric(as.vector(VarcompProt$biorep))
VarcompProt$residual <- as.numeric(as.vector(VarcompProt$residual))

#add ppMajority IDs and copynumber information to VarcompProt
VarcompProt$ppMajorityIDs <- multExpanded1_withDE[multExpanded1_withDE$idmult %in% row.names(VarcompProt), "ppMajorityProteinIDs"]

##merge by medCNDFvc rownames
#rowNames <- row.names(VarcompProt)
VarcompProtMerged <- merge(VarcompProt,MedCNDFvc, by.x = "ppMajorityIDs", by.y = "row.names")
# row.names(VarcompProt) <- rowNames#lost rownames somehow. A few missing here.

#using red gray blue color gradient to represent copynumber relative rank
rgbpal <- colorRampPalette(c("Red","Gray","Blue"))

#plot final choice. RGB using relcnrank using 1 row and two columns
par(mfcol = c(1,2))
index <- length(unique(VarcompProtMerged$ppMajorityIDs))
plot(VarcompProtMerged$RelRank, VarcompProtMerged$copynumber, 
     col = rgbpal(index)[findInterval(VarcompProtMerged$RelRank,
                                      seq(range(VarcompProtMerged$RelRank)[1], 
                                          range(VarcompProtMerged$RelRank)[2], 
                                          length.out = index), rightmost.closed = T)], 
     pch = 19, ylab = "log10(copynumber)", xlab = "Relative Rank")

#quadrant plot with copynumber RelRank heatmap overlay
plot(log10(VarcompProtMerged$individual),log10(VarcompProtMerged$biorep),
     col = rgbpal(index)[findInterval(VarcompProtMerged$RelRank,
                                      seq(range(VarcompProtMerged$RelRank)[1], 
                                          range(VarcompProtMerged$RelRank)[2], 
                                          length.out = index), rightmost.closed = T)], 
     pch = 19, ylab = "log10 BioRep VarComp", xlab = "log10 Ind VarComp")

##################
##correlation between phosphoworkup copynumber estimates and protein workup copynumber estimates. 
keepers <- grep("copy.", names(Zia), value = T, ignore.case = T)
keepers <- keepers[c(1,4,7)]# only copynumbers kept
copynumbers <- Zia[,names(Zia) %in% keepers]
names(copynumbers) <- c("cn18486","cn18862","cn19160")
copynumbers <- log10(copynumbers)
row.names(copynumbers) <- Zia$Majority.protein.IDs
boxplot(copynumbers, ylab = "log10(copynumber)")#looks better than cns1


keepers <- grep("copy.", names(Brett), value = T, ignore.case = T)
keepers <- keepers[c(1,4,7)]# only copynumbers kept
copynumbersPhos <- Brett[,names(Brett) %in% keepers]
names(copynumbersPhos) <- c("cn18486Phos","cn18862Phos","cn19160Phos")
copynumbersPhos <- log10(copynumbersPhos)
row.names(copynumbersPhos) <- Brett$Majority.protein.IDs
boxplot(copynumbersPhos, ylab = "log10(copynumber)")#looks better than cns1


#subset dataframes
copynumbers <- copynumbers[row.names(copynumbers) %in% row.names(copynumbersPhos),]
copynumbersPhos <- copynumbersPhos[row.names(copynumbersPhos) %in% row.names(copynumbers),]

#merge by rownames
VarcompProtMerged <- merge(VarcompProt,MedCNDFvc, by.x = "ppMajorityIDs", by.y = "row.names")

combined <- merge(copynumbers,copynumbersPhos, by.x = 0, by.y = 0)
row.names(combined) <- combined$Row.names
combined <- combined[,2:length(combined)]
boxplot(combined)

#calculate spearman correlations and do the-damn-thing (it's the hyphens that make it funny). The protein workup values agree more closely with one another across samples. This could be due to an increased number of peptides contributing to the protein group concentration estimates in the protein workup vs the phospho. It also could be due to batch effects obscuring the measurements.  
pairs(combined)
cn.corr <- rcorr(as.matrix(combined), type = "spearman")
heatmap.2(
  cn.corr$r,
  key = TRUE,
  density.info = "none",#no histogram in key
  key.xlab = "Spearman R", key.ylab=NULL, key.title = "",
  dendrogram = "none",
  trace = "none",
  cellnote = round(cn.corr$r,3),
  notecex = 1.25,
  notecol = "black",
  col=rev(rgbpal(20)),
  cexCol=1.254,
  cexRow=1.254,
  margins = c(10,10)
)





