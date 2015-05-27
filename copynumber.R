####################perseus
#output data matrix for perseus absolute quantification
PerseusOut(directory = "D:/EnsemblDBProteome/txt/", dataset = "Zia")
PerseusOut(directory = "D:/EnsemblDBPhospho/PilotPhosphoensemblDB/combined/txt/", dataset = "Brett")

#input perseus estimates (some info about perseus plugin)

#used the 'raw' light intensity values from MQ for each of the three Zia samples.

#USED Same normalization for all columns: The plugin will calculate the total protein mass per cell for all columns and then use the median thereof as scaling factor for all columns. (alternatively I could quantile normalize the raw output)

#plugin assumes linearity between the Intensity values and the cumulative mass of each protein. In other words, the molecular mass of each protein serves as the correction factor if one wants to calculate copy numbers. HERE I USED THE MEDIAN MW OF EACH OF THE MAJORITY PROTEINS. 'AVERAGE MOLECULR MASS' THE MEDIAN THEORETICAL PEPTIDES IS ALSO REPORTED. 

#NOTE THIS IS DIFFERENT THAN HOW MQ DOES IBAQ. IT USES MASS FROM LEADING PROTEIN/NUMBER OF THEORETICAL PEPTIDES FROM LEADING PROTEIN AS WELL

# The total protein amount per cell was used as a scaling factor. In this case the total intensity will be considered proportional to the total protein amount per cell.

#input and explore copynumber data
Zianonorm <- read.table(file="./Perseus/TPcolumnind.txt", sep = "\t", header=T, fill = T, quote = "")#column independent (no normalization)
keepers <- grep("copy.", names(Zianonorm), value = T, ignore.case = T)
copynumbersall <- Zianonorm[,names(Zianonorm) %in% keepers]
copynumbers <- copynumbersall[,c(1,4,7)]
names(copynumbers) <- c("cn19160","cn18862","cn18486")
copynumbers <- log10(copynumbers)
boxplot(copynumbers, ylab = "log10(copynumber)")#looks better than cns1


# The spearman correlation across samples is good. Scatterplots reveal (as expected) most of variation at low intensity end of distribution.

require(gplots)
require(Hmisc)
pairs(copynumbers)
cn.corr <- rcorr(as.matrix(copynumbers), type = "spearman")
heatmap.2(
  cn.corr$r,
  key = FALSE,
  dendrogram = "none",
  trace = "none",
  cellnote = round(cn.corr$r,2),
  notecex = 1.5,
  notecol = "black",
  col=bluered(25),
  cexCol=1.5,
  cexRow=1.5,
  margins = c(7,7)
)

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



# 
# #using red to yellow "heatmap style colors" with absolute copynumber estimates
# plot(MedCNDFvc$RelRank, MedCNDFvc$copynumber, 
#      col = rev(heat.colors(1152))[findInterval(MedCNDFvc$copynumber,seq(range(MedCNDFvc$copynumber)[1], range(MedCNDFvc$copynumber)[2], length.out = 1152), rightmost.closed = T)], pch = 19)#could get fancy with the range dealy
# 
# #using red to yellow "heatmap style colors" with RelRank copynumber estimates
# plot(MedCNDFvc$RelRank, MedCNDFvc$copynumber, 
#      col = heat.colors(1152)[findInterval(MedCNDFvc$RelRank,seq(range(MedCNDFvc$RelRank)[1], range(MedCNDFvc$RelRank)[2], length.out = 1152), rightmost.closed = T)], pch = 19)#could get fancy with the range dealy




#using red gray blue
rgbpal <- colorRampPalette(c("Red","Gray","Blue"))

plot(MedCNDFvc$RelRank, MedCNDFvc$copynumber, 
     col = rev(rgb(1000))[findInterval(MedCNDFvc$copynumber,seq(range(MedCNDFvc$copynumber)[1], range(MedCNDFvc$copynumber)[2], length.out = 1000), rightmost.closed = T)], pch = 19)#could get fancy with the range dealy


#using red to yellow "heatmap style colors" with RelRank copynumber estimates
plot(MedCNDFvc$RelRank, MedCNDFvc$copynumber, 
     col = rgb(1152)[findInterval(MedCNDFvc$RelRank,seq(range(MedCNDFvc$RelRank)[1], range(MedCNDFvc$RelRank)[2], length.out = 1152), rightmost.closed = T)], pch = 19)#could get fancy with the range dealy




plot(MedCNDFvc$RelRank, MedCNDFvc$copynumber, 
     col = rev(rgb(1000))[findInterval(MedCNDFvc$copynumber,seq(range(MedCNDFvc$copynumber)[1], range(MedCNDFvc$copynumber)[2], length.out = 1000), rightmost.closed = T)], pch = 19)#could get fancy with the range dealy



col=bluefunc(5)[findInterval(dat$value, seq(2:9))] )

#plot with a color gradient. This color gradient will then be overlaid on the quadrant plot. 

VarcompProt$individual <- as.numeric(as.vector(VarcompProt$individual))
VarcompProt$biorep <- as.numeric(as.vector(VarcompProt$biorep))
VarcompProt$residual <- as.numeric(as.vector(VarcompProt$residual))

#add ppMajority IDs and copynumber information to VarcompProt
VarcompProt$ppMajorityIDs <- multExpanded1_withDE[multExpanded1_withDE$idmult %in% row.names(VarcompProt), "ppMajorityProteinIDs"]

##merge by medCNDFvc rownames
rowNames <- row.names(VarcompProt)
VarcompProt <- merge(VarcompProt,MedCNDFvc, by.x = "ppMajorityIDs", by.y = "row.names")
row.names(VarcompProt) <- rowNames#lost rownames somehow. A few missing here.



#quadrant plot with copynumber heatmap overlay
plot(log10(VarcompProt$individual),log10(VarcompProt$biorep),
     col = rev(heat.colors(1152))[findInterval(VarcompProt$copynumber,
                                               seq(range(VarcompProt$copynumber)[1], 
                                                   range(VarcompProt$copynumber)[2], 
                                                   length.out = 1152), rightmost.closed = T)], pch = 19)#could get fancy with the range dealy


#quadrant plot with copynumber RelRank heatmap overlay
plot(log10(VarcompProt$individual),log10(VarcompProt$biorep),
     col = heat.colors(1152)[findInterval(VarcompProt$RelRank,
                                               seq(range(VarcompProt$RelRank)[1], 
                                                   range(VarcompProt$RelRank)[2], 
                                                   length.out = 1152), rightmost.closed = T)], pch = 19)#could get fancy with the range dealy


#same using rgb
plot(log10(VarcompProt$individual),log10(VarcompProt$biorep),
     col = rev(rgb(1152))[findInterval(VarcompProt$copynumber,
                                               seq(range(VarcompProt$copynumber)[1], 
                                                   range(VarcompProt$copynumber)[2], 
                                                   length.out = 1152), rightmost.closed = T)], pch = 19)#could get fancy with the range dealy


#quadrant plot with copynumber RelRank heatmap overlay
plot(log10(VarcompProt$individual),log10(VarcompProt$biorep),
     col = rgb(1152)[findInterval(VarcompProt$RelRank,
                                          seq(range(VarcompProt$RelRank)[1], 
                                              range(VarcompProt$RelRank)[2], 
                                              length.out = 1152), rightmost.closed = T)], pch = 19)#could get fancy with the range dealy


#plot final choice. RGB using relcnrank using 1 row and two columns
par(mfcol = c(1,2))

plot(MedCNDFvc$RelRank, MedCNDFvc$copynumber, 
     col = rgb(1152)[findInterval(MedCNDFvc$RelRank,seq(range(MedCNDFvc$RelRank)[1], range(MedCNDFvc$RelRank)[2], length.out = 1152), rightmost.closed = T)], pch = 19, ylab = "log10(copynumber)", xlab = "Relative Rank")

#quadrant plot with copynumber RelRank heatmap overlay
plot(log10(VarcompProt$individual),log10(VarcompProt$biorep),
     col = rgb(1152)[findInterval(VarcompProt$RelRank,
                                  seq(range(VarcompProt$RelRank)[1], 
                                      range(VarcompProt$RelRank)[2], 
                                      length.out = 1152), rightmost.closed = T)], pch = 19,
     ylab = "log10 BioRep VarComp", xlab = "log10 Ind VarComp")





#add varcomp columns and color by quadrant
varcomp <- multExpanded1_withDE[multExpanded1_withDE$ppSubtoVarcomp == "+", c(77,101:104)]

#I need to add the copynumber and rank information to this DF I would think.
test <- merge(MedCNDFvc,varcomp, 
              
              
              ##comparisons between this approach to copy number estimation and ibaq 
              
              ##S plot of Rank v logarithmized copy number with colored overlay of Varcomp 'quandrants' etc to demonstrate bias/lack of bias in the expression profile
              