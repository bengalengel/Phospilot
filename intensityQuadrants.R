#is there a relationship between peptide intensity (ionization efficiency and concentration dependent) and bimodal distribution of varcomp?
#there really shouldn't be. and THERE IS NOT!

#is there a relationship between median H/L ratio magnitude (perhaps issues with incorrect silac pair assignment are more likely when the )
##############

require(gplots)
require(Hmisc)


#variance component estimates from phosphopeptide data (normalized/batch corrected) without protein normalization. 
Varcomp <- varcomp

high_ind_var <- ifelse(log10(Varcomp[,1]) >= -5, "+", "-")
low_ind_var <- ifelse(log10(Varcomp[,1]) < -5, "+", "-")
high_bio_var <- ifelse(log10(Varcomp[,2]) >= -6, "+", "-")
low_bio_var <- ifelse(log10(Varcomp[,2]) < -6, "+", "-")
plot(log10(Varcomp[,1]),log10(Varcomp[,2]), xlab = "individual variance", ylab = "biological variance")
abline(v = -5)
abline(h = -6)
text(1, 1, sum(high_ind_var == "+" & high_bio_var == "+"), col = "darkred", cex = 2, xpd = T)
text(0,-15, sum(low_bio_var == "+" & high_ind_var == "+"), col = "darkred",cex = 2)
text(-15,0, sum(high_bio_var == "+" & low_ind_var == "+"), col = "darkred",cex = 2)
text(-15,-15, sum(low_bio_var == "+" & low_ind_var == "+"), col = "darkred",cex = 2)

#combine into one dataframe
Varcomp <- as.data.frame(Varcomp)
Varcomp <- cbind(Varcomp,high_ind_var,low_ind_var,high_bio_var,low_bio_var)


#subset MEw_DE to those phospho sub to varcomp with intensity information. intensities are the sum of heavy and light intensities per peptide
intensities <- grep("intensity(.*)[0-9]$", names(multExpanded1_withDE), ignore.case = T, value = T)
ratios <- grep("^HL(.*)[0-9]$", names(multExpanded1_withDE), ignore.case = T, value = T)

STvcmp <- multExpanded1_withDE[multExpanded1_withDE$SubtoVarcomp == "+", c("idmult",intensities,ratios)]

#note that many intensities are not multiplicity specific, they are summed over all multiplicity states.
#only mult 1 will be used for this analysis
 
STvcmpsingle <- STvcmp[grepl("[0-9]+_1", STvcmp$idmult),]

#extract the intensities, log transform (otherwise zeros will fuck things up!!!) and replace the infinities
STvcmpsingle2 <- STvcmpsingle[,2:13]
STvcmpsingle2 <- as.matrix(STvcmpsingle2)
STvcmpsingle2 <- log10(STvcmpsingle2)
STvcmpsingle2[!is.finite(STvcmpsingle2)] <- NA



# median of absolute HL ratio magnitude -----------------------------------

#I will use the batch corrected and quantile normalized values 'adata' dataframe which contains one sample in each biological replicate. First I will examine if there is a bias with the direction of the difference. Then I will examine if there is a bias with the magnitude.

#get median intensity and median rank, calculate relative rank.
RatioMedian <- apply(adata, 1, function(x) median(x, na.rm = T))#some duplicates here!!

Rank <- rank(RatioMedian)#ranks smallest to largest
RelRank <- sapply(Rank,function(x) x/max(Rank))
#combine to create DF
MedRatioDF <- data.frame(Ratio = RatioMedian, Rank = Rank, RelRank = RelRank)
par(mfrow = c(1,1))
plot(MedRatioDF$RelRank, MedRatioDF$Ratio) #Median Ratios span 4 OOM in log2 units.

#merge with Varcomp. 
MedRatioVarcomp <- merge(Varcomp,MedRatioDF, by = "row.names")

#using red gray blue color gradient to represent copynumber relative rank
rgbpal <- colorRampPalette(c("Red","Gray","Blue"))


#plot final choice. RGB using relcnrank using 1 row and two columns. There may be some enrichment for low intensity in the low ind & bio var category
par(mfcol = c(1,2))
index <- length(unique(row.names(MedRatioVarcomp)))
plot(MedRatioVarcomp$RelRank, MedRatioVarcomp$Ratio, 
     col = rgbpal(index)[findInterval(MedRatioVarcomp$RelRank,
                                      seq(range(MedRatioVarcomp$RelRank)[1], 
                                          range(MedRatioVarcomp$RelRank)[2], 
                                          length.out = index), rightmost.closed = T)], 
     pch = 19, ylab = "log10(Median Ratio)", xlab = "Relative Rank")

#quadrant plot with copynumber RelRank heatmap overlay
plot(log10(MedRatioVarcomp$individual),log10(MedRatioVarcomp$biorep),
     col = rgbpal(index)[findInterval(MedRatioVarcomp$RelRank,
                                      seq(range(MedRatioVarcomp$RelRank)[1], 
                                          range(MedRatioVarcomp$RelRank)[2], 
                                          length.out = index), rightmost.closed = T)], 
     pch = 19, ylab = "log10 BioRep VarComp", xlab = "log10 Ind VarComp")



#perform the same analysis using ratio magnitudes

RatioMedian <- abs(RatioMedian)
Rank <- rank(1/RatioMedian)#ranks smallest to largest
RelRank <- sapply(Rank,function(x) x/max(Rank))
#combine to create DF
MedRatioDF <- data.frame(Ratio = RatioMedian, Rank = Rank, RelRank = RelRank)
par(mfrow = c(1,1))
plot(MedRatioDF$RelRank, MedRatioDF$Ratio) #Median Ratios span 4 OOM in log2 units.

#merge with Varcomp. 
MedRatioVarcomp <- merge(Varcomp,MedRatioDF, by = "row.names")

#run the above scripts
par(mfcol = c(1,2))
index <- length(unique(row.names(MedRatioVarcomp)))
plot(MedRatioVarcomp$RelRank, MedRatioVarcomp$Ratio, 
     col = rgbpal(index)[findInterval(MedRatioVarcomp$RelRank,
                                      seq(range(MedRatioVarcomp$RelRank)[1], 
                                          range(MedRatioVarcomp$RelRank)[2], 
                                          length.out = index), rightmost.closed = T)], 
     pch = 19, ylab = "log10(Median Ratio)", xlab = "Relative Rank")

#quadrant plot with copynumber RelRank heatmap overlay
plot(log10(MedRatioVarcomp$individual),log10(MedRatioVarcomp$biorep),
     col = rgbpal(index)[findInterval(MedRatioVarcomp$RelRank,
                                      seq(range(MedRatioVarcomp$RelRank)[1], 
                                          range(MedRatioVarcomp$RelRank)[2], 
                                          length.out = index), rightmost.closed = T)], 
     pch = 19, ylab = "log10 BioRep VarComp", xlab = "log10 Ind VarComp")








#what about the histograms of biological, individual, technical and cumulative? Cut median Ratio into deciles and overlay on histograms
#use the hmisc package 'cut2' function to cut dataframe into deciles and quantiles
MedRatioVarcomp$Ratiodecile <- cut2(MedRatioVarcomp$Ratio, g=10)
MedRatioVarcomp$Ratioquartile <- cut2(MedRatioVarcomp$Ratio, g=4)

#quick redo of above using ggplot2
library(ggplot2)
par(mfcol = c(1,1))

#CONCLUSION
#slightly larger overall variance for larger ratios. This seems to be patitioned into the individual variance category

#using color
qplot(log10(individual), log10(biorep), data = MedRatioVarcomp, color = Ratiodecile)
qplot(log10(individual), log10(biorep), data = MedRatioVarcomp, color = Ratioquartile)

#using facets
qplot(log10(individual), log10(biorep), data = MedRatioVarcomp, facets = .~Ratiodecile)
qplot(log10(individual), log10(biorep), data = MedRatioVarcomp, facets = .~Ratioquartile)

#looking at histograms. It's a double rainbow! (except residual)
qplot(log10(individual), data = MedRatioVarcomp, fill = Ratiodecile)
qplot(log10(individual), data = MedRatioVarcomp, fill = Ratioquartile)
qplot(log10(individual), data = MedRatioVarcomp, geom = "density", color = Ratiodecile)
qplot(log10(individual), data = MedRatioVarcomp, geom = "density", color = Ratioquartile)

qplot(log10(biorep), data = MedRatioVarcomp, fill = Ratiodecile)
qplot(log10(biorep), data = MedRatioVarcomp, fill = Ratioquartile)
qplot(log10(biorep), data = MedRatioVarcomp, geom = "density", color = Ratiodecile)
qplot(log10(biorep), data = MedRatioVarcomp, geom = "density", color = Ratioquartile)

#residual has quite small shift towards larger variance with with lower intensities. This can be expected to produce M/A plots with balooning spread at lower intensities. This was not observed so this effect must be (and clearly seems to be) very small. 
qplot(log10(residual), data = MedRatioVarcomp, fill = Ratiodecile)
qplot(log10(residual), data = MedRatioVarcomp, fill = Ratioquartile)
qplot(log10(residual), data = MedRatioVarcomp, geom = "density", color = Ratiodecile)
qplot(log10(residual), data = MedRatioVarcomp, geom = "density", color = Ratioquartile)

#slightly larger overall variance for larger ratios. This seems to be patitioned into the individual variance category
MedRatioVarcomp$cumulative <- rowSums(MedRatioVarcomp[,2:4])
qplot(log10(cumulative), data = MedRatioVarcomp, fill = Ratiodecile)
qplot(log10(cumulative), data = MedRatioVarcomp, fill = Ratioquartile)
qplot(log10(cumulative), data = MedRatioVarcomp, geom = "density", color = Ratiodecile)
qplot(log10(cumulative), data = MedRatioVarcomp, geom = "density", color = Ratioquartile)















#####
#get median intensity and median rank, calculate relative rank.
intMedian <- apply(STvcmpsingle2, 1, function(x) median(x, na.rm = T))#some duplicates here!!

Rank <- rank(1/intMedian)#rank largest to smallest
RelRank <- sapply(Rank,function(x) x/max(Rank))
#combine to create DF
MedIntDF <- data.frame(Intensity = intMedian, Rank = Rank, RelRank = RelRank)
row.names(MedIntDF) <- STvcmpsingle$idmult
par(mfrow = c(1,1))
plot(MedIntDF$RelRank, MedIntDF$Intensity) #intensities span 5 orders of magnitude

     
     
#merge with Varcomp. 
MedIntVarcomp <- merge(Varcomp,MedIntDF, by = "row.names")

#using red gray blue color gradient to represent copynumber relative rank
rgbpal <- colorRampPalette(c("Red","Gray","Blue"))


#plot final choice. RGB using relcnrank using 1 row and two columns. There may be some enrichment for low intensity in the low ind & bio var category
par(mfcol = c(1,2))
index <- length(unique(row.names(MedIntVarcomp)))
plot(MedIntVarcomp$RelRank, MedIntVarcomp$Intensity, 
     col = rgbpal(index)[findInterval(MedIntVarcomp$RelRank,
                                      seq(range(MedIntVarcomp$RelRank)[1], 
                                          range(MedIntVarcomp$RelRank)[2], 
                                          length.out = index), rightmost.closed = T)], 
     pch = 19, ylab = "log10(Median Intensity)", xlab = "Relative Rank")

#quadrant plot with copynumber RelRank heatmap overlay
plot(log10(MedIntVarcomp$individual),log10(MedIntVarcomp$biorep),
     col = rgbpal(index)[findInterval(MedIntVarcomp$RelRank,
                                      seq(range(MedIntVarcomp$RelRank)[1], 
                                          range(MedIntVarcomp$RelRank)[2], 
                                          length.out = index), rightmost.closed = T)], 
     pch = 19, ylab = "log10 BioRep VarComp", xlab = "log10 Ind VarComp")


#what about the histograms of biological, individual, technical and cumulative? Cut median intensity into deciles and overlay on histograms
#use the hmisc package 'cut2' function to cut dataframe into deciles and quantiles
MedIntVarcomp$Intdecile <- cut2(MedIntVarcomp$Intensity, g=10)
MedIntVarcomp$Intquartile <- cut2(MedIntVarcomp$Intensity, g=4)

#quick redo of above using ggplot2
library(ggplot2)
par(mfcol = c(1,1))

#using color
qplot(log10(individual), log10(biorep), data = MedIntVarcomp, color = Intquartile)
qplot(log10(individual), log10(biorep), data = MedIntVarcomp, color = Intdecile)

#using facets
qplot(log10(individual), log10(biorep), data = MedIntVarcomp, facets = .~Intquartile)
qplot(log10(individual), log10(biorep), data = MedIntVarcomp, facets = .~Intdecile)

#looking at histograms. It's a double rainbow! (except residual)
qplot(log10(individual), data = MedIntVarcomp, fill = Intquartile)
qplot(log10(individual), data = MedIntVarcomp, fill = Intdecile)
qplot(log10(individual), data = MedIntVarcomp, geom = "density", color = Intdecile)
qplot(log10(individual), data = MedIntVarcomp, geom = "density", color = Intquartile)

qplot(log10(biorep), data = MedIntVarcomp, fill = Intquartile)
qplot(log10(biorep), data = MedIntVarcomp, fill = Intdecile)
qplot(log10(biorep), data = MedIntVarcomp, geom = "density", color = Intdecile)
qplot(log10(biorep), data = MedIntVarcomp, geom = "density", color = Intquartile)

#residual has quite small shift towards larger variance with with lower intensities. This can be expected to produce M/A plots with balooning spread at lower intensities. This was not observed so this effect must be (and clearly seems to be) very small. 
qplot(log10(residual), data = MedIntVarcomp, fill = Intquartile)
qplot(log10(residual), data = MedIntVarcomp, fill = Intdecile)
qplot(log10(residual), data = MedIntVarcomp, geom = "density", color = Intdecile)
qplot(log10(residual), data = MedIntVarcomp, geom = "density", color = Intquartile)

#lets look at the cumulative now. One can make out the very slight bias again but all of this seems to have been partitioned into the technical variance. This is a good thing!
MedIntVarcomp$cumulative <- rowSums(MedIntVarcomp[,2:4])
qplot(log10(cumulative), data = MedIntVarcomp, fill = Intquartile)
qplot(log10(cumulative), data = MedIntVarcomp, fill = Intdecile)
qplot(log10(cumulative), data = MedIntVarcomp, geom = "density", color = Intdecile)
qplot(log10(cumulative), data = MedIntVarcomp, geom = "density", color = Intquartile)

#######################################

##ibaq now to see if median relative expression of proteins influences the distribution of phosphopeptides to the quadrants
##ibaq values will be taken from PhosPrep and GelPrep
#################################

#subset MEw_DE to those phospho sub to varcomp with intensity information. intensities are the sum of heavy and light intensities per peptide
ibaq <- grep("ibaq", names(multExpanded1_withDE), ignore.case = T, value = T)

STvcmp <- multExpanded1_withDE[multExpanded1_withDE$SubtoVarcomp == "+", c("idmult","Leading.proteins",ibaq)]




#gelprep
ppnames <- grep("ppiBAQ.L", names(STvcmp), value = T)
GelPrepibaq <- STvcmp[,ppnames]
GelPrepibaq <- as.matrix(GelPrepibaq)
GelPrepibaq <- log10(GelPrepibaq)
GelPrepibaq[!is.finite(GelPrepibaq)] <- NA
row.names(GelPrepibaq) <- STvcmp$idmult

#keep cases where each of the replicates are represented
PhosPrepibaq <- PhosPrepibaq[rowSums(is.na(PhosPrepibaq[, 1:4])) < 4 & rowSums(is.na(PhosPrepibaq[, 5:8])) < 4 &
                               rowSums(is.na(PhosPrepibaq[, 9:12])) < 4, ]


GelPrepibaq <- na.omit(GelPrepibaq)



#note the batch effect signature for phosprep. I may correct for this later.
boxplot(PhosPrepibaq)
boxplot(GelPrepibaq)



# batch correction of phosprepibaq --------------------------------------

#here i will remove duplicates, add the proteins to row.names, batch correct and then remerge with PhosPrepibaq using rownames and leading proteins column.

#extract the phospprep and gelprep ibaq measurements and subset to unique 'leading protein entries'
PhosPrepibaq <- STvcmp[,c(2,4:15)]

PhosPrepibaqUnique <- PhosPrepibaq[!duplicated(PhosPrepibaq$Leading.proteins),]
protnames <- PhosPrepibaqUnique$Leading.proteins
PhosPrepibaqUnique <- PhosPrepibaqUnique[,2:length(PhosPrepibaqUnique)]

# log transform (otherwise zeros will fuck things up!!!) and replace the infinities
PhosPrepibaqUnique <- as.matrix(PhosPrepibaqUnique)
PhosPrepibaqUnique <- log10(PhosPrepibaqUnique)
PhosPrepibaqUnique[!is.finite(PhosPrepibaqUnique)] <- NA
row.names(PhosPrepibaqUnique) <- protnames


##EDA before batch correction reveals as expected a batch effect
cdata <- na.omit(PhosPrepibaqUnique)#559 observations
dataZ <- scale(cdata)##Z-scored column wise

# dendogram using euclidian distance (default) and ward or complete agglomeration
dend.ward<- as.dendrogram(hclust(dist(t(dataZ)),method="ward"))
dend.complete<- as.dendrogram(hclust(dist(t(dataZ))))

ward.o<- order.dendrogram(dend.ward)
complete.o<- order.dendrogram(dend.complete)

plot(dend.complete,ylab="height", main = "Euclidian/Complete")
plot(dend.ward, leaflab = "perpendicular", ylab = "height", main = "Euclidian/Ward")

##normalization and batch correction


quantiledibaq <- normalizeQuantiles(PhosPrepibaqUnique)
colnames(quantiledibaq) <- gsub("PhosPrepProteiniBAQ.", "", colnames(quantiledibaq))

# remove exp obs if not observed two or more times in each batch to ensure a variance measurement
quantiledibaq <- quantiledibaq[rowSums(is.na(quantiledibaq[ , c("H.18486_1_1", "H.18486_1_2", "H.18862_1_1", "H.18862_1_2", "H.19160_1_1", "H.19160_1_2")])) < 5 & rowSums(is.na(quantiledibaq[, c("H.18486_2_1", "H.18486_2_2", "H.18862_2_1", "H.18862_2_2", "H.19160_2_1", "H.19160_2_2")])) < 5,] 



##MUST BE FIXED!! I NEED TO DO THIS WITHOUT DUPLICATES AND THEN INSERT BACK IN          
#batch effect identification and adjustment using swamp/combat*************************
##################
swamp <- as.matrix(cdata)
swamp <- swamp[,1:12]
##### sample annotations (data.frame)
set.seed(50)
o1<-data.frame(Factor1=factor(rep(c("A","A","B","B"),3)),
               Numeric1=rnorm(12),row.names=colnames(swamp))

# PCA analysis
res1<-prince(swamp,o1,top=10,permute=T)
str(res1)
a <- res1$linp#plot p values
b <- res1$linpperm#plot p values for permuted data
prince.plot(prince=res1)

#There is a batch effect associated with the process date.
# I must combat this
##batch adjustment using PhosPrepibaq
swamp <- as.matrix(quantiledibaq)
##### sample annotations (data.frame)
set.seed(50)
o1<-data.frame(Factor1=factor(rep(c("A","A","B","B"),3)),
               Numeric1=rnorm(12),row.names=colnames(swamp))


ccquantiledibaq <- combat(swamp, o1$Factor1, batchcolumn=1) #WORKS AFTER ENSURING AT LEAST TWO IN A BATCH.

par(mfrow = c(1,1))
boxplot(ccquantiledibaq)

# dendrograms and PCA of batch corrected data. STILL not totally corrected
#############
# now for the full dataset n=8560
cdata <- na.omit(ccquantiledibaq)#also used below
prince.plot(prince(cdata,o1,top=10)) #huzzah!

# PCA analysis and dendrograms on corrected data
x <- t(cdata)#samples are the rows of the column matrix
pc <- prcomp(x)#scale = T, center = T) as of now I am not scaling

cols <- as.factor(substr(colnames(cdata), 3, 7))##use 5 digit exp name.
plot(pc$x[, 1], pc$x[, 2], col=as.numeric(cols), main = "PCA", xlab = "PC1", ylab = "PC2")
legend("bottomleft", levels(cols), col = seq(along=levels(cols)), pch = 1)

#dendrograms
dataZ <- scale(cdata)##Z-scored column wise

# dendogram using euclidian distance (default) and ward or complete agglomeration
dend.ward<- as.dendrogram(hclust(dist(t(dataZ)),method="ward"))
dend.complete<- as.dendrogram(hclust(dist(t(dataZ))))

ward.o<- order.dendrogram(dend.ward)
complete.o<- order.dendrogram(dend.complete)

plot(dend.complete,ylab="height", main = "Euclidian/Complete")
plot(dend.ward, leaflab = "perpendicular", ylab = "height", main = "Euclidian/Ward")
############################
# OK while batch correction works on ratios it does not work for ibaq, at least using the subset of proteins corresponding to the subset of peptides subjected to varcomp analysis.

# This time I will calculate the medians and rank based on the unique data frame and then merge to create duplicates. (avoids averaging ties  when ranking medians directly within the phospho data frame)
##################

#get median intensity and median rank, calculate relative rank for gel and protein prep ibaq values.
ibaqMedianPhos <- apply(ccquantiledibaq, 1, function(x) median(x, na.rm = T))#some duplicates here!!
ibaqMedianGel <- apply(GelPrepibaq, 1, function(x) median(x, na.rm = T))#some duplicates here!!

#ranks and combine to create DF
Rank <- rank(1/ibaqMedianPhos)#rank largest to smallest
RelRank <- sapply(Rank,function(x) x/max(Rank))
MedPhosPrepibaq <- data.frame(Ibaq = ibaqMedianPhos, Rank = Rank, RelRank = RelRank)

#ranks and combine to create DF
Rank <- rank(1/ibaqMedianGel)#rank largest to smallest
RelRank <- sapply(Rank,function(x) x/max(Rank))
MedGelPrepibaq <- data.frame(Ibaq = ibaqMedianGel, Rank = Rank, RelRank = RelRank)


#s curve plots
plot(MedPhosPrepibaq$RelRank, MedPhosPrepibaq$Ibaq) #intensities span 5 orders of magnitude
plot(MedGelPrepibaq$RelRank, MedGelPrepibaq$Ibaq) #intensities span 5 orders of magnitude



#merge with Varcomp. 
#first add leading proteins from Subtovcmp to varcomp so I can merge MedPhosPrepibaq with varcomp
test <- ifelse(STvcmp$idmult == row.names(Varcomp), STvcmp$Leading.proteins, "")

MedIbaqPhosPrepVarcomp <- merge(Varcomp,MedPhosPrepibaq, by = "row.names")
MedIbaqGelPrepVarcomp <- merge(Varcomp,MedGelPrepibaq, by = "row.names")


#quadrant plot shows little enrichment for either
par(mfcol = c(1,2))
index <- length(unique(row.names(MedIbaqPhosPrepVarcomp)))
plot(MedIbaqPhosPrepVarcomp$RelRank, MedIbaqPhosPrepVarcomp$Ibaq, 
     col = rgbpal(index)[findInterval(MedIbaqPhosPrepVarcomp$RelRank,
                                      seq(range(MedIbaqPhosPrepVarcomp$RelRank)[1], 
                                          range(MedIbaqPhosPrepVarcomp$RelRank)[2], 
                                          length.out = index), rightmost.closed = T)], 
     pch = 19, ylab = "log10(Median Intensity)", xlab = "Relative Rank")

#quadrant plot with copynumber RelRank heatmap overlay
plot(log10(MedIbaqPhosPrepVarcomp$individual),log10(MedIbaqPhosPrepVarcomp$biorep),
     col = rgbpal(index)[findInterval(MedIbaqPhosPrepVarcomp$RelRank,
                                      seq(range(MedIbaqPhosPrepVarcomp$RelRank)[1], 
                                          range(MedIbaqPhosPrepVarcomp$RelRank)[2], 
                                          length.out = index), rightmost.closed = T)], 
     pch = 19, ylab = "log10 BioRep VarComp", xlab = "log10 Ind VarComp")

index <- length(unique(row.names(MedIbaqGelPrepVarcomp)))
plot(MedIbaqGelPrepVarcomp$RelRank, MedIbaqGelPrepVarcomp$Ibaq, 
     col = rgbpal(index)[findInterval(MedIbaqGelPrepVarcomp$RelRank,
                                      seq(range(MedIbaqGelPrepVarcomp$RelRank)[1], 
                                          range(MedIbaqGelPrepVarcomp$RelRank)[2], 
                                          length.out = index), rightmost.closed = T)], 
     pch = 19, ylab = "log10(Median Intensity)", xlab = "Relative Rank")

#quadrant plot with copynumber RelRank heatmap overlay
plot(log10(MedIbaqGelPrepVarcomp$individual),log10(MedIbaqGelPrepVarcomp$biorep),
     col = rgbpal(index)[findInterval(MedIbaqGelPrepVarcomp$RelRank,
                                      seq(range(MedIbaqGelPrepVarcomp$RelRank)[1], 
                                          range(MedIbaqGelPrepVarcomp$RelRank)[2], 
                                          length.out = index), rightmost.closed = T)], 
     pch = 19, ylab = "log10 BioRep VarComp", xlab = "log10 Ind VarComp")



#what about the histograms of biological, individual, technical and cumulative? Cut median intensity into deciles and overlay on histograms
#use the hmisc package 'cut2' function to cut dataframe into deciles and quantiles
MedIbaqPhosPrepVarcomp$Ibaqdecile <- cut2(MedIbaqPhosPrepVarcomp$Ibaq, g=10)
MedIbaqPhosPrepVarcomp$Ibaqquartile <- cut2(MedIbaqPhosPrepVarcomp$Ibaq, g=4)
MedIbaqGelPrepVarcomp$Ibaqdecile <- cut2(MedIbaqGelPrepVarcomp$Ibaq, g=10)
MedIbaqGelPrepVarcomp$Ibaqquartile <- cut2(MedIbaqGelPrepVarcomp$Ibaq, g=4)

#quick redo of above using ggplot2
library(ggplot2)
par(mfcol = c(1,1))

#using color
qplot(log10(individual), log10(biorep), data = MedIbaqPhosPrepVarcomp, color = Ibaqquartile)
qplot(log10(individual), log10(biorep), data = MedIbaqPhosPrepVarcomp, color = Ibaqdecile)
qplot(log10(individual), log10(biorep), data = MedIbaqGelPrepVarcomp, color = Ibaqquartile)
qplot(log10(individual), log10(biorep), data = MedIbaqGelPrepVarcomp, color = Ibaqdecile)

#using facets
qplot(log10(individual), log10(biorep), data = MedIbaqPhosPrepVarcomp, facets = .~Ibaqquartile)
qplot(log10(individual), log10(biorep), data = MedIbaqPhosPrepVarcomp, facets = .~Ibaqdecile)
qplot(log10(individual), log10(biorep), data = MedIbaqGelPrepVarcomp, facets = .~Ibaqquartile)
qplot(log10(individual), log10(biorep), data = MedIbaqGelPrepVarcomp, facets = .~Ibaqdecile)

#looking at histograms of bio and ind. higher expressed proteins are actually not over represented in higher ind variation as opposed to lower. However the shape of the distributions appears to be different, with a clear trend of a 'tighter' distribution for the highly expressed proteins. There is no difference in the biological histograms. Perhaps these are low stoicheometry sites varying across individuals?
qplot(log10(individual), data = MedIbaqPhosPrepVarcomp, fill = Ibaqquartile)
qplot(log10(individual), data = MedIbaqPhosPrepVarcomp, fill = Ibaqdecile)
qplot(log10(individual), data = MedIbaqPhosPrepVarcomp, geom = "density", color = Ibaqdecile)
qplot(log10(individual), data = MedIbaqPhosPrepVarcomp, geom = "density", color = Ibaqquartile)
qplot(log10(individual), data = MedIbaqGelPrepVarcomp, fill = Ibaqquartile)
qplot(log10(individual), data = MedIbaqGelPrepVarcomp, fill = Ibaqdecile)
qplot(log10(individual), data = MedIbaqGelPrepVarcomp, geom = "density", color = Ibaqdecile)
qplot(log10(individual), data = MedIbaqGelPrepVarcomp, geom = "density", color = Ibaqquartile)

#What fraction of high individual data are represented by each quartile/decile. The change in frequency is muted when looking at decile
table(MedIbaqGelPrepVarcomp$high_ind_var,MedIbaqGelPrepVarcomp$Ibaqquartile)
table(MedIbaqGelPrepVarcomp$high_ind_var,MedIbaqGelPrepVarcomp$Ibaqdecile)
table(MedIbaqPhosPrepVarcomp$high_ind_var,MedIbaqPhosPrepVarcomp$Ibaqquartile)
table(MedIbaqPhosPrepVarcomp$high_ind_var,MedIbaqPhosPrepVarcomp$Ibaqdecile)


qplot(log10(biorep), data = MedIbaqPhosPrepVarcomp, fill = Ibaqquartile)
qplot(log10(biorep), data = MedIbaqPhosPrepVarcomp, fill = Ibaqdecile)
qplot(log10(biorep), data = MedIbaqPhosPrepVarcomp, geom = "density", color = Ibaqdecile)
qplot(log10(biorep), data = MedIbaqPhosPrepVarcomp, geom = "density", color = Ibaqquartile)
qplot(log10(biorep), data = MedIbaqGelPrepVarcomp, fill = Ibaqquartile)
qplot(log10(biorep), data = MedIbaqGelPrepVarcomp, fill = Ibaqdecile)
qplot(log10(biorep), data = MedIbaqGelPrepVarcomp, geom = "density", color = Ibaqdecile)
qplot(log10(biorep), data = MedIbaqGelPrepVarcomp, geom = "density", color = Ibaqquartile)


#hist of residual displays no bias in protein expression level
qplot(log10(residual), data = MedIbaqPhosPrepVarcomp, fill = Ibaqquartile)
qplot(log10(residual), data = MedIbaqPhosPrepVarcomp, fill = Ibaqdecile)
qplot(log10(residual), data = MedIbaqPhosPrepVarcomp, geom = "density", color = Ibaqdecile)
qplot(log10(residual), data = MedIbaqPhosPrepVarcomp, geom = "density", color = Ibaqquartile)
qplot(log10(residual), data = MedIbaqGelPrepVarcomp, fill = Ibaqquartile)
qplot(log10(residual), data = MedIbaqGelPrepVarcomp, fill = Ibaqdecile)
qplot(log10(residual), data = MedIbaqGelPrepVarcomp, geom = "density", color = Ibaqdecile)
qplot(log10(residual), data = MedIbaqGelPrepVarcomp, geom = "density", color = Ibaqquartile)


#lets look at the cumulative now and no bias is present
MedIbaqPhosPrepVarcomp$cumulative <- rowSums(MedIbaqPhosPrepVarcomp[,2:4])
qplot(log10(cumulative), data = MedIbaqPhosPrepVarcomp, fill = Ibaqquartile)
qplot(log10(cumulative), data = MedIbaqPhosPrepVarcomp, fill = Ibaqdecile)
qplot(log10(cumulative), data = MedIbaqPhosPrepVarcomp, geom = "density", color = Ibaqdecile)
qplot(log10(cumulative), data = MedIbaqPhosPrepVarcomp, geom = "density", color = Ibaqquartile)

MedIbaqGelPrepVarcomp$cumulative <- rowSums(MedIbaqGelPrepVarcomp[,2:4])
qplot(log10(cumulative), data = MedIbaqGelPrepVarcomp, fill = Ibaqquartile)
qplot(log10(cumulative), data = MedIbaqGelPrepVarcomp, fill = Ibaqdecile)
qplot(log10(cumulative), data = MedIbaqGelPrepVarcomp, geom = "density", color = Ibaqdecile)
qplot(log10(cumulative), data = MedIbaqGelPrepVarcomp, geom = "density", color = Ibaqquartile)

#using median ibaq values I am assuming normalization and batch effect correction are not having a huge impact. I will address this after the run


     
