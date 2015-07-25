
#This function will accept confounded data that is median/quantile normalized (See BMC bioinformatics for an alternative approach (Joyce sent link): Normalize within batches perhaps?). It may work with technical replicates to estimate and correct for batch effects. It is unclear how to incorporate technical replication that contains missing values!!

#Limma will be used and will account for the batch effect as a covariate in the analysis. This will be compared with using combat to correct for the batch effect.


# packages ----------------------------------------------------------------
require(limma)
require(sva)
require(statmod)

#limma help and guides
library(help = limma)
limmaUsersGuide()
?limma #choose chapters from here
# input data --------------------------------------------------------------
quantiled <- CorrectedData[[3]]#median/quantile normalized
quantiledBio <- quantiled[rowSums(is.na(quantiled[ , 1:2])) < 2 & rowSums(is.na(quantiled[ , 3:4])) < 2 & rowSums(is.na(quantiled[ , 5:6])) < 2 
              & rowSums(is.na(quantiled[ , 7:8])) < 2 & rowSums(is.na(quantiled[ , 9:10])) < 2 & rowSums(is.na(quantiled[ , 11:12])) < 2,] #med/quantile with at least one measurement per biological replicate
str(adata)#normalized/combat batch corrected data frame with at lesat 1 obs in each bio rep
str(pilot)#combat corrected median and quantile normalized dataframe with averages for each biological replicate



# Test 1, include technical replicates with missing values. Can I get a fit? -----------------------------------------------


# First convert the adata matrix to add the meta-data information to the matrix (protein level information can be added here later as well). I should use a melted data format similar to the one used in the random effects model. each row is an observation. Each column specifies the phospho (confounded with protein), individual, 

melted <- melt(adata, measure.vars = names(adata))#using batch corrected/normalized DF with at least one measurement per biological rep

#identify individual name and add it to the table
matches <- gregexpr("[0-9]{5}", melted$Var2, perl=T)
individual <- regmatches(melted$Var2,matches)
individual <- as.character(individual)
individual <- as.factor(individual)
melted$individual <- individual

#identify the biological replicate
biorep <- rm_between(melted$Var2, "_", "_", extract=TRUE)
biorep <- as.character(biorep)
biorep <- as.factor(biorep)
melted$biorep <- biorep

#identify the technical replicate
matches <- gregexpr("[0-9]$", melted$Var2, perl=T)
techrep <- regmatches(melted$Var2,matches)
techrep <- as.character(techrep)
techrep <- as.factor(techrep)
melted$techrep <- techrep

#Create a small model matrix of variables of interest and blocking variables (replicates)
SingleCase <- melted[melted$Var1 %in% levels(melted$Var1)[1],]
row.names(SingleCase) <- SingleCase$Var2
SingleCase <- SingleCase[,4:6]

# make the simple model matrix
design_base <- model.matrix(~individual, data = SingleCase)#with constant intercept term
design_base <- model.matrix(~0 + individual, data = SingleCase)#with more explicit values

#calculate the correlation between technical replicates from a series of arrays.

# From dupcor help: 
# If block is not null, this function estimates the correlation between repeated observations on the blocking variable. Typically the blocks are biological replicates and the repeated observations are technical replicates. In either case, the correlation is estimated by fitting a mixed linear model by REML individually for each gene. The function also returns a consensus correlation, which is a robust average of the individual correlations, which can be used as input for functions lmFit or gls.series.

#My case is exactly that described above. 
block = c(1,1,2,2,3,3,4,4,5,5,6,6)
dupcor <- duplicateCorrelation(adata,design_base,block=block)
dupcor$consensus.correlation
fit <- lmFit(adata,design_base,block=block,correlation=dupcor$consensus)

#lets try the contrasts now. Also works
contrast.matrix <- makeContrasts(individual18862-individual18486, individual19160-individual18862, 
                                 individual19160-individual18486, levels = design_base)
fit2 <- contrasts.fit(fit, contrast.matrix)

#now for ebays. Also works. I have gone throught the rest of the pipeline and the results are very similar to those from before.
fit2 <- eBayes(fit2)


# Trying approach with batch modeled explicitly ---------------------------

#calculate the duplicate correlation before combat batch effect removal
quantiledBio <- as.matrix(quantiledBio)
dupcor <- duplicateCorrelation(quantiledBio,design_base,block=block)#note the warning message about too much damping - convergence tol not achievable
dupcor$consensus.correlation #.843

#I need to change the model matrix first. To do this I need to add batch to 'SingleCase'. biorep = batch so simply add that term to the model.
design_batch <- model.matrix(~0 + individual + biorep, data = SingleCase)
dupcor <- duplicateCorrelation(quantiledBio,design_batch,block=block)#the warning message is gone!
dupcor$consensus.correlation# a bit lower than above .608 here, but was half this after combat was applied
fit <- lmFit(quantiledBio,design_batch,block=block,correlation=dupcor$consensus)

#lets try the contrasts using the new design matrix.
contrast.matrix <- makeContrasts(individual18862-individual18486, individual19160-individual18862, 
                                 individual19160-individual18486, levels = design_batch)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
#so how does this compare with the above?....Significantly less DiffPhos!!

#Remove batch effects corrects for the 
batchcorrected <- removeBatchEffect(x = quantiledBio, batch = SingleCase$biorep,
                               design = design_base)
#EDA on batchcorrected
#PCA analysis.  
cdata <- na.omit(batchcorrected)
x <- t(cdata)#samples are the rows of the column matrix
pc <- prcomp(x, retx = T, scale = T, center = T) #With scaling the segregation is less robust
cols <- as.factor(substr(colnames(cdata), 3, 7))##5 digit sample  name.
plot(pc$x[, 1], pc$x[, 2], col=as.numeric(cols), main = "PCA", xlab = "PC1", ylab = "PC2")
legend("bottomleft", levels(cols), col = seq(along=levels(cols)), pch = 1)
#dendrograms
dataZ <- scale(cdata)##Z-scored column wise

# now all data excepting complete cases (note that the sample dendograms look the same)
#hist(dataZ[,6], breaks = 100)

# dendogram using euclidian distance (default) and ward or complete agglomeration
dend.ward<- as.dendrogram(hclust(dist(t(dataZ)),method="ward.D"))
dend.ward2<- as.dendrogram(hclust(dist(t(dataZ)),method="ward.D2"))
dend.complete<- as.dendrogram(hclust(dist(t(dataZ))))

ward.o<- order.dendrogram(dend.ward)
ward2.o<- order.dendrogram(dend.ward2)
complete.o<- order.dendrogram(dend.complete)

plot(dend.complete,ylab="height", main = "Euclidian/Complete")
plot(dend.ward, leaflab = "perpendicular", ylab = "height", main = "Euclidian/Ward")
plot(dend.ward2, leaflab = "perpendicular", ylab = "height", main = "Euclidian/Ward2")

#how many less DiffPhos using limma covariate vs combat batch correction?...(later)

# model fits with PhosPrep protein as a covariate

#retrive data with phosphopeptide and phosprep protein quantification in at least one biological replicate representation
str(quantiledBioProt)
#add 'protein' to colnames
names(quantiledBioProt) <- paste("Protein",names(quantiledBioProt), sep = "")

#merge phospho and protein estimates. Each has at least one measurement in all bio replicates
PhosProt <- merge(quantiledBio, quantiledBioProt, by = "row.names")
row.names(PhosProt) <- PhosProt$Row.names
PhosProt <- as.matrix(PhosProt)
PhosProt <- PhosProt[,2:length(PhosProt)]
PhosProt <- as.matrix(PhosProt)

#create the model matrix... Here I think SingleCase needs to be doubled in size with the last 12 rows corresponding to protein estimates.
SingleCase2 <- SingleCase
row.names(SingleCase2) <- names(quantiledBioProt)
PhosProtSingleCase <- rbind(SingleCase,SingleCase2)
#now add protein binary flag
PhosProtSingleCase$Protein <- rep(c(0,1), each = 12)

design_base <- model.matrix(~0 + individual + Protein, data = PhosProtSingleCase)#with more explicit values


#calculate the duplicate correlation before batch effect removal
block = rep(seq(1:12), each = 2)

dupcor <- duplicateCorrelation(PhosProt,design_base,block=block)#note the warning message about too much damping - convergence tol not achievable
dupcor$consensus.correlation #.865 

#This isn't correct. I need to account for 'batch' while calculating the duplicate correlation.

#I need to change the model matrix first. To do this I need to add batch to 'SingleCase'. biorep = batch so simply add that term to the model.
design_batch <- model.matrix(~0 + individual + biorep + Protein, data = PhosProtSingleCase)
dupcor <- duplicateCorrelation(PhosProt,design_batch,block=block)#the warning message is gone!
dupcor$consensus.correlation# here the consensus correlation is .687.
fit <- lmFit(PhosProt, design_batch, block=block, correlation=dupcor$consensus)#Is this right?...I obviously think so.

#lets try the contrasts using the new design matrix.
contrast.matrix <- makeContrasts(individual18862-individual18486, individual19160-individual18862, 
                                 individual19160-individual18486, levels = design_batch)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
#so how does this compare with the combat approach? I don't know because I haven't done it....
# Well what about diffphos. how does it look? Seems reasonable! 

# Now for the batch effect correction. Just to check for proper coding....


#Remove batch effects corrects for the 
batchcorrectedProtcovar <- removeBatchEffect(x = PhosProt, batch = PhosProtSingleCase$biorep,
                                    design = design_base)
#EDA on batchcorrectedProtcovar
batchcorPhos <- batchcorrectedProtcovar[,1:12]
batchcorProt <- batchcorrectedProtcovar[,13:24]
#PCA analysis.  
cdata <- na.omit(batchcorProt)#only 271 complete cases here!!!
x <- t(cdata)#samples are the rows of the column matrix
pc <- prcomp(x, retx = T, scale = T, center = T) #With scaling the segregation is less robust
cols <- as.factor(substr(colnames(cdata), 10, 14))##5 digit sample  name.
plot(pc$x[, 1], pc$x[, 2], col=as.numeric(cols), main = "PCA", xlab = "PC1", ylab = "PC2")
legend("bottomleft", levels(cols), col = seq(along=levels(cols)), pch = 1)

#dendrograms
dataZ <- scale(cdata)##Z-scored column wise

# now all data excepting complete cases (note that the sample dendograms look the same)
#hist(dataZ[,6], breaks = 100)

# dendogram using euclidian distance (default) and ward or complete agglomeration
dend.ward<- as.dendrogram(hclust(dist(t(dataZ)),method="ward.D"))
dend.ward2<- as.dendrogram(hclust(dist(t(dataZ)),method="ward.D2"))
dend.complete<- as.dendrogram(hclust(dist(t(dataZ))))

ward.o<- order.dendrogram(dend.ward)
ward2.o<- order.dendrogram(dend.ward2)
complete.o<- order.dendrogram(dend.complete)

plot(dend.complete,ylab="height", main = "Euclidian/Complete")
plot(dend.ward, leaflab = "perpendicular", ylab = "height", main = "Euclidian/Ward")
plot(dend.ward2, leaflab = "perpendicular", ylab = "height", main = "Euclidian/Ward2")


# I can also check for proper coding by not including protein and seeing if the results change... (dupcor values will be different I think)
#actually I think the above fitting isn't correct because duplicate corr is on two separate "molecular phenotypes". The best approach it seems is to average biological replicates and fit.

#PhosProt with biorep averaging
HL18486_1 <- rowMeans(PhosProt[,1:2], na.rm = T)
HL18486_2 <- rowMeans(PhosProt[,3:4], na.rm = T)
HL18862_1 <- rowMeans(PhosProt[,5:6], na.rm = T)
HL18862_2 <- rowMeans(PhosProt[,7:8], na.rm = T)
HL19160_1 <- rowMeans(PhosProt[,9:10], na.rm = T)
HL19160_2 <- rowMeans(PhosProt[,11:12], na.rm = T)
ProtHL18486_1 <- rowMeans(PhosProt[,13:14], na.rm = T)
ProtHL18486_2 <- rowMeans(PhosProt[,15:16], na.rm = T)
ProtHL18862_1 <- rowMeans(PhosProt[,17:18], na.rm = T)
ProtHL18862_2 <- rowMeans(PhosProt[,19:20], na.rm = T)
ProtHL19160_1 <- rowMeans(PhosProt[,21:22], na.rm = T)
ProtHL19160_2 <- rowMeans(PhosProt[,23:24], na.rm = T)


PhosProtBio <- cbind(HL18486_1, HL18486_2, HL18862_1, HL18862_2, HL19160_1, HL19160_2,
                     ProtHL18486_1, ProtHL18486_2, ProtHL18862_1, ProtHL18862_2, ProtHL19160_1, ProtHL19160_2)#1308 class 1 measurements with at least one quant in each biologic

#a new singlecase/design matrix
matches <- gregexpr("[0-9]{5}", colnames(PhosProtBio), perl=T)
individual <- regmatches(colnames(PhosProtBio),matches)
individual <- as.character(individual)
individual <- as.factor(individual)
SingleCase2 <- data.frame(individual = individual, biorep = rep(c(1,2), 6))
SingleCase2$Protein <- rep(c(0,1), each = 6)
row.names(SingleCase2) <- colnames(PhosProtBio)

design <- model.matrix(~0 + individual + biorep + Protein, data = SingleCase2)#design matrix with batch explicitly modeled and protein as a covariate

fit <- lmFit(PhosProtBio, design)
contrast.matrix <- makeContrasts(individual18862-individual18486, individual19160-individual18862, 
                                 individual19160-individual18486, levels = design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)




