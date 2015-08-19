##PhosPrep diffphos. All matrices combat corrected with a common 'concensus correlation' calculated from each block. How similar are they?
1) diffphos with protein as covariate
2) diffphos with proein levels substracted from phospho
3) overlap of results from 1 and the common confounded data

# 1)
# model fits with PhosPrep protein as a covariate

#retrive data with phosphopeptide and phosprep protein quantification in at least one biological replicate representation
str(PhosPrepCombatBio)#combat correccted PhosPrep with at least one obs in each biorep
#add 'protein' to colnames
colnames(PhosPrepCombatBio) <- paste("Protein",colnames(PhosPrepCombatBio), sep = "")#4122 obs

#this will be merged with combat corrected phospho data with at least one measurement in each biorep (adata) #4738 measurements
str(adata)

#merge phospho and protein estimates. Each has at least one measurement in all bio replicates
PhosProt <- merge(adata, PhosPrepCombatBio, by = "row.names")#1308 observations
row.names(PhosProt) <- PhosProt$Row.names
PhosProt <- PhosProt[,2:length(PhosProt)]
PhosProt <- as.matrix(PhosProt)

#create the model matrix... Here I think SingleCase needs to be doubled in size with the last 12 rows corresponding to protein estimates.
SingleCase2 <- SingleCase
row.names(SingleCase2) <- colnames(PhosPrepCombatBio)
PhosProtSingleCase <- rbind(SingleCase,SingleCase2)
#now add protein binary flag
PhosProtSingleCase$Protein <- as.factor(rep(c(0,1), each = 12))

design_base <- model.matrix(~0 + individual + Protein, data = PhosProtSingleCase)#with more explicit values


#calculate the duplicate correlation before batch effect removal. ASSUMES DUPLICATE CORRELATION IS SIMILAR FOR BOTH MOLECULAR PHENOTYPES!!
block = rep(seq(1:12), each = 2)

dupcorCombo <- duplicateCorrelation(PhosProt,design_base,block=block)#note the warning message about too much damping - convergence tol not achievable
dupcorCombo$consensus.correlation #.5012503

#what is the concensus correlation for each molecular phenotype?
# 1st the phospho data
design_base <- model.matrix(~0 + individual, data = SingleCase)#with more explicit values
block = rep(seq(1:6), each = 2)
dupcorPhos <- duplicateCorrelation(adata,design_base,block=block)#note the warning message about too much damping - convergence tol not achievable
dupcorPhos$consensus.correlation
[1] 0.3701672

#2nd is protein data
design_base <- model.matrix(~0 + individual, data = SingleCase2)#with more explicit values
block = rep(seq(1:6), each = 2)
dupcorProt <- duplicateCorrelation(PhosPrepCombatBio,design_base,block=block)#note the warning message about too much damping - convergence tol not achievable
dupcorProt$consensus.correlation
[1] 0.3316078

#boxplots of genewise correlations for protein and phospho
PhosCorrelations <- tanh(dupcorPhos$atanh.correlations)
ProtCorrelations <- tanh(dupcorProt$atanh.correlations)
boxplot(PhosCorrelations, main = "PhosCorr")
boxplot(ProtCorrelations, main = "ProtCorr")

#THE AVE OF THESE NUMBERS IS USED FOR THE FITTING
CorForFitting <- mean(c(.37, .33))


#the separate consensus correlations are very similar! 
# If this is the case, where does the .5 number come from? Is this because I am only looking at common sites? NO the common sites give similar corr estimates

#Molecular phenotype consensus correlation for common sites
ComPhos <- PhosProt[,1:12]
ComProt <- PhosProt[,13:24]

#comphos
design_base <- model.matrix(~0 + individual, data = SingleCase)#with more explicit values
dupcorComPhos <- duplicateCorrelation(ComPhos,design_base,block=block)#note the warning message about too much damping - convergence tol not achievable
dupcorComPhos$consensus.correlation
[1] 0.3826011

#comprot
design_base <- model.matrix(~0 + individual, data = SingleCase2)#with more explicit values
dupcorComProt <- duplicateCorrelation(ComProt,design_base,block=block)#note the warning message about too much damping - convergence tol not achievable
dupcorComProt$consensus.correlation
[1] 0.3093484


###fitting with CorForFitting function and blocking on biological replication
design_base <- model.matrix(~0 + individual + Protein, data = PhosProtSingleCase)#with more explicit values


#calculate the duplicate correlation before batch effect removal. ASSUMES DUPLICATE CORRELATION IS SIMILAR FOR BOTH MOLECULAR PHENOTYPES!!
block = rep(seq(1:12), each = 2)

fit <- lmFit(PhosProt, design_base, block=block, correlation=CorForFitting)

#lets try the contrasts using the new design matrix.
contrast.matrix <- makeContrasts(individual18862-individual18486, individual19160-individual18862, 
                                 individual19160-individual18486, levels = design_base)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
#####################################
3)
#fitting confounded phospho levels from the phosprot matrix to compare novelty
PhosphoCommon <- PhosProt[,1:12]
design_base <- model.matrix(~0 + individual, data = SingleCase)#with more explicit values
block = c(1,1,2,2,3,3,4,4,5,5,6,6)
dupcor <- duplicateCorrelation(PhosphoCommon,design_base,block=block)
all.correlations <- tanh(dupcor$atanh.correlations)
boxplot(all.correlations)
dupcor$consensus.correlation
fit <- lmFit(PhosphoCommon,design_base,block=block,correlation=dupcor$consensus)

#lets try the contrasts now. Also works
contrast.matrix <- makeContrasts(individual18862-individual18486, individual19160-individual18862, 
                                 individual19160-individual18486, levels = design_base)
fit2 <- contrasts.fit(fit, contrast.matrix)

#now for ebays. Also works. I have gone throught the rest of the pipeline and the results are very similar to those from before.
fit2 <- eBayes(fit2)

###########################################

##GelPrep diffphos. All matrices combat corrected
1) diffphos with protein as covariate
2) diffphos with proein levels substracted from phospho
3) overlap of results from 1 and the common confounded data


# 1)
# model fits with GelPrep protein as a covariate. Here I am averaging the techreps and fitting without a consensus correlation


##Combine the 'pilot' phospho data and the GelPrep protein levels in a single dataframe
GelPrep <- multExpanded1[,c("idmult", "LH18862", "LH18486", "LH19160")]
row.names(GelPrep) <- GelPrep$idmult
GelPrep <- GelPrep[,2:4]
GelPrep <- na.omit(GelPrep)#11671
GelPrep <- as.matrix(GelPrep)
PhosProt2 <- merge(pilot,GelPrep, by = "row.names") #3488
row.names(PhosProt2) <- PhosProt2$Row.names
PhosProt2 <- PhosProt2[,2:10]
PhosProt2 <- as.matrix(PhosProt2)

##run the diffphos with averaged bioreps and GelPrep protein as a covariate

#create the design matrix
individual <- as.factor(c("18486", "18486", "18862", "18862", "19160", "19160", "18862", "18486", "19160"))  
Protein <- as.factor(c(rep(0,times = 6),rep(1, times = 3)))
SingleCaseGel <- data.frame(individual = individual, Protein = Protein)
row.names(SingleCaseGel) <- (c("18486_1", "18486_2", "18862_1", "18862_2", "19160_1", "19160_2", "Prot18862", "Prot18486", "Prot19160"))



design_GelPrep <- model.matrix(~0 + individual + Protein, data = SingleCaseGel)

#fit the data
fit <- lmFit(PhosProt2,design_GelPrep)

#lets try the contrasts now. CHECK THIS FOR AN ERROR!!!
contrast.matrix <- makeContrasts(individual18862-individual18486, individual19160-individual18862, 
                                 individual19160-individual18486, levels = design_GelPrep)
fit2 <- contrasts.fit(fit, contrast.matrix)

#now for ebays. Also works. I have gone throught the rest of the pipeline and the results are very similar to those from before.
fit2 <- eBayes(fit2)

##run through the entire worklow to get F stats


3) comparison with the confounded data not using correlation structure



#fitting confounded phospho levels from the phosprot matrix to compare novelty
PhosphoCommonGel <- PhosProt2[,1:6]

individual <- as.factor(c("18486", "18486", "18862", "18862", "19160", "19160"))
SingleCase <- data.frame(individual = individual)            

design_GelPrep2 <- model.matrix(~0 + individual, data = SingleCase)

fit <- lmFit(PhosphoCommonGel, design_GelPrep2)

#lets try the contrasts now. Also works
contrast.matrix <- makeContrasts(individual18862-individual18486, individual19160-individual18862, 
                                 individual19160-individual18486, levels = design_GelPrep2)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

##############

























#retrive data with phosphopeptide and phosprep protein quantification in at least one biological replicate representation
str(PhosPrepCombatBio)#combat correccted PhosPrep with at least one obs in each biorep
#add 'protein' to colnames
colnames(PhosPrepCombatBio) <- paste("Protein",colnames(PhosPrepCombatBio), sep = "")#4122 obs

#this will be merged with combat corrected phospho data with at least one measurement in each biorep (adata) #4738 measurements
str(adata)

#merge phospho and protein estimates. Each has at least one measurement in all bio replicates
PhosProt <- merge(adata, PhosPrepCombatBio, by = "row.names")#1308 observations
row.names(PhosProt) <- PhosProt$Row.names
PhosProt <- PhosProt[,2:length(PhosProt)]
PhosProt <- as.matrix(PhosProt)

#create the model matrix... Here I think SingleCase needs to be doubled in size with the last 12 rows corresponding to protein estimates.
SingleCase2 <- SingleCase
row.names(SingleCase2) <- colnames(PhosPrepCombatBio)
PhosProtSingleCase <- rbind(SingleCase,SingleCase2)
#now add protein binary flag
PhosProtSingleCase$Protein <- as.factor(rep(c(0,1), each = 12))

design_base <- model.matrix(~0 + individual + Protein, data = PhosProtSingleCase)#with more explicit values













