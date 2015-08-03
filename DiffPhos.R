DiffPhos <- function(phosdata, PhosPrep, GelPrep, multExpanded1){
  #this function accepts phospho and protein matrices runs diffphos analysis using limma. Columns are appended
  #to the multexpanded1 file according to the presence of DE or not. As of now combat correction is being performed upstream and protein levels are being fitted as a covariate as opposed to normalizing phosdata and fitting withiout the covariate. It is unclear how normalization would affect fitting biorep as defacto RE.
  
  require(limma)
  require(sva)
  require(statmod)
  require(qdapRegex)
  require(plyr)
  require(reshape2)
  
  ##############################  Confounded data -----
  # First convert the adata matrix to add the meta-data information to the matrix (protein level information can be added here later as well). I should use a melted data format similar to the one used in the random effects model. each row is an observation. Each column specifies the phospho (confounded with protein), individual, 
  
  melted <- melt(phosdata, measure.vars = names(phosdata))#using batch corrected/normalized DF with at least one measurement per biological rep
  
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
  design_base <- model.matrix(~0 + individual, data = SingleCase)#with more explicit values
  
  #calculate the correlation between technical replicates from a series of arrays.
  block = rep(seq(1:6), each = 2)
  dupcor <- duplicateCorrelation(phosdata,design_base,block=block)
  all.correlations <- tanh(dupcor$atanh.correlations)
  boxplot(all.correlations)
  fit <- lmFit(phosdata,design_base,block=block,correlation=dupcor$consensus)
  
  #construct the contrast matrix
  contrast.matrix <- makeContrasts(individual18862-individual18486, individual19160-individual18862, 
                                   individual19160-individual18486, levels = design_base)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  
  #eBayes
  ConfoundedFit <- eBayes(fit2)
  
  #################### PhosPrep as a covariate ---------
  
  #add 'protein' to colnames of PhosPrep
  colnames(PhosPrep) <- paste("Protein",colnames(PhosPrep), sep = "")#4122 obs
    
  #merge phospho and protein estimates. Each has at least one measurement in all bio replicates
  PhosProt <- merge(phosdata, PhosPrep, by = "row.names")#1308 observations
  row.names(PhosProt) <- PhosProt$Row.names
  PhosProt <- PhosProt[,2:length(PhosProt)]
  PhosProt <- as.matrix(PhosProt)
  
  #create the model matrix. Here SingleCase needs to be doubled in size with the last 12 rows corresponding to protein estimates.
  SingleCase2 <- SingleCase
  row.names(SingleCase2) <- colnames(PhosPrep)
  PhosProtSingleCase <- rbind(SingleCase,SingleCase2)
  #now add protein binary flag
  PhosProtSingleCase$Protein <- as.factor(rep(c(0,1), each = 12))
  design_PhosProt <- model.matrix(~0 + individual + Protein, data = PhosProtSingleCase)#with more explicit values
  
  #calculate the duplicate correlation before fitting. DupCor is similar for each molecular phenotype
  dupcorPhos <- dupcor #from confounded fitting
  
  #calculate dupcor for protein
  design_Prot <- model.matrix(~0 + individual, data = SingleCase2)#with more explicit values
  block = rep(seq(1:6), each = 2)
  dupcorProt <- duplicateCorrelation(PhosPrep,design_Prot,block=block)
  dupcorProt$consensus.correlation
  
  #boxplots of genewise correlations for protein and phospho
  PhosCorrelations <- tanh(dupcorPhos$atanh.correlations)
  ProtCorrelations <- tanh(dupcorProt$atanh.correlations)
  boxplot(PhosCorrelations, main = "PhosCorr")
  boxplot(ProtCorrelations, main = "ProtCorr")
  
  #THE AVE OF THESE NUMBERS IS USED FOR THE FITTING
  CorForFitting <- mean(c(dupcorProt$consensus.correlation, dupcorPhos$consensus.correlation))
  
  #Fit
  block = rep(seq(1:12), each = 2)
  fit <- lmFit(PhosProt, design_PhosProt, block=block, correlation=CorForFitting)
  
  #lets try the contrasts using the new design matrix.
  contrast.matrix <- makeContrasts(individual18862-individual18486, individual19160-individual18862, 
                                   individual19160-individual18486, levels = design_PhosProt)
  
  fit2 <- contrasts.fit(fit, contrast.matrix)
  PhosPrepCovFit <- eBayes(fit2)
  
  #################### GelPrep as a covariate!! ---------
  
#   RIGHT NOW THIS IS DONE INCORRECTLY! HOW TO FIT WITH GELPREP AS A COVARIATE AND WITH A DIFFERENT BLOCKING STRUCTURE?
  
  #if the phosprep above with the values removed is the same then just fit with protein estimate regressed and not more design issues!
  #create data matrix with replicated phospho data and single point estimate protein data
  PhosProt2 <- merge(phosdata, GelPrep, by = "row.names")
  row.names(PhosProt2) <- PhosProt2$Row.names
  PhosProt2 <- PhosProt2[,2:16]
  PhosProt2 <- as.matrix(PhosProt2)

#   #The design matrix
#   individual <- as.factor(c(rep("18486", times = 4), rep("18862", times = 4), rep("19160", times = 4), "18862", "18486", "19160"))  
#   Protein <- as.factor(c(rep(0,times = 12), rep(1, times = 3)))
#   SingleCaseGel <- data.frame(individual = individual, Protein = Protein)
#   row.names(SingleCaseGel) <- colnames(PhosProt2)
#   design_GelPrep <- model.matrix(~0 + individual + Protein, data = SingleCaseGel)
# 
#   #Blocking and Fitting 
#   block = c(1,1,2,2,3,3,4,4,5,5,6,6,0,0,0)
#   dupcor <- duplicateCorrelation(PhosProt3, design_GelPrep, block=block)
#   dupcor$consensus.correlation
#   fit <- lmFit(PhosProt2,design_GelPrep,block=block,correlation=dupcor$consensus)
#   contrast.matrix <- makeContrasts(individual18862-individual18486, individual19160-individual18862, 
#                                    individual19160-individual18486, levels = design_GelPrep)
#   fit2 <- contrasts.fit(fit, contrast.matrix)
#   GelPrepCovFit <- eBayes(fit2)

  #An alternative fitting approach with gelprep estimates subtracted away
  
  #normalize each line (using mapply and apply here for each set of replicates to return dataframe of four columns)
  PhosProt2 <- as.data.frame(PhosProt2)
  ind18486Norm <- PhosProt2[,1:4] - PhosProt2$LH18862
  ind18862Norm <- PhosProt2[,5:8] - PhosProt2$LH18486
  ind19160Norm <- PhosProt2[,9:12] - PhosProt2$LH19160
  PhosGelnorm <- cbind(ind18486Norm, ind18862Norm, ind19160Norm)
  PhosGelnorm <- as.matrix(PhosGelnorm)

#calculate the correlation between technical replicates from a series of arrays.
block = rep(seq(1:6), each = 2)
dupcor <- duplicateCorrelation(PhosGelnorm,design_base,block=block)
all.correlations <- tanh(dupcor$atanh.correlations)
boxplot(all.correlations)
fit <- lmFit(PhosGelnorm,design_base,block=block,correlation=dupcor$consensus)#consensus here is .393

#construct the contrast matrix
contrast.matrix <- makeContrasts(individual18862-individual18486, individual19160-individual18862, 
                                 individual19160-individual18486, levels = design_base)
fit2 <- contrasts.fit(fit, contrast.matrix)

#eBayes
GelPrepNormFit <- eBayes(fit2)
  
  ############# Collect the sig hits and annotate ME df --------------

ProcessFit <- function(fit2, header, FitData, multExpanded1){
  #this function accepts the ebays moderated fit data from a given processing choice, returns charts and annotates the ME dataframe.
#   It requires the ebays modifed contrast fits, a header string (such as "confounded") and the matrix passed to limma.
  
  
#Look at pairwise DE using toptable and the coef parameter to id which genes you are interested in 
sig1 <- topTable(fit2, coef = 1, adjust = "BH", n=Inf, sort="p", p=.05)#sorts by adjusted p up to the threshold of .05, which is the default FDR chosen for differential expression ("results" function). This actually seems a conservative way to sort.
sig2 <- topTable(fit2, coef = 2, adjust = "BH", n=Inf, sort="p", p=.05)
sig3 <- topTable(fit2, coef = 3, adjust = "BH", n=Inf, sort="p", p=.05)

# sig1 - 18862-18486
# sig2 - 19160-18862
# sig3 - 19160-18486

c1up  <- sig1[sig1$logFC > 0,]
c1down <- sig1[sig1$logFC < 0,]
c2up <- sig2[sig2$logFC > 0,]
c2down <- sig2[sig2$logFC < 0,]
c3up <- sig3[sig3$logFC > 0,]
c3down <- sig3[sig3$logFC < 0,]


tt1 <- topTable(fit2, coef = 1, adjust = "BH", n=Inf)#sorts by adjusted p up to the threshold of .
tt2 <- topTable(fit2, coef = 2, adjust = "BH", n=Inf)#sorts by adjusted p up to the threshold of .
tt3 <- topTable(fit2, coef = 3, adjust = "BH", n=Inf)#sorts by adjusted p up to the threshold of .

hist(tt1$P.Value, nc=40, xlab="P values", main = colnames(contrast.matrix)[1])
hist(tt2$P.Value, nc=40, xlab="P values", main = colnames(contrast.matrix)[2])
hist(tt3$P.Value, nc=40, xlab="P values", main = colnames(contrast.matrix)[3])

plot(tt1$logFC,-log10(tt1$P.Value), xlab = colnames(contrast.matrix)[1], pch = 20, ylab = "-log10(P)",xlim = c(-5, 5))
#sites with sig difference in comparison 1
names <- row.names(sig1)
names2 <- row.names(tt1)
index <- which(names2 %in% names)
points(tt1$logFC[index],-log10(tt1$P.Value)[index], col="red3", pch = 20)


plot(tt2$logFC,-log10(tt2$P.Value), xlab = colnames(contrast.matrix)[2], pch = 20, ylab = "-log10(P)",xlim = c(-5, 5))
#sites with sig difference in comparison 1
names <- row.names(sig2)
names2 <- row.names(tt2)
index <- which(names2 %in% names)
points(tt2$logFC[index],-log10(tt2$P.Value)[index], col="red3", pch = 20)


plot(tt3$logFC,-log10(tt3$P.Value), xlab = colnames(contrast.matrix)[3], pch = 20, ylab = "-log10(P)",xlim = c(-5, 5))
#sites with sig difference in comparison 1
names <- row.names(sig3)
names2 <- row.names(tt3)
index <- which(names2 %in% names)
points(tt3$logFC[index],-log10(tt3$P.Value)[index], col="red3", pch = 20)


results <- decideTests(fit2, adjust.method = "BH", method = "separate")#results is a 'TestResults' matrix
#separate compares each sample individually and is the default approach
colnames(results) <- c("18862 - 18486", "19160 - 18862", "19160 - 18486")
summary(results)
vennDiagram(results, cex=c(1.2,1,0.7)) #good DE across conditions
vennDiagram(results, cex=c(1.2,1,0.7), include = "up") #good DE across conditions
vennDiagram(results, cex=c(1.2,1,0.7), include = "down") #good DE across conditions
vennDiagram(results, cex=c(1.2,1,0.7), include = c("up", "down")) #good DE across conditions


table("18862-18486" =results[,1],"19160-18862"=results[,2])


#DE sites by contrast type
#DE sites using 'separate' contrasts
DE <- results[results[,1] != 0 | results[,2] != 0 | results[,3] != 0,]
#sites only DE in exactly one contrast
absDE <- abs(DE)
DE1 <- absDE[rowSums(absDE)==1,]
#sites DE in exactly two contrast
DE2 <- absDE[rowSums(absDE)==2,]
#site DE in all 3 contrasts
DE3 <- results[results[,1] != 0 & results[,2] != 0 & results[,3] != 0,]

#F stats by contrast type

# #Sorting F Values
# test <- topTable(fit2, coef = c(1,2,3), adjust = "BH", n=Inf, sort.by="F", p=.05)#equivalent to below
# test2 <- topTableF(fit2, adjust = "BH", n=Inf, sort.by="F", p=.05)

Fvals <- topTableF(fit2, adjust = "BH", n=Inf, sort.by="F")#all F values
sigFvals <- topTableF(fit2, adjust = "BH", n=Inf, sort.by="F", p=.05)#gives 1355 compared to 1549 DE total for separate comparisons

#subsets by contrast specific DE
FDE1 <- Fvals[which(row.names(DE1)%in%row.names(Fvals)),5:6]
FDE2 <- Fvals[which(row.names(DE2)%in%row.names(Fvals)),5:6]
FDE3 <- Fvals[which(row.names(DE3)%in%row.names(Fvals)),5:6]

#below gives adjusted pvalue
FDE1 <- Fvals[match(row.names(DE1), row.names(Fvals), nomatch = F),5:7]
FDE2 <- Fvals[match(row.names(DE2), row.names(Fvals), nomatch = F),5:7]
FDE3 <- Fvals[match(row.names(DE3), row.names(Fvals), nomatch = F),5:7]

#boxplot(FDE1$F,FDE2$F,FDE3$F)
boxplot(log10(FDE1$F),log10(FDE2$F),log10(FDE3$F))
summary(FDE1$F)
summary(FDE2$F)
summary(FDE3$F)

plot(density(log10(FDE1$F)),xlim = c(0,3), main = "Sig F Stats Cut by Number of DiffPhos Contrasts")
lines(density(log10(FDE2$F)), col = 2)
lines(density(log10(FDE3$F)), col = 3)

#add annotation to multexpanded DF
multExpanded1$SubtoDE = ifelse(multExpanded1$idmult %in% row.names(FitData),"+","-")

#add F test values to the table
multExpanded1$globalFsig = ifelse(multExpanded1$idmult %in% row.names(sigFvals),"+","-")


#add DE to table
multExpanded1$DEcont1 = ifelse(multExpanded1$idmult %in% row.names(sig1),"+","-")
multExpanded1$DEcont2 = ifelse(multExpanded1$idmult %in% row.names(sig2),"+","-")
multExpanded1$DEcont3 = ifelse(multExpanded1$idmult %in% row.names(sig3),"+","-")

#add DE direction to table
multExpanded1$cont1up = ifelse(multExpanded1$idmult %in% row.names(c1up),"+","-")
multExpanded1$cont1down = ifelse(multExpanded1$idmult %in% row.names(c1down),"+","-")
multExpanded1$cont2up = ifelse(multExpanded1$idmult %in% row.names(c2up),"+","-")
multExpanded1$cont2down = ifelse(multExpanded1$idmult %in% row.names(c2down),"+","-")
multExpanded1$cont3up = ifelse(multExpanded1$idmult %in% row.names(c3up),"+","-")
multExpanded1$cont3down = ifelse(multExpanded1$idmult %in% row.names(c3down),"+","-")

#replace the names
propernames <- c(paste(header, "SubtoDE", sep = ""), paste(header, "globalFsig", sep = ""),paste(header, "DEcont1", sep = ""),
                 paste(header, "DEcont2", sep = ""), paste(header, "DEcont3", sep = ""), paste(header, "cont1up", sep = ""),
                 paste(header, "cont1down", sep = ""), paste(header, "cont2up", sep = ""), paste(header, "cont2down", sep = ""),
                 paste(header, "cont3up", sep = ""), paste(header, "cont3down", sep = ""))
names(multExpanded1)[(length(names(multExpanded1))-(length(propernames) - 1)):length(names(multExpanded1))] <- propernames
return(multExpanded1)
}

#process fits
multExpanded1 <- ProcessFit(fit2 = ConfoundedFit, header = "Confounded", FitData = phosdata, multExpanded1)
multExpanded1 <- ProcessFit(fit2 = PhosPrepCovFit, header = "PhosPrepCov", FitData = PhosProt, multExpanded1)
multExpanded1 <- ProcessFit(fit2 = GelPrepNormFit, header = "GelPrepNorm", FitData = PhosGelnorm, multExpanded1)

  ############# Annotate the ME dataframe --------------
  
  #
  # NOTE THE FOLLOWING WHEN COMPARING ALL THE CONTRASTS AT ONCE PG 62 IN LIMMA USER GUIDE - decideTests method global should be used here.
  
  # method="global" is recommended when a set of closely related contrasts are being tested. This
  # method simply appends all the tests together into one long vector of tests, i.e., it treats all the tests
  # as equivalent regardless of which probe or contrast they relate to. An advantage is that the raw
  # p-value cutoff is consistent across all contrasts. For this reason, method="global" is recommended if
  # you want to compare the number of DE genes found for different contrasts, for example interpreting
  # the number of DE genes as representing the strength of the contrast. However users need to be aware
  # that the number of DE genes for any particular contrasts will depend on which other contrasts are
  # tested at the same time. Hence one should include only those contrasts which are closely related to
  # the question at hand. Unnecessary contrasts should be excluded as these would affect the results for
  # the contrasts of interest. Another more theoretical issue is that there is no theorem which proves that
  # adjust.method="BH" in combination with method="global" will correctly control the false discovery
  # rate for combinations of negatively correlated contrasts, however simulations, experience and some
  # theory suggest that the method is safe in practice.
  
 
  
  #write the multExpanded table with DE information
  write.table(multExpanded1,"multExpanded1_withDE1.csv", sep=',',col.names=T, row.names=F)
  return(multExpanded1)
}