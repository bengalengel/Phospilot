#Do alternative solvers with nlme produce better results?
#Do alternatives to nested random effects design produce biomodal distributions?

#all tests with balanced design on normalized and batch effect corrected confounded phospho data (com3) with BALANCED DESIGN

#basic linear model with compete cases.


# The input matrix is converted into short form for linear modeling
require(plyr)
require(reshape2)
require(qdapRegex)
require(nlme)

ratios <- com3
ratios <- as.matrix(ratios)
ratios <- na.omit(ratios)#balanced data
dim(ratios)
boxplot(ratios)

#melting the dataframe and adding 'individual' 'biorep' and 'techrep' factors works properly.
melted <- melt(ratios, measure.vars = names(ratios))

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


  #all random model with a fixed intercept.technical replicate not explicitly modeled. 
  sites <- c()
  Varcomp <- c()
  Expindvar <- c()
  Expindvarmeans <- c()
  Expbiovar <- c()
  Expbiovarmeans <- c()
  Exptechvar <- c()
  #new loop that gives a model for each unique phosphosite. Fits biorep as a character! also added an error catch
  for(id in levels(melted$Var1)){
    test <- melted[melted$Var1 %in% id,]
    test1 <- test[3:6]
    test1$biorep <- as.factor(test1$biorep)
    pos_err <- tryCatch(lme(value~1, data=test1, random =~1|individual/biorep, na.action = na.exclude),error=function(e) e)
    if(!inherits(pos_err, "error")){
      lmemodel <- lme(value~1, data=test1, random =~1|individual/biorep, na.action = na.exclude)
      temp <- as.numeric(VarCorr(lmemodel)[,1])
      temp <- temp[!is.na(temp)]
      temp <- na.omit(temp)
      Expindvartmp <- 4*temp[1]+2*temp[2]+temp[3]
      Expindvar <- c(Expindvar,Expindvartmp)
      Expindvarmeanstmp <- Expindvartmp/4
      Expindvarmeans <- c(Expindvarmeans,Expindvarmeanstmp)
      Expbiovartmp <- 2*temp[2]+temp[3]
      Expbiovar <- c(Expbiovar,Expbiovartmp)
      Expbiovarmeanstmp <- Expbiovartmp/2
      Expbiovarmeans <- c(Expbiovarmeans,Expbiovarmeanstmp)
      Exptechvartmp <- temp[3]
      Exptechvar <- c(Exptechvar,Exptechvartmp)
      Varcomp <- cbind(Varcomp,temp)
      sites <- c(sites,as.character(unique(test$Var1)))
    }
  }
  
#warning NAs introduced by coercion is not from the model fitting but from retrieving the variance component estimates


  #plot the variance component distributions
  colnames(Varcomp) <- sites
  row.names(Varcomp) <- c("individual","biorep","residual")
  dim(Varcomp)
  Varcomp <- t(Varcomp)
  boxplot(log10(Varcomp), ylab = "log10 variance component")
  summary(Varcomp)
  head(Varcomp)
  
  ##cumulative Varcomp
  Varcomp2 <- Varcomp
  IndCum <- rowSums(Varcomp2)
  Varcomp2[,1] <- 0
  BioCum <- rowSums(Varcomp2)
  Varcomp2[,2] <- 0
  techCum <- rowSums(Varcomp2) #for symmetry!
  CumMat <- cbind(IndCum,BioCum,techCum)
  boxplot(log10(CumMat), ylab = "log10 cumulative variance components")
  
  
  #Expected Variance components
  ExpVars <- cbind(Expindvar,Expbiovar,Exptechvar)
  boxplot(log10(ExpVars), ylab = "log10 Expected variance")
  
  ##Ratio of component estimates relative to technical
  Ind_tech <- Varcomp[,1]/Varcomp[,3]
  Bio_tech <- Varcomp[,2]/Varcomp[,3]
  
  ratios_tech <- cbind(Ind_tech,Bio_tech)
  boxplot(log10(ratios_tech), ylab = "log10 ratios of variance components")
  
  ##scatter plots show discrete groupings?
  plot(log10(Varcomp[,1]),log10(Varcomp[,3]), xlab = "individual variance", ylab = "technical variance")
  plot(log10(Varcomp[,1]),log10(Varcomp[,2]), xlab = "individual variance", ylab = "biological variance")
  plot(log10(Varcomp[,2]),log10(Varcomp[,3]), xlab = "biological variance", ylab = "technical variance")
  
  #histograms of log10 variance are bimodal for biological and individual variance.
  plot(density(log10(Varcomp)), xlab = "variance", main = "hist of variance components")
  plot(density(log10(Varcomp[,1])), xlab = "variance", main = "hist of individual variance components")
  plot(density(log10(Varcomp[,2])), xlab = "variance", main = "hist of biological variance components")
  plot(density(log10(Varcomp[,3])), xlab = "variance", main = "hist of technical variance components")
  


# Do fits to a simple linear model also have bimodal varcomp distributions? --------


#Test one below is incorrect because it does not model the nested structure of the data.
test <- melted[melted$Var1 %in% id,]
test1 <- test[3:6]

#test 1 a straight up linear model

#ok what about a linear model
linear <- lm(value ~ individual + biorep + techrep, data = test1)
aov(linear)
linear <- lm(value ~ individual + biorep, data = test1)
aov(linear)

#lets try fitting it all with a linear model with residual equal to the mass spec workup
sites <- c()
sumsq <- c()
#new loop that gives a model for each unique phosphosite. Fits biorep as a character! also added an error catch
for(id in levels(melted$Var1)){
  test <- melted[melted$Var1 %in% id,]
  test1 <- test[3:6]
  pos_err <- tryCatch(lm(value ~ individual + biorep, data = test1),error=function(e) e)
  if(!inherits(pos_err, "error")){
    linear <- lm(value ~ individual + biorep, data = test1)
    temp <- aov(linear)
    temp <- summary(temp)
    temp <- temp[[1]][,2]
    sumsq <- cbind(sumsq,temp)
    sites <- c(sites,as.character(unique(test$Var1)))
  }
}

#plot the variance component distributions
colnames(sumsq) <- sites
row.names(sumsq) <- c("individual","biorep","residual")
dim(sumsq)
sumsq <- t(sumsq)
boxplot(log10(sumsq), ylab = "log10 sumsq")
summary(sumsq)
head(sumsq)

##scatter plots show discrete groupings?
plot(log10(sumsq[,1]),log10(sumsq[,3]), xlab = "individual variance", ylab = "technical variance")
plot(log10(sumsq[,1]),log10(sumsq[,2]), xlab = "individual variance", ylab = "biological variance")
plot(log10(sumsq[,2]),log10(sumsq[,3]), xlab = "biological variance", ylab = "technical variance")

#histograms of log10 variance are bimodal for biological and individual variance.
plot(density(log10(sumsq)), xlab = "variance", main = "hist of variance components")
plot(density(log10(sumsq[,1])), xlab = "variance", main = "hist of individual variance components")
plot(density(log10(sumsq[,2])), xlab = "variance", main = "hist of biological variance components")
plot(density(log10(sumsq[,3])), xlab = "variance", main = "hist of technical variance components")

#The above does not take into account the correlated nature of the replicates (nested design)

#test 2. a nested design. note that using error doesn't make a difference on the estimates it seems. That is error(ind/bio/tech) gives the same df/ms estimates as ind/bio/tech without specifiying error. However it will not calculate F-ratios 'properly' (the residual MS will always be the denominator).
# a nested design is what is appropriate

#example script using method of moments to estimate variance components
#first note that the following give equivalent estimates of SS/MS
t <- aov(value ~ Error(individual/biorep/techrep), data=test1)
summary(t)
t <- aov(value ~ individual/biorep/techrep, data=test1)
summary(t)
#therefore the latter is used to estimate variance components
t <- aov(value ~ individual/biorep/techrep, data=test1)
t <- summary(t)
IndVarComp <- (t[[1]][1,3] - t[[1]][2,3])/4
BioVarComp <- (t[[1]][2,3] - t[[1]][3,3])/2
TechVarComp <- t[[1]][3,3]

#all 'balanced' data run with this variance component estimation approach. Note there are some negative estimates
sites <- c()
IndVarComp <- c()
BioVarComp <- c()
TechVarComp <- c()

for(id in levels(melted$Var1)){
  test <- melted[melted$Var1 %in% id,]
  test1 <- test[3:6]
  pos_err <- tryCatch(aov(value ~ individual/biorep/techrep, data=test1),error=function(e) e)
  if(!inherits(pos_err, "error")){
    linear <- aov(value ~ individual/biorep/techrep, data=test1)
    linear <- summary(linear)
    IndVarComptmp <- (linear[[1]][1,3] - linear[[1]][2,3])/4
    BioVarComptmp <- (linear[[1]][2,3] - linear[[1]][3,3])/2
    TechVarComptmp <- linear[[1]][3,3]
    IndVarComp <- c(IndVarComp, IndVarComptmp)
    BioVarComp <- c(BioVarComp, BioVarComptmp)
    TechVarComp <- c(TechVarComp, TechVarComptmp)
    sites <- c(sites,as.character(unique(test$Var1)))
  }
}

# Variance component estimates through 'aov' function and method of moments. Note that negative estimates are found with the mse is greater for the lower level relative to the higher level. These are removed by logging the matrix (produces the 'NaNs produced' warning messages)

mmVarcomp <- cbind(IndVarComp,BioVarComp,TechVarComp)
row.names(mmVarcomp) <- sites
boxplot(log10(mmVarcomp), ylab = "log10 mmVarcomp")

##scatter plots no longer show discrete groupings. quadrant plot is no longer present
plot(log10(mmVarcomp[,1]),log10(mmVarcomp[,3]), xlab = "individual variance", ylab = "technical variance")
plot(log10(mmVarcomp[,1]),log10(mmVarcomp[,2]), xlab = "individual variance", ylab = "biological variance")
plot(log10(mmVarcomp[,2]),log10(mmVarcomp[,3]), xlab = "biological variance", ylab = "technical variance")


#histograms of log10 variance are no longer bimodal
plot(density(log10(mmVarcomp), na.rm = T), xlab = "variance", main = "hist of variance components")
plot(density(log10(mmVarcomp[,1]), na.rm = T), xlab = "variance", main = "hist of individual variance components")
plot(density(log10(mmVarcomp[,2]), na.rm = T), xlab = "variance", main = "hist of biological variance components")
plot(density(log10(mmVarcomp[,3]), na.rm = T), xlab = "variance", main = "hist of technical variance components")



###The distributions fitted via this ANOVA approach match the nested random effet model fits sans the ind/bio components with very low values 
par(mfrow = c(2,3))
plot(density(log10(mmVarcomp[,1]), na.rm = T), xlab = "variance", main = "hist of ANOVA derived individual variance components")
plot(density(log10(mmVarcomp[,2]), na.rm = T), xlab = "variance", main = "hist of ANOVA derived biological variance components")
plot(density(log10(mmVarcomp[,3]), na.rm = T), xlab = "variance", main = "hist of ANOVA derived technical variance components")
plot(density(log10(Varcomp[,1])), xlab = "variance", main = "hist of lme derived individual variance components")
plot(density(log10(Varcomp[,2])), xlab = "variance", main = "hist of lme derived biological variance components")
plot(density(log10(Varcomp[,3])), xlab = "variance", main = "hist of lme derived technical variance components")


#### The standardized variance components display roughly the same trend as the lme model?

#replace negative varcomp estimates with zero (tis a bit convoluted, sorry future self but I'm hungry)
test <- log10(mmVarcomp)
test[is.na(test)] <- 0
test2 <- 10^test
test2[test2==1] <- 0

mmvarprop = test2/rowSums(test2)
mmvarprop[mmvarprop==0] <- NA

labs = c("individual","biorep","tech")
boxplot(mmvarprop,axes=F)
axis(1,at=c(1,2,3),labels=labs,col="white");axis(2)

#Central tendencies remain the same for for complete cases (no NAs) however the bio/tech quartiles are affected.
complete <- na.omit(mmvarprop)
boxplot(complete,axes=F)
axis(1,at=c(1,2,3),labels=labs,col="white");axis(2)

# Heatmap representation of the standardized VCs. 
require(gplots)
require(RColorBrewer)
colnames(mmvarprop) = c("individual","bio","tech")
heatmap.2(as.matrix(mmvarprop),
          col=brewer.pal(9,"YlGnBu"),
          Colv=F,
          labRow="",
          trace="none",
          srtCol=45,  ,adjCol = c(1,1),
          margins = c(6,5),
          cexCol=1.5,
          key.xlab = "Standardized VC", key.ylab=NULL, key.title = "",
)

#complete cases
colnames(complete) = c("individual","bio","tech")
heatmap.2(as.matrix(complete),
          col=brewer.pal(9,"YlGnBu"),
          Colv=F,
          labRow="",
          trace="none",
          srtCol=45,  ,adjCol = c(1,1),
          margins = c(6,5),
          cexCol=1.5,
          key.xlab = "Standardized VC", key.ylab=NULL, key.title = "",
)



# Are the low valued variance components (RMLE via lme model) assigned a negative value in the ANOVA approach (method of moments approach) --------

#ind/bio var categorization by REML estimation (lme)
Varcomp2 <- as.data.frame(Varcomp)
Varcomp2$high_ind_var <- ifelse(log10(Varcomp[,1]) >= -5, "TRUE", "FALSE")
Varcomp2$low_ind_var <- ifelse(log10(Varcomp[,1]) < -5, "TRUE", "FALSE")
Varcomp2$high_bio_var <- ifelse(log10(Varcomp[,2]) >= -6, "TRUE", "FALSE")
Varcomp2$low_bio_var <- ifelse(log10(Varcomp[,2]) < -6, "TRUE", "FALSE")

#fnd the low variance sites
RMLELowIndVar <- row.names(Varcomp2[Varcomp2$low_ind_var == "TRUE",])
RMLELowBioVar <- row.names(Varcomp2[Varcomp2$low_bio_var == "TRUE",])
#either
RMLELowVar <- c(RMLELowIndVar,RMLELowBioVar)
RMLELowVar <- unique(RMLELowVar)#681

#what are the zero variance sits by ANOVA method of moments?
tmp <- log10(mmVarcomp)
zeros <- tmp[rowSums(is.na(tmp[ , 1:3])) >= 1,]
zeros <- row.names(zeros)

#how many of these zeros are in the RMLE low variance category?
length(intersect(zeros,RMLELowVar))
length(zeros)
#ALL of them!! the only one missing was because the estimate did not converge for this site in lme






