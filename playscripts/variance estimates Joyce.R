# pulling variance component estimates from lme model using only the paired design. Some convergence issues.
library(nlme)

#all random model with a fixed intercept.technical replicate not explicitly modeled. 
sites <- c()
Varcomp <- c()
Expindvar <- c()
Expindvarmeans <- c()
Expbiovar <- c()
Expbiovarmeans <- c()
Exptechvar <- c()

for(i in 1:1562){
  test <- melted[(melted$Var1 %in% melted$Var1[i]),]
  test1 <- test[3:6]
  pos_err <- tryCatch(lme(value~1, data=test1, random =~1|individual/biorep),error=function(e) e)
  if(!inherits(pos_err, "error")){
    lmemodel <- lme(value~1, data=test1, random =~1|individual/biorep)
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

#****OR***********

#new loop that gives a model for each unique phosphosite. Fits biorep as a character! also added an error catch
for(id in levels(melted$Var1)){
  test <- melted[melted$Var1 %in% id,]
  test1 <- test[3:6]
  test1$biorep <- as.factor(test1$biorep)
  pos_err <- tryCatch(lme(value~1, data=test1, random =~1|individual/biorep),error=function(e) e)
  if(!inherits(pos_err, "error")){
    lmemodel <- lme(value~1, data=test1, random =~1|individual/biorep)
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
boxplot(Expindvar,Expbiovar)
ExpVars <- cbind(Expindvar,Expbiovar,Exptechvar)
boxplot(log10(ExpVars), ylab = "log10 Expected variance")

boxplot(Expindvarmeans,Expbiovarmeans)
boxplot(log10(Expindvarmeans),log10(Expbiovarmeans))


##Ratio of component estimates relative to technical
Ind_tech <- Varcomp[,1]/Varcomp[,3]
Bio_tech <- Varcomp[,2]/Varcomp[,3]

ratios_tech <- cbind(Ind_tech,Bio_tech)
ratios_tech_inv <- ratios_tech^-1

boxplot(log10(ratios_tech), ylab = "log10 ratios of variance components")
boxplot(log10(ratios_tech_inv), ylab = "log10 inverse ratios of variance components")

##scatter plots show discrete groupings?
plot(log10(Varcomp[,1]),log10(Varcomp[,3]), xlab = "individual variance", ylab = "technical variance")
plot(log10(Varcomp[,1]),log10(Varcomp[,2]), xlab = "individual variance", ylab = "biological variance")
plot(log10(Varcomp[,2]),log10(Varcomp[,3]), xlab = "biological variance", ylab = "technical variance")

#histograms of log10 variance are bimodal for biological and individual variance.
hist(Varcomp[,1])
hist(Varcomp[,2])
hist(Varcomp[,3])
hist(log10(Varcomp[,1]))
hist(log10(Varcomp[,2]))
hist(log10(Varcomp[,3]))
plot(density(log10(Varcomp)))
plot(density(log10(Varcomp[,1])), main = "hist of individual variance components")
plot(density(log10(Varcomp[,2])),main = "hist of biological variance components")
plot(density(log10(Varcomp[,3])), main = "hist of technical variance components")


#is this true for the stdDevs of the intercepts? yes
stdDevs <- sqrt(Varcomp)
boxplot(stdDevs)
hist(stdDevs[,1])
hist(stdDevs[,2])
hist(stdDevs[,3])
hist(log10(stdDevs[,1]))
hist(log10(stdDevs[,2]))
hist(log10(stdDevs[,3]))
plot(density(log10(stdDevs)))


#does this have something to do with combat or the normalization scheme?...
# quantiled5 is median and quantile normalized with no missing values

require(plyr)
require(reshape2)
library(qdapRegex)


#no missing values version
quantiled6 <- as.matrix(quantiled5)
melted <- melt(quantiled6, measure.vars = names(quantiled6))


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


#after running above it seems that there is the same signature in the non batch corrected data but quantile normalized data.

#now I will run the same model using non-normalized log2 ratio data. Irene thought that any normalization with 0s can result in small variances apprearing...

#no missing values version
data <- as.matrix(data)
melted <- melt(data, measure.vars = names(data))


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

#run workflow above. The same signature is present even in the non-normalized data.

#what about the very presense of the batches? This bimodal distribution is striking. Cutting non normalized data into two batches
data <- as.data.frame(data)
batch1 <- data[,c(1:2,5:6,9:10)]
batch2 <- data[,c(3:4,7:8,11:12)]

batch1 <- na.omit(batch1)
batch2 <- na.omit(batch2)

#no missing values version
batch2 <- as.matrix(batch2)
melted <- melt(batch2, measure.vars = names(batch2))


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



#now fitted to a random effects model with a fixed intercept. Only fitted variable is individual. Does it still have biomodal distribution?
YES!

  #however if we fit individual as both a fixed and a random effect I get unimodal distributions. However the variance components are 
  
  
#all random model with a fixed intercept.technical replicate not explicitly modeled. 
sites <- c()
Varcomp <- c()
Expindvar <- c()
Expindvarmeans <- c()
Expbiovar <- c()
Expbiovarmeans <- c()
Exptechvar <- c()
#ctrl <- lmeControl(opt='optim') with this arguement nothing converges? throws an error
for(i in 1:1562){
  test <- melted[(melted$Var1 %in% melted$Var1[i]),]
  test1 <- test[3:6]
  lmemodel <- lme(value~1, data=test1, random =~1|individual)
  temp <- as.numeric(VarCorr(lmemodel)[,1])
  temp <- temp[!is.na(temp)]
  temp <- na.omit(temp)
  Varcomp <- cbind(Varcomp,temp)
  sites <- c(sites,as.character(unique(test$Var1)))
}
#plot the variance component distributions
colnames(Varcomp) <- sites
row.names(Varcomp) <- c("individual","residual")
dim(Varcomp)
Varcomp <- t(Varcomp)
boxplot(log10(Varcomp), ylab = "log10 variance component")
summary(Varcomp)
head(Varcomp)

##scatter plots show discrete groupings. Yes
plot(log10(Varcomp[,1]),log10(Varcomp[,2]))

#histograms of log10 variance are bimodal for biological and individual variance.
hist(Varcomp[,1])
hist(Varcomp[,2])
hist(log10(Varcomp[,1]))
hist(log10(Varcomp[,2]))
plot(density(log10(Varcomp)))
plot(density(log10(Varcomp[,1])))
plot(density(log10(Varcomp[,2])))

# plot.design(test1)

#OK now how about the non-normalized ratios from MQ? Do they show this bimodal distribution?
#worked up the the data and the non-normalized has the outlier as negative as opposed to positive?
#this is removed

#workflow modified with an error catch. 1560 models converge. The bimodal variance distribution is present, but mostly captured at the individual level as opposed to the biological level. 

sites <- c()
Varcomp <- c()
Expindvar <- c()
Expindvarmeans <- c()
Expbiovar <- c()
Expbiovarmeans <- c()
Exptechvar <- c()
#ctrl <- lmeControl(opt='optim') with this arguement nothing converges? throws an error
for(i in 1:1178){
  #for(i in levels(melted$Var1)){
  test <- melted[(melted$Var1 %in% melted$Var1[i]),]
  test1 <- test[3:6]
  pos_err <- tryCatch(lme(value~1, data=test1, random =~1|individual/biorep),error=function(e) e)
  if(!inherits(pos_err, "error")){
  lmemodel <- lme(value~1, data=test1, random =~1|individual/biorep)
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


#is the multiplicity a factor? restrict analysis to non-normalized, class1 data with single multiplicites
melted <- melted[ grepl("[0-9]+_1", melted$Var1),]

##the structure is still there but there are duplicates. After removing duplicates the structure is still there. batch correction, normalization (mine or MQs), and multiplicity does not seem to impact the bimodal distribution of the variances. 


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
boxplot(Expindvar,Expbiovar)
ExpVars <- cbind(Expindvar,Expbiovar,Exptechvar)
boxplot(log10(ExpVars), ylab = "log10 Expected variance")

boxplot(Expindvarmeans,Expbiovarmeans)
boxplot(log10(Expindvarmeans),log10(Expbiovarmeans))


##Ratio of component estimates relative to technical
Ind_tech <- Varcomp[,1]/Varcomp[,3]
Bio_tech <- Varcomp[,2]/Varcomp[,3]

ratios_tech <- cbind(Ind_tech,Bio_tech)
ratios_tech_inv <- ratios_tech^-1

boxplot(log10(ratios_tech), ylab = "log10 ratios of variance components")
boxplot(log10(ratios_tech_inv), ylab = "log10 inverse ratios of variance components")

##scatter plots show discrete groupings?
plot(log10(Varcomp[,1]),log10(Varcomp[,3]))
plot(log10(Varcomp[,1]),log10(Varcomp[,2]))
plot(log10(Varcomp[,2]),log10(Varcomp[,3]))

#histograms of log10 variance are bimodal for biological and individual variance.
hist(Varcomp[,1])
hist(Varcomp[,2])
hist(Varcomp[,3])
hist(log10(Varcomp[,1]))
hist(log10(Varcomp[,2]))
hist(log10(Varcomp[,3]))
plot(density(log10(Varcomp)))
plot(density(log10(Varcomp[,1])))
plot(density(log10(Varcomp[,2])))
plot(density(log10(Varcomp[,3])))


