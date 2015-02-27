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
#ctrl <- lmeControl(opt='optim') with this arguement nothing converges? throws an error
for(i in 1:1562){
  test <- melted[(melted$Var1 %in% melted$Var1[i]),]
  test1 <- test[3:6]
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






