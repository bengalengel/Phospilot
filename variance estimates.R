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
colnames(Varcomp) <- sites
row.names(Varcomp) <- c("individual","biorep","residual")
dim(Varcomp)
boxplot(log10(t(Varcomp)), ylab = "log10 variance component")
summary(t(Varcomp))
head(t(Varcomp))

boxplot(Expindvar,Expbiovar)
boxplot(log10(Expindvar),log10(Expbiovar), log10(Exptechvar), ylab = "log10 Expected variance")

boxplot(Expindvarmeans,Expbiovarmeans)
boxplot(log10(Expindvarmeans),log10(Expbiovarmeans))








#random nested with fixed individual. with optim arguement low variances for biorep are gone.
sites <- c()
Varcomp <- c()
ctrl <- lmeControl(opt='optim')
for(i in 1:1562){
  test <- melted[(melted$Var1 %in% melted$Var1[i]),]
  test1 <- test[3:6]
  lmemodel <- lme(value~individual, control = ctrl, data=test1, random =~1|individual/biorep)
  temp <- as.numeric(VarCorr(lmemodel)[,1])
  temp <- temp[!is.na(temp)]
  temp <- na.omit(temp)
  Varcomp <- cbind(Varcomp,temp)
  sites <- c(sites,as.character(unique(test$Var1)))
}
colnames(Varcomp) <- sites
row.names(Varcomp) <- c("individual","biorep","residual")
dim(Varcomp)
boxplot(log10(t(Varcomp)), ylab = "log10 variance component")
summary(t(Varcomp))
head(t(Varcomp))


#random nested with fixed individual.leave optim out
sites <- c()
Varcomp <- c()
ctrl <- lmeControl(opt='optim')
for(i in 1:1562){
  test <- melted[(melted$Var1 %in% melted$Var1[i]),]
  test1 <- test[3:6]
  lmemodel <- lme(value~individual, data=test1, random =~1|individual/biorep)
  temp <- as.numeric(VarCorr(lmemodel)[,1])
  temp <- temp[!is.na(temp)]
  temp <- na.omit(temp)
  Varcomp <- cbind(Varcomp,temp)
  sites <- c(sites,as.character(unique(test$Var1)))
}
colnames(Varcomp) <- sites
row.names(Varcomp) <- c("individual","biorep","residual")
dim(Varcomp)
boxplot(log10(t(Varcomp)))
summary(t(Varcomp))
head(t(Varcomp))



#random nested with fixed individual/biorep. no optim arguement because with it things fail to converge.HIGHER BIOREP VARIANCE IN THIS SITUATION
sites <- c()
Varcomp <- c()
#ctrl <- lmeControl(opt='optim')
for(i in 1:1562){
  test <- melted[(melted$Var1 %in% melted$Var1[i]),]
  test1 <- test[3:6]
  lmemodel <- lme(value~individual/biorep, data=test1, random =~1|individual/biorep)
  temp <- as.numeric(VarCorr(lmemodel)[,1])
  temp <- temp[!is.na(temp)]
  temp <- na.omit(temp)
  Varcomp <- cbind(Varcomp,temp)
  sites <- c(sites,as.character(unique(test$Var1)))
}
colnames(Varcomp) <- sites
row.names(Varcomp) <- c("individual","biorep","residual")
dim(Varcomp)
boxplot(log10(t(Varcomp)))
summary(t(Varcomp))
head(t(Varcomp))




#averages at the different levels to compare globally for if I should use a fixed factor for individual or biorep. measurements at each tech level.
sites <- c()
indmeans <- c()
biomeans <- c()
techvals <- c()
for(i in 1:1562){
  test <- melted[(melted$Var1 %in% melted$Var1[i]),]
  sites <- c(sites,as.character(unique(test$Var1)))
  values <- test$value
  indmeanTmp <- c(mean(values[1:4]),mean(values[5:8]),mean(values[9:12]))
  indmeans <- cbind(indmeans,indmeanTmp)
  biomeanTmp <- c(mean(values[1:2]),mean(values[3:4]),mean(values[5:6]),mean(values[7:8]),mean(values[9:10]),mean(values[11:12]))
  biomeans <- cbind(biomeans,biomeanTmp)
  techvalTmp <- values
  techvals <- cbind(techvals,techvalTmp)
}
colnames(indmeans) <- sites
colnames(biomeans) <- sites
colnames(techvals) <- sites
boxplot(t(indmeans))
boxplot(t(biomeans))
boxplot(t(techvals))
##my conclusion is that this isn't helpful. everything looks pretty much the same. Perhaps I can do ANOVA on the different levels. 

# LRT tests to decide on the appropriate model at a per peptide level?

lmemodel <- lme(value~individual/biorep, data=test1, random = ~1 | individual/biorep, method = "ML")
lmemodel2 <- lme(value~individual, data=test1, random = ~1 | individual/biorep, method = "ML")
lmemodel3 <- lme(value~1, data=test1, random = ~1 | individual/biorep, method = "ML")
anova.lme(lmemodel3,lmemodel2,lmemodel)##same as anova gives LRT p value

# pairwise
anova(lmemodel3,lmemodel2)

#if better than perform this comparison

#the null hypothesis is that the simpler model is better. the lower the p the more likely you are to reject the null and go with the better model.
#but can my model include both a fixed and random 'individual' term. Here I'm done.

pchisq(logLik(lmemodel3) - logLik(lmemodel), 1)

anova(lmemodel2,lmemodel)

anova(lmemodel,lmemodel3)


anova(lmemodel3,lmemodel2,lmemodel)





















#nesting amounts to adding one main effect and one interaction
aov(value~individual/biorep,data=test1)
anova(lm(value~individual/biorep,data=test1))
# same as
aov(value~individual + individual:biorep, data = test1)

#how about nested random?
aov(value~ Error(individual/biorep), data = test1)
anova(lm(value~ Error(individual/biorep), data = test1))##doesn'twork



