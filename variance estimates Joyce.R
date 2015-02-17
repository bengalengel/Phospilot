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
