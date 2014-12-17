# pulling variance component estimates from lme model using only the paired design. Some convergence issues.
library(nlme)

#all random model with a fixed intercept.technical replicate not explicitly modeled. 
sites <- c()
Varcomp <- c()
#ctrl <- lmeControl(opt='optim') with this arguement nothing converges? throws an error
for(i in 1:1562){
  test <- melted[(melted$Var1 %in% melted$Var1[i]),]
  test1 <- test[3:6]
  lmemodel <- lme(value~1, data=test1, random =~1|individual/biorep)
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
boxplot(log10(t(Varcomp)))
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

#random nested with fixed individual/biorep. no optim arguement because with it things fail to converge.
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



#all random model with a fixed intercept.technical replicate explicitly modeled. NO DF FOR TECH REP ESTIMATE!
sites <- c()
Varcomp <- c()
ctrl <- lmeControl(opt='optim')
for(i in 1:1562){
  test <- melted[(melted$Var1 %in% melted$Var1[i]),]
  test1 <- test[3:6]
  lmemodel <- lme(value~1, control = ctrl, data=test1, random =~1|individual/biorep/techrep)
  temp <- as.numeric(VarCorr(lmemodel)[,1])
  temp <- temp[!is.na(temp)]
  temp <- na.omit(temp)
  Varcomp <- cbind(Varcomp,temp)
  sites <- c(sites,as.character(unique(test$Var1)))
}
colnames(Varcomp) <- sites
row.names(Varcomp) <- c("individual","biorep","techrep", "residual")
boxplot(log10(t(Varcomp)))
summary(t(Varcomp))








Fvalues <- as.data.frame(Fvalues)
row.names(Fvalues) <- sites
boxplot(log10(Fvalues))
hist(log10(as.matrix(Fvalues)))
plot(density(log10(as.matrix(Fvalues))))
plot(density(as.matrix(Fvalues)))
density(as.matrix(Fvalues))