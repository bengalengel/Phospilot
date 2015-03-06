require(nlme)

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
  pos_err <- tryCatch(lme(value~individual, data=test1, random =~1|biorep),error=function(e) e)
  if(!inherits(pos_err, "error")){
    lmemodel <- lme(value~individual, data=test1, random =~1|biorep)
    temp <- as.numeric(VarCorr(lmemodel)[,1])
    Varcomp <- cbind(Varcomp,temp)
    sites <- c(sites,as.character(unique(test$Var1)))
  }
}
#plot the variance component distributions
colnames(Varcomp) <- sites
row.names(Varcomp) <- c("biorep","residual")
Varcomp <- t(Varcomp)
dim(Varcomp)
boxplot(log10(Varcomp), ylab = "log10 variance component")
summary(Varcomp)
head(Varcomp)

##scatter plots show discrete groupings?
plot(log10(Varcomp[,1]),log10(Varcomp[,2]), xlab = "biological variance", ylab = "technical variance")

#histograms of log10 variance are bimodal for biological and individual variance.
hist(log10(Varcomp[,1]))
hist(log10(Varcomp[,2]))
plot(density(log10(Varcomp[,1])),main = "hist of biological variance components", xlab = "variance")
plot(density(log10(Varcomp[,2])), main = "hist of technical variance components", xlab = "variance")






