# combat melting for anova working with com2

require(plyr)
require(reshape2)
library(qdapRegex)

melted <- melt(com2, measure.vars = names(com2))

#no missing values version
com3 <- na.omit(com2)
melted <- melt(com3, measure.vars = names(com3))


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


# for each unique var1 (phosphosite) make a list of values/individual/biorep/techrep
# for i in levels for factorized for looping
sites <- c()
biovar <- c()
techvar <- c()
indvar <- c()
for(i in 1:1562){
  test <- melted[(melted$Var1 %in% melted$Var1[i]),]
  sites <- c(sites,as.character(unique(test$Var1)))
  values <- test$value
  indvar <- c(indvar,var(values))
  biovarTmp <- c(var(values[1:4]),var(values[5:8]),var(values[9:12]))
  biovar <- c(biovar,mean(biovarTmp))
  techvarTmp <- c(var(values[1:2]),var(values[3:4]),var(values[5:6]),var(values[7:8]),var(values[9:10]),var(values[11:12]))
  techvar <- c(techvar,mean(techvarTmp))
}
var_breakdown <- as.data.frame(cbind(techvar,biovar,indvar))
row.names(var_breakdown) <- sites      
var_breakdown <- log10(var_breakdown)
boxplot(var_breakdown, main = "log10 variance per phosphomeasurement (n=1562)")
summary(var_breakdown)



# here I am getting the F stat from the anova table froma simple  model <- lm(value~individual, data=test1)

sites <- c()
Fvalues <- c()
for(i in 1:1562){
  test <- melted[(melted$Var1 %in% melted$Var1[i]),]
  test1 <- test[3:6]
  model <- lm(value~individual, data=test1)
  Fval <- anova(model)[1,4]
  Fvalues <- c(Fvalues,Fval)
  sites <- c(sites,as.character(unique(test$Var1)))
}
Fvalues <- as.data.frame(Fvalues)
row.names(Fvalues) <- sites
boxplot(log10(Fvalues))
hist(log10(as.matrix(Fvalues)))
plot(density(log10(as.matrix(Fvalues))))
plot(density(as.matrix(Fvalues)))
density(as.matrix(Fvalues))



lot.new()
par(mfrow = c(1, 1))
for (i in 1:(ncol(data)-1)){
  if(i==1) plot(density(data[, i], na.rm=T), col = i, ylim = c(0,.7))
  else lines(density(data[, i], na.rm=T), col = i)
}






#variance calculation for test data
indvar <- var(test$value)
biovarTmp <- c(var(test$value[1:4]),var(test$value[5:8]),var(test$value[9:12]))
biovar <- mean(biovarTmp)
techvarTmp <- c(var(test$value[1:2]),var(test$value[3:4]),var(test$value[5:6]),var(test$value[7:8]),var(test$value[9:10]),var(test$value[11:12]))
techvar <- mean(techvarTmp)
#note how the variance is different than means and medians of a vector. very interesting. I think it has something to do with the degrees of freedom
# > data
# [1] 1.1416511 1.0039299 1.0957643 0.9887157 1.0085880 1.0422567 1.2912567 0.9889490 0.9749079 0.9291389
# > var(data)
# [1] 0.01113381
# > var(data[1:5])
# [1] 0.004515427
# > var(data[6:10])
# [1] 0.02053197
# > data[6:10]
# [1] 1.0422567 1.2912567 0.9889490 0.9749079 0.9291389
# > c(var(data[1:5]),var(data[6:10]))
# [1] 0.004515427 0.020531965
# > mean(c(var(data[1:5]),var(data[6:10])))
# [1] 0.0125237
# > var(data)
# [1] 0.01113381
# > mean(data)
# [1] 1.046516
# > c(mean(data[1:5]),mean(data[6:10]))
# [1] 1.047730 1.045302
# > mean(c(mean(data[1:5]),mean(data[6:10])))
# [1] 1.046516


# Per peptide one-way ANOVA with hierarchical/nested measurement is the way to go I believe. I think it should be all random effects if I am interested in how much variability to expect in phosphorylation amount at the individual population level.

#open questions: 

# 1) what does Yoav want?; my thinking is that he wants to calculate F statistics (essentially variance) for each peptide to see if peptides that are differentially expressed between individuals tend to be more variable in general;  Doesn't the omnibus F test already do this? Well all of the variability could come from one sample. It may be more appropriate to stay within limma to calculate F statistics! 'fit2$F'



# 2)Am I estimating variance appropriately above?

library(nlme)
#as of now unclear if intercept should be fixed or what. results are different for below
lmemodel <- lme(value~individual, data=test1, random =~1|individual/biorep/techrep)
lmemodel <- lme(value~1, data=test1, random =~1|individual/biorep/techrep)
lmemodel <- lme(value~1, data=test1, random =~1|individual/biorep)


lmemodel <- lme(value~individual, data=test1, random =~1|individual)
##lme model and linear models for F stats
lmemodel <- lme(value~individual, data=test1, random =~1|biorep/techrep)#with the 1 random intercepts are fitted. Without the 1 int and slopes are fitted
lmemodel <- lme(value~individual, data=test1, random =~1|individual/biorep/techrep)#individual must be treated as a random effect in order to get a different standard error for the slope estimates and to get variable intercept estimates depending on the individual/biological replicate level
anova.lme(lmemodel)
model <- lm(value~individual, data=test1)
anova(model)
anova(lmemodel)
#anova table reports the same result for F stat
summary(lmemodel)
summary(model)
coef(model)
coef(lmemodel)
#naive model
model1 <- lm(yield~Variety*nitro, data=Oats)
summary(model1)#intercept is the yield of variety golden rain
#the nitrogen effect is given by slopes...

model2 <- lme(yield~Variety*nitro, data=Oats, random=~1|Block/Variety/nitro)#model now saturated with random effects (no df for within group estimates)

model2 <- lme(yield~Variety*nitro, data=Oats, random=~1|Block/Variety)#This version gives the appropriate model formalism

#much to learn here
summary(model2)

# there are 72 observations and
# there are six blocks and three varieties in each block
# restricted ml fit (REML)
#note that the std errors and p values have changed

coef(model1)#gives the intercept (effect for varieties) terms and
#also gives the slope terms for nitrogen given a particular variety. Note these values are given
# as differencs between the first value

coef(model2)

colMeans(model2)#gives model1 coefficients
plot(ranef(model2))#random effects should so random spread of parameter estimate. Here it does
plot(model2)#standardized residuals to try and identify heteroskedasticity


Variogram.lme(model2)


s <- anova(model2)
summary(s)




#test$biorep <- as.factor(c(1,1,2,2,3,3,4,4,5,5,6,6))
#test$techrep <- as.factor(c(1:12))
# test$techrep <- as.factor(c(1:12))

test$biorep <- as.factor(c(1,1,2,2,1,1,2,2,1,1,2,2))
test$techrep <- as.factor(c(1,2,1,2,3,4,3,4,5,6,5,6))

#nested/hierarchical anova
res1 <- lm(test$value ~ test$individual/test$biorep/test$techrep)
anova(res1)

res1 <- lm(test$value ~ test$individual + test$individual/test$biorep/test$techrep)
anova(res1)

res1 <- lm(test$value ~ test$individual + test$individual/test$biorep + test$individual/test$biorep/test$techrep)
anova(res1)

res1 <- lm(test$value ~ test$individual + test$individual/test$biorep )
anova(res1)


res1 <- lm(test$value ~ test$individual + test$biorep + test$techrep)
res1 <- lm(test$value ~ test$techrep + test$individual)

res1 <- lm(test$value ~ test$individual + test$techrep)

res1 <- lme(test$value ~ test$individual + test$techrep)
and nlme too!!

anova(res1)


aov(test$value ~ test$individual + test$biorep + test$techrep)





