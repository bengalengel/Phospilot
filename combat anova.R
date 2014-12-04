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





