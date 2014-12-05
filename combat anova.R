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

# for(i in 1:1562)
test <- melted[(melted$Var1 %in% melted$Var1[3]),]

test$biorep <- as.factor(c(1,1,2,2,3,3,4,4,5,5,6,6))
test$techrep <- as.factor(c(1:12))
# test$techrep <- as.factor(c(1:12))

test$biorep <- as.factor(c(1,1,2,2,1,1,2,2,1,1,2,2))
test$techrep <- as.factor(c(1,2,1,2,3,4,3,4,5,6,5,6))



library(nlme)
lme.mod1 <- lme(value ~ individual, random=~1|individual/biorep/techrep, test)

rsquared.glmm(lme.mod1)

#Change to old optimizer to solve convergence issueslme4.models2[2:4, -(1:3)]
lme.mod1.1 <- lme(y ~ fixed1, random=~fixed1|rand2/rand1, control = lmeControl(opt = "optim"), data) 
lme.mod2 <- lme(y ~ fixed1 + fixed2, random=~1|rand2/rand1, data)
(lme.models <- rsquared.glmm(list(lme.mod1, lme.mod1.1, lme.mod2)))
# Compare to lme4 models, minor differences
all.equal(lme4.models2[2:4, -(1:3)], lme.models[,-(1:3)], tol = 1e-4)


library(lme4)
value <- test$value
individual <- test$individual
biorep <- test$biorep
techrep <- test$techrep

model<-lmer(value~individual +(1|individual/biorep)) 

model<-lmer(y~treatment+(1|treatment/mouse),family=binomial) 

test$value ~ test$individual/test$biorep/test$techrep





library(nlme)
lme.mod1 <- lme(y ~ fixed1, random=~1|rand2/rand1, data)
#Change to old optimizer to solve convergence issueslme4.models2[2:4, -(1:3)]
lme.mod1.1 <- lme(y ~ fixed1, random=~fixed1|rand2/rand1, control = lmeControl(opt = "optim"), data) 
lme.mod2 <- lme(y ~ fixed1 + fixed2, random=~1|rand2/rand1, data)
(lme.models <- rsquared.glmm(list(lme.mod1, lme.mod1.1, lme.mod2)))
# Compare to lme4 models, minor differences
all.equal(lme4.models2[2:4, -(1:3)], lme.models[,-(1:3)], tol = 1e-4)
















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





