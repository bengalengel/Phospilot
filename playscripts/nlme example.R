#nlme model example using the lme function and oats dataset
# https://www.youtube.com/watch?v=VhMWPkTbXoY

#how to model a nested random effects design
library(nlme)
data(Oats)
str(Oats)
# $ Block  : Ord.factor w/ 6 levels "VI"<"V"<"III"<..: 6 6 6 6 6 6 6 6 6 6 ...
# $ Variety: Factor w/ 3 levels "Golden Rain",..: 3 3 3 3 1 1 1 1 2 2 ...
# $ nitro  : num  0 0.2 0.4 0.6 0 0.2 0.4 0.6 0 0.2 ...
# $ yield  : num  111 130 157 174 117 114 161 141 105 140 ...

plot(Oats)#"possible only with a grouped dataset" shows the yield for for each plot (1-6)
#as a function of nitrogen (4 values). For each nitrogen amount each crop variety gives a response.

#naive model
model1 <- lm(yield~Variety*nitro, data=Oats)
summary(model1)#intercept is the yield of variety golden rain and the other estimates are average differences
#the nitrogen effect is given by slopes (this is an interaction effect)




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
