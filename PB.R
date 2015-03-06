#Pinheiro and Bates examples
library(nlme)

##first is the Rail data
data(Rail)
plot(Rail)#nice plots for 'grouped datasets' see ch3 for how to do this. 

#Here I want to know the average travel time for a 'typical' rail. That is I am interested in the population. 

#the plot shows considerable between rail variability with little within rail variability

#basic lm fit to estimate population average
fm1Rail.lm <- lm( travel ~ 1, data = Rail )
fm1Rail.lm

#average is sound but doesn't take into account the correlated errors within the groups of Rail. That is, the error is not independent and these 'goup effects' are inserted into the residuals, leading to an increased estimate of the within rails variability. This would also be considered pseudoreplication to estimate population average. However, this model does provide an estimate of between rail variability
x <- resid(fm1Rail.lm)
y <- Rail$Rail
y <- as.vector(y)
plot(x,y)

#one can account for this by using the individual estimates to calculate per rail averages.  
fm2Rail.lm <- lm( travel ~ Rail -1, data = Rail )#removes intercept term
fm2Rail.lm
#note the residual standard error is lower, and similar boxplots show the absolute residuals centered closer to 0
x <- resid(fm2Rail.lm)
y <- Rail$Rail
y <- as.vector(y)
plot(x,y)

#Here we are accounting for the correlated measurements within the different rails but this does not give an estimate of between rail variability?...(isn't this the Mean Sq?) Also the number of terms in the model increase linearly with the number of rails, thats a bad thing for sure.

#The above FE model is equivalent to
obsij = Railavg + (Raili-Railavg) + errorij
#in fact this is what you would get returned if I didn't put the intercept correction (-1) term in.

#A random effects model replaces Railavg with the mean travel time over the population of rails and (Raili-Railavg) deviations are replaced with random variables whos DISTRIBUTION IS TO BE ESTIMATED. That is,
obsij = Railpopmean + bi + eij
#where Railpopmean is the mean travel time for the population being sampled, bi is a random variable representing the deviation of the mean travel time for the ith rail from the population, and eij is a random variable representing the deviation of the travel time for observation j on rail i from the mean travel time for all i.
#each random variable is assumed to be normally distributed, independant with constant variances. 

#this model is nested because it has to levels of random variables. They are random becauuse they are selected from a population and they are effects because they represent the shift from the population mean. 

#the parameters of this model are Railpopmean, the variance of bi and the variance of eij. we do produce predictions of bi, bihat but they are used to estimate the variance

#here is the random effects model
fm1Rail.lme <- lme(travel ~ 1, data = Rail, random = ~ 1 | Rail)

summary(fm1Rail.lme)

#note that the random variable standard deviations are close to those for the residuals in the linear models (see text).

#if I want to change the optimization approach
fm1Rail.lmeML <- update( fm1Rail.lme, method = "ML" )
summary(fm1Rail.lmeML)

#assess constant variance in eij (heteroskedasticity)
plot(fm1Rail.lme)

#shown are the residuals eij = obsij-Railpopmean(hat)-bi(hat)
#heteroskedasticity may be apparent over larger ranges of continuous variables

#95% confidence intervals
intervals(fm1Rail.lme)
#note the large intervals for the random effect and the fixed effect (intercept). This will likely be the case for me as well

#for fixed effects
anova(fm1Rail.lme)
#the null here is that the intercept equals 0. Which is not relevant

#Randomized block design 1.2 **************************************************************8

#two factors. A fixed experimental and a random 'blocking' factor/effect.

data(ergoStool)
plot.design(ergoStool)#see the magnitude of each effect
ergoStool

#model with fixed type of stool effect and random subject effect
obsij = Typej + Subjecti + errorij

#Each subject tries each stool once. The above is equivalent to a matrix formulation, where for each subject, four stool fits, four subject specific random effects, and four subject specific errors are estimated. See book. 

#which stool is easier to get up from, knowing that subjects have varying difficulties arising
#I have four stools and 9 subjects. I need a constrast matrix for the four stools.

contrasts(ergoStool$Type)
#a baseline for T1 is produced and three rows which are the contrasts. This is referred to the 'contrasts' matrix
model.matrix(effort ~ Type, ergoStool)
#note how there is a model matrix (or Xi matrix) for each individual with a fixed intercept.

#for just one contrast
ergoStool1 <- ergoStool[ ergoStool$Subject == "1", ]
model.matrix(effort ~ Type, ergoStool1)

#fit the model
fm2Stool <- lme(effort ~ Type, data = ergoStool, random = ~ 1 | Subject)
summary(fm2Stool)
#anova
anova(fm2Stool)
#some confustion about the anova output here.Testing B2-B4 = 0.

#fitting using the cell means parameterization. The column of 1s is removed
model.matrix( effort ~ Type - 1, ergoStool1 )
fm3Stool <- lme(effort ~ Type - 1, data = ergoStool, random = ~ 1 | Subject)
summary(fm3Stool)
anova(fm3Stool)
#here I am testing B1-B4 = 0. This is not an appropriate test. Here's what to know:
# 1) The overall effect of the factor should be assessed with anova, not
# by examining the t-value’s or p-value’s associated with the fixedeffects parameters. The anova output does not depend on the choice of contrasts as long as the intercept term is retained in the model.

# 2)Interpretation of the parameter estimates for a fixed-effects term depends on the contrasts being

# 3) # For REML estimation, likelihood-ratio tests or comparisons of AIC
# or BIC require the same fixed-effects structure and the same choice
# of contrasts in all models.

# 4) The “cell means” parameters can be estimated by adding -1 to a model formula but this will usually make the results of anova meaningless.

#assessing model fit
intervals(fm2Stool)
plot(fm2Stool)#residuals do not show heteroskedasticity 

plot( fm2Stool, form = resid(., type = "p") ~ fitted(.) | Subject, abline = 0 )#residuals cut by subject.



































