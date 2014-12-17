zia <- read.table("C:/Users/Brett/Documents/Dropbox/", sep = "\t", header=T, fill = T)

zia <- read.table("E:/My Documents/Dropbox/Postdoc-Gilad/Zia/AllProFixed.txt", sep = "\t", header=T, fill = T)

proteincomp <- na.omit(zia)
pro <- proteincomp[2:63]
pro <- as.matrix(pro)
protbiovar <- apply(pro,1,var)
summary(protbiovar)
protbiovar <- log10(protbiovar)

#melted phospho dataset
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


##combine with other data plotting interindividual variance estimate of phospho vs zia's protein variance estimates
boxplot(var_breakdown$indvar,protbiovar)

# t-test with unequal sample number?

qqnorm(var_breakdown$indvar)#normal
qqnorm(protbiovar)#normal

t.test(var_breakdown$indvar,protbiovar)

# Welch Two Sample t-test
# 
# data:  var_breakdown$indvar and protbiovar
# t = 14.7946, df = 2999.016, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.1781315 0.2325614
# sample estimates:
#   mean of x  mean of y 
# -0.9786413 -1.1839878 



