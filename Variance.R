gene.clean = read.table('Pilot_geneexpr_finalclean.txt', as.is=T)
names = read.csv('SampleKey_PrunednoiPSC_clean_techavg.csv', header=T, as.is=T)
colnames(gene.clean) = names$Name

Cont = c(1,2,4,9,16,17,20,23,26,28)
T40 = c(3,8,13,14,15,19,24,29,32)
eCont<-gene.clean[,Cont]
eT40<-gene.clean[,T40]

colnames(gene.clean[,Cont])
colnames(gene.clean[,T40])

C40 = cbind(eCont, eT40)

thislist = c("16A", "16B", "1EA", "1EB", "18", "17A", "17B", "24A", "24B")
output.2 = matrix(data=NA, ncol=9, nrow= 14625)
for(i in 1:length(thislist)){
  samp1 = paste(thislist[i],"0",sep="")
  samp2 = paste(thislist[i],"40",sep="")
  ind1 = match(samp1,colnames(C40))
  ind2 = match(samp2,colnames(C40))
  output.2[,i] = abs(C40[,ind1] - C40[,ind2])
}

colnames(output.2) = thislist
thislist.2 = c("16", "1E", "18", "17", "24")

#varlist.2= list()
#for(i in 1:length(thislist.2)){
#  currind = grep(thislist.2[i], colnames(output.2))
#  if(length(currind)<2){
#    next
#  }
#  currvar = apply(output[,currind],1, var)
#  varlist.2[[1]]=c(varlist.2[[1]], currvar)
#}

varlist.2= list()
varlist.2[[1]] = NA
varlist.2[[2]] = NA
for(i in 1:length(thislist.2)){
  currind = grep(thislist.2[i], colnames(output.2))
  for(j in 1:length(currind)){
    currcol = currind[j]
    currmatrix = output.2[,-currind]
    for(k in 1:dim(currmatrix)[2]){
      currvar = apply(cbind(output.2[,currcol], currmatrix[,k]),1,var)
      varlist.2[[2]] = c(varlist.2[[2]], currvar)
    }
  }
  if(length(currind)<2){
    next
  }
  currvar = apply(output.2[,currind],1, var)
  varlist.2[[1]]=c(varlist.2[[1]], currvar)
}
varlist.2[[1]] = varlist.2[[1]][-1]
varlist.2[[2]] = varlist.2[[2]][-1]

boxplot(log(varlist.2[[1]]), log(varlist.2[[2]]), ylab = "Log(Variance)",main = "Variance in Fold Change of all Genes")

thislist.3 = c("16", "16", "1E", "1E", "18", "17", "17", "24", "24")
thislist.3.f = as.factor(thislist.3)

ouranova = function(x){
  anova(lm(x~thislist.3.f))$'Pr(>F)'[1]
}
pval.var = apply(output.2,1,ouranova)
sum(pval.var <= .05)
which.min(pval.var)
hist(pval.var)
library(qvalue)
qvals.var <- qvalue(pval.var)$qvalues
hist(qvals.var)
length(which(qvals.var <= .05))


##Variance and test of only DE genes
thislist = c("16A", "16B", "1EA", "1EB", "18", "17A", "17B", "24A", "24B")
output = matrix(data=NA, ncol=9, nrow= 5729)
for(i in 1:length(thislist)){
  samp1 = paste(thislist[i],"0",sep="")
  samp2 = paste(thislist[i],"40",sep="")
  ind1 = match(samp1,colnames(C40.DE))
  ind2 = match(samp2,colnames(C40.DE))
  output[,i] = abs(C40.DE[,ind1] - C40.DE[,ind2])
}

colnames(output) = thislist
thislist.2 = c("16", "1E", "18", "17", "24")

#varlist= list()
#for(i in 1:length(thislist.2)){
#  currind = grep(thislist.2[i], colnames(output))
#  if(length(currind)<2){
#    next
#  }
#  currvar = apply(output[,currind],1, var)
#  varlist[[1]]=c(varlist[[1]], currvar)
#}

varlist= list()
varlist[[1]] = NA
varlist[[2]] = NA
for(i in 1:length(thislist.2)){
  currind = grep(thislist.2[i], colnames(output))
  for(j in 1:length(currind)){
    currcol = currind[j]
    currmatrix = output[,-currind]
    for(k in 1:dim(currmatrix)[2]){
      currvar = apply(cbind(output[,currcol], currmatrix[,k]),1,var)
      varlist[[2]] = c(varlist[[2]], currvar)
    }
  }
  if(length(currind)<2){
    next
  }
  currvar = apply(output[,currind],1, var)
  varlist[[1]]=c(varlist[[1]], currvar)
}
varlist[[1]] = varlist[[1]][-1]
varlist[[2]] = varlist[[2]][-1]


mtext("Within Individuals                             Across Individuals", side =1, line = 0, cex=.7, padj = 2, pos=2) 
boxplot(log(varlist[[1]]), log(varlist[[2]]),  ylab = "Log(Variance)", main = "Variance in Fold Change of \nDifferentially Expressed Genes")

thislist.3 = c("16", "16", "1E", "1E", "18", "17", "17", "24", "24")
thislist.3.f = as.factor(thislist.3)

ouranova = function(x){
  anova(lm(x~thislist.3.f))$'Pr(>F)'[1]
}
pval.var.DE = apply(output,1,ouranova)
hist(pval.var.DE, main = "Histogram of pvals")
length(which(pval.var.DE <= .05))
library(qvalue)
qvals.var.DE <- qvalue(pval.var.DE)$qvalues
hist(qvals.var.DE)
length(which(qvals.var.DE <= .05))
which.min(pval.var.DE)
C40.DE[2721,]
output[2721,]

var.test(varlist[[1]], varlist[[2]])