### R code from vignette source 'sva.Rnw'

###################################################
### code chunk number 1: sva.Rnw:5-6
###################################################
options(width=65)


###################################################
### code chunk number 2: style-Sweave
###################################################
BiocStyle::latex()


###################################################
### code chunk number 3: input
###################################################
library(sva)
library(bladderbatch)
data(bladderdata)
library(pamr)
library(limma)


###################################################
### code chunk number 4: input
###################################################
pheno = pData(bladderEset)


###################################################
### code chunk number 5: input
###################################################
edata = exprs(bladderEset)


###################################################
### code chunk number 6: input
###################################################
mod = model.matrix(~as.factor(cancer), data=pheno)


###################################################
### code chunk number 7: input
###################################################
mod0 = model.matrix(~1,data=pheno)


###################################################
### code chunk number 8: input
###################################################
n.sv = num.sv(edata,mod,method="leek")
n.sv


###################################################
### code chunk number 9: input
###################################################
svobj = sva(edata,mod,mod0,n.sv=n.sv)


###################################################
### code chunk number 10: input
###################################################
pValues = f.pvalue(edata,mod,mod0)
qValues = p.adjust(pValues,method="BH")


###################################################
### code chunk number 11: input
###################################################
modSv = cbind(mod,svobj$sv)
mod0Sv = cbind(mod0,svobj$sv)

pValuesSv = f.pvalue(edata,modSv,mod0Sv)
qValuesSv = p.adjust(pValuesSv,method="BH")


###################################################
### code chunk number 12: input
###################################################
fit = lmFit(edata,modSv)


###################################################
### code chunk number 13: input
###################################################
contrast.matrix <- cbind("C1"=c(-1,1,0,rep(0,svobj$n.sv)),"C2"=c(0,-1,1,rep(0,svobj$n.sv)),"C3"=c(-1,0,1,rep(0,svobj$n.sv)))
fitContrasts = contrasts.fit(fit,contrast.matrix)


###################################################
### code chunk number 14: input
###################################################
eb = eBayes(fitContrasts)
topTableF(eb, adjust="BH")


###################################################
### code chunk number 15: input
###################################################
batch = pheno$batch


###################################################
### code chunk number 16: input
###################################################
modcombat = model.matrix(~1, data=pheno)


###################################################
### code chunk number 17: input
###################################################
combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, numCovs=NULL, par.prior=TRUE, prior.plots=T)
combat_edata2 = ComBat(dat=edata, batch=batch, mod=mod, numCovs=NULL, par.prior=TRUE, prior.plots=T)




###################################################
### code chunk number 18: input
###################################################
pValuesComBat = f.pvalue(combat_edata,mod,mod0)
qValuesComBat = p.adjust(pValuesComBat,method="BH")

pValuesComBat2 = f.pvalue(combat_edata2,mod,mod0)
qValuesComBat2 = p.adjust(pValuesComBat2,method="BH")

data <- cbind(qValuesComBat,qValuesComBat2)
# Density plot overlays by color show a difference in the distributions
plot.new()
par(mfrow = c(1, 1))
for (i in 1:(ncol(data))){
  if(i==1) plot(density(data[, i], na.rm=T), col = i, ylim = c(0,2))
  else lines(density(data[, i], na.rm=T), col = i)
}


###################################################
### code chunk number 19: input
###################################################
modBatch = model.matrix(~as.factor(cancer) + as.factor(batch),data=pheno)
mod0Batch = model.matrix(~as.factor(batch),data=pheno)
pValuesBatch = f.pvalue(edata,modBatch,mod0Batch)
qValuesBatch = p.adjust(pValuesBatch,method="BH")


###################################################
### code chunk number 20: input2
###################################################
n.sv = num.sv(edata,mod,vfilter=2000,method="leek")
svobj = sva(edata,mod,mod0,n.sv=n.sv,vfilter=2000)


###################################################
### code chunk number 21: input
###################################################
set.seed(12354)
trainIndicator = sample(1:57,size=30,replace=F)
testIndicator = (1:57)[-trainIndicator]

trainData = edata[,trainIndicator]
testData = edata[,testIndicator]

trainPheno = pheno[trainIndicator,]
testPheno = pheno[testIndicator,]


###################################################
### code chunk number 22: input
###################################################
mydata = list(x=trainData,y=trainPheno$cancer)
mytrain = pamr.train(mydata)
table(pamr.predict(mytrain,testData,threshold=2),testPheno$cancer)


###################################################
### code chunk number 23: input
###################################################
trainMod = model.matrix(~cancer,data=trainPheno)
trainMod0 = model.matrix(~1,data=trainPheno)
trainSv = sva(trainData,trainMod,trainMod0)


###################################################
### code chunk number 24: input
###################################################
fsvaobj = fsva(trainData,trainMod,trainSv,testData)
mydataSv = list(x=fsvaobj$db,y=trainPheno$cancer)
mytrainSv = pamr.train(mydataSv)
table(pamr.predict(mytrainSv,fsvaobj$new,threshold=1),testPheno$cancer)


###################################################
### code chunk number 25: input
###################################################
library(zebrafishRNASeq)
library(genefilter)
data(zfGenes)
filter = apply(zfGenes, 1, function(x) length(x[x>5])>=2)
filtered = zfGenes[filter,]
genes = rownames(filtered)[grep("^ENS", rownames(filtered))]
controls = grepl("^ERCC", rownames(filtered))
group = as.factor(rep(c("Ctl", "Trt"), each=3))
dat0 = as.matrix(filtered)


###################################################
### code chunk number 26: input3
###################################################
## Set null and alternative models (ignore batch)
mod1 = model.matrix(~group)
mod0 = cbind(mod1[,1])

svseq = svaseq(dat0,mod1,mod0,n.sv=1)$sv
plot(svseq,pch=19,col="blue")


###################################################
### code chunk number 27: input4
###################################################
sup_svseq = svaseq(dat0,mod1,mod0,controls=controls,n.sv=1)$sv
plot(sup_svseq, svseq,pch=19,col="blue")

