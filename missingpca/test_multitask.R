

source("multitask_nuc.R")

# generate some synthetic data
P=200
D=5
N=100

beta=matrix(rnorm(P*D),P)
beta[,1]=beta[,2] # low rank!
X=matrix(rnorm(N*P),N)
Y=X %*% beta + .1 * matrix(rnorm(N*D),N)

hatbeta=multitaskNuclearNorm(X,Y,10,verbose=T)

X=x
Y=y
N=nrow(Y)
D=ncol(Y)
P=ncol(X)

#rs=10^seq(0,3,length.out = 10)
rs=10^seq(-2,0,length.out = 5)
nfolds=10
foldid = sample(rep(seq(nfolds), length = N))
res=list()
for (i in 1:nfolds) {
  trainErrors=numeric(length(rs))
  testErrors=numeric(length(rs))
  trainY=Y[foldid!=i,]
  trainX=X[foldid!=i,]
  testY=Y[foldid==i,]
  testX=X[foldid==i,]
  times=numeric(length(rs))
  for (ri in 1:length(rs)){
    r=rs[ri]
    cat("------- r=",r,"--------\n")
    times[ri]=system.time( betahat<-multitaskNuclearNorm(trainX,trainY,r,verbose=T,warmStart=if (ri>1) betahat else NULL ) )[1]
    trainYrecon=trainX %*% betahat
    testYrecon=testX %*% betahat
    trainErrors[ri]=sqrt(mean( (trainY-trainYrecon)^2,na.rm=T))
    testErrors[ri]=sqrt(mean( (testY[!is.na(testY)]-testYrecon[!is.na(testY)])^2))
  }
  res[[i]]=list(trainRMSE=trainErrors,testRMSE=testErrors)
}

trainErrs=do.call(rbind,lapply(res,function(g) g[[1]]))
testErrs=do.call(rbind,lapply(res,function(g) g[[2]]))
colnames(trainErrs)=format(rs,digits=2)
colnames(testErrs)=format(rs,digits=2)
boxplot(trainErrs)
boxplot(testErrs)

y=scale(sens)

yNoNa=y
yNoNa[is.na(y)]=0.0

alpha=0.8
print("Running cross validation to choose lambda")
cv.res=cv.glmnet(x,yNoNa,alpha=alpha,family="mgaussian")



