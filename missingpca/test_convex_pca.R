source("./missingpca/convex_pca.R")

# generate some synthetic data
P=100
N=200
K=min(P,N)
Y=matrix(rnorm(P*K),P) %*% diag(rgamma(K,shape=.1)) %*% matrix(rnorm(K*N),K)
Y=scale(Y)

# get the SVD on the full data
sFull=svd(Y)

# break the matrix into training, test and validation sets, equally and at random
rand=matrix( runif(nrow(Y)*ncol(Y)),nrow(Y))
train=rand < 1/3
test=(1/3 <= rand) & (rand < 2/3)
validation=2/3 <= rand 

# training data
trainY=Y
trainY[!train]=NA

w=wrapperPCA(trainY)

reconstructedY = w$u %*% diag(w$d) %*% t(w$v)

sqrt(mean((reconstructedY-trainY)^2,na.rm=T))

# try a range of regularisation parameters
rs=10^seq(0,log10(sum(trainY^2,na.rm=T)),length.out = 10)
trainErrors=numeric(length(rs))
testErrors=numeric(length(rs))
times=numeric(length(rs))
for (ri in 1:length(rs)){
  r=rs[ri]
  cat("------- r=",r,"--------\n")
  times[ri]=system.time( X<-convexPCA(trainY,r,warmStart=if (ri>1) X else NULL ) )[1]
  trainErrors[ri]=sqrt(mean( (X-trainY)^2,na.rm=T))
  testErrors[ri]=sqrt(mean( (X[test]-Y[test])^2))
  cat("\nError is",trainErrors[ri],"\n")
}

# plot training and test error
plot(rs,trainErrors,log="x",type="b",xlab="regularization parameter",ylab="reconstruction error")
lines(rs,testErrors,type="b",col="red")
legend("topright",legend=c("train error","test error"),col=c("black","red"),lwd=1)
optimalIndex=which.min(testErrors)
rOptimal=rs[optimalIndex]
X=convexPCA(trainY,rOptimal)
validationError=sqrt(mean( (X[validation]-Y[validation])^2))
cat("Optimal r:",rOptimal," train RMSE:",trainErrors[optimalIndex]," test RMSE:",testErrors[optimalIndex]," validation RMSE:", validationError,"\n")

# get the SVD on the reconstruction
svdPartial=svd(X)

# plot true vs inferred eigenvalues
plot(sFull$d,pch=16,col="red",ylim=c(0,max(sFull$d)),ylab="eigenvalues",type="b",log="x")
points(svdPartial$d,col="black",type="b")
legend("topright",legend=c("true","inferred"),col=c("black","red"),lwd=1)