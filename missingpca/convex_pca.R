# code to solve the optimization problem
# || X - Y ||_obs ^ 2 subject to ||X||_* < r
# where _obs denotes that we only care about the non-missing (NA) values in Y and
# ||X||_* is the nuclear (trace) norm of X

require(irlba)

convexPCA=function(Y,r,Ytest=NULL,its=1000,rmseTol=1e-3,warmStart=NULL,verbose=F){
  P=nrow(Y)
  N=ncol(Y)
  X=if (is.null(warmStart)) matrix(0,P,N) else warmStart
  ertest=NA
  oldRmse=Inf
  for (it in 0:its){
    tol=max( 1e-1/(it+1)^2, 1e-6 )
    g=X-Y
    g[is.na(Y)]=0.0
    svdTemp=irlba(-g,1,1,tol=tol)
    ruv=(r * svdTemp$u) %*% t(svdTemp$v)
    er=g[!is.na(Y)]
    erv=X[!is.na(Y)]-ruv[!is.na(Y)]
    stepSize=sum(er*erv)/sum(erv*erv)
    if (stepSize<0)
      cat('Warning: step size is',stepSize,'\n')
    stepSize=min(stepSize,1)
    X = (1.0-stepSize)*X + stepSize* ruv
    er=X[!is.na(Y)]-Y[!is.na(Y)]
    if (!is.null(Ytest)) 
      ertest=X[is.na(Y)]-Ytest[is.na(Y)]
    rmse=sqrt(mean(er^2))
    rmseDelta=abs(rmse-oldRmse)
    if (verbose) 
      cat(it,rmse,sqrt(mean(ertest^2)),stepSize,rmseDelta,'\n')
    if ( rmseDelta < rmseTol)
      break
    oldRmse=rmse
  }
  X
}

cv.convexPCA=function(Y,its=1000,rmseTol=1e-3) {
  
  # break the matrix into training, test and validation sets, equally and at random
  rand=matrix( runif(nrow(Y)*ncol(Y)),nrow(Y))
  train=rand < 2/3
  test= !train & !is.na(Y)
  
  # training data
  trainY=Y
  trainY[!train]=NA
  
  # try a range of regularisation parameters
  rs=10^seq(0,log10(sum(trainY^2,na.rm=T)),length.out = 10)
  trainErrors=numeric(length(rs))
  testErrors=numeric(length(rs))
  times=numeric(length(rs))
  for (ri in 1:length(rs)){
    r=rs[ri]
    cat("Trying regularisation parameter r=",r,"...\n")
    times[ri]=system.time( X<-convexPCA(trainY,r,its=its,rmseTol=rmseTol,warmStart=if (ri>1) X else NULL ) )[1]
    trainErrors[ri]=sqrt(mean( (X-trainY)^2,na.rm=T))
    testErrors[ri]=sqrt(mean( (X[test]-Y[test])^2))
    cat("Error is",testErrors[ri],"\n")
  }
  rs[which.min(testErrors)]
}

wrapperPCA=function(Y,pcsRequired=5) {
  optimalR=cv.convexPCA(Y)
  cat("Picked regularisation parameter r=",optimalR,"\n")
  Yfull=convexPCA(Y,optimalR)
  irlba(Yfull,pcsRequired,pcsRequired)
}