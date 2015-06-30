# code to solve the optimization problem
# || X - Y ||_obs ^ 2 subject to ||X||_* < r
# where _obs denotes that we only care about the non-missing (NA) values in Y and
# ||X||_* is the nuclear (trace) norm of X

require(irlba)
multitaskNuclearNorm=function(X,Y,r,its=1000,rmseTol=1e-3,warmStart=NULL,verbose=F){
  N=nrow(Y)
  D=ncol(Y)
  P=ncol(X)
  beta=if (is.null(warmStart)) matrix(0,P,D) else warmStart
  ertest=NA
  oldRmse=Inf
  for (it in 0:its){
    tol=max( 1e-1/(it+1)^2, 1e-6 )
    resid=Y-X %*% beta
    resid[is.na(Y)]=0.0 # copes with missingness
    g=-2*t(X) %*% resid
    svdTemp=irlba(-g,1,1,tol=tol)
    ruv=(r * svdTemp$u) %*% t(svdTemp$v)
    W=ruv-beta
    XWM=X %*% W
    XWM[is.na(Y)]=0.0
    stepSize=sum( (-.5*g) * W ) / sum( (t(X) %*% XWM) * W )
    if (stepSize<0)
      cat('Warning: step size is',stepSize,'\n')
    stepSize=min(stepSize,1)
    beta = (1.0-stepSize)*beta + stepSize * ruv
    resid=Y-X %*% beta
    resid[is.na(Y)]=0.0
    rmse=sqrt(mean(resid^2))
    rmseDelta=abs(rmse-oldRmse)
    if (verbose) 
      cat(it,rmse,sqrt(mean(ertest^2)),stepSize,rmseDelta,'\n')
    if ( rmseDelta < rmseTol)
      break
    oldRmse=rmse
  }
  beta
}

