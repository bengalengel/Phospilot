# you need the irlba package installed! Run install.packages('irlba') if you don't have it. 

source("./missingpca/convex_pca.R")

# generate some synthetic data
P=100
N=200
K=min(P,N)
Y=matrix(rnorm(P*K),P) %*% diag(rgamma(K,shape=.1)) %*% matrix(rnorm(K*N),K)
Y=scale(Y)
Y[ runif(P*N)<.2 ]=NA # add some missing entries

mysvd=wrapperPCA(trainY)

reconstructedY = mysvd$u %*% diag(mysvd$d) %*% t(mysvd$v)

sqrt(mean((reconstructedY-trainY)^2,na.rm=T))