# 
# Batch adjustment
# 
# To illustrate how we can adjust for batch effects using statistcal methods, we will create a data example in which the outcome of interest is confounded with batch but not completely. We will also select a outcome for which we have an expectation of what genes should be differentially expressed. Namely, we make sex the outcome of interest and expect genes on the Y chromosome to be differentially expressed. Note that we may also see genes from the X chromosome as differentially expressed as some escape X inactivation.
# 
# We start by finding the genes on the Y chromosome.


# install bioconductor stuff if not previously installed
# source("http://bioconductor.org/biocLite.R")
# library(dagdata)##Rafa script
# biocLite()

# installing packages from github!
install_github('rafalib','ririzarr')
install_github('dagdata','genomicsclass')

# library(RColorBrewer)
library(Biobase)
library(genefilter)
library(rafalib)
library(dagdata)##Rafa script
data(GSE5859)

biocLite("hgfocus.db")
library(hgfocus.db) ##get the gene chromosome
biocLite("sva")
library(sva)


chr<-mget(featureNames(e),hgfocusCHRLOC)
chr <- sapply(chr,function(x){ tmp<-names(x[1]); ifelse(is.null(tmp),NA,paste0("chr",tmp))})
y<- colMeans(exprs(e)[which(chr=="chrY"),])
sex <- ifelse(y<4.5,"F","M")

# Now we select samples so that sex and month of hybridization are confounded. 


batch <- format(pData(e)$date,"%y%m")
ind<-which(batch%in%c("0506","0510"))
set.seed(1)
N <- 12; N1 <-3; M<-12; M1<-9
ind <- c(sample(which(batch=="0506" & sex=="F"),N1),
         sample(which(batch=="0510" & sex=="F"),N-N1),
         sample(which(batch=="0506" & sex=="M"),M1),
         sample(which(batch=="0510" & sex=="M"),M-M1))
table(batch[ind],sex[ind])

#To illustrate the confounding we will pick some genes to show in a heatmap plot. We pick all Y chromosome genes, some genes that we see correlate with batch, and then some randomly selected genes. 


set.seed(1)
tt<-genefilter::rowttests(exprs(e)[,ind],factor(batch[ind]))
ind1 <- which(chr=="chrY") ##real differences
ind2 <- setdiff(c(order(tt$dm)[1:25],order(-tt$dm)[1:25]),ind1)
ind0 <- setdiff(sample(seq(along=tt$dm),50),c(ind2,ind1))
geneindex<-c(ind2,ind0,ind1)
mat<-exprs(e)[geneindex,ind]
mat <- mat-rowMeans(mat)#;mat[mat>3]<-3;mat[mat< -3]<- -3
icolors <- rev(brewer.pal(11,"RdYlBu"))
mypar(1,1)
image(t(mat),xaxt="n",yaxt="n",col=icolors)

# So what follows is like the analysis we would do in practice. We don't know there is a batch and we are interested in finding genes that are different between males and females. We start by computing t-statistics and p-values comparing males and females. We use histograms to notice the problem introduced by the batch.
# 
# The batch effect adjustment methods are best described with the linear models so we start by writing down the linear more for this particular case:


dat <- exprs(e)[,ind]
X <- sex[ind] ## the covariate
Z <- batch[ind]
tt<-genefilter::rowttests(dat,factor(X))
HLIM<-c(0,1500)
mypar(1,2)
hist(tt$p[!chr%in%c("chrX","chrY")],nc=20,xlab="p-value",ylim=HLIM,main="")
hist(tt$p[chr%in%c("chrY")],nc=20,xlab="p-value",ylim=c(0,9),main="")

# Combat
# 
# Here we show how to implement Combat. 


mod<-model.matrix(~X)

# mod<-model.matrix(~0+X) doesn't work

cleandat <- ComBat(dat,Z,mod)

tt<-genefilter::rowttests(cleandat,factor(X))
mypar(1,1)
hist(tt$p[!chr%in%c("chrX","chrY")],nc=20,xlab="p-value",ylim=HLIM,main="")

hist(tt$p[chr%in%c("chrY")],nc=20,xlab="p-value",ylim=c(0,9),main="")


# But what exactly is a batch?

times <- (pData(e)$date)
mypar(1,2)
o=order(times)
plot(times[o],pch=21,bg=as.fumeric(batch)[o],ylab="date")
o=order(times[ind])
plot(times[ind][o],pch=21,bg=as.fumeric(batch)[ind][o],ylab="date")





