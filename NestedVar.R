NestedVar <- function(ratios, batch = F, balanced = T){
  #this function accepts a matrix of H/L ratios and calculates variance components within a nested design. A random effects model from the nlme package is used. The batch flag indicates if only one batch is passed to this function, requiring a model that does not include 'biological' or tissue culture replication. The nlme package from Piniero and Bates is used. A matrix variance components is returned.

  # The input matrix is converted into short form for linear modeling
  require(plyr)
  require(reshape2)
  require(qdapRegex)
  require(nlme)
  
  ratios <- as.matrix(ratios)
  if(balanced){
    ratios <- na.omit(ratios)#balanced data
  }
  #no missing values version
  melted <- melt(ratios, measure.vars = names(ratios))
  
  #identify individual name and add it to the table
  matches <- gregexpr("[0-9]{5}", melted$Var2, perl=T)
  individual <- regmatches(melted$Var2,matches)
  individual <- as.character(individual)
  individual <- as.factor(individual)
  melted$individual <- individual
  
  #identify the biological replicate
  biorep <- rm_between(melted$Var2, "_", "_", extract=TRUE)
  biorep <- as.character(biorep)
  biorep <- as.factor(biorep)
  melted$biorep <- biorep
  
  #identify the technical replicate
  matches <- gregexpr("[0-9]$", melted$Var2, perl=T)
  techrep <- regmatches(melted$Var2,matches)
  techrep <- as.character(techrep)
  techrep <- as.factor(techrep)
  melted$techrep <- techrep
  
  if(!batch){
  #all random model with a fixed intercept.technical replicate not explicitly modeled. 
  sites <- c()
  Varcomp <- c()
  Expindvar <- c()
  Expindvarmeans <- c()
  Expbiovar <- c()
  Expbiovarmeans <- c()
  Exptechvar <- c()
  #new loop that gives a model for each unique phosphosite. Fits biorep as a character! also added an error catch
  for(id in levels(melted$Var1)){
    test <- melted[melted$Var1 %in% id,]
    test1 <- test[3:6]
    test1$biorep <- as.factor(test1$biorep)
    pos_err <- tryCatch(lme(value~1, data=test1, random =~1|individual/biorep, na.action = na.exclude),error=function(e) e)
    if(!inherits(pos_err, "error")){
      lmemodel <- lme(value~1, data=test1, random =~1|individual/biorep, na.action = na.exclude)
      temp <- as.numeric(VarCorr(lmemodel)[,1])
      temp <- temp[!is.na(temp)]
      temp <- na.omit(temp)
      Expindvartmp <- 4*temp[1]+2*temp[2]+temp[3]
      Expindvar <- c(Expindvar,Expindvartmp)
      Expindvarmeanstmp <- Expindvartmp/4
      Expindvarmeans <- c(Expindvarmeans,Expindvarmeanstmp)
      Expbiovartmp <- 2*temp[2]+temp[3]
      Expbiovar <- c(Expbiovar,Expbiovartmp)
      Expbiovarmeanstmp <- Expbiovartmp/2
      Expbiovarmeans <- c(Expbiovarmeans,Expbiovarmeanstmp)
      Exptechvartmp <- temp[3]
      Exptechvar <- c(Exptechvar,Exptechvartmp)
      Varcomp <- cbind(Varcomp,temp)
      sites <- c(sites,as.character(unique(test$Var1)))
    }
  }
  
  #plot the variance component distributions
  colnames(Varcomp) <- sites
  row.names(Varcomp) <- c("individual","biorep","residual")
  dim(Varcomp)
  Varcomp <- t(Varcomp)
  boxplot(log10(Varcomp), ylab = "log10 variance component")
  summary(Varcomp)
  head(Varcomp)
  
  ##cumulative Varcomp
  Varcomp2 <- Varcomp
  IndCum <- rowSums(Varcomp2)
  Varcomp2[,1] <- 0
  BioCum <- rowSums(Varcomp2)
  Varcomp2[,2] <- 0
  techCum <- rowSums(Varcomp2) #for symmetry!
  CumMat <- cbind(IndCum,BioCum,techCum)
  boxplot(log10(CumMat), ylab = "log10 cumulative variance components")
  
  
  #Expected Variance components
  ExpVars <- cbind(Expindvar,Expbiovar,Exptechvar)
  boxplot(log10(ExpVars), ylab = "log10 Expected variance")
  
  ##Ratio of component estimates relative to technical
  Ind_tech <- Varcomp[,1]/Varcomp[,3]
  Bio_tech <- Varcomp[,2]/Varcomp[,3]
  
  ratios_tech <- cbind(Ind_tech,Bio_tech)
  boxplot(log10(ratios_tech), ylab = "log10 ratios of variance components")
  
  ##scatter plots show discrete groupings?
  plot(log10(Varcomp[,1]),log10(Varcomp[,3]), xlab = "individual variance", ylab = "technical variance")
  plot(log10(Varcomp[,1]),log10(Varcomp[,2]), xlab = "individual variance", ylab = "biological variance")
  plot(log10(Varcomp[,2]),log10(Varcomp[,3]), xlab = "biological variance", ylab = "technical variance")
  
  #histograms of log10 variance are bimodal for biological and individual variance.
  plot(density(log10(Varcomp)), xlab = "variance", main = "hist of variance components")
  plot(density(log10(Varcomp[,1])), xlab = "variance", main = "hist of individual variance components")
  plot(density(log10(Varcomp[,2])), xlab = "variance", main = "hist of biological variance components")
  plot(density(log10(Varcomp[,3])), xlab = "variance", main = "hist of technical variance components")
  
  
  ##Some Joyce additions below
  hist(log10(rowSums(Varcomp)))
  
#   As expected, larger sample variability (total VC) permits variability between individual samples, between biological replicates, and betweeen techincial replicates. 
#   
#   Note that individual VC and biological VC are extremely small in a good amount phosphopeptides. More importantly, these phosphopeptides are not limited in the range of sample variability and span the entire range of total VC.
par(mfrow=c(2,2))
ylims=c(-20,0);xlims=c(-3,1)
plot(log10(rowSums(Varcomp)),log10(Varcomp[,1]),xlim=xlims,ylim=ylims,
     xlab="log10 total VC",ylab="log10 individual VC",axes=F)
axis(1);axis(2)
plot(log10(rowSums(Varcomp)),log10(Varcomp[,2]),xlim=xlims,ylim=ylims,
     xlab="log10 total VC",ylab="log10 biorep VC",axes=F)
axis(1);axis(2)
plot(log10(rowSums(Varcomp)),log10(Varcomp[,3]),xlim=xlims,ylim=ylims,
     xlab="log10 total VC",ylab="log10 tech VC",axes=F)
axis(1);axis(2)
# Here we'd like to identify phosphpeptides with little or no variability at the individual level and at the biological replicate level. To do so, we standardized the values of the variance components for each phosphopeptides with respect to its sum of variance components. The standardized variance components are the proportion of the total variation in each phosphopeptides attributed to individuals, biological replicates, and technical replicates. 


# Boxplots of the standardized VCs confirm our observations from the raw VC values. Proportion of variability attributed to biological replicates is the smallest, followed by technical replicates, with individaul samples contributing the largest portion of variabilty in expression levels. 
par(mfrow = c(1,1))
varprop = Varcomp/rowSums(Varcomp)

labs = c("individual","biorep","tech")
boxplot((varprop),axes=F)
axis(1,at=c(1,2,3),labels=labs,col="white");axis(2)


# Heatmap representation of the standardized VCs. 
require(gplots)
require(RColorBrewer)
colnames(varprop) = c("individual","bio","tech")
heatmap.2(as.matrix(varprop),
          col=brewer.pal(9,"YlGnBu"),
          Colv=F,
          labRow="",
          trace="none",
          srtCol=45,  ,adjCol = c(1,1),
          margins = c(6,5),
          cexCol=1.5,
          key.xlab = "Standardized VC", key.ylab=NULL, key.title = "",
          )




#```{r}
# Some phospeptides with large biological variability. 
# melted$techrep = as.factor(unlist(melted$techrep))
# ii = rownames(varprop)[which(rank(varprop[,2])<10)]
# i=1
# par(mfrow=c(2,2))
# # 
# foo = melted[melted$Var1==ii[i],]
# grand = mean(foo$value)
# bio=aggregate(value ~ biorep,data=foo,FUN=mean)
# tech=aggregate(value ~ techrep,data=foo,FUN=mean)
# plot(rep(1,5),c(grand,bio$value,tech$value))
# # 
# boxplot(foo$value~foo$biorep)


#assign flag to varcomp file according to four categories. high ind/high biological. high ind/low bio. low ind/high bio and low ind/low bio
#cutoff for high/low individual is log10 = -5. cutoff for high/low biological variance is log10 = -6.
high_ind_var <- ifelse(log10(Varcomp[,1]) >= -5, "+", "-")
low_ind_var <- ifelse(log10(Varcomp[,1]) < -5, "+", "-")
high_bio_var <- ifelse(log10(Varcomp[,2]) >= -6, "+", "-")
low_bio_var <- ifelse(log10(Varcomp[,2]) < -6, "+", "-")
Varcomp <- cbind(Varcomp, high_ind_var, low_ind_var, high_bio_var, low_bio_var)
  
  }
  
  else{
    #fit a nested random effects model with a single random effect term
    sites <- c()
    Varcomp <- c()
    Expindvar <- c()
    Expindvarmeans <- c()
    Expbiovar <- c()
    Expbiovarmeans <- c()
    Exptechvar <- c()
    #new loop that gives a model for each unique phosphosite. Fits biorep as a character! also added an error catch
    for(id in levels(melted$Var1)){
      test <- melted[melted$Var1 %in% id,]
      test1 <- test[3:6]
      test1$biorep <- as.factor(test1$biorep)
      pos_err <- tryCatch(lme(value~1, data=test1, random =~1|individual),error=function(e) e)
      if(!inherits(pos_err, "error")){
        lmemodel <- lme(value~1, data=test1, random =~1|individual)
        temp <- as.numeric(VarCorr(lmemodel)[,1])
        temp <- temp[!is.na(temp)]
        temp <- na.omit(temp)
        Varcomp <- cbind(Varcomp,temp)
        sites <- c(sites,as.character(unique(test$Var1)))
      }
    }
    #plot the variance component distributions
    colnames(Varcomp) <- sites
    row.names(Varcomp) <- c("individual","residual")
    dim(Varcomp)
    Varcomp <- t(Varcomp)
    boxplot(log10(Varcomp), ylab = "log10 variance component")
    summary(Varcomp)
    head(Varcomp)
    
    ##cumulative Varcomp
    Varcomp2 <- Varcomp
    IndCum <- rowSums(Varcomp2)
    Varcomp2[,1] <- 0
    tech <- rowSums(Varcomp2)
    CumMat <- cbind(IndCum, tech)
    boxplot(log10(CumMat), ylab = "log10 cumulative variance components")
    
    ##Ratio of component estimates relative to technical
    Ind_tech <- Varcomp[,1]/Varcomp[,2]
    boxplot(log10(Ind_tech), ylab = "log10 ratios of variance components")
    
    ##scatter plots show discrete groupings?
    plot(log10(Varcomp[,1]),log10(Varcomp[,2]), xlab = "individual variance", ylab = "technical variance")
    
    #histograms of log10 variance are bimodal for biological and individual variance.
    plot(density(log10(Varcomp)), xlab = "variance", main = "hist of variance components")
    plot(density(log10(Varcomp[,1])), xlab = "variance", main = "hist of individual variance components")
    plot(density(log10(Varcomp[,2])), xlab = "variance", main = "hist of technical variance components")    
  }
  return(Varcomp)
}
  
  
  
  
  