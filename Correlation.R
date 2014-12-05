expr_gene = read.table('Pilot_geneexpr_noiPSC.txt', header=T, as.is=T)
expr_quant.all = read.table('Pilot_probeexpr_noiPSC.txt', header=T, as.is=F)

library(gplots)

#To get the correlation between replicates vs non-replicates
cor.g <- cor(gene,method="spearman", use="complete.obs")
replicates=c()
non.replicates=c()
samplenames = names[,9]
for(i in 1:35){
  for(j in (i+1):36){
    if(samplenames[i]==samplenames[j]){
      replicates=c(replicates,cor.g[i,j])
    }
    else{
      non.replicates=c(non.replicates,cor.g[i,j])
    }
  }
}

boxplot.n(replicates,non.replicates, main = "Correlation of Samples", ylab = "Correlation", xlab = 'Sample Relation')          

boxplot.n(replicates.tech, replicates.bio, replicates.indiv.card, non.replicates, main = "Correlation of Samples", ylab = "Correlation", xlab = 'Sample Relation')          

replicates.tech = replicates
nonreplicates.tech =non.replicates
replicates.bio =replicates
non.replicatesbio = non.replicates
replicates.indiv.card = replicates
non.replicates.indiv.card = non.replicates
replicates.indiv = replicates
non.replicates.indiv = non.replicates
boxplot.n(replicates.tech, replicates.bio, replicates.indiv.card, replicates.indiv,non.replicatesbio, non.replicates.indiv, main = "Correlation of Samples", ylab = "Correlation", xlab = 'Sample Relation')          

cor.bc <- cor(expr_gene[,c(1,9)],method="spearman", use="complete.obs")
samplenames = read.csv('SampleKey_PrunednoiPSC.csv', as.is=T)

##After averaging tech replicates
cor.g.t <- cor(gene,method="spearman", use="complete.obs")
replicates=c()
non.replicates=c()
samplenames = names[,10]
for(i in 1:32){
  for(j in (i+1):33){
    if(samplenames[i]==samplenames[j]){
      replicates=c(replicates,cor.g.t[i,j])
    }
    else{
      non.replicates=c(non.replicates,cor.g.t[i,j])
    }
  }
}

##Across Indiv Cont + 40
cor.g <- cor(gene.clean,method="spearman", use="complete.obs")
replicates=c()
samplenames = names[,9]
naming = c()
thislist = c(16,24,10,17,18,25)
z=1
for(g in 1:length(thislist)){
  indiv = thislist[g]
  for(i in 1:33){
    for(j in 1:33){
      if(samplenames[i]== paste(indiv, "0", sep="") & (as.numeric(samplenames[j]) < 1000) & samplenames[j] != paste(indiv,"0",sep="")){
        replicates=c(replicates,cor.g[i,j])
        naming[z] = paste(samplenames[i],samplenames[j],sep="&")
        z = z+1
      }
      if(samplenames[i]== paste(indiv, "40", sep="") & (as.numeric(samplenames[j]) > 1000) & ((as.numeric(samplenames[j])-40) %% 100 == 0) & samplenames[j] != paste(indiv,"40",sep="")){
        replicates=c(replicates,cor.g[i,j])
        naming[z] = paste(samplenames[i],samplenames[j],sep="&")
        z = z+1
      }
    }
  }
}
reps.across.indiv.treat =replicates

##Across Indiv across treat (cont + 40)
cor.g <- cor(gene.clean,method="spearman", use="complete.obs")
replicates=c()
samplenames = names[,9]
naming = c()
thislist = c(16,24,10,17,18,25)
z=1
for(g in 1:length(thislist)){
  indiv = thislist[g]
   for(i in 1:33){
     for(j in 1:33){
       if(samplenames[i]== paste(indiv, "0", sep="") & (as.numeric(samplenames[j]) > 1000) & ((as.numeric(samplenames[j])-40) %% 100 == 0) & samplenames[j] != paste(indiv,"40",sep="")){
         replicates=c(replicates,cor.g[i,j])
         naming[z] = paste(samplenames[i],samplenames[j],sep="&")
         z = z+1
      }
  
    }
  }
}

reps.across.indiv =replicates
boxplot.n(replicates,main = "Correlation of Samples", ylab = "Correlation", xlab = 'Sample Relation')          

##Within Indiv
cor.g <- cor(gene.clean,method="spearman", use="complete.obs")
replicates=c()
samplenames = names[,9]
  for(g in 1:5){
      for(i in 1:33){
        for(j in 1:33){
          indiv = list[g]
            if(samplenames[i]== paste(indiv, "0", sep="") & samplenames[j] == paste(indiv, "40", sep="")){
              replicates=c(replicates,cor.g[i,j])
              print(paste(samplenames[i],samplenames[j],sep="&"))
            }
          
        }
      }
    }
reps.within.indiv =replicates

boxplot.n(replicates.tech, replicates.bio, reps.across.indiv.treat, reps.within.indiv, reps.across.indiv, main = "Correlation of Samples", ylab = "Correlation", xlab = 'Sample Relation')          