GOpath <- function(multExpandedwithDE){
  #I will need to decide on the imput matrix. This function recieves Uniprot annotations and significance binary indicators.
#   It then performs enrichment analysis using the zia proteome (I think) as background. 
#   GO biological process and reactome pathway enrichmet is performed.
#   Open question is whether or not subtoDE should be the background or ziaproteins should be the background.
#   I have already performed a selection with the subtoDE subset that is driven by the data so I think enrichment should be recognized within
#   this subset because enrichment should be analyzed relative to what 'could be' enriched.
#   However if there is no bias in the selection of the subtoDE set then this should be moot. however there will be bias due to concentration dependent identification.Whatever the bias (ID or quantification), if there is no bias in the expression level and length of the protein in determining differential phosphorylation (biologically small and large proteins are equally probable of being differentially phosphorylated) However this must take into account the length of the protein. 
  Did Banovich do this with the methylation data? any sort of enrichment analysis?
  
  
  The background needs to be the subtoDE cohort.
#will a multiply phosphorylated protein be more likely to be represented in diffphos than singly phosphorylated proteins? Is the degree of 

#converting uniprotIDs to entrezIDs
# biocLite("mygene")
# library(mygene)

source("uniprot to entrezID.R")
# source("http://bioconductor.org/biocLite.R")
# biocLite("biomaRt")
require("biomaRt")
# ##subset of DE in any contrast
# DE <- multExpanded1[(multExpanded1$DEcont1 == "+"| multExpanded1$DEcont2 == "+"| multExpanded1$DEcont3 == "+"),]
# #convert row of uniprot ids to ensemble ids
# test <- strsplit(as.character(DE$Proteins), ";")
# test <- as.character(unlist(test))
# test <- unique(test)
# out <- Uniprot2EG(test)
# DEentrez <- as.character(out$entrezgene)
# DEentrez <- unique(DEentrez)

##subset of global DE using F statistics
DE <- multExpanded1[multExpanded1$globalFsig == "+",]
#convert row of uniprot ids to ensemble ids
MajorityProteins <- strsplit(as.character(DE$Proteins), ";")
MajorityProteins <- as.character(unlist(test))
#test <- unique(test)
out <- Uniprot2EG(MajorityProteins)
##add (potentially redundant)


DEentrezwithdup <- as.character(out$entrezgene)
#DEentrez <- unique(DEentrez)


##now to add the enrezid next to the 'test' vector if it is found in the enrezid list.
poopy <- out[which(duplicated(out$uniprot)),1]
write.table(poopy, "duplicateduniprot.csv", sep=',',row.names=FALSE, col.names=F)


##subset of DEcont1
DE1 <- multExpanded1[multExpanded1$DEcont1 == "+",]
#convert row of uniprot ids to ensemble ids
test <- strsplit(as.character(DE1$Proteins), ";")
test <- as.character(unlist(test))
test <- unique(test)
out <- Uniprot2EG(test)
DE1entrez <- as.character(out$entrezgene)
DE1entrez <- unique(DE1entrez)

##subset of DEcont2
DE2 <- multExpanded1[multExpanded1$DEcont2 == "+",]
#convert row of uniprot ids to ensemble ids
test <- strsplit(as.character(DE2$Proteins), ";")
test <- as.character(unlist(test))
test <- unique(test)
out <- Uniprot2EG(test)
DE2entrez <- as.character(out$entrezgene)
DE2entrez <- unique(DE2entrez)

##subset of DEcont3
DE3 <- multExpanded1[multExpanded1$DEcont3 == "+",]
#convert row of uniprot ids to ensemble ids
test <- strsplit(as.character(DE3$Proteins), ";")
test <- as.character(unlist(test))
test <- unique(test)
out <- Uniprot2EG(test)
DE3entrez <- as.character(out$entrezgene)
DE3entrez <- unique(DE3entrez)


# sig1 - 18862-18486
# sig2 - 19160-18862
# sig3 - 19160-18486


DEcontrasts <- list(DE1entrez,DE2entrez,DE3entrez)
names(DEcontrasts) <- c("18862-18486","19160-18862","19160-18486")

##directional subsetting**************************************************************

c1up <- multExpanded1[multExpanded1$cont1up == "+",]
#convert row of uniprot ids to ensemble ids
test <- strsplit(as.character(c1up$Proteins), ";")
test <- as.character(unlist(test))
test <- unique(test)
out <- Uniprot2EG(test)
c1upentrez <- as.character(out$entrezgene)
c1upentrez <- unique(c1upentrez)

c1down <- multExpanded1[multExpanded1$cont1down == "+",]
#convert row of uniprot ids to ensemble ids
test <- strsplit(as.character(c1down$Proteins), ";")
test <- as.character(unlist(test))
test <- unique(test)
out <- Uniprot2EG(test)
c1downentrez <- as.character(out$entrezgene)
c1downentrez <- unique(c1downentrez)

c2up <- multExpanded1[multExpanded1$cont2up == "+",]
#convert row of uniprot ids to ensemble ids
test <- strsplit(as.character(c2up$Proteins), ";")
test <- as.character(unlist(test))
test <- unique(test)
out <- Uniprot2EG(test)
c2upentrez <- as.character(out$entrezgene)
c2upentrez <- unique(c2upentrez)

c2down <- multExpanded1[multExpanded1$cont2down == "+",]
#convert row of uniprot ids to ensemble ids
test <- strsplit(as.character(c2down$Proteins), ";")
test <- as.character(unlist(test))
test <- unique(test)
out <- Uniprot2EG(test)
c2downentrez <- as.character(out$entrezgene)
c2downentrez <- unique(c2downentrez)

c3up <- multExpanded1[multExpanded1$cont3up == "+",]
#convert row of uniprot ids to ensemble ids
test <- strsplit(as.character(c3up$Proteins), ";")
test <- as.character(unlist(test))
test <- unique(test)
out <- Uniprot2EG(test)
c3upentrez <- as.character(out$entrezgene)
c3upentrez <- unique(c3upentrez)

c3down <- multExpanded1[multExpanded1$cont3down == "+",]
#convert row of uniprot ids to ensemble ids
test <- strsplit(as.character(c3down$Proteins), ";")
test <- as.character(unlist(test))
test <- unique(test)
out <- Uniprot2EG(test)
c3downentrez <- as.character(out$entrezgene)
c3downentrez <- unique(c3downentrez)

# sig1 - 18862-18486
# sig2 - 19160-18862
# sig3 - 19160-18486

"18862-18486","19160-18862","19160-18486"

Directioncontrasts <- list(c1upentrez,c1downentrez,c2upentrez,c2downentrez,c3upentrez,c3downentrez)
names(Directioncontrasts) <- c("18862-18486up","18862-18486down","19160-18862up","19160-18862down","19160-18486up","19160-18486down")
#names(Directioncontrasts) <- c("C1up","C1down","C2up","C2down","C3up","C3down")



#a background set for enrichment analysis the yoruba LCLs using zia's data. Not limiting background set to those data subjected to DE, assume this sampling to be random? But it isn't I don't think...
#used bioDBnet to convert ensemble to entrez gene ids and kept all unique (n=10827)
ziac <- read.table("E:/My Documents/Dropbox/Postdoc-Gilad/Zia/ziaConversion.txt", sep = "\t", header=T, fill = T)
entrezbg <- strsplit(as.character(ziac$Gene.ID), ";")
entrezbgdups <- as.character(unlist(entrezbg))
entrezbg <- unique(entrezbgdups)
####################################################################################

require(ReactomePA)

##all DE pathway enrichment using the F-statistics. pvalue cutoff is the adjusted p value using BH! minimum geneset and minimum geneset size of 5. uses the union of annotated genes and the genes expressed in the LCLs.
Reactome <- enrichPathway(gene=DEentrez,pvalueCutoff=0.05, readable=T, universe = entrezbg)#see help file. uses hypergeometric distribution.
ReactomeDEenrich <- as.data.frame(summary(Reactome))

##some tests here. First with no duplication
Reactome <- enrichPathway(gene=DEentrez,pvalueCutoff=0.05, readable=T)#see help file. uses hypergeometric distribution.
ReactomeDEenrich <- as.data.frame(summary(Reactome))

##some tests here. With duplication
Reactomedups <- enrichPathway(gene=DEentrezwithdup,pvalueCutoff=0.05, readable=T)#see help file. uses hypergeometric distribution.
ReactomeDEenrichdups <- as.data.frame(summary(Reactome))

#seems to return the same result


#lets look at the source
capture.output(getAnywhere('enrichPathway'), file='source_enrichPathway.r')


capture.output(getAnywhere('enrich.internal'), file='enrich.internal.r')

EXTID2TERMID
capture.output(getAnywhere('EXTID2TERMID'), file='enrich.internal.r')



#now for each contrast
ReactomeDE1 <- enrichPathway(gene=DE1entrez,pvalueCutoff=0.05, readable=T, universe = entrezbg)#see help file. uses hypergeometric distribution.
ReactomeDEenrichDE1 <- as.data.frame(summary(ReactomeDE1))



#MF molecular function BP biological process CC cellular component
GOMF <- enrichGO(gene=DEentrez, pvalueCutoff=0.05, readable=T, ont = "MF", universe = entrezbg)
GOMFDE <- as.data.frame(summary(GOMF))
GOBP <- enrichGO(gene=DEentrez, pvalueCutoff=0.05, readable=T, ont = "BP", universe = entrezbg)
GOBPDE <- as.data.frame(summary(GOBP))
GOCC <- enrichGO(gene=DEentrez, pvalueCutoff=0.05, readable=T, ont = "CC", universe = entrezbg)
GOCCDE <- as.data.frame(summary(GOCC))


# head(summary(Reactome))#note that the number of entrez ids has been cut down to 273

# barplot(Reactome,showCategory=8)#showCategory is the number of categories to show

#output to data frame for further analysis
ReactomeDEenrich <- as.data.frame(summary(Reactome))

enrichMap(Reactome, fixed =F)#not sure what this 
enrichMap(GOMF, fixed=F)


#compare enrichment across clusters###############################################

I think it makes the most sense to stick with just the global DE for now
require(clusterProfiler)
#first global DE. Can't use showCategory
res <- compareCluster(DEcontrasts, organism="human", fun = "enrichPathway", pvalueCutoff=0.05,universe = entrezbg)
plot(res,showCategory=NULL)#saved
plot(res,showCategory=NULL, by="rowPercentage")
plot(res,showCategory=10, by="rowPercentage")
plot(res,showCategory=10, by="count")
plot(res,showCategory=3, by="geneRatio")#default

res <- compareCluster(DEcontrasts, organism="human", fun = "enrichGO", pvalueCutoff=0.05,universe = entrezbg)
plot(res,showCategory=NULL)#saved


#res <- compareCluster(DEcontrasts, organism="human", fun = "enrichPathway", pvalueCutoff=0.05)
#plot(res)
res <- compareCluster(DEcontrasts, organism="human", fun = "enrichGO", pvalueCutoff=0.05, universe = entrezbg)
plot(res)
res <- compareCluster(DEcontrasts, organism="human", fun = "enrichGO", ont="BP",pvalueCutoff=0.05, universe = entrezbg)
plot(res)

#res <- compareCluster(DEcontrasts, organism="human", fun = "groupGO")
#plot(res)
res <- compareCluster(DEcontrasts, organism="human", fun = "enrichKEGG", universe = entrezbg)
plot(res)




#now directional DE
res <- compareCluster(Directioncontrasts, organism="human", fun = "enrichPathway", pvalueCutoff=0.05,universe = entrezbg)
plot(res)
#res <- compareCluster(DEcontrasts, organism="human", fun = "enrichPathway", pvalueCutoff=0.05)
#plot(res)
res <- compareCluster(Directioncontrasts, organism="human", fun = "enrichGO", pvalueCutoff=0.05, universe = entrezbg)
plot(res)
#res <- compareCluster(DEcontrasts, organism="human", fun = "groupGO")
#plot(res)
res <- compareCluster(Directioncontrasts, organism="human", fun = "enrichKEGG", universe = entrezbg)
plot(res)







