#diffphos ProtLevel enrichment analysis - (all enrichment analysis at the protein group level)

#diffphos motif enrichment analysis in separate script. Contrast specific (volcanos, overlap venns, hprd enrichment wordclouds, plogo charts and frequency motifs)

#snp/indel enrichments in separate script


#packages
require(plyr)
require(dplyr)
require(ggplot)
library(gplots)

##########Enrichment of highly expressed proteins in diffphos? ----

##ibaq v p.value density scatter with regression line + significance overlay

#get median ibaq and median rank, calculate relative rank.
ibaq <- multExpanded1_withDE_annotated[multExpanded1_withDE_annotated$GelPrepNormSubtoDE == "+",
                                       c("ppiBAQ.L.18486", "ppiBAQ.L.18862", "ppiBAQ.L.19160", "GelPrepNormFPval")] #note the NAs
ibaq$ibaq.median <- apply(as.matrix(ibaq[,1:3]), 1, median)
ibaq <- lapply(ibaq, as.numeric)

y <- -log10(ibaq$GelPrepNormFPval)
x <- log10(ibaq$ibaq.median)

R <- cor(x,y, use = "complete.obs")
R
cor.test(x,y)$p.value

#make and save plot
pdf("ibaq_pvalue_density.pdf", 7, 5)
smoothScatter(x,y, nbin = 150, bandwidth = 0.1,
              cex = .3,
              xlab = expression(log[10](iBAQ~abundance~estimate)),
              ylab = expression(-log[10](P~value)), lwd = 10
              )
reg.line <- lm(y~x, na.action = "na.omit")
abline(reg.line, lwd = 1.5, lty = 2)
text(8.8, 10.2, expression(R^2 == -.03), col = "darkred", cex = 1) # rsquared and pvalue
text(8.8, 9.4, expression(p == .1), col = "darkred", cex = 1)
dev.off()





##reactome and GO enrichments ----
source("Enrichment.R")

enrichment_tables <- Enrichment(multExpanded1_withDE_annotated)
dir.create("./Enrichments") 
for(i in 1:length(enrichment_tables)){
  write.table(enrichment_tables[[i]], file.path(getwd(), paste0("Enrichments/", names(enrichment_tables)[i], ".txt")),
              sep = "\t", row.names = F)
}

#output for cytoscape


##barplot for connectivity, % Disorder, extent of PTM modifications ----
# FE test 'with domains', with phosphorelevant domains. These are binary assignments and the latter should be tested with using a glm
# and pfam domain positive? as a covariate

#three barplots produced at equal size to be joined with legend in illustrator.
# 
####counfounded bars ----
# confounded.data <- multExpanded1_withDE_annotated[, c("ConfoundedSubtoDE", "ConfoundedglobalFsig", "ConfoundedFAdjPval",
#                                                       "confoundedPFamIDs", "confoundedPFamIDPhospho",
#                                                        "confoundedInteractCount", "confoundedPercentDisorder", "total.mod.count.confounded")]
# #note the factors. Revert. remember a dataframe is a list of vectors.
# str(confounded.data)
# i <- sapply(confounded.data, is.factor)
# confounded.data[i] <- lapply(confounded.data[i], as.character)
# 
# #subset to those phosphopeptides subjected to diffphos analysis
# confounded.data <- confounded.data[confounded.data$ConfoundedSubtoDE == "+",]
# 
# # plot(-log10(as.numeric(confounded.data$confoundedInteractCount)), log10(as.numeric(confounded.data$ConfoundedFAdjPval)))
# 
# sig.index <- confounded.data$ConfoundedglobalFsig == "+"
# 
# #function for mean and se calculations
# mean.se.calc <- function(pos.set, neg.set){
#   m1 <- mean(pos.set[is.finite(pos.set)], na.rm = T)
#   sem1 <- sd(pos.set[is.finite(pos.set)])/ sqrt(sum(is.finite(pos.set)))
#   m2 <- mean(neg.set[is.finite(neg.set)], na.rm = T)
#   sem2 <- sd(neg.set[is.finite(neg.set)])/ sqrt(sum(is.finite(neg.set)))
#   list(means = c(m1, m2), stderr = c(sem1, sem2))
# }
# 
# ##connectivity with biogrid
# 
# #these data are clearly not normally distributed (shapiro.test). Test for differences
# # For multiple groups Note there is the kruskal.test and the anderson darling test with package "KSamples"
# 
# 
# interactions.pos <- as.numeric(confounded.data[sig.index, "confoundedInteractCount"])
# interactions.neg <- as.numeric(confounded.data[!sig.index, "confoundedInteractCount"])
# 
# #this difference is significant. proteins with less interactions (according to biogrid) are slightly de-enriched.
# wilcox.test(interactions.pos, interactions.neg, alternative="two.sided")$p.value#performs mann-whitney test
# # ks.test(interactions.pos, interactions.neg, alternative="two.sided")
# 
# data <- unlist(mean.se.calc(interactions.pos, interactions.neg))
# 
# 
# pdf("confoundedInteractions.pdf", 5, 7)
# barplot2(data[1:2], plot.ci=T,
#          ci.l=c(data[1] - data[3], data[2] - data[3]), ci.u=c(data[1] + data[3], data[2] + data[3]),
#          lwd=3, ci.lwd=3, col = c("Red", "Blue"), axisnames = F,
#          main = "Number of protein-protein interactions")
# dev.off()
# 
# 
# 
# ##Disorder
# confounded.data$FinalPercentDisorder <- sapply(confounded.data$confoundedPercentDisorder, function(x){
#   mean(as.numeric(unlist(strsplit(x, ";"))))
# })
# 
# disorder.pos <- as.numeric(confounded.data[sig.index, "FinalPercentDisorder"])
# disorder.neg <- as.numeric(confounded.data[!sig.index, "FinalPercentDisorder"])
# wilcox.test(disorder.pos, disorder.neg, alternative="two.sided")$p.value#performs mann-whitney test
# 
# 
# data <- unlist(mean.se.calc(disorder.pos,disorder.neg))
# pdf("confoundedDisorder.pdf", 5, 7)
# barplot2(data[1:2], plot.ci=T,
#          ci.l=c(data[1] - data[3], data[2] - data[3]), ci.u=c(data[1] + data[3], data[2] + data[3]),
#          lwd=3, ci.lwd=3, col = c("Red", "Blue"), axisnames = F,
#          main = "Percent disorder")
# dev.off()
# 
# ##modification count
# 
# modcount.pos <- as.numeric(confounded.data[sig.index, "total.mod.count.confounded"])
# modcount.neg <- as.numeric(confounded.data[!sig.index, "total.mod.count.confounded"])
# wilcox.test(modcount.pos[is.finite(modcount.pos)], modcount.neg[is.finite(modcount.neg)], alternative="two.sided")$p.value#performs mann-whitney test
# 
# 
# data <- unlist(mean.se.calc(modcount.pos, modcount.neg))
# pdf("confoundedModcount.pdf", 5, 7)
# barplot2(data[1:2], plot.ci=T,
#          ci.l=c(data[1] - data[3], data[2] - data[3]), ci.u=c(data[1] + data[3], data[2] + data[3]),
#          lwd=3, ci.lwd=3, col = c("Red", "Blue"), axisnames = F,
#          main = "Total number of PTMs")
# dev.off()
# 
# 
# #all of these results have been the opposite of my expectation...
# 
# ##FE test for pfam domains and phospho domains (separately for now)
# 
# #1st enrichment of proteins containing pfam domains
# confounded.data$pfam.pos <- sapply(confounded.data$confoundedPFamIDs, function(x){
#   any(unlist(strsplit(x, ";")) != "NA")
# })
# 
# #row 1 diffphos/withpfam AND diffphos/withoutpfam
# #row 2 !diffphos/withpfam AND !diffphos/withoutpfam
# 
# row1 <- c(nrow(confounded.data[confounded.data$ConfoundedglobalFsig == "+" & confounded.data$pfam.pos == T, ]),
#           nrow(confounded.data[confounded.data$ConfoundedglobalFsig == "+" & confounded.data$pfam.pos == F, ])
# )
# row2 <- c(nrow(confounded.data[confounded.data$ConfoundedglobalFsig == "-" & confounded.data$pfam.pos == T, ]),
#           nrow(confounded.data[confounded.data$ConfoundedglobalFsig == "-" & confounded.data$pfam.pos == F, ])
# )
# 
# #FEtest. no sig enrich or depletion
# contmatrix <- rbind(row1,row2)
# result <- fisher.test(contmatrix, alternative = "g")
# result$p.value
# row1 
# 0.8355088 
# 
# #second is enrichment for proteins containing phospho relevant domains
# row1 <- c(nrow(confounded.data[confounded.data$ConfoundedglobalFsig == "+" & confounded.data$confoundedPFamIDPhospho == "yes", ]),
#          nrow(confounded.data[confounded.data$ConfoundedglobalFsig == "+" & confounded.data$confoundedPFamIDPhospho == "no", ])
# )
# row2 <- c(nrow(confounded.data[confounded.data$ConfoundedglobalFsig == "-" & confounded.data$confoundedPFamIDPhospho == "yes", ]),
#           nrow(confounded.data[confounded.data$ConfoundedglobalFsig == "-" & confounded.data$confoundedPFamIDPhospho == "no", ])
# )
# 
# #FEtest. Again no sig enrich or depletion
# contmatrix <- rbind(row1,row2)
# result <- fisher.test(contmatrix, alternative = "g")
# result$p.value
# 0.4332596



#### Gelprot normalized bars ----

GelPrep.data <- multExpanded1_withDE_annotated[, c("GelPrepNormSubtoDE", "GelPrepNormglobalFsig", "GelPrepNormFAdjPval",
                                                      "GelPrepPFamIDs", "GelPrepPFamIDPhospho",
                                                      "GelPrepInteractCount", "GelPrepPercentDisorder", "total.mod.count.GelPrep")]
#note the factors. Revert. remember a dataframe is a list of vectors.
str(GelPrep.data)
i <- sapply(GelPrep.data, is.factor)
GelPrep.data[i] <- lapply(GelPrep.data[i], as.character)

#subset to those phosphopeptides subjected to diffphos analysis (n = 3257)
GelPrep.data <- GelPrep.data[GelPrep.data$GelPrepNormSubtoDE == "+",]

sig.index <- GelPrep.data$GelPrepNormglobalFsig == "+"


##connectivity with biogrid

#these data are clearly not normally distributed (shapiro.test). Test for differences
# For multiple groups Note there is the kruskal.test and the anderson darling test with package "KSamples"


interactions.pos <- as.numeric(GelPrep.data[sig.index, "GelPrepInteractCount"])
interactions.neg <- as.numeric(GelPrep.data[!sig.index, "GelPrepInteractCount"])

#this difference is significant. proteins with less interactions (according to biogrid) are slightly de-enriched.
wilcox.test(interactions.pos, interactions.neg, alternative="two.sided")$p.value#performs mann-whitney test
# ks.test(interactions.pos, interactions.neg, alternative="two.sided")

data <- unlist(mean.se.calc(interactions.pos, interactions.neg))
pdf("GelPrepInteractions.pdf", 5, 7)
barplot2(data[1:2], plot.ci=T,
         ci.l=c(data[1] - data[3], data[2] - data[3]), ci.u=c(data[1] + data[3], data[2] + data[3]),
         lwd=3, ci.lwd=3, col = c("Red", "Blue"), axisnames = F,
         main = "Number of protein-protein interactions")
dev.off()



##Disorder
GelPrep.data$FinalPercentDisorder <- sapply(GelPrep.data$GelPrepPercentDisorder, function(x){
  mean(as.numeric(unlist(strsplit(x, ";"))))
})

disorder.pos <- as.numeric(GelPrep.data[sig.index, "FinalPercentDisorder"])
disorder.neg <- as.numeric(GelPrep.data[!sig.index, "FinalPercentDisorder"])
wilcox.test(disorder.pos, disorder.neg, alternative="two.sided")$p.value#performs mann-whitney test


data <- unlist(mean.se.calc(disorder.pos,disorder.neg))
pdf("GelPrepDisorder.pdf", 5, 7)
barplot2(data[1:2], plot.ci=T,
         ci.l=c(data[1] - data[3], data[2] - data[3]), ci.u=c(data[1] + data[3], data[2] + data[3]),
         lwd=3, ci.lwd=3, col = c("Red", "Blue"), axisnames = F,
         main = "Percent disorder")
dev.off()

##modification count
modcount.pos <- as.numeric(GelPrep.data[sig.index, "total.mod.count.GelPrep"])
modcount.neg <- as.numeric(GelPrep.data[!sig.index, "total.mod.count.GelPrep"])
wilcox.test(modcount.pos[is.finite(modcount.pos)], modcount.neg[is.finite(modcount.neg)], alternative="two.sided")$p.value#performs mann-whitney test


data <- unlist(mean.se.calc(modcount.pos, modcount.neg))
pdf("GelPrepModcount.pdf", 5, 7)
barplot2(data[1:2], plot.ci=T,
         ci.l=c(data[1] - data[3], data[2] - data[3]), ci.u=c(data[1] + data[3], data[2] + data[3]),
         lwd=3, ci.lwd=3, col = c("Red", "Blue"), axisnames = F,
         main = "Total number of PTMs")
dev.off()

#all non significant changes

##FE test for pfam domains and phospho domains (separately for now)

#1st enrichment of proteins containing pfam domains
GelPrep.data$pfam.pos <- sapply(GelPrep.data$GelPrepPFamIDs, function(x){
  any(unlist(strsplit(x, ";")) != "NA")
})

#row 1 diffphos/withpfam AND diffphos/withoutpfam
#row 2 !diffphos/withpfam AND !diffphos/withoutpfam

row1 <- c(nrow(GelPrep.data[GelPrep.data$GelPrepNormglobalFsig == "+" & GelPrep.data$pfam.pos == T, ]),
          nrow(GelPrep.data[GelPrep.data$GelPrepNormglobalFsig == "+" & GelPrep.data$pfam.pos == F, ])
)
row2 <- c(nrow(GelPrep.data[GelPrep.data$GelPrepNormglobalFsig == "-" & GelPrep.data$pfam.pos == T, ]),
          nrow(GelPrep.data[GelPrep.data$GelPrepNormglobalFsig == "-" & GelPrep.data$pfam.pos == F, ])
)

#FEtest. A significant *depletion* of domain containing proteins in the diffphos list
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix)
result$p.value
[1] 0.002440458

#second is enrichment for proteins containing phospho relevant domains
row1 <- c(nrow(GelPrep.data[GelPrep.data$GelPrepNormglobalFsig == "+" & GelPrep.data$GelPrepPFamIDPhospho == "yes", ]),
          nrow(GelPrep.data[GelPrep.data$GelPrepNormglobalFsig == "+" & GelPrep.data$GelPrepPFamIDPhospho == "no", ])
)
row2 <- c(nrow(GelPrep.data[GelPrep.data$GelPrepNormglobalFsig == "-" & GelPrep.data$GelPrepPFamIDPhospho == "yes", ]),
          nrow(GelPrep.data[GelPrep.data$GelPrepNormglobalFsig == "-" & GelPrep.data$GelPrepPFamIDPhospho == "no", ])
)

#FEtest. Here there is a *slightly* significant enrichment of phospho relevant domains in the diffphos dataset!
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g")
result$p.value
row1 
0.03180943 

















