#diffphos ProtLevel enrichment analysis - (all enrichment analysis at the protein group level)

#diffphos motif enrichment analysis in separate script. Contrast specific (volcanos, overlap venns, hprd enrichment wordclouds, plogo charts and frequency motifs)

#snp/indel enrichments in separate script


#packages
require(plyr)
require(dplyr)
require(gplots)



### relationship between phosphoPvalues and protein length


GelPrep.length <- multExpanded1_withDE_annotated[multExpanded1_withDE_annotated$GelPrepCovSubtoDE == "+",
                                                   c("GelPrepCovglobalFsig", "GelPrepCovFAdjPval", "ppSequence.length")]
y <- -log10(as.numeric(GelPrep.length$GelPrepCovFAdjPval))
x <- log10(GelPrep.length$ppSequence.length)

plot(x,y)
R <- cor(x, y, method = "spearman", use = "complete.obs")
R
# [1] -0.07458249
cor.test(x,y, method = "spearman", exact = F)$p.value
# [1] 2.037569e-05

#make and save plot if necessary
# pdf("nscount_pvalue_density.pdf", 7, 5)
smoothScatter(x,y, nbin = 150, bandwidth = 0.1,
              cex = .3,
              pch = 19, nrpoints = .15*length(x),
              colramp = colorRampPalette(c("white", "light gray", "dark gray", "red")),
              xlab = expression(log[10](Protein~Length)),
              ylab = expression(-log[10](P~value)), lwd = 10,
              family = "serif"
)
reg.line <- lm(y~x, na.action = "na.omit")
abline(reg.line, lwd = 2, lty = 2)
# text(3, 7.25, expression(R == -.12), col = "darkred", cex = 1, family = "serif") # rsquared and pvalue
# text(3, 6.85, expression(p == 9.90e-06), col = "darkred", cex = 1, family = "serif")
# dev.off()


#relationship between protein length and germline variance component estimates

# Here I use the variance componnent estimates derived from "mcmcVarcomp.proteinCov". (batch corrected phospho and pQTL regressed gelprot protein estimates)

tmp <- multExpanded1_withDE_annotated[, c("idmult", "ppSequence.length")]
length.variance <- merge(tmp, mcmcVarcomp.proteinCov, by.x = "idmult", by.y = "row.names")

y <- log10(length.variance$individual)
x <- log10(length.variance$ppSequence.length)

plot(x,y)
R <- cor(x, y, method = "pearson", use = "complete.obs")
R
cor.test(x,y)$p.value




##########Enrichment of highly expressed proteins in diffphos? ----

##ibaq v p.value density scatter with regression line + significance overlay

#get median ibaq and median rank, calculate relative rank.
ibaq <- multExpanded1_withDE_annotated[multExpanded1_withDE_annotated$GelPrepCovSubtoDE == "+",
                                       c("ppiBAQ.L.18486", "ppiBAQ.L.18862", "ppiBAQ.L.19160", "GelPrepCovFPval")] #note the NAs
ibaq$ibaq.median <- apply(as.matrix(ibaq[,1:3]), 1, median)
ibaq[] <- lapply(ibaq, as.numeric)

y <- -log10(ibaq$GelPrepCovFPval)
x <- log10(ibaq$ibaq.median)
cor(x,y, use = "complete.obs", method = "spearman")
# [1] 0.07929577
cor.test(x,y, alternative = "two", method = "spearman", exact = F)$p.value
# [1] 5.875555e-06


#make and save plot
pdf("ibaq_pvalue_density.pdf", 7, 5)
smoothScatter(x,y, nbin = 150, bandwidth = 0.1,
              cex = .3,
              pch = 19, nrpoints = .15*length(ibaq$GelPrepCovFPval),
              colramp = colorRampPalette(c("white", "light gray", "dark gray", "red")),
              xlab = expression(log[10](iBAQ~abundance~estimate)),
              ylab = expression(-log[10](P~value)), lwd = 10,
              family = "serif"
              )
reg.line <- lm(y~x, na.action = "na.omit")
abline(reg.line, lwd = 2, lty = 2)
text(8.8, 10.2, expression(R == .053), col = "darkred", cex = 1, family = "serif") # rsquared and pvalue
text(8.8, 9.4, expression(p == .003), col = "darkred", cex = 1, family = "serif")
dev.off()


#protein expression level is negatively correlated with the number of sites identified and quantified!
ibaq.sites <- multExpanded1_withDE_annotated[ , c("ppiBAQ.L.18486", "ppiBAQ.L.18862", "ppiBAQ.L.19160", "GelPrepCovFPval", "GelPrepCovSubtoDE", 
                                         "ppMajorityProteinIDs", "ppProteinIDs", "ppSequence.length", "id")] #note the NAs
ibaq.sites$ibaq.median <- apply(as.matrix(ibaq.sites[,1:3]), 1, median)

#x = # of sites per protein
#y = expression level

#each element of x is the number of unique 'id'/unique ppMajorityProteinID (x)

ProtID.sites.expression <- ibaq.sites %>% group_by(ppProteinIDs) %>% summarise(sites = length(unique(id)), expression.level = unique(ibaq.median))



y <- log2(ProtID.sites.expression$sites)
x <- log10(ProtID.sites.expression$expression.level)
cor(x,y, use = "comp")
# -0.1277885
cor.test(x,y, alternative = "two", method = "pearson")$p.value
# [1] 1.337768e-09

pdf("ibaq_sitecount_density.pdf", 7, 5)
smoothScatter(x,y, nbin = 150, bandwidth = 0.1,
              cex = .3,
              pch = 19, nrpoints = .15*length(ibaq$GelPrepCovFPval),
              colramp = colorRampPalette(c("white", "light gray", "dark gray", "red")),
              xlab = expression(log[10](iBAQ~abundance~estimate)),
              ylab = expression(log[2](sites~per~protein)), lwd = 10,
              family = "serif"
)
reg.line <- lm(y~x, na.action = "na.omit")
abline(reg.line, lwd = 2, lty = 2)
text(8.6, 7, expression(R == -0.13), col = "darkred", cex = 1, family = "serif") # rsquared and pvalue
text(8.7, 6.5, expression(p == 1.34e-09), col = "darkred", cex = 1, family = "serif")
dev.off()




#majority protein ids produce the same plot
#remove those identifications that have a ppProteinID assignment but not a ppMajProteinID assignment
index <- which(ibaq.sites$ppMajorityProteinIDs != "")
ibaq.sites <- ibaq.sites[index,]

# MajProtID.sites.expression <- ibaq.sites %>% group_by(ppMajorityProteinIDs) %>% summarise(sites = length(unique(id)), expression.level = unique(ibaq.median))
# 
# x <- log2(MajProtID.sites.expression$sites)
# y <- log10(MajProtID.sites.expression$expression.level)
# plot(x,y)
# reg.line <- lm(y~x, na.action = "na.omit")
# abline(reg.line, lwd = 2, lty = 2)
# cor.test(x,y, alternative = "two", method = "pearson")$p.value


#number of sites identified strongly correlated with protein length 
MajProtID.sites.length <- ibaq.sites %>% group_by(ppMajorityProteinIDs) %>% summarise(sites = length(unique(id)), length = unique(ppSequence.length))

x <- log2(MajProtID.sites.length$length)
y <- log2(MajProtID.sites.length$sites)
cor.test(x,y, alternative = "two", method = "pearson")$p.value
cor(x,y, use = "complete.obs")
# [1] 0.3397038

pdf("length_sites_density.pdf", 7, 5)
smoothScatter(x,y, nbin = 150, bandwidth = 0.1,
              cex = .3,
              pch = 19, nrpoints = .15*length(ibaq$GelPrepCovFPval),
              colramp = colorRampPalette(c("white", "light gray", "dark gray", "red")),
              xlab = expression(log[2](protein~length)),
              ylab = expression(log[2](sites~per~protein)), lwd = 10,
              family = "serif"
)
reg.line <- lm(y~x, na.action = "na.omit")
abline(reg.line, lwd = 2, lty = 2)
text(6.3, 7, expression(R == 0.34), col = "darkred", cex = 1, family = "serif") # rsquared and pvalue
# text(8.7, 6.5, expression(p == 1.34e-09), col = "darkred", cex = 1, family = "serif")
dev.off()



#number of phosphosites is negatively correlated with phosphopeptide variability?

# note there is a splicing factor protein (Serine/arginine repetitive matrix protein 2) with over 200 identified phosphorylation sites!

#add a column with the number of phosphosites (repeated for each protein using Majority Protein ID)
ibaq.sites <- ibaq.sites %>% group_by(ppMajorityProteinIDs) %>% mutate(sites = length(unique(id)))
ibaq.sites.subtoDE <- ibaq.sites[ibaq.sites$GelPrepCovSubtoDE == "+",]
ibaq.sites.subtoDE["GelPrepCovFPval"] <- as.numeric(ibaq.sites.subtoDE$GelPrepCovFPval)

y <- -log10(ibaq.sites.subtoDE$GelPrepCovFPval)
x <- log2(ibaq.sites.subtoDE$sites)
cor.test(x,y, alternative = "two", method = "spearman", exact = F)$p.value
# [1] 4.270045e-06
cor(x,y, method = "spearman", use = "complete.obs")
# [1] -0.08053805


pdf("sites_pval_density.pdf", 7, 5)
smoothScatter(x,y, nbin = 150, bandwidth = 0.1,
              cex = .3,
              pch = 19, nrpoints = .15*length(ibaq$GelPrepCovFPval),
              colramp = colorRampPalette(c("white", "light gray", "dark gray", "red")),
              xlab = expression(log[2](sites~per~protein)),
              ylab = expression(-log[10](P~value)), lwd = 10,
              family = "serif"
)
reg.line <- lm(y~x, na.action = "na.omit")
abline(reg.line, lwd = 2, lty = 2)
text(6.3, 10, expression(R == -0.09), col = "darkred", cex = 1, family = "serif") # rsquared and pvalue
text(6.5, 9.25, expression(p == 6.97e-07), col = "darkred", cex = 1, family = "serif")
dev.off()



#Sequence length is negatively correlated with variability
y <- -log10(ibaq.sites.subtoDE$GelPrepCovFPval)
x <- log10(ibaq.sites.subtoDE$ppSequence.length)
cor.test(x,y, alternative = "two", method = "pearson")$p.value
# [1] 0.007607161
cor(x,y, use = "complete.obs")
# [1] -0.04680295

pdf("length_pval_density.pdf", 7, 5)
smoothScatter(x,y, nbin = 150, bandwidth = 0.1,
              cex = .3,
              pch = 19, nrpoints = .15*length(ibaq$GelPrepCovFPval),
              colramp = colorRampPalette(c("white", "light gray", "dark gray", "red")),
              xlab = expression(log[10](protein~length)),
              ylab = expression(-log[10](P~value)), lwd = 10,
              family = "serif"
)
reg.line <- lm(y~x, na.action = "na.omit")
abline(reg.line, lwd = 2, lty = 2)
text(3.7, 10.75, expression(R == -0.05), col = "darkred", cex = 1, family = "serif") # rsquared and pvalue
text(3.7, 10, expression(p == 7.6e-03), col = "darkred", cex = 1, family = "serif")
dev.off()




#Sequence length is negatively correlated with expression, which explains why number of identifications per protein is negatively correlated with expression
ibaq.sites["ibaq.median"] <- as.numeric(ibaq.sites$ibaq.median)

#phosphosite level plot
x <- log10(ibaq.sites$ppSequence.length)
y <- log10(ibaq.sites$ibaq.median)

#protein level plot
MajProtID.length.expression <- ibaq.sites %>% group_by(ppMajorityProteinIDs) %>% summarise(length = unique(ppSequence.length), expression = unique(ibaq.median))

x <- log10(MajProtID.length.expression$length)
y <- log10(MajProtID.length.expression$expression)
cor.test(x,y, alternative = "two", method = "pearson")$p.value
[1] 2.818012e-183
cor(x,y, use = "complete.obs")
[1] -0.5620399

pdf("length_expression_density.pdf", 7, 5)
smoothScatter(x,y, nbin = 150, bandwidth = 0.1,
              cex = .3,
              pch = 19, nrpoints = .15*length(ibaq$GelPrepCovFPval),
              colramp = colorRampPalette(c("white", "light gray", "dark gray", "red")),
              xlab = expression(log[10](protein~length)),
              ylab = expression(log[10](iBAQ~abundance~estimate)), lwd = 10,
              family = "serif"
)
reg.line <- lm(y~x, na.action = "na.omit")
abline(reg.line, lwd = 2, lty = 2)
text(3.7, 8.5, expression(R == -0.56), col = "darkred", cex = 1, family = "serif") # rsquared and pvalue
# text(3.7, 10, expression(p == 7.6e-03), col = "darkred", cex = 1, family = "serif")
dev.off()



### Number of phosphosites, cut by expression decile and correlation with phosphopeptide variability
require(ggplot2)
require("Hmisc")

ibaq.sites.subtoDE$Ibaqdecile <- cut2(ibaq.sites.subtoDE$ibaq.median, g=10)
ibaq.sites.subtoDE$Ibaqquintile <- cut2(ibaq.sites.subtoDE$ibaq.median, g=5)


y <- -log10(ibaq.sites.subtoDE$GelPrepCovFPval)
x <- log2(ibaq.sites.subtoDE$sites)



qplot(log2(sites), -log10(GelPrepCovFPval), data = ibaq.sites.subtoDE, color = Ibaqquintile)
qplot(log2(sites), -log10(GelPrepCovFPval), data = ibaq.sites.subtoDE, facets = .~Ibaqquintile)


site.variation.quintile <- ggplot(ibaq.sites.subtoDE, aes(x = log2(sites),
                                                          y = -log10(GelPrepCovFPval),
                                                          color = Ibaqquintile)) +
  geom_point() + 
#   facet_grid(~.Ibaqquintile) + 
  geom_smooth(method = "lm")

site.variation.quintile

site.variation.decile <- ggplot(ibaq.sites.subtoDE, aes(x = log2(sites),
                                                          y = -log10(GelPrepCovFPval),
                                                          color = Ibaqdecile)) +
  geom_point() + 
  #   facet_grid(~.Ibaqquintile) + 
  geom_smooth(method = "lm")

site.variation.decile


#quantile with facet plots. At the higher concentrations there is a connection. This effect may also be driven by length.
site.variation.quintile <- ggplot(ibaq.sites.subtoDE, aes(x = log2(sites),
                                                          y = -log10(GelPrepCovFPval))) +
  geom_point() + 
  facet_grid(.~Ibaqquintile) + 
  geom_smooth(method = "lm")

site.variation.quintile

site.variation.decile <- ggplot(ibaq.sites.subtoDE, aes(x = log2(sites),
                                                          y = -log10(GelPrepCovFPval))) +
  geom_point() + 
  facet_grid(.~Ibaqdecile) + 
  geom_smooth(method = "lm")

site.variation.decile






#Normalize the number of sites by protein length. I hypothesize that when the number of sites is normalized by protein length there is relationship
# between number of sites per protein and phosphopeptide variability. This relationship should be positive and should become stronger with increasing concentration.

ibaq.sites.subtoDE$site.normalized <- ibaq.sites.subtoDE$sites/ibaq.sites.subtoDE$ppSequence.length
ibaq.sites.subtoDE$site.normalized2 <- ibaq.sites.subtoDE$ppSequence.length / ibaq.sites.subtoDE$sites

site.variation.normalized <- ggplot(ibaq.sites.subtoDE, aes(x = log2(site.normalized),
                                                            y = -log10(GelPrepCovFPval))) + 
  geom_point() + 
  geom_smooth(method = "lm")
site.variation.normalized

#The results confirm the hypothesis and the trend is significant.

x = log2(ibaq.sites.subtoDE$site.normalized)
y = -log10(ibaq.sites.subtoDE$GelPrepCovFPval)

cor.test(x,y, alternative = "two", method = "pearson")$p.value
[1] 0.0004833163
cor(x,y, use = "complete.obs")
[1] 0.06117257










#Evidence of an interaction with concentration? I hypothesize that concentration has an additive affect on the relationship.
require(MASS)#for rlm function

site.variation.normalized.concentration <- ggplot(ibaq.sites.subtoDE, aes(x = log2(site.normalized),
                                                            y = -log10(GelPrepCovFPval))) + 
  geom_point() + 
  facet_grid(.~Ibaqquintile) +
  geom_smooth(method = "rlm")
site.variation.normalized.concentration


site.variation.normalized.concentration.decile <- ggplot(ibaq.sites.subtoDE, aes(x = log2(site.normalized),
                                                                          y = -log10(GelPrepCovFPval))) + 
  geom_point() + 
  facet_grid(.~Ibaqdecile) +
  geom_smooth(method = "rlm")
site.variation.normalized.concentration.decile



#these results imply that the number of sites/unit length is a more relevant factor for the lowest expressed protien complement, even after correcting for protein length. I was expecting something of an additive effect. HERE THERE SEEMS TO BE A NEGATIVE RELATIONSHIP BETWEEN THE NUMBER OF SITES/UNIT LENGTH AND VARIABILITY WHEN CUT BY EXPRESSION LEVEL.



# what about a lm with Fstat pvalue as resposnse variable with sites + protein.length + expression.level as terms. These are all correlated with each other. What about sites/length and expression level?
site.expression <- lm(-log10(GelPrepCovFPval) ~ log2(sites) + log10(ibaq.median), data = ibaq.sites.subtoDE)
sites <- lm(-log10(GelPrepCovFPval) ~ log2(sites), data = ibaq.sites.subtoDE)
expression <- lm(-log10(GelPrepCovFPval) ~ log10(ibaq.median), data = ibaq.sites.subtoDE)
site.expression.interaction <- lm(-log10(GelPrepCovFPval) ~ log2(sites) * log10(ibaq.median), data = ibaq.sites.subtoDE)
site.expression.length.interactions <- lm(-log10(GelPrepCovFPval) ~ log2(sites) * log10(ibaq.median) * log2(ppSequence.length)
                                         , data = ibaq.sites.subtoDE)



#summaries
summary(site.expression.interaction)
summary(site.expression.length.interactions)



#is it better to take the length normalized number of sites to make the term independent of expression level?... I think yes
normalizedSites.expression <- lm(-log10(GelPrepCovFPval) ~ log2(site.normalized) + log10(ibaq.median), data = ibaq.sites.subtoDE)
normalizedSites.expression.interaction <- lm(-log10(GelPrepCovFPval) ~ log2(site.normalized) * log10(ibaq.median), data = ibaq.sites.subtoDE)
normalizedSites.expression.interaction.robust <- rlm(-log10(GelPrepCovFPval) ~ log2(site.normalized) * log10(ibaq.median), data = ibaq.sites.subtoDE)

#summaries
summary(normalizedSites.expression)
summary(normalizedSites.expression.interaction)#interaction term is not significant


#robust summary
robust.info <- data.frame(summary(normalizedSites.expression.interaction.robust)$coefficients)
robust.info$pval <- 2*pt(abs(robust.info$t.value), summary(normalizedSites.expression.interaction.robust)$df[2], lower.tail=FALSE)
robust.info



# sites per unit length vs. variability. Relationship is still significant after normalization. I will produce this plot and note that an interaction term is not significant in the text.

y <- -log10(ibaq.sites.subtoDE$GelPrepCovFPval)
x <- log2(ibaq.sites.subtoDE$site.normalized)
cor.test(x,y, alternative = "two", method = "pearson")$p.value
# [1] 0.0004833163
cor(x,y, use = "complete.obs")
# [1] -0.06117257


pdf("normalized_sites_pval_density.pdf", 7, 5)
smoothScatter(x,y, nbin = 150, bandwidth = 0.1,
              cex = .3,
              pch = 19, nrpoints = .15*length(ibaq$GelPrepCovFPval),
              colramp = colorRampPalette(c("white", "light gray", "dark gray", "red")),
              xlab = expression(log[2](sites~per~amino~acid)),
              ylab = expression(-log[10](P~value)), lwd = 10,
              family = "serif"
)
reg.line <- lm(y~x, na.action = "na.omit")
abline(reg.line, lwd = 2, lty = 2)
text(-11.5, 10, expression(R == -0.06), col = "darkred", cex = 1, family = "serif") # rsquared and pvalue
text(-11.5, 9.25, expression(p == 4.83e-04), col = "darkred", cex = 1, family = "serif")
dev.off()












##reactome and GO enrichments ----
source("Enrichment.R")

enrichment_tables <- Enrichment(multExpanded1_withDE_annotated)
dir.create("./Enrichments") 
for(i in 1:length(enrichment_tables)){
  write.table(enrichment_tables[[i]], file.path(getwd(), paste0("Enrichments/", names(enrichment_tables)[i], ".txt")),
              sep = "\t", row.names = F)
}

## output for cytoscape ----


#alter output of each enrichment table for cytoscape
for(i in 1:length(enrichment_tables)){
  enrich.table <- enrichment_tables[[i]]
  enrich.table <- enrich.table[, c(1, 6, 3, 4)]
  names(enrich.table) <- c("ID", "DESCRIPTION", "P-VALUE", "FDR")
  write.table(enrich.table, file.path(getwd(), paste0("Enrichments/modified/", names(enrichment_tables)[i], "modified", ".txt")),
              sep = "\t", row.names = F, qmethod = "d")
}



# Require that I create my own 'gmt' file..

require(biomaRt)
ensembl_75 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", path="/biomart/martservice", 
                     dataset="hsapiens_gene_ensembl")

##create the required file:
# goid tab go discription tab gen1 tab gene2 tab ....

ensembl_hgncid_go <- getBM(attributes = c("go_id", "hgnc_symbol"), filters = 'with_go_id', values = T, mart = ensembl_75)


#For each goid, collect all the matching hgncid. originall tried cast but this took to long and there are differing numbers of genes per go_id, requiring a list output
ensembl_hgncid_go$go_id <- as.factor(ensembl_hgncid_go$go_id)

require(foreach)
require(iterators)
require(GO.db)

goid.list <- foreach(go.id = levels(ensembl_hgncid_go$go_id)) %do% {
  #the biomart db and the go.db may not match
  pos_err <- tryCatch(select(GO.db, keys = as.character(go.id), keytype = "GOID", columns = c("TERM")), error = function(e) e)
  if(!inherits(pos_err, "error")){
    info <- select(GO.db, keys = as.character(go.id), keytype = "GOID", columns = c("TERM"))
    genes <- c(ensembl_hgncid_go[ensembl_hgncid_go$go_id == go.id, "hgnc_symbol"])
    genes <- genes[sapply(genes, function(x) x != "")]
    c(as.character(info), genes)
  }
}

#remove null elements from the list
goid.list <- goid.list[lapply(goid.list, length) > 0]


#write out in 'gmt' format
lapply(goid.list, write, "GO.gmt", append = TRUE, ncolumns = 500, sep = "\t")





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

GelPrep.data <- multExpanded1_withDE_annotated[, c("GelPrepCovSubtoDE", "GelPrepCovglobalFsig", "GelPrepCovFAdjPval",
                                                   "GelPrepPFamIDPhosphoST", "GelPrepPFamIDs", "GelPrepPFamIDPhospho", "GelPrepCovFPval",
                                                      "GelPrepInteractCount", "GelPrepPercentDisorder", "total.mod.count.GelPrep",
                                                   "ppSequence.length", "GelPrep.Disorder.MaxLength", "ppMajorityProteinIDs"
                                                   )]
#note the factors. Revert. remember a dataframe is a list of vectors.
str(GelPrep.data)
i <- sapply(GelPrep.data, is.factor)
GelPrep.data[i] <- lapply(GelPrep.data[i], as.character)


#quick test of interactivity proportional to %disorder. HIGHLY SIGNIFICANT 
y <- as.numeric(GelPrep.data$GelPrepInteractCount)
x <- as.numeric(GelPrep.data$GelPrepPercentDisorder)
plot(x,y)
cor(x,y, use = "complete.obs")
[1] 0.1846606
cor.test(x,y, alternative = "two", method = "pearson")$p.value #spearman also highly sig
[1] 0


#length vs % disorder. pearson sig but spearman not. VERY WEAK
GelPrep.data.uniqueProteins <- GelPrep.data[!duplicated(GelPrep.data$ppMajorityProteinIDs), names(GelPrep.data) %in% c("ppSequence.length", "GelPrepPercentDisorder", "GelPrep.Disorder.MaxLength")]

y <- as.numeric(GelPrep.data.uniqueProteins$ppSequence.length)
x <- as.numeric(GelPrep.data.uniqueProteins$GelPrepPercentDisorder)
plot(x,y)
cor(x,y, use = "complete.obs", method = "spearman")
[1] 0.02124543
cor.test(x,y, alternative = "two", method = "spearman", exact = F)$p.value 
[1] 0.1027654


#length vs longest disorder stretch. strong positive association
y <- as.numeric(GelPrep.data$ppSequence.length)
x <- as.numeric(GelPrep.data$GelPrep.Disorder.MaxLength)
plot(x,y)
cor(x,y, use = "complete.obs", method = "spearman")
[1] 0.4925966
cor.test(x,y, alternative = "two", method = "spearman", exact = F)$p.value 
[1] 1.015572e-198



#length vs interactivity. Sig positive association
y <- as.numeric(GelPrep.data$ppSequence.length)
x <- as.numeric(GelPrep.data$GelPrepInteractCount)
plot(x,y)
cor(x,y, use = "complete.obs", method = "spearman")
[1] 0.06451147
cor.test(x,y, alternative = "two", method = "spearman", exact = F)$p.value 
[1] 1.480578e-11

#length vs ptmcount. Obviously very large positive association
y <- as.numeric(GelPrep.data$total.mod.count.GelPrep)
x <- as.numeric(GelPrep.data$ppSequence.length)
plot(x,y)
cor(x,y, use = "complete.obs", method = "spearman")
[1] 0.6640799
cor.test(x,y, alternative = "two", method = "spearman", exact = F)$p.value 
[1] 0

#ptmcount vs interactivity. Large positive association.
y <- as.numeric(GelPrep.data$total.mod.count.GelPrep)
x <- as.numeric(GelPrep.data$GelPrepInteractCount)
plot(x,y)
cor(x,y, use = "complete.obs", method = "spearman")
[1] 0.384984
cor.test(x,y, alternative = "two", method = "spearman", exact = F)$p.value 
[1] 0






#subset to those phosphopeptides subjected to diffphos analysis (n = 3257)
GelPrep.data <- GelPrep.data[GelPrep.data$GelPrepCovSubtoDE == "+",]

sig.index <- GelPrep.data$GelPrepCovglobalFsig == "+"

#function for mean and se calculations
mean.se.calc <- function(pos.set, neg.set){
  m1 <- mean(pos.set[is.finite(pos.set)], na.rm = T)
  sem1 <- sd(pos.set[is.finite(pos.set)])/ sqrt(sum(is.finite(pos.set)))
  m2 <- mean(neg.set[is.finite(neg.set)], na.rm = T)
  sem2 <- sd(neg.set[is.finite(neg.set)])/ sqrt(sum(is.finite(neg.set)))
  list(means = c(m1, m2), stderr = c(sem1, sem2))
}

##connectivity with biogrid

#these data are clearly not normally distributed (shapiro.test). Test for differences
# For multiple groups Note there is the kruskal.test and the anderson darling test with package "KSamples"


interactions.pos <- as.numeric(GelPrep.data[sig.index, "GelPrepInteractCount"])
interactions.neg <- as.numeric(GelPrep.data[!sig.index, "GelPrepInteractCount"])

#this difference is significant. proteins with less interactions (according to biogrid) are slightly de-enriched.
wilcox.test(interactions.pos, interactions.neg, alternative="two.sided")$p.value#performs mann-whitney test
# ks.test(interactions.pos, interactions.neg, alternative="two.sided")
[1] 8.337766e-06


data <- unlist(mean.se.calc(interactions.pos, interactions.neg))
pdf("GelPrepInteractions.pdf", 5, 7)
barplot2(data[1:2], plot.ci=T,
         ci.l=c(data[1] - data[3], data[2] - data[3]), ci.u=c(data[1] + data[3], data[2] + data[3]),
         lwd=3, ci.lwd=3, col = c("Red", "Blue"), axisnames = F,
         main = "Number of protein-protein interactions")
dev.off()

## threshold independent spearman correlation test of interaction count and F-test p-values
y <- as.numeric(GelPrep.data$GelPrepInteractCount)
x <- -log10(as.numeric(GelPrep.data$GelPrepCovFPval))
plot(x,y)
cor(x,y, use = "complete.obs", method = "spearman")
[1] -0.1009617
cor.test(x,y, alternative = "two", method = "spearman", exact = F)$p.value 
[1] 7.968663e-09


##Disorder

#take the mean disorder percentage from all majority protein ids
GelPrep.data$FinalPercentDisorder <- sapply(GelPrep.data$GelPrepPercentDisorder, function(x){
  mean(as.numeric(unlist(strsplit(x, ";"))))
})

disorder.pos <- as.numeric(GelPrep.data[sig.index, "FinalPercentDisorder"])
disorder.neg <- as.numeric(GelPrep.data[!sig.index, "FinalPercentDisorder"])
wilcox.test(disorder.pos, disorder.neg, alternative="two.sided")$p.value#performs mann-whitney test
[1] 0.09880896


data <- unlist(mean.se.calc(disorder.pos,disorder.neg))
pdf("GelPrepDisorder.pdf", 5, 7)
barplot2(data[1:2], plot.ci=T,
         ci.l=c(data[1] - data[3], data[2] - data[3]), ci.u=c(data[1] + data[3], data[2] + data[3]),
         lwd=3, ci.lwd=3, col = c("Red", "Blue"), axisnames = F,
         main = "Percent disorder")
dev.off()


## threshold independent spearman correlation test of percent disordered residues and F-test p-values
y <- as.numeric(GelPrep.data$FinalPercentDisorder)
x <- -log10(as.numeric(GelPrep.data$GelPrepCovFPval))
plot(x,y)
cor(x,y, use = "complete.obs", method = "spearman")
[1] -0.03758456
cor.test(x,y, alternative = "two", method = "spearman", exact = F)$p.value 
[1] 0.03212017




##Disorder max run
disorder.run.pos <- as.numeric(GelPrep.data[sig.index, "GelPrep.Disorder.MaxLength"])
disorder.run.neg <- as.numeric(GelPrep.data[!sig.index, "GelPrep.Disorder.MaxLength"])
wilcox.test(disorder.run.pos, disorder.run.neg, alternative="two.sided")$p.value#performs mann-whitney test
[1] 1.250861e-06


data <- unlist(mean.se.calc(disorder.run.pos,disorder.run.neg))
pdf("GelPrepDisorderRun.pdf", 5, 7)
barplot2(data[1:2], plot.ci=T,
         ci.l=c(data[1] - data[3], data[2] - data[3]), ci.u=c(data[1] + data[3], data[2] + data[3]),
         lwd=3, ci.lwd=3, col = c("Red", "Blue"), axisnames = F,
         main = "Disorder Run Length")
dev.off()


## threshold independent spearman correlation test of disordered max-run-length and F-test p-values
y <- as.numeric(GelPrep.data$GelPrep.Disorder.MaxLength)
x <- -log10(as.numeric(GelPrep.data$GelPrepCovFPval))
plot(x,y)
cor(x,y, use = "complete.obs", method = "spearman")
[1] -0.08446736
cor.test(x,y, alternative = "two", method = "spearman", exact = F)$p.value 
[1] 1.383675e-06


##modification count
modcount.pos <- as.numeric(GelPrep.data[sig.index, "total.mod.count.GelPrep"])
modcount.neg <- as.numeric(GelPrep.data[!sig.index, "total.mod.count.GelPrep"])
wilcox.test(modcount.pos[is.finite(modcount.pos)], modcount.neg[is.finite(modcount.neg)], alternative="two.sided")$p.value#performs mann-whitney test
[1] 3.187623e-05


data <- unlist(mean.se.calc(modcount.pos, modcount.neg))
pdf("GelPrepModcount.pdf", 5, 7)
barplot2(data[1:2], plot.ci=T,
         ci.l=c(data[1] - data[3], data[2] - data[3]), ci.u=c(data[1] + data[3], data[2] + data[3]),
         lwd=3, ci.lwd=3, col = c("Red", "Blue"), axisnames = F,
         main = "Total number of PTMs")
dev.off()

## threshold independent spearman correlation test of PTM count and F-test p-values
y <- as.numeric(GelPrep.data$total.mod.count.GelPrep)
x <- -log10(as.numeric(GelPrep.data$GelPrepCovFPval))
plot(x,y)
cor(x,y, use = "complete.obs", method = "spearman")
[1] -0.07493412
cor.test(x,y, alternative = "two", method = "spearman", exact = F)$p.value 
[1] 1.861572e-05







#all non significant changes

##FE test for pfam domains and phospho domains (separately for now)

#1st enrichment of proteins containing pfam domains
GelPrep.data$pfam.pos <- sapply(GelPrep.data$GelPrepPFamIDs, function(x){
  any(unlist(strsplit(x, ";")) != "NA")
})

#row 1 diffphos/withpfam AND diffphos/withoutpfam
#row 2 !diffphos/withpfam AND !diffphos/withoutpfam

row1 <- c(nrow(GelPrep.data[GelPrep.data$GelPrepCovglobalFsig == "+" & GelPrep.data$pfam.pos == T, ]),
          nrow(GelPrep.data[GelPrep.data$GelPrepCovglobalFsig == "+" & GelPrep.data$pfam.pos == F, ])
)
row2 <- c(nrow(GelPrep.data[GelPrep.data$GelPrepCovglobalFsig == "-" & GelPrep.data$pfam.pos == T, ]),
          nrow(GelPrep.data[GelPrep.data$GelPrepCovglobalFsig == "-" & GelPrep.data$pfam.pos == F, ])
)

#FEtest. A significant *depletion* of domain containing proteins in the diffphos list
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix)
result$p.value
# [1] 0.09047412
# p-value = 0.04572 (when alternative is "l")



# Threshold independent rho test

# Threshold independent test of association using spearman rank cor coef. Here using nominal ps in the event I want to produce a qq plot
domain.matrix <- GelPrep.data[ , c("GelPrepPFamIDs", "GelPrepCovFPval")]

#switch to 0/1 designation. for now the NAs are a bug
domain.matrix$GelPrepPFamIDs <- ifelse(domain.matrix$GelPrepPFamIDs != "NA", 1, 0)

domain.matrix[] <- lapply(domain.matrix, as.numeric)

plot(domain.matrix[[1]], -log10(domain.matrix[[2]]))
plot(-log10(domain.matrix[[2]]), domain.matrix[[1]])

# Working with negative transform where a positive association indicates enrichment. Negative depletion. Here enrichment!
cor(domain.matrix[[1]], -log10(domain.matrix[[2]]), method = "spearman")
cor(-log10(domain.matrix[[2]]), domain.matrix[[1]], method = "spearman")
[1] -0.03411577

#the correlation IS significant.
cor.test(domain.matrix[[1]], domain.matrix[[2]], method = "spearman", exact = F)$p.value
[1] 0.05155744






#second is enrichment for proteins containing phospho relevant domains
row1 <- c(nrow(GelPrep.data[GelPrep.data$GelPrepCovglobalFsig == "+" & GelPrep.data$GelPrepPFamIDPhospho == "yes", ]),
          nrow(GelPrep.data[GelPrep.data$GelPrepCovglobalFsig == "+" & GelPrep.data$GelPrepPFamIDPhospho == "no", ])
)
row2 <- c(nrow(GelPrep.data[GelPrep.data$GelPrepCovglobalFsig == "-" & GelPrep.data$GelPrepPFamIDPhospho == "yes", ]),
          nrow(GelPrep.data[GelPrep.data$GelPrepCovglobalFsig == "-" & GelPrep.data$GelPrepPFamIDPhospho == "no", ])
)

#FEtest.
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix, alternative = "g")
result$p.value
# row1 
# 0.2740768 



# Threshold independent rho test
# ST specific
domain.matrix <- GelPrep.data[ , c("GelPrepPFamIDPhosphoST", "GelPrepCovFPval")]

#switch to 0/1 designation. for now the NAs are a bug
domain.matrix$GelPrepPFamIDPhosphoST <- ifelse(domain.matrix$GelPrepPFamIDPhosphoST == "yes", 1, 0)

domain.matrix[] <- lapply(domain.matrix, as.numeric)

plot(domain.matrix[[1]], -log10(domain.matrix[[2]]))
plot(-log10(domain.matrix[[2]]), domain.matrix[[1]])

# Working with negative transform where a positive association indicates enrichment. Negative depletion. Here enrichment!
cor(domain.matrix[[1]], -log10(domain.matrix[[2]]), method = "spearman")
cor(-log10(domain.matrix[[2]]), domain.matrix[[1]], method = "spearman")
[1] 0.009986746

#the correlation is not significant.
cor.test(domain.matrix[[1]], domain.matrix[[2]], method = "spearman", exact = F)$p.value
[1] 0.5688538



# Threshold independent rho test
# Threshold independent test of association using spearman rank cor coef. Here using nominal ps in the event I want to produce a qq plot
domain.matrix <- GelPrep.data[ , c("GelPrepPFamIDPhospho", "GelPrepCovFPval")]

#switch to 0/1 designation. for now the NAs are a bug
domain.matrix$GelPrepPFamIDPhospho <- ifelse(domain.matrix$GelPrepPFamIDPhospho == "yes", 1, 0)

domain.matrix[] <- lapply(domain.matrix, as.numeric)

plot(domain.matrix[[1]], -log10(domain.matrix[[2]]))
plot(-log10(domain.matrix[[2]]), domain.matrix[[1]])

# Working with negative transform where a positive association indicates enrichment. Negative depletion. Here enrichment!
cor(domain.matrix[[1]], -log10(domain.matrix[[2]]), method = "spearman")
cor(-log10(domain.matrix[[2]]), domain.matrix[[1]], method = "spearman")
[1] 0.007596505

#the correlation is not significant.
cor.test(domain.matrix[[1]], domain.matrix[[2]], method = "spearman", exact = F)$p.value
[1] 0.6647436




# Threshold independent rho test on proteins with a domain as background

# Threshold independent test of association using spearman rank cor coef. Here using nominal ps in the event I want to produce a qq plot
domain.matrix <- GelPrep.data[GelPrep.data$GelPrepPFamIDs != "NA", c("GelPrepPFamIDPhospho", "GelPrepCovFPval")]

#switch to 0/1 designation. for now the NAs are a bug
domain.matrix$GelPrepPFamIDPhospho <- ifelse(domain.matrix$GelPrepPFamIDPhospho == "yes", 1, 0)

domain.matrix[] <- lapply(domain.matrix, as.numeric)

plot(domain.matrix[[1]], -log10(domain.matrix[[2]]))
plot(-log10(domain.matrix[[2]]), domain.matrix[[1]])

# Working with negative transform where a positive association indicates enrichment. Negative depletion. Here enrichment!
cor(domain.matrix[[1]], -log10(domain.matrix[[2]]), method = "spearman")
cor(-log10(domain.matrix[[2]]), domain.matrix[[1]], method = "spearman")
[1] 0.01275158

#the correlation IS significant.
cor.test(domain.matrix[[1]], domain.matrix[[2]], method = "spearman", exact = F)$p.value
[1] 0.4987089



# I hypothesize that proteins with a higher fraction of phosphosites in disordered regions are depleted within differential phosphorylation subset

# per protein count phosphosites in disorder and phospphosites in ordered. fraction assignment per protein. a single fractional assignment is projected to each phosphosite and a linear regression is performed.

# This will require a dplyr like approach

disorder.fraction <- multExpanded1_withDE_annotated[ , c("GelPrepCovFPval", "GelPrepCovSubtoDE", "GelPrep.Pos.Disorder",
                                                  "ppMajorityProteinIDs", "ppProteinIDs", "ppSequence.length", "idmult")] #note the NAs

# I need two vectors, one of fraction of sites that map to disordered regions and another with pvalues. pvalues will be unique but proteins ids will be repeated. The depletion is significant but not very strong.

disorder.fraction <- disorder.fraction %>% group_by(ppProteinIDs) %>% mutate(fraction.disordered = sum(GelPrep.Pos.Disorder / n(), na.rm = T), sites = n())

disorder.fraction <- arrange(disorder.fraction, ppProteinIDs)
disorder.fraction <- filter(disorder.fraction, GelPrepCovSubtoDE == "+") %>% arrange(ppProteinIDs)

y <- -log10(as.numeric(disorder.fraction$GelPrepCovFPval))
x <- disorder.fraction$fraction.disordered
plot(x,y)
R <- cor(x, y, method = "spearman", use = "complete.obs")
R
[1] -0.05551729
cor.test(x,y, method = "spearman", exact = F)$p.value
[1] 0.001526451

R <- cor(x, y, method = "pearson", use = "complete.obs")
R
[1] -0.05198285
cor.test(x,y, method = "pearson", exact = F)$p.value
[1] 0.003001972


#protein length





# A threshold independent test of a correlation between percent disorder and pvalues using spearman is barely sig while pearson is not.

# Threshold independent test of association using spearman rank cor coef. Here using nominal ps in the event I want to produce a qq plot
disorder.matrix <- GelPrep.data[, c("FinalPercentDisorder", "GelPrepCovFPval")]

# remove reverted identifications
disorder.matrix <- disorder.matrix[complete.cases(disorder.matrix),]

# y <- -log10(as.numeric(disorder.matrix$GelPrepCovFPval))
# x <- log10(disorder.matrix$FinalPercentDisorder) 

y <- -log10(as.numeric(disorder.matrix$GelPrepCovFPval))
x <- disorder.matrix$FinalPercentDisorder 

plot(x,y)
R <- cor(x, y, method = "spearman", use = "complete.obs")
R
[1] -0.03758456
cor.test(x,y, method = "spearman", exact = F)$p.value
[1] 0.03212017



#Hypothesis is that proteins with more motifs are depleted if information nodes are buffered from differential phosphorylation. Motif counts is very sparse.


motifCounts <- multExpanded1_withDE_annotated[multExpanded1_withDE_annotated$GelPrepCovSubtoDE == "+", 
                                                 c("GelPrepCovFPval", "GelPrepMotifCount"
                                                 )]
y <- -log10(as.numeric(motifCounts$GelPrepCovFPval))
x <- log10(motifCounts$GelPrepMotifCount + 1) 
plot(x,y)

R <- cor(x, y, method = "spearman", use = "complete.obs")
R
[1] 0.09044962
cor.test(x,y, method = "spearman", exact = F)$p.value
[1] 2.333546e-07

#how many have >0?. This is very small sample space.
table(motifCounts$GelPrepMotifCount)





## proteome wide correlation between protein length and disordered residue percentage












