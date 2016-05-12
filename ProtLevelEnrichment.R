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
R <- cor(x, y, method = "pearson", use = "complete.obs")
R
cor.test(x,y)$p.value

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
cor(x,y, use = "complete.obs")
# [1] 0.05184077
cor.test(x,y, alternative = "two", method = "pearson")$p.value
# [1] 0.003082246

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
cor.test(x,y, alternative = "two", method = "pearson")$p.value
# [1] 6.969347e-07
cor(x,y, use = "complete.obs")
# [1] -0.08689839


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
                                                      "GelPrepInteractCount", "GelPrepPercentDisorder", "total.mod.count.GelPrep"
                                                   )]
#note the factors. Revert. remember a dataframe is a list of vectors.
str(GelPrep.data)
i <- sapply(GelPrep.data, is.factor)
GelPrep.data[i] <- lapply(GelPrep.data[i], as.character)

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



##Disorder
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












