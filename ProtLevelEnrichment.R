#diffphos ProtLevel enrichment analysis - (all enrichment analysis at the protein group level)

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


#protein expression level is negatively correlated with the number of sites identified and quantified
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


##barplot for connectivity, % Disorder, extent of PTM modifications ----
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



# Threshold independent rho test of presence of pfam id and phosphopeptide variation

# Threshold independent test of association using spearman rank cor coef. Here using nominal ps in the event I want to produce a qq plot
domain.matrix <- GelPrep.data[ , c("GelPrepPFamIDs", "GelPrepCovFPval")]

#switch to 0/1 designation. for now the NAs are a bug
domain.matrix$GelPrepPFamIDs <- ifelse(domain.matrix$GelPrepPFamIDs != "NA", 1, 0)

domain.matrix[] <- lapply(domain.matrix, as.numeric)

plot(domain.matrix[[1]], -log10(domain.matrix[[2]]))
plot(-log10(domain.matrix[[2]]), domain.matrix[[1]])

# Working with negative transform where a positive association indicates enrichment. Negative depletion.
cor(domain.matrix[[1]], -log10(domain.matrix[[2]]), method = "spearman")
cor(-log10(domain.matrix[[2]]), domain.matrix[[1]], method = "spearman")
[1] -0.03411577

#the correlation is not significant.
cor.test(domain.matrix[[1]], domain.matrix[[2]], method = "spearman", exact = F)$p.value
[1] 0.05155744


# Threshold independent rho test of phosphorylation regulation association and phosphopeptide variation (ST specific)
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



# Threshold independent rho test phosphorylation regulation association and phosphopeptide variation
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

cor.test(domain.matrix[[1]], domain.matrix[[2]], method = "spearman", exact = F)$p.value
[1] 0.4987089
