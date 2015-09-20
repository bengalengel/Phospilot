# This funcion applies diffphos nrichment tests at the local (peptide) level. The tests are that diff phosphorylated sites are:
#   1) enriched within disordered regions
#   2) enriched within protein domains
#   3) enriched within phospho relevant domains
#   4) enriched with specific kinase motifs (HPRD binary assignment)
# 
# Future: Enriched within motifs, enriched within motifs relevant to phosphorylation
# 
# 

# Load Data ----
GelPrep.data <- multExpanded1_withDE_annotated[, c("GelPrepNormSubtoDE", "GelPrepNormglobalFsig", "GelPrepNormFAdjPval", "GelPrepNormFPval",
                                                   "GelPrep.Pos.Disorder", "site.in.domain.GelPrep", "phospho.relevant.GelPrep")]


#note their may be factors. Revert. remember a dataframe is a list of vectors.
str(GelPrep.data)
i <- sapply(GelPrep.data, is.factor)
GelPrep.data[i] <- lapply(GelPrep.data[i], as.character)

#subset to those phosphopeptides subjected to diffphos analysis (n = 3257)
GelPrep.data <- GelPrep.data[GelPrep.data$GelPrepNormSubtoDE == "+",]


# Diffphos Enrichement in Disordered Regions ----
row1 <- c(nrow(GelPrep.data[GelPrep.data$GelPrepNormglobalFsig == "+" & GelPrep.data$GelPrep.Pos.Disorder == T, ]),
          nrow(GelPrep.data[GelPrep.data$GelPrepNormglobalFsig == "+" & GelPrep.data$GelPrep.Pos.Disorder == F, ])
)
row2 <- c(nrow(GelPrep.data[GelPrep.data$GelPrepNormglobalFsig == "-" & GelPrep.data$GelPrep.Pos.Disorder == T, ]),
          nrow(GelPrep.data[GelPrep.data$GelPrepNormglobalFsig == "-" & GelPrep.data$GelPrep.Pos.Disorder == F, ])
)


#FEtest. No enrichment
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix)
result$p.value
0.3882814

# Threshold independent rho test

# Threshold independent test of association using spearman rank cor coef. Here using nominal ps in the event I want to produce a qq plot
disorder.matrix <- GelPrep.data[, c("GelPrep.Pos.Disorder", "GelPrepNormFPval")]

#switch to 0/1 designation
disorder.matrix$GelPrep.Pos.Disorder <- ifelse(disorder.matrix$GelPrep.Pos.Disorder == TRUE, 1, 0)
disorder.matrix[] <- lapply(disorder.matrix, as.numeric)

plot(disorder.matrix[[1]], -log10(disorder.matrix[[2]]))
plot(-log10(disorder.matrix[[2]]), disorder.matrix[[1]])

# Working with negative transform where a positive association indicates enrichment. Here rho is negative
disorder.matrix <- na.omit(disorder.matrix)#omit 6 NA values...why am I getting NA values?...
cor(disorder.matrix[[1]], -log10(disorder.matrix[[2]]), method = "spearman")
cor(-log10(disorder.matrix[[2]]), disorder.matrix[[1]], method = "spearman")
-0.01276385 


#the correlation is not significant
cor.test(disorder.matrix[[1]], -log10(disorder.matrix[[2]]), method = "spearman", exact = F)$p.value
[1] 0.4669117




# Diffphos Enrichment for sites within a domain ----

#issues with background. It should be only those sites

row1 <- c(nrow(GelPrep.data[GelPrep.data$GelPrepNormglobalFsig == "+" & GelPrep.data$site.in.domain.GelPrep == T, ]),
          nrow(GelPrep.data[GelPrep.data$GelPrepNormglobalFsig == "+" & GelPrep.data$site.in.domain.GelPrep == F, ])
)
row2 <- c(nrow(GelPrep.data[GelPrep.data$GelPrepNormglobalFsig == "-" & GelPrep.data$site.in.domain.GelPrep == T, ]),
          nrow(GelPrep.data[GelPrep.data$GelPrepNormglobalFsig == "-" & GelPrep.data$site.in.domain.GelPrep == F, ])
)


#FEtest. No enrichment
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix)
result$p.value
0.2770237



# Threshold independent rho test

# Threshold independent test of association using spearman rank cor coef. Here using nominal ps in the event I want to produce a qq plot
domain.matrix <- GelPrep.data[, c("site.in.domain.GelPrep", "GelPrepNormFPval")]

#switch to 0/1 designation
domain.matrix$site.in.domain.GelPrep <- ifelse(domain.matrix$site.in.domain.GelPrep == TRUE, 1, 0)
domain.matrix[] <- lapply(domain.matrix, as.numeric)

plot(domain.matrix[[1]], -log10(domain.matrix[[2]]))
plot(-log10(domain.matrix[[2]]), domain.matrix[[1]])

# Working with negative transform where a positive association indicates enrichment. Here rho is negative
domain.matrix <- na.omit(domain.matrix)#omit NA values
cor(domain.matrix[[1]], -log10(domain.matrix[[2]]), method = "spearman")
cor(-log10(domain.matrix[[2]]), domain.matrix[[1]], method = "spearman")
[1] 0.008545924

#the correlation is not significant
cor.test(domain.matrix[[1]], -log10(domain.matrix[[2]]), method = "spearman", exact = F)$p.value
[1] 0.6371396


# Diffphos Enrichment for sites within a phospho relevant domain ----

row1 <- c(nrow(GelPrep.data[GelPrep.data$GelPrepNormglobalFsig == "+" & GelPrep.data$phospho.relevant.GelPrep == T, ]),
          nrow(GelPrep.data[GelPrep.data$GelPrepNormglobalFsig == "+" & GelPrep.data$phospho.relevant.GelPrep == F, ])
)
row2 <- c(nrow(GelPrep.data[GelPrep.data$GelPrepNormglobalFsig == "-" & GelPrep.data$phospho.relevant.GelPrep == T, ]),
          nrow(GelPrep.data[GelPrep.data$GelPrepNormglobalFsig == "-" & GelPrep.data$phospho.relevant.GelPrep == F, ])
)


#FEtest. No enrichment
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix)
result$p.value
[1] 0.64664



