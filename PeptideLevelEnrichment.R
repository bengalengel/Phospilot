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
GelPrep.data <- multExpanded1_withDE_annotated[, c("GelPrepNormSubtoDE", "GelPrepNormglobalFsig", "GelPrepNormFAdjPval",
                                                   "GelPrep.Pos.Disorder", "site.in.domain.GelPrep", "phospho.relevant.GelPrep")]


#note their may be factors. Revert. remember a dataframe is a list of vectors.
str(GelPrep.data)
i <- sapply(GelPrep.data, is.factor)
GelPrep.data[i] <- lapply(GelPrep.data[i], as.character)

#subset to those phosphopeptides subjected to diffphos analysis (n = 3257)
GelPrep.data <- GelPrep.data[GelPrep.data$GelPrepNormSubtoDE == "+",]

sig.index <- GelPrep.data$GelPrepNormglobalFsig == "+"

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



