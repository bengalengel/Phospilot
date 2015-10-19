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
GelPrep.data <- multExpanded1_withDE_annotated[, c("GelPrepCovSubtoDE", "GelPrepCovglobalFsig", "GelPrepCovFAdjPval", "GelPrepCovFPval",
                                                   "GelPrep.Pos.Disorder", "site.in.domain.GelPrep", "phospho.relevant.GelPrep")]


#note their may be factors. Revert. remember a dataframe is a list of vectors.
str(GelPrep.data)
i <- sapply(GelPrep.data, is.factor)
GelPrep.data[i] <- lapply(GelPrep.data[i], as.character)

#subset to those phosphopeptides subjected to diffphos analysis (n = 3257)
GelPrep.data <- GelPrep.data[GelPrep.data$GelPrepCovSubtoDE == "+",]


# Diffphos Enrichement in Disordered Regions ----
row1 <- c(nrow(GelPrep.data[GelPrep.data$GelPrepCovglobalFsig == "+" & GelPrep.data$GelPrep.Pos.Disorder == T, ]),
          nrow(GelPrep.data[GelPrep.data$GelPrepCovglobalFsig == "+" & GelPrep.data$GelPrep.Pos.Disorder == F, ])
)
row2 <- c(nrow(GelPrep.data[GelPrep.data$GelPrepCovglobalFsig == "-" & GelPrep.data$GelPrep.Pos.Disorder == T, ]),
          nrow(GelPrep.data[GelPrep.data$GelPrepCovglobalFsig == "-" & GelPrep.data$GelPrep.Pos.Disorder == F, ])
)


#FEtest. No enrichment
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix)
result$p.value
0.02103861

# Threshold independent rho test

# Threshold independent test of association using spearman rank cor coef. Here using nominal ps in the event I want to produce a qq plot
disorder.matrix <- GelPrep.data[, c("GelPrep.Pos.Disorder", "GelPrepCovFPval")]

#switch to 0/1 designation
disorder.matrix$GelPrep.Pos.Disorder <- ifelse(disorder.matrix$GelPrep.Pos.Disorder == TRUE, 1, 0)
disorder.matrix[] <- lapply(disorder.matrix, as.numeric)

plot(disorder.matrix[[1]], -log10(disorder.matrix[[2]]))
plot(-log10(disorder.matrix[[2]]), disorder.matrix[[1]])

# Working with negative transform where a positive association indicates enrichment. Here rho is negative
disorder.matrix <- na.omit(disorder.matrix)#omit 6 NA values...why am I getting NA values?...
cor(disorder.matrix[[1]], -log10(disorder.matrix[[2]]), method = "spearman")
cor(-log10(disorder.matrix[[2]]), disorder.matrix[[1]], method = "spearman")
-0.04761224 

#the correlation is significant
cor.test(disorder.matrix[[1]], -log10(disorder.matrix[[2]]), method = "spearman", exact = F)$p.value
[1] 0.006623004




# Diffphos Enrichment for sites within a domain ----

#issues with background. It should be only those sites

row1 <- c(nrow(GelPrep.data[GelPrep.data$GelPrepCovglobalFsig == "+" & GelPrep.data$site.in.domain.GelPrep == T, ]),
          nrow(GelPrep.data[GelPrep.data$GelPrepCovglobalFsig == "+" & GelPrep.data$site.in.domain.GelPrep == F, ])
)
row2 <- c(nrow(GelPrep.data[GelPrep.data$GelPrepCovglobalFsig == "-" & GelPrep.data$site.in.domain.GelPrep == T, ]),
          nrow(GelPrep.data[GelPrep.data$GelPrepCovglobalFsig == "-" & GelPrep.data$site.in.domain.GelPrep == F, ])
)


#FEtest. Enrichment
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix)
result$p.value
0.02972399



# Threshold independent rho test

# Threshold independent test of association using spearman rank cor coef. Here using nominal ps in the event I want to produce a qq plot
domain.matrix <- GelPrep.data[, c("site.in.domain.GelPrep", "GelPrepCovFPval")]

#switch to 0/1 designation
domain.matrix$site.in.domain.GelPrep <- ifelse(domain.matrix$site.in.domain.GelPrep == TRUE, 1, 0)
domain.matrix[] <- lapply(domain.matrix, as.numeric)

plot(domain.matrix[[1]], -log10(domain.matrix[[2]]))
plot(-log10(domain.matrix[[2]]), domain.matrix[[1]])

# Working with negative transform where a positive association indicates enrichment. Here rho is negative
domain.matrix <- na.omit(domain.matrix)#omit NA values
cor(domain.matrix[[1]], -log10(domain.matrix[[2]]), method = "spearman")
cor(-log10(domain.matrix[[2]]), domain.matrix[[1]], method = "spearman")
[1] 0.06167596

#the correlation is not significant
cor.test(domain.matrix[[1]], -log10(domain.matrix[[2]]), method = "spearman", exact = F)$p.value
[1] 0.0006557162



# Diffphos Enrichment for sites within a phospho relevant domain ----

row1 <- c(nrow(GelPrep.data[GelPrep.data$GelPrepCovglobalFsig == "+" & GelPrep.data$phospho.relevant.GelPrep == T, ]),
          nrow(GelPrep.data[GelPrep.data$GelPrepCovglobalFsig == "+" & GelPrep.data$phospho.relevant.GelPrep == F, ])
)
row2 <- c(nrow(GelPrep.data[GelPrep.data$GelPrepCovglobalFsig == "-" & GelPrep.data$phospho.relevant.GelPrep == T, ]),
          nrow(GelPrep.data[GelPrep.data$GelPrepCovglobalFsig == "-" & GelPrep.data$phospho.relevant.GelPrep == F, ])
)


#FEtest. No enrichment
contmatrix <- rbind(row1,row2)
result <- fisher.test(contmatrix)
result$p.value
[1] 0.7033173


# Threshold independent test of association using spearman rank cor coef. Here using nominal ps in the event I want to produce a qq plot
phospho.domain.matrix <- GelPrep.data[, c("phospho.relevant.GelPrep", "GelPrepCovFPval")]

#switch to 0/1 designation
phospho.domain.matrix$phospho.relevant.GelPrep <- ifelse(phospho.domain.matrix$phospho.relevant.GelPrep == TRUE, 1, 0)
phospho.domain.matrix[] <- lapply(phospho.domain.matrix, as.numeric)

plot(phospho.domain.matrix[[1]], -log10(phospho.domain.matrix[[2]]))
plot(-log10(phospho.domain.matrix[[2]]), phospho.domain.matrix[[1]])

# Working with negative transform where a positive association indicates enrichment. Here rho is negative
phospho.domain.matrix <- na.omit(phospho.domain.matrix)#omit NA values
cor(phospho.domain.matrix[[1]], -log10(phospho.domain.matrix[[2]]), method = "spearman")
cor(-log10(phospho.domain.matrix[[2]]), phospho.domain.matrix[[1]], method = "spearman")
[1] -0.008525573

#the correlation is not significant
cor.test(phospho.domain.matrix[[1]], -log10(phospho.domain.matrix[[2]]), method = "spearman", exact = F)$p.value
[1] 0.6386046


# Threshold independent rho test on proteins with a domain as background

# Threshold independent test of association using spearman rank cor coef. Here using nominal ps in the event I want to produce a qq plot
phospho.domain.matrix <- GelPrep.data[GelPrep.data$site.in.domain.GelPrep == TRUE, c("phospho.relevant.GelPrep", "GelPrepCovFPval")]

#switch to 0/1 designation
phospho.domain.matrix$phospho.relevant.GelPrep <- ifelse(phospho.domain.matrix$phospho.relevant.GelPrep == TRUE, 1, 0)
phospho.domain.matrix[] <- lapply(phospho.domain.matrix, as.numeric)

plot(phospho.domain.matrix[[1]], -log10(phospho.domain.matrix[[2]]))
plot(-log10(phospho.domain.matrix[[2]]), phospho.domain.matrix[[1]])

# Working with negative transform where a positive association indicates enrichment. Here rho is negative
phospho.domain.matrix <- na.omit(phospho.domain.matrix)#omit NA values
cor(phospho.domain.matrix[[1]], -log10(phospho.domain.matrix[[2]]), method = "spearman")
cor(-log10(phospho.domain.matrix[[2]]), phospho.domain.matrix[[1]], method = "spearman")
[1] -0.0389816

#the correlation is not significant
cor.test(phospho.domain.matrix[[1]], -log10(phospho.domain.matrix[[2]]), method = "spearman", exact = F)$p.value
[1] 0.2759511













# HPRD motif enrichment analysis ----

#no significant motif enrichments

Enrich <- function(x,y){
  # This function accepts character vectors of a selected subset and background and returns a DF of adjusted pvalues for categorical enrichment using a one sided fisher's exact test.
  # x=background and y=enriched. 
  require(plyr)
  BGtable <- as.matrix(table(x))
  #remove entries with 0
  BGtable <- BGtable[BGtable!=0,,drop=F]
  DEtable <- as.matrix(table(y))
  DEtable <- as.matrix(DEtable[row.names(DEtable) %in% row.names(BGtable),])#removing zeros and all factors not present in BG data. 
  #subset the background table in a similar way to ensure we are making the proper comparisons
  BGtable <- as.matrix(BGtable[row.names(BGtable) %in% row.names(DEtable),])
  NotDE <- BGtable-DEtable
  facttemp <- as.factor(row.names(DEtable))
  pvals <- c()
  frequency <- c()
  ids <- c()
  #make contingency table and perform FE test one sided enrichment for each term
  for(i in levels(facttemp)){
    if(DEtable[as.character(i),] & NotDE[as.character(i),] >= 0){#this condition should always be true
      DErow <- c(DEtable[as.character(i),],(sum(DEtable)-DEtable[as.character(i),]))
      NotDErow <- c(NotDE[as.character(i),],(sum(NotDE)-NotDE[as.character(i),]))
      contmatrix <- rbind(DErow,NotDErow)
      tmp <- fisher.test(contmatrix, alternative = "g")
      pvals <- c(pvals,tmp$p.value)
      ids <- c(ids,i)
      tmp2 <- paste(as.character(contmatrix[1,1]), as.character(colSums(contmatrix)[1]), sep = "/")
      frequency <- c(frequency, tmp2)
    }
  }
  ##multiple testing correction for pvals
  adjps <- p.adjust(pvals,method="BH")#none pass significance in this instance
  #put into a dataframe object, sort by adjusted p values
  DF <- data.frame('ID' = ids, 'Frequency' = frequency, 'pvalue' = pvals, 'adjpvalue' = adjps)
  DF <- arrange(DF, adjpvalue)#arrange
  DF <- DF[DF$adjpvalue <= .05, ]
  if(nrow(DF) < 1){
    DF <- "No Enrichment"
  }
  return(DF)
}

##########################

#function for annotation
SplitNClean <- function(x){
  #SplitNClean takes a character vector of annotations separated by a colon and returns a character vector where each term is an element, and no element of the character vector is empty/NA
  y <- sapply(x, function(x) strsplit(x, ";"))
  y <- unlist(y)
  index <- which(y == "")
  y[index] <- NA
  y <- na.omit(y)
}


#Gelprep normalized analysis
GelPrep.BGmotifs <- SplitNClean(as.character(multExpanded1_withDE_annotated[multExpanded1_withDE_annotated$GelPrepCovSubtoDE == "+",
                                                                                  "motifs"]))
GelPrep.DEmotifs.omnibus <- SplitNClean(as.character(multExpanded1_withDE_annotated[multExpanded1_withDE_annotated$GelPrepCovglobalFsig == "+",
                                                                                  "motifs"]))  
GelPrep.DEmotifs.cont1 <- SplitNClean(as.character(multExpanded1_withDE_annotated[multExpanded1_withDE_annotated$GelPrepCovDEcont1 == "+",
                                                                                  "motifs"]))  
GelPrep.DEmotifs.cont2 <- SplitNClean(as.character(multExpanded1_withDE_annotated[multExpanded1_withDE_annotated$GelPrepCovDEcont2 == "+",
                                                                                  "motifs"]))  
GelPrep.DEmotifs.cont3 <- SplitNClean(as.character(multExpanded1_withDE_annotated[multExpanded1_withDE_annotated$GelPrepCovDEcont3 == "+",
                                                                                  "motifs"]))  
#directional
GelPrep.DEmotifs.cont1up <- SplitNClean(as.character(multExpanded1_withDE_annotated[multExpanded1_withDE_annotated$GelPrepCovcont1up == "+",
                                                                                  "motifs"]))  
GelPrep.DEmotifs.cont1down <- SplitNClean(as.character(multExpanded1_withDE_annotated[multExpanded1_withDE_annotated$GelPrepCovcont1down == "+",
                                                                                    "motifs"]))  

GelPrep.DEmotifs.cont2up <- SplitNClean(as.character(multExpanded1_withDE_annotated[multExpanded1_withDE_annotated$GelPrepCovcont2up == "+",
                                                                                    "motifs"]))  
GelPrep.DEmotifs.cont2down <- SplitNClean(as.character(multExpanded1_withDE_annotated[multExpanded1_withDE_annotated$GelPrepCovcont2down == "+",
                                                                                      "motifs"]))  

GelPrep.DEmotifs.cont3up <- SplitNClean(as.character(multExpanded1_withDE_annotated[multExpanded1_withDE_annotated$GelPrepCovcont3up == "+",
                                                                                    "motifs"]))  
GelPrep.DEmotifs.cont3down <- SplitNClean(as.character(multExpanded1_withDE_annotated[multExpanded1_withDE_annotated$GelPrepCovcont3down == "+",
                                                                                      "motifs"]))  



#GelPrep normalized contrast enrichments. No enrichments at FDR 5%
Enrich.motifs.omnibus <- Enrich(GelPrep.BGmotifs, GelPrep.DEmotifs.omnibus)
Enrich.motifs.1 <- Enrich(GelPrep.BGmotifs, GelPrep.DEmotifs.cont1)
Enrich.motifs.2 <- Enrich(GelPrep.BGmotifs, GelPrep.DEmotifs.cont2)
Enrich.motifs.3 <- Enrich(GelPrep.BGmotifs, GelPrep.DEmotifs.cont3)

# Directional enrichments. Very weak enrichments and nothing is popping in terms of motifs
Enrich.motifs.1up <- Enrich(GelPrep.BGmotifs, GelPrep.DEmotifs.cont1up) # 1 enrichment
Enrich.motifs.1down <- Enrich(GelPrep.BGmotifs, GelPrep.DEmotifs.cont1down)
Enrich.motifs.2up <- Enrich(GelPrep.BGmotifs, GelPrep.DEmotifs.cont2up) # 1 enrichment
Enrich.motifs.2down <- Enrich(GelPrep.BGmotifs, GelPrep.DEmotifs.cont2down)
Enrich.motifs.3up <- Enrich(GelPrep.BGmotifs, GelPrep.DEmotifs.cont3up) # 1 enrichment
Enrich.motifs.3down <- Enrich(GelPrep.BGmotifs, GelPrep.DEmotifs.cont3down)


# motif-x heatmap ----

# why? does the enrichment of certain motifs buttress our function claim? shows that something non-random is happening




