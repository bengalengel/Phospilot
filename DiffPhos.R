DiffPhos <- function(pilot){
  #this function accepts the pilot dataframe and runs DE analysis using limma. Columns are appended
  #to the multexpanded1 file according to the presence of DE or not.
  
  require(limma)
  require(statmod)
  #Produce the design matrix
  
  fac <- factor(c(1,1,2,2,3,3))##codes the grouping for the ttests
  design <- model.matrix(~0 + fac)
  dnames <- levels(as.factor(substr(colnames(pilot), 1, 7))) ##check me out. use 5 digit exp name.
  colnames(design) <- dnames
  
  #limma fit using all common for now.
  # The philosophy of the approach is as follows. You have to start by fitting a linear model to
  # your data which fully models the systematic part of your data. The model is specified by the design
  # matrix. Each row of the design matrix corresponds to an array in your experiment and each column
  # corresponds to a coefficient that is used to describe the RNA sources in your experiment.
  # The main purpose of this step is to estimate the variability in the data, hence the systematic part needs to be modelled so it can be distinguished from random variation.
  
  #Perhaps I can add replication correlation information. As of now sparseness throws an error. Example below:
  
  # 17.3.6 Within-patient correlations
  # The study involves multiple cell types from the same patient. Arrays from the same donor are not
  # independent, so we need to estimate the within-dinor correlation:
  #   > ct <- factor(targets$CellType)
  # > design <- model.matrix(~0+ct)
  # > colnames(design) <- levels(ct)
  # > dupcor <- duplicateCorrelation(y,design,block=targets$Donor)
  # > dupcor$consensus.correlation
  # [1] 0.134
  # As expected, the within-donor correlation is small but positive.
  
  fit <- lmFit(pilot, design)
  
  # In practice the requirement to have exactly as many coefficients as RNA sources is too restrictive
  # in terms of questions you might want to answer. You might be interested in more or fewer comparisons
  # between the RNA source. Hence the contrasts step is provided so that you can take the initial
  # coefficients and compare them in as many ways as you want to answer any questions you might have,
  # regardless of how many or how few these might be.
  
  #Now to make all pairwise comparisons (group2-1, group3-2, group3-1)
  contrast.matrix <- makeContrasts(HL18862-HL18486, HL19160-HL18862, HL19160-HL18486, levels = design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  
  # For statistical analysis and assessing differential expression, limma uses an empirical Bayes method
  # to mo derate the standard errors of the estimated log-fold changes. This results in more stable
  # inference and improved p ower, esp ecially for exp eriments with small numb ers of arrays
  fit2 <- eBayes(fit2)
  
  
  
  #Look at pairwise DE using toptable and the coef parameter to id which genes you are interested in 
  sig1 <- topTable(fit2, coef = 1, adjust = "BH", n=Inf, sort="p", p=.05)#sorts by adjusted p up to the threshold of .05, which is the default FDR chosen for differential expression ("results" function). This actually seems a conservative way to sort.
  sig2 <- topTable(fit2, coef = 2, adjust = "BH", n=Inf, sort="p", p=.05)
  sig3 <- topTable(fit2, coef = 3, adjust = "BH", n=Inf, sort="p", p=.05)
  
  # sig1 - 18862-18486
  # sig2 - 19160-18862
  # sig3 - 19160-18486
  
  c1up  <- sig1[sig1$logFC > 0,]
  c1down <- sig1[sig1$logFC < 0,]
  c2up <- sig2[sig2$logFC > 0,]
  c2down <- sig2[sig2$logFC < 0,]
  c3up <- sig3[sig3$logFC > 0,]
  c3down <- sig3[sig3$logFC < 0,]
  
  
  tt1 <- topTable(fit2, coef = 1, adjust = "BH", n=Inf)#sorts by adjusted p up to the threshold of .
  tt2 <- topTable(fit2, coef = 2, adjust = "BH", n=Inf)#sorts by adjusted p up to the threshold of .
  tt3 <- topTable(fit2, coef = 3, adjust = "BH", n=Inf)#sorts by adjusted p up to the threshold of .
  
  hist(tt1$P.Value, nc=40, xlab="P values", main = colnames(contrast.matrix)[1])
  hist(tt2$P.Value, nc=40, xlab="P values", main = colnames(contrast.matrix)[2])
  hist(tt3$P.Value, nc=40, xlab="P values", main = colnames(contrast.matrix)[3])
  
  plot(tt1$logFC,-log10(tt1$P.Value), xlab = colnames(contrast.matrix)[1], pch = 20, ylab = "-log10(P)",xlim = c(-5, 5))
  #sites with sig difference in comparison 1
  names <- row.names(sig1)
  names2 <- row.names(tt1)
  index <- which(names2 %in% names)
  points(tt1$logFC[index],-log10(tt1$P.Value)[index], col="red3", pch = 20)
  
  
  plot(tt2$logFC,-log10(tt2$P.Value), xlab = colnames(contrast.matrix)[2], pch = 20, ylab = "-log10(P)",xlim = c(-5, 5))
  #sites with sig difference in comparison 1
  names <- row.names(sig2)
  names2 <- row.names(tt2)
  index <- which(names2 %in% names)
  points(tt2$logFC[index],-log10(tt2$P.Value)[index], col="red3", pch = 20)
  
  
  plot(tt3$logFC,-log10(tt3$P.Value), xlab = colnames(contrast.matrix)[3], pch = 20, ylab = "-log10(P)",xlim = c(-5, 5))
  #sites with sig difference in comparison 1
  names <- row.names(sig3)
  names2 <- row.names(tt3)
  index <- which(names2 %in% names)
  points(tt3$logFC[index],-log10(tt3$P.Value)[index], col="red3", pch = 20)
  
  
  
  
  results <- decideTests(fit2, adjust.method = "BH", method = "separate")#results is a 'TestResults' matrix
  #separate compares each sample individually and is the default approach
  summary(results)
  
  
  vennDiagram(results, cex=c(1.2,1,0.7)) #good DE across conditions
  vennDiagram(results, cex=c(1.2,1,0.7), include = "up") #good DE across conditions
  vennDiagram(results, cex=c(1.2,1,0.7), include = "down") #good DE across conditions
  vennDiagram(results, cex=c(1.2,1,0.7), include = c("up", "down")) #good DE across conditions
  
  
  table("18862-18486" =results[,1],"19160-18862"=results[,2])
  
  
  # volcanoplot(fit2, coef=1, main = colnames(contrast.matrix)[1], highlight = 900)
  # abline(v=1))
  # abline(v=-1)
  # volcanoplot(fit2, coef=2, main = colnames(contrast.matrix)[2])
  # abline(v=1)
  # abline(v=-1)
  # volcanoplot(fit2, coef=3, main = colnames(contrast.matrix)[3])
  # abline(v=1)
  # abline(v=-1)
  
  # F statistic distributions and cuts by DE across individuals
  # Fvalues <- as.data.frame(Fvalues)
  # row.names(Fvalues) <- sites
  # boxplot(log10(Fvalues))
  # hist(log10(as.matrix(Fvalues)))
  # plot(density(log10(as.matrix(Fvalues))))
  # plot(density(as.matrix(Fvalues)))
  # density(as.matrix(Fvalues))
  
#   plot(density(fit2$F))
#   plot(density(log10(fit2$F)))
  
  # NOTE THE FOLLOWING WHEN COMPARING ALL THE CONTRASTS AT ONCE PG 62 IN LIMMA USER GUIDE - decideTests method global should be used here.
  
  # method="global" is recommended when a set of closely related contrasts are being tested. This
  # method simply appends all the tests together into one long vector of tests, i.e., it treats all the tests
  # as equivalent regardless of which probe or contrast they relate to. An advantage is that the raw
  # p-value cutoff is consistent across all contrasts. For this reason, method="global" is recommended if
  # you want to compare the number of DE genes found for different contrasts, for example interpreting
  # the number of DE genes as representing the strength of the contrast. However users need to be aware
  # that the number of DE genes for any particular contrasts will depend on which other contrasts are
  # tested at the same time. Hence one should include only those contrasts which are closely related to
  # the question at hand. Unnecessary contrasts should be excluded as these would affect the results for
  # the contrasts of interest. Another more theoretical issue is that there is no theorem which proves that
  # adjust.method="BH" in combination with method="global" will correctly control the false discovery
  # rate for combinations of negatively correlated contrasts, however simulations, experience and some
  # theory suggest that the method is safe in practice.
  
  #DE sites by contrast type
  #DE sites using 'separate' contrasts
  DE <- results[results[,1] != 0 | results[,2] != 0 | results[,3] != 0,]
  #sites only DE in exactly one contrast
  absDE <- abs(DE)
  DE1 <- absDE[rowSums(absDE)==1,]
  #sites DE in exactly two contrast
  DE2 <- absDE[rowSums(absDE)==2,]
  #site DE in all 3 contrasts
  DE3 <- results[results[,1] != 0 & results[,2] != 0 & results[,3] != 0,]
  
  #F stats by contrast type
  
  # #Sorting F Values
  # test <- topTable(fit2, coef = c(1,2,3), adjust = "BH", n=Inf, sort.by="F", p=.05)#equivalent to below
  # test2 <- topTableF(fit2, adjust = "BH", n=Inf, sort.by="F", p=.05)
  
  Fvals <- topTableF(fit2, adjust = "BH", n=Inf, sort.by="F")#all F values
  sigFvals <- topTableF(fit2, adjust = "BH", n=Inf, sort.by="F", p=.05)#gives 1355 compared to 1549 DE total for separate comparisons
  
  #subsets by contrast specific DE
  FDE1 <- Fvals[which(row.names(DE1)%in%row.names(Fvals)),5:6]
  FDE2 <- Fvals[which(row.names(DE2)%in%row.names(Fvals)),5:6]
  FDE3 <- Fvals[which(row.names(DE3)%in%row.names(Fvals)),5:6]
  
  #below gives adjusted pvalue
  FDE1 <- Fvals[match(row.names(DE1), row.names(Fvals), nomatch = F),5:7]
  FDE2 <- Fvals[match(row.names(DE2), row.names(Fvals), nomatch = F),5:7]
  FDE3 <- Fvals[match(row.names(DE3), row.names(Fvals), nomatch = F),5:7]
  
  #boxplot(FDE1$F,FDE2$F,FDE3$F)
  boxplot(log10(FDE1$F),log10(FDE2$F),log10(FDE3$F))
  summary(FDE1$F)
  summary(FDE2$F)
  summary(FDE3$F)
  
  plot(density(log10(FDE1$F)),xlim = c(0,3))
  lines(density(log10(FDE2$F)), col = 2)
  lines(density(log10(FDE3$F)), col = 3)
  
  
  
  
  # GSEA of differentially expressed lists across contrasts
  
  #add annotation to multexpanded DF
  head(row.names(pilot))
  
  #add F test values to the table
  multExpanded1$globalFsig = ifelse(multExpanded1$idmult %in% row.names(sigFvals),"+","-")
  
  
  #add DE to table
  multExpanded1$SubtoDE = ifelse(multExpanded1$idmult %in% row.names(pilot),"+","-")
  multExpanded1$DEcont1 = ifelse(multExpanded1$idmult %in% row.names(sig1),"+","-")
  multExpanded1$DEcont2 = ifelse(multExpanded1$idmult %in% row.names(sig2),"+","-")
  multExpanded1$DEcont3 = ifelse(multExpanded1$idmult %in% row.names(sig3),"+","-")
  
  #add DE direction to table
  multExpanded1$cont1up = ifelse(multExpanded1$idmult %in% row.names(c1up),"+","-")
  multExpanded1$cont1down = ifelse(multExpanded1$idmult %in% row.names(c1down),"+","-")
  multExpanded1$cont2up = ifelse(multExpanded1$idmult %in% row.names(c2up),"+","-")
  multExpanded1$cont2down = ifelse(multExpanded1$idmult %in% row.names(c2down),"+","-")
  multExpanded1$cont3up = ifelse(multExpanded1$idmult %in% row.names(c3up),"+","-")
  multExpanded1$cont3down = ifelse(multExpanded1$idmult %in% row.names(c3down),"+","-")
  
  #write the multExpanded table with DE information
  write.table(multExpanded1,"multExpanded1_withDE.csv", sep=',',col.names=T, row.names=F)
  
  return(multExpanded1)
  
  
  
  
  
  #add F test values to the table loess
  #multExpanded1$globalFsigloess = ifelse(multExpanded1$idmult %in% row.names(sigFvals),"+","-")
  
  
  #add DE to table loess
  # multExpanded1$SubtoDEloess = ifelse(multExpanded1$idmult %in% row.names(pilot),"+","-")
  # multExpanded1$DEcont1loess = ifelse(multExpanded1$idmult %in% row.names(sig1),"+","-")
  # multExpanded1$DEcont2loess = ifelse(multExpanded1$idmult %in% row.names(sig2),"+","-")
  # multExpanded1$DEcont3loess = ifelse(multExpanded1$idmult %in% row.names(sig3),"+","-")
  # 
  # #add DE direction to table loess
  # multExpanded1$cont1uploess = ifelse(multExpanded1$idmult %in% row.names(c1up),"+","-")
  # multExpanded1$cont1downloess = ifelse(multExpanded1$idmult %in% row.names(c1down),"+","-")
  # multExpanded1$cont2uploess = ifelse(multExpanded1$idmult %in% row.names(c2up),"+","-")
  # multExpanded1$cont2downloess = ifelse(multExpanded1$idmult %in% row.names(c2down),"+","-")
  # multExpanded1$cont3uploess = ifelse(multExpanded1$idmult %in% row.names(c3up),"+","-")
  # multExpanded1$cont3downloess = ifelse(multExpanded1$idmult %in% row.names(c3down),"+","-")
  
  
  
  
  
  
  # write output table to perform enrichment analysis in perseus
  #write.table(multExpanded1,"multExpanded1_withDE and loess.csv",sep=',',col.names=T,row.names=F)
  
  
}