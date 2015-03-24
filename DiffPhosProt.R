DiffPhosProt <- function(multExpanded1_withDE, phosphonorm){
##For now, this program takes normalized, batch corrected, and biological replicate averaged confounded phospho data (4991x6) and subtracts protein estimates from this dataframe. A limma diffphos pipeline is run on this corrected data. A limma diffphos pipeline is also run with same settings on the confounded data, and a Venn diagram depicting overlap in diffphos calls is produced. 
  require(limma)
  require(statmod)
  require(VennDiagram)
  require(gridExtra)
  
  #subset the protein data
  protein <- multExpanded1_withDE[multExpanded1_withDE$idmult %in% row.names(phosphonorm),]
  row.names(protein) <- protein$idmult
  protein <- protein[,(length(protein)-2):length(protein)]
  protein <- na.omit(protein)
  
  #subset the phospho data
  phosphonorm <- phosphonorm[row.names(phosphonorm) %in% row.names(protein),]
  
  #combined dataframe
  combined <- cbind(phosphonorm,protein)
  
#protein normalized. It is okay to normalize the average in this instance because the sets (techreps) have the same cardinality (1). This needs to be updated s.t. the entire matrix is normalized if technical replicate co-variance is used for lmfit. However the idmults should be subsetted s.t. there is information in at least one biological replicate for each condition like 'adata'.
  HL18486_1norm <- combined$HL18486_1-combined$LH18486
  HL18486_2norm <- combined$HL18486_2-combined$LH18486
  HL18862_1norm <- combined$HL18862_1-combined$LH18862
  HL18862_2norm <- combined$HL18862_2-combined$LH18862
  HL19160_1norm <- combined$HL19160_1-combined$LH19160
  HL19160_2norm <- combined$HL19160_2-combined$LH19160
  
  ProtNormalized <- cbind(HL18486_1norm, HL18486_2norm, HL18862_1norm, HL18862_2norm, HL19160_1norm, HL19160_2norm)
  row.names(ProtNormalized) <- row.names(combined)

  
  
  #***Differential PHOS******************************************************
  fac <- factor(c(1,1,2,2,3,3))##codes the grouping for the ttests
  design <- model.matrix(~0 + fac)
  dnames <- levels(as.factor(substr(colnames(ProtNormalized), 1, 7))) ##use 5 digit exp name.
  colnames(design) <- dnames
  
  #here is the lmFit result. Needed are technical replicate correlation and protein as a blocking factor as opposed to normalizing the data.
  fit <- lmFit(ProtNormalized, design)

  # Now to make all pairwise comparisons (group2-1, group3-2, group3-1)
  contrast.matrix <- makeContrasts(HL18862-HL18486, HL19160-HL18862, HL19160-HL18486, levels = design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
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
  
  plot(tt1$logFC,-log10(tt1$P.Value), xlab = colnames(contrast.matrix)[1], pch = 20, ylab = "-log10(P)",xlim = c(-6, 6))
  #sites with sig difference in comparison 1
  names <- row.names(sig1)
  names2 <- row.names(tt1)
  index <- which(names2 %in% names)
  points(tt1$logFC[index],-log10(tt1$P.Value)[index], col="red3", pch = 20)
  
  
  plot(tt2$logFC,-log10(tt2$P.Value), xlab = colnames(contrast.matrix)[2], pch = 20, ylab = "-log10(P)",xlim = c(-6, 6))
  #sites with sig difference in comparison 1
  names <- row.names(sig2)
  names2 <- row.names(tt2)
  index <- which(names2 %in% names)
  points(tt2$logFC[index],-log10(tt2$P.Value)[index], col="red3", pch = 20)
  
  
  plot(tt3$logFC,-log10(tt3$P.Value), xlab = colnames(contrast.matrix)[3], pch = 20, ylab = "-log10(P)",xlim = c(-6, 6))
  #sites with sig difference in comparison 1
  names <- row.names(sig3)
  names2 <- row.names(tt3)
  index <- which(names2 %in% names)
  points(tt3$logFC[index],-log10(tt3$P.Value)[index], col="red3", pch = 20)
  
  
  results <- decideTests(fit2, adjust.method = "BH", method = "separate")#results is a 'TestResults' matrix
  #separate compares each sample individually and is the default approach
  summary(results)
  
  
  vennDiagram(results, cex=c(1.2,1,0.7)) #good DE across conditions
  vennDiagram(results, cex=c(1.2,1,0.7), include = "up")
  vennDiagram(results, cex=c(1.2,1,0.7), include = "down") 
  vennDiagram(results, cex=c(1.2,1,0.7), include = c("up", "down")) #reproduces 'signature'
  
  
  #F statistics
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
  
  #add annotation to multexpanded DF
  multExpanded1_withDE$SubtoDEpn = ifelse(multExpanded1_withDE$idmult %in% row.names(ProtNormalized),"+","-")
  
  #add F test values to the table
  multExpanded1_withDE$globalFsigpn = ifelse(multExpanded1_withDE$idmult %in% row.names(sigFvals),"+","-")
  
  #add DE to table
  multExpanded1_withDE$DEcont1pn = ifelse(multExpanded1_withDE$idmult %in% row.names(sig1),"+","-")
  multExpanded1_withDE$DEcont2pn = ifelse(multExpanded1_withDE$idmult %in% row.names(sig2),"+","-")
  multExpanded1_withDE$DEcont3pn = ifelse(multExpanded1_withDE$idmult %in% row.names(sig3),"+","-")
  
  #add DE direction to table
  multExpanded1_withDE$cont1uppn = ifelse(multExpanded1_withDE$idmult %in% row.names(c1up),"+","-")
  multExpanded1_withDE$cont1downpn = ifelse(multExpanded1_withDE$idmult %in% row.names(c1down),"+","-")
  multExpanded1_withDE$cont2uppn = ifelse(multExpanded1_withDE$idmult %in% row.names(c2up),"+","-")
  multExpanded1_withDE$cont2downpn = ifelse(multExpanded1_withDE$idmult %in% row.names(c2down),"+","-")
  multExpanded1_withDE$cont3uppn = ifelse(multExpanded1_withDE$idmult %in% row.names(c3up),"+","-")
  multExpanded1_withDE$cont3downpn = ifelse(multExpanded1_withDE$idmult %in% row.names(c3down),"+","-")
  
  #write the multExpanded table with DE information
  write.table(multExpanded1_withDE,"multExpanded1_withDE_protein.csv", sep=',',col.names=T, row.names=F)
  
  
#######################################

##check for overlap between protein normalized and confounded diffphos using omnibus F values from limma. Use only sites sub to diffphos in each condition. I have also run the diffphos on the common subset of confounded phosphosites and had essentially the same result.
###############################
  
  diffphos <- nrow(multExpanded1_withDE[multExpanded1_withDE$globalFsig == "+" & multExpanded1_withDE$SubtoDEpn == "+",])
  diffphosnorm <- nrow(multExpanded1_withDE[multExpanded1_withDE$globalFsigpn == "+",])
  
  
  intersection <- nrow(multExpanded1_withDE[multExpanded1_withDE$globalFsigpn == "+" & multExpanded1_withDE$globalFsig == "+" & 
                                              multExpanded1_withDE$SubtoDEpn == "+",])
  
  #make a double venn
  plot.new()
  venn.plot <- draw.pairwise.venn(
    area1 = diffphos,
    area2 = diffphosnorm,
    cross.area = intersection,
    category = c("PhosDE", "normPhosDE"),
    fill = c("green", "blue"),
    lty = "blank",
    cex = 2,
    cat.cex = 2,
    cat.col = c("green", "blue"), 
    margin = .1,
    main="test"
  )
  plot.new()
  grid.arrange(gTree(children=venn.plot), main="Differential Phosphorylation Overlap")
  
  
  ##almost all of the phosDE is picked up. But there is just as many new DE in the non-confounded data!...
  return(multExpanded1_withDE)
}