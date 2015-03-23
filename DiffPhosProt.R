DiffPhosProt <- function(proteinfull, proteinnorm, multExpanded1_withDE, phosphonorm){
##This program assigns protein groups from previous proteomic analysis to phosphosites from the SCX-TiO2 workflow for normalization. 
  
  
  
  
  #first loops for phosphopeptides subjected to diffphos and the second for all phos.
  
  #loop1 - for every protein assigned to a phosphopeptide ask is this protein is found in the zia dataset? 
  nPTable <- length(SubtoDEtable[,1])
  proteinindex <- integer(nPTable)
  #morethan1 <- c()#need to come back to this because this isn't work properly
  
  for(i in seq_along(SubtoDEtable[,1])){#for every protein assigned to a phosphopeptide ask; Is this protein found in the zia dataset? 
    tmp <- grep(SubtoDEtable$Protein[i], Ziaproteins$Majority.protein.IDs)#Searches all IDs in comma separated majority protein list.
    #more than one value? This can happen with isoforms
    if(length(tmp)>1){
      #compare the razor plus unique count across the two matches and choose the one with the most   matches
      counts <- Ziaproteins$Razor...unique.peptides[tmp]
      proteinindex[i] <- tmp[which.max(counts)]
      #morethan1 <- c(morethan1,tmp)
    }
    #if length of tmp >0 add to index
    if(length(tmp) == 1){
      proteinindex[i] <- tmp
    }
    else{
      proteinindex[i] <- NA
    }
  }
  
  SubtoDE$ziaindex <- proteinindex#add the index to the parent dataframe
  
  ##for each proteinindex value add the three protein quantification (H/L) values and the majority ids.
  protein_norm <- data.frame()
  for(i in seq_along(proteinindex)){
    if(!is.na(proteinindex[i])){
      tmp <- datacomp[c("Majority.protein.IDs","HL18862", "HL18486", "HL19160")][proteinindex[i],]#adds the quantile normalized values
      protein_norm <- rbind(protein_norm,tmp)
    }
    if(is.na(proteinindex[i])){
      tmp <- rep(NA,4)
      protein_norm <- rbind(protein_norm,tmp)
      if(dim(protein_norm)[1]==1){#reset names of dataframe in the event the first loop is an NA
        names(protein_norm) <- c("Majority.protein.IDs","HL18862", "HL18486", "HL19160")
        protein_norm$Majority.protein.IDs <- as.factor(protein_norm$Majority.protein.IDs)
        protein_norm$HL18862 <- as.numeric(protein_norm$HL18862)
        protein_norm$HL18486 <- as.numeric(protein_norm$HL18486)
        protein_norm$HL19160 <- as.numeric(protein_norm$HL19160)
      }
    }
  }    
  ####################################################################
  
  
  
  
  #***Differential Expression******************************************************
  
  fac <- factor(c(1,1,2,2,3,3))##codes the grouping for the ttests
  design <- model.matrix(~0 + fac)
  dnames <- levels(as.factor(substr(colnames(ProtNormalized), 1, 7))) ##check me out. use 5 digit exp name.
  colnames(design) <- dnames
  
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
  
  
  
  
  
  ##overlap with phospho dataset
  #add annotation to multExpanded withDE and save the new file
  
  
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
  
  
  ##check for overlap using omnibus F
  
  #subset
  SubtoDEpn <- multExpanded1_withDE[multExpanded1_withDE$SubtoDEpn == "+",]
  
  #for an apples to apples comparison I will need to subject the same subset to DE in both cases.
  
  #limma on common subset for phosphodata
  
  #subset the phospho data by id_mult found in the protein data
  phosdata <- pilot[row.names(pilot) %in% SubtoDEpn$idmult,]
  
  require(limma)
  require(statmod)
  #Produce the design matrix
  
  fac <- factor(c(1,1,2,2,3,3))##codes the grouping for the ttests
  design <- model.matrix(~0 + fac)
  dnames <- levels(as.factor(substr(colnames(phosdata), 1, 7))) ##check me out. use 5 digit exp name.
  colnames(design) <- dnames
  
  fit <- lmFit(phosdata, design)
  
  #Now to make all pairwise comparisons (group2-1, group3-2, group3-1)
  contrast.matrix <- makeContrasts(HL18862-HL18486, HL19160-HL18862, HL19160-HL18486, levels = design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  Fvals <- topTableF(fit2, adjust = "BH", n=Inf, sort.by="F")#all F values
  sigFvals <- topTableF(fit2, adjust = "BH", n=Inf, sort.by="F", p=.05)#gives 1355 compared to 1549 DE total for separate comparisons
  
  diffphos <- nrow(phosdata[row.names(phosdata) %in% row.names(sigFvals),])
  diffphosnorm <- nrow(SubtoDEpn[SubtoDEpn$globalFsigpn == "+",])
  
  ##append to SubtoDEpn
  SubtoDEpn$diffphos = ifelse(SubtoDEpn$idmult %in% row.names(sigFvals),"+","-")
  
  intersection <- nrow(SubtoDEpn[SubtoDEpn$globalFsigpn == "+" & SubtoDEpn$diffphos == "+",])
  
  
  #make a double venn
  require(VennDiagram)
  require(gridExtra)
  plot.new()
  venn.plot <- draw.pairwise.venn(
    area1 = nrow(SubtoDEpn[SubtoDEpn$globalFsig == "+",]),
    area2 = nrow(SubtoDEpn[SubtoDEpn$globalFsigpn == "+",]),
    cross.area = nrow(SubtoDEpn[SubtoDEpn$globalFsigpn == "+" & SubtoDEpn$globalFsig == "+",]),
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
  ##better now that the inverse is used!!!
  
  DFs <- list(multExpanded1_withDE, ProtNormalized)
  return(DFs)
}

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ##proteinnorm is the quantile normalized protein data. proteinfull is the 60 estimates from zias work. pilot is the normalized/batch corrected confounded phospho data from my experimental runs.multExpanded1_withDE is the parent dataframe for the class one sites with DiffExp annotation.  

#subset proteinfull by protein group id found in the three normalized samples of interest. 
proteinfull <- proteinfull[proteinfull$id %in% row.names(proteinnorm),]

#revert names of proteinnorm back to 'HL' for continuity with the program structure below
names(proteinnorm) <- gsub(names(proteinnorm), pattern = "LH", replacement = "HL")

##combine the normalized protein information and the annotation data into a common dataframe for matching to the phospho data
datacomp <- cbind(proteinnorm,proteinfull[c("Protein.IDs","Majority.protein.IDs","Protein.names","Gene.names","Sequence.coverage....",
                                  "Number.of.proteins", "Peptides", "Razor...unique.peptides", "Unique.peptides",
                                  "Razor...unique.peptides.18862", "Razor...unique.peptides.18486", 
                                  "Razor...unique.peptides.19160")])

#how many unique proteins in the phospho data are identified and subjected to DE analysis? Below are counts with and without isoform designation.
######################
SubtoDE <- multExpanded1_withDE[as.character(multExpanded1_withDE$SubtoDE) == "+",]

#remove the reverse for now and add 'AllPhos', a data frame for all phosphopeptides.
SubtoDE <- SubtoDE[!grepl(SubtoDE$Protein, pattern = "REV"),]#the logical subsets the DF by row. 4991 observations.
AllPhos <- multExpanded1_withDE[!grepl(multExpanded1_withDE$Protein, pattern = "REV"),] #18238
  
#How many 'protiens' subject to DE? These are majority protein ids from the protein group assigned to the phosphopeptide.
SubtoDEtable <- as.matrix(table(SubtoDE$Protein))
SubtoDEtable <- SubtoDEtable[SubtoDEtable!=0,,drop=F]#1991 proteins

#How many 'proteins' in the entire class1 dataset?
AllPhostable <- as.matrix(table(AllPhos$Protein))
AllPhostable <- AllPhostable[AllPhostable!=0,,drop=F]#4184 proteins

#Remove isoform designation (this should be noted in the discussion of the paper. Unless there is a specific sequence aligning the peptide to a particular isoform the peptide could reasonably belong to any of the isoforms. This workflow (like those that use only the annotated uniprot databases) assumes this is not the case. Only a minor issue using this approach anyway 1991 to 1968 and 4184 to 4127.

row.names(SubtoDEtable) <- substr(row.names(SubtoDEtable),1,6)#eventually 
SubtoDEproteins <- row.names(SubtoDEtable)#1991
SubtoDEproteins <- unique(SubtoDEproteins)#1968 unique proteins (excluding isoforms) subjected to DE

row.names(AllPhostable) <- substr(row.names(AllPhostable),1,6)#eventually 
AllPhosproteins <- row.names(AllPhostable)#4184
AllPhosproteins <- unique(AllPhosproteins)#4127 unique proteins (excluding isoforms) subjected to DE

###############################################
#How many proteins that are subjected to DE analysis are also IDd and quantified by proteomic analysis (Zia)? 
#Here I will use the majority protein IDs. I will use the majority protein IDs for each protein group quantification (these proteins have at least half the peptides of the leading protein within the group)
Ziaproteins <- datacomp$Majority.protein.IDs

#unparse and turn into a long string for comparison
#now I need to design a loop to dig into the phospho file 
pcount <- strsplit(as.character(Ziaproteins), ";")
pcount <- as.character(unlist(pcount))
pcount <- unique(pcount)
any(duplicated(pcount))#12642 proteins identified with at least 1/2 of the peptides of the majority protein(s) in Zias work within 4270 protein groups

pcount <- substr(pcount,1,6)#now I have some duplicates due to isoform designation (3500 duplicates)
pcount <- unique(pcount) #Now I have 8885 unique protein with at least 1/2 of the peptides of the majority protein(s) within 4270 groups.

table(SubtoDEproteins%in%pcount)#isoform free
#1206 of 1968  sub to diffphos (61.2% are quantified in all three samples in Zia's work)

table(AllPhosproteins%in%pcount)
#2117 of 4127 (51% of proteins)

#####################################################################
#Assigning protein group quantifications to phosphopeptide quantifications
Ziaproteins <- datacomp[c("Majority.protein.IDs","Razor...unique.peptides", "Unique.peptides", "Razor...unique.peptides.18862", "Razor...unique.peptides.18486", "Razor...unique.peptides.19160")]
Ziaproteins$id <- row.names(datacomp)#this is a bit redundant
Ziaproteins$Majority.protein.IDs <- gsub("-.", "", Ziaproteins$Majority.protein.IDs)#removes the isoform indicator

#tables for comparison with protein with protein data. 
SubtoDEtable <- SubtoDE[c("Protein","idmult")]#note that reverse entries have already been removed above
SubtoDEtable$Protein <- substr(SubtoDEtable$Protein,1,6)#remove isoform designation
AllPhostable <- AllPhos[c("Protein","idmult")]#note that reverse entries have already been removed above
AllPhostable$Protein <- substr(AllPhostable$Protein,1,6)#remove isoform designation

#For every protein assigned to an id_mult from the phosphotable, a paired protein group from the ziaproteins table is found (if present) using any of the majority protein ids within that group. If the phospho id maps to multiple protein groups, the one with the most peptides is used.

#first loops for phosphopeptides subjected to diffphos and the second for all phos.

#loop1 - for every protein assigned to a phosphopeptide ask is this protein is found in the zia dataset? 
nPTable <- length(SubtoDEtable[,1])
proteinindex <- integer(nPTable)
#morethan1 <- c()#need to come back to this because this isn't work properly

for(i in seq_along(SubtoDEtable[,1])){#for every protein assigned to a phosphopeptide ask; Is this protein found in the zia dataset? 
  tmp <- grep(SubtoDEtable$Protein[i], Ziaproteins$Majority.protein.IDs)#Searches all IDs in comma separated majority protein list.
  #more than one value? This can happen with isoforms
  if(length(tmp)>1){
    #compare the razor plus unique count across the two matches and choose the one with the most   matches
    counts <- Ziaproteins$Razor...unique.peptides[tmp]
    proteinindex[i] <- tmp[which.max(counts)]
    #morethan1 <- c(morethan1,tmp)
  }
  #if length of tmp >0 add to index
  if(length(tmp) == 1){
    proteinindex[i] <- tmp
  }
  else{
    proteinindex[i] <- NA
  }
}

SubtoDE$ziaindex <- proteinindex#add the index to the parent dataframe

##for each proteinindex value add the three protein quantification (H/L) values and the majority ids.
protein_norm <- data.frame()
for(i in seq_along(proteinindex)){
  if(!is.na(proteinindex[i])){
    tmp <- datacomp[c("Majority.protein.IDs","HL18862", "HL18486", "HL19160")][proteinindex[i],]#adds the quantile normalized values
    protein_norm <- rbind(protein_norm,tmp)
  }
  if(is.na(proteinindex[i])){
    tmp <- rep(NA,4)
    protein_norm <- rbind(protein_norm,tmp)
    if(dim(protein_norm)[1]==1){#reset names of dataframe in the event the first loop is an NA
      names(protein_norm) <- c("Majority.protein.IDs","HL18862", "HL18486", "HL19160")
      protein_norm$Majority.protein.IDs <- as.factor(protein_norm$Majority.protein.IDs)
      protein_norm$HL18862 <- as.numeric(protein_norm$HL18862)
      protein_norm$HL18486 <- as.numeric(protein_norm$HL18486)
      protein_norm$HL19160 <- as.numeric(protein_norm$HL19160)
    }
  }
}    
####################################################################

#all phos loops. I can pre-allocate memory in this loop by defining the length of the poroteinindex!
nPTable <- length(AllPhostable[,1])
proteinindex <- integer(nPTable)
#morethan1 <- c()#need to come back to this because this isn't work properly
for(i in seq_along(AllPhostable[,1])){#for every protein assigned to a phosphopeptide ask is this protein found in the zia dataset? 
  tmp <- grep(AllPhostable$Protein[i], Ziaproteins$Majority.protein.IDs)#Searches all IDs in comma separated majority protein list.
  #more than one value? This can happen with isoforms
  if(length(tmp)>1){
    #compare the razor plus unique count across the two matches and choose the one with the most   matches
    counts <- Ziaproteins$Razor...unique.peptides[tmp]
    proteinindex[i] <- tmp[which.max(counts)]
    #morethan1 <- c(morethan1,tmp)
  }
  #if length of tmp >0 add to index
  if(length(tmp) == 1){
    proteinindex[i] <- tmp
  }
  else{
    proteinindex[i] <- NA
  }
}

AllPhostable$ziaindex <- proteinindex#add the index to the parent dataframe

##for each proteinindex value add the three protein quantification (H/L) values and the majority ids.
protein_norm2 <- data.frame()
for(i in seq_along(proteinindex)){
  if(!is.na(proteinindex[i])){
    tmp <- datacomp[c("Majority.protein.IDs","HL18862", "HL18486", "HL19160")][proteinindex[i],]#uses the quantile normalized values
    protein_norm2 <- rbind(protein_norm2,tmp)
  }
  if(is.na(proteinindex[i])){
    tmp <- rep(NA,4)
    protein_norm2 <- rbind(protein_norm2,tmp)
    if(dim(protein_norm2)[1]==1){#reset names of dataframe in the event the first loop is an NA
      names(protein_norm2) <- c("Majority.protein.IDs","HL18862", "HL18486", "HL19160")
      protein_norm2$Majority.protein.IDs <- as.factor(protein_norm2$Majority.protein.IDs)
      protein_norm2$HL18862 <- as.numeric(protein_norm2$HL18862)
      protein_norm2$HL18486 <- as.numeric(protein_norm2$HL18486)
      protein_norm2$HL19160 <- as.numeric(protein_norm2$HL19160)
    }
  }
}    


################################################################
#link the protein quants to the phospho ids to make a dataframe with normalized protein quants appended. Note "REV_" entries are removed.
AllPhos <- cbind(AllPhos,protein_norm2)
colnames(AllPhos)[87] <- "ProtPrep Majority Protein IDs"

#this is to be returned


#what follows should be within the combined workflow script here using adata but this REV stuff should be omitted because the new multExpanded 
# data frame will have the REV entries removed from diff phos analysis. Perhpas I can add the normalized values to the DF as well. 

#subset adata s.t. peptides were mapped to a protein in Zia's data

#first generate the ids for subsetting
normphos2 <- AllPhos[,c("idmult", "Protein", "Leading.proteins")]
normphos2 <- cbind(normphos2,protein_norm2)
colnames(normphos2)[4] <- "ProtPrep Majority Protein IDs"
MappedPhos <- na.omit(normphos2)

#subset adata
adata2 <- adata[row.names(adata) %in% MappedPhos$idmult , ]

#now subset the MappedPhospho
MappedPhos2 <- MappedPhos[MappedPhos$idmult %in% row.names(adata), ]

#combine
adataprot <- cbind(adata2,MappedPhos2)

#normalize and subject to EDA and varcomp. protein will be run as a covariate for DiffPhos.

##now subset and normalize to make a final table for limma DE and all other downstream analysis
expCol <- grep("HL(.*)", colnames(Phos_Protein))
data <- Phos_Protein[,expCol]
row.names(data) <- Phos_Protein$idmult
data <- na.omit(data)

#perform the normalization
HL18486_1norm <- data$HL18486_1-data$HL18486
HL18486_2norm <- data$HL18486_2-data$HL18486
HL18862_1norm <- data$HL18862_1-data$HL18862
HL18862_2norm <- data$HL18862_2-data$HL18862
HL19160_1norm <- data$HL19160_1-data$HL19160
HL19160_2norm <- data$HL19160_2-data$HL19160

ProtNormalized <- cbind(HL18486_1norm, HL18486_2norm, HL18862_1norm, HL18862_2norm, HL19160_1norm, HL19160_2norm)
row.names(ProtNormalized) <- row.names(data)

boxplot(ProtNormalized)
summary(ProtNormalized)
par(mfrow = c(1, 1))
for (i in 1:(ncol(ProtNormalized))){
  if(i==1) plot(density(ProtNormalized[, i], na.rm=T), col = i, ylim = c(0,1.5))
  else lines(density(ProtNormalized[, i], na.rm=T), col = i)
}





#phospho information. id, id_mult, protein and "pilot dataframe"
normphos <- SubtoDE[,c("id","idmult")]
normphos$PhosphoProtein <- SubtoDE$Protein

# #add the pilot data after removing the REV proteins. (these are not removed upstream because at least one of the majority protein ids is NOT a reverse protein hit)
# RevProtMaps <- which(grepl(multExpanded1_withDE$Protein, pattern = "REV") & multExpanded1_withDE$SubtoDE == "+")
# idmults_rev <- as.character(multExpanded1_withDE$idmult[RevProtMaps])
# 
# #here pilot is all those phosphopeptides that have been subject to diffphos
# #subset out reverse protein maps from normalized batch corrected phoshporylation dataframe
# pilot2 <- pilot[!rownames(pilot) %in% idmults_rev,]


#combine the two
normphos <- cbind(normphos,pilot2)

#now add the protein data to make a table of normalized/processed data from both molecular phenotypes:
Phos_Protein <- cbind(normphos,protein_norm)

#write out this table
write.csv(Phos_Protein, "Phos_Protein.csv", row.names=F)

##now subset and normalize to make a final table for limma DE and all other downstream analysis
expCol <- grep("HL(.*)", colnames(Phos_Protein))
data <- Phos_Protein[,expCol]
row.names(data) <- Phos_Protein$idmult
data <- na.omit(data)

#perform the normalization
HL18486_1norm <- data$HL18486_1-data$HL18486
HL18486_2norm <- data$HL18486_2-data$HL18486
HL18862_1norm <- data$HL18862_1-data$HL18862
HL18862_2norm <- data$HL18862_2-data$HL18862
HL19160_1norm <- data$HL19160_1-data$HL19160
HL19160_2norm <- data$HL19160_2-data$HL19160

ProtNormalized <- cbind(HL18486_1norm, HL18486_2norm, HL18862_1norm, HL18862_2norm, HL19160_1norm, HL19160_2norm)
row.names(ProtNormalized) <- row.names(data)

boxplot(ProtNormalized)
summary(ProtNormalized)
par(mfrow = c(1, 1))
for (i in 1:(ncol(ProtNormalized))){
  if(i==1) plot(density(ProtNormalized[, i], na.rm=T), col = i, ylim = c(0,1.5))
  else lines(density(ProtNormalized[, i], na.rm=T), col = i)
}

#do these still cluster?..scary YES!!!

dataZ <- scale(ProtNormalized)##Z-scored column wise the complete data matrix

# now all data excepting complete cases (note that the sample dendograms look the same)
#hist(dataZ[,6], breaks = 100)

# dendogram using euclidian distance (default) and ward or complete agglomeration
dend.ward<- as.dendrogram(hclust(dist(t(dataZ)),method="ward"))
dend.complete<- as.dendrogram(hclust(dist(t(dataZ))))

ward.o<- order.dendrogram(dend.ward)
complete.o<- order.dendrogram(dend.complete)

plot(dend.complete,ylab="height", main = "Euclidian/Complete")
plot(dend.ward, leaflab = "perpendicular", ylab = "height", main = "Euclidian/Ward")



# row scaled
r <- t(scale(t(ProtNormalized)))#transpose to zscale the rows then transpose back to original format

# sample scaled
c <- scale(ProtNormalized)


# install heatmap.2 package
# install.packages("gplots")
library(gplots)

# Create dendrogram using the data without NAs
feature.dend<- as.dendrogram(hclust(dist(r),method="ward"))
sample.dend<- as.dendrogram(hclust(dist(t(c)),method="ward"))##note that dist caclculates distance between rows by default


##produce the heatmap. Note that the help page has a nice section on identifying subregions by color. Although I will likely have to cut the dendogram to id clusters of interest

heatmap.2(
  r,#row Z scores
  Colv=sample.dend,
  Rowv=feature.dend,
  col=bluered(25),
  scale="none",
  trace="none",
  density.info="none",
  key.xlab = "Row Z scores", key.ylab=NULL, key.title = "",
  srtCol=45,  ,adjCol = c(1,1),
  margins = c(6,5),
  cexCol=1,
  labRow = NA#remove row labels
)


# plot.new()

#PCA analysis 
x <- t(ProtNormalized)#samples are the rows of the column matrix
pc <- prcomp(x)#scale = T, center = T) as of now I am not scaling

cols <- as.factor(substr(colnames(ProtNormalized), 3, 7))##check me out. use 5 digit exp name.
plot(pc$x[, 1], pc$x[, 2], col=as.numeric(cols), main = "PCA", xlab = "PC1", ylab = "PC2")
legend("bottomleft", levels(cols), col = seq(along=levels(cols)), pch = 1)


summary(pc)

#SVD for calculating variance explained; see Rafa's notes for an explaination
cx <- sweep(x, 2, colMeans(x), "-")
sv <- svd(cx)
names(sv)
plot(sv$u[, 1], sv$u[, 2], col = as.numeric(cols), main = "SVD", xlab = "U1", ylab = "U2")


plot(sv$d^2/sum(sv$d^2), xlim = c(1, 12), type = "b", pch = 16, xlab = "principal components", 
     ylab = "variance explained")


#***Differential Expression******************************************************

fac <- factor(c(1,1,2,2,3,3))##codes the grouping for the ttests
design <- model.matrix(~0 + fac)
dnames <- levels(as.factor(substr(colnames(ProtNormalized), 1, 7))) ##check me out. use 5 digit exp name.
colnames(design) <- dnames

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





##overlap with phospho dataset
#add annotation to multExpanded withDE and save the new file


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


##check for overlap using omnibus F

#subset
SubtoDEpn <- multExpanded1_withDE[multExpanded1_withDE$SubtoDEpn == "+",]

#for an apples to apples comparison I will need to subject the same subset to DE in both cases.

#limma on common subset for phosphodata

#subset the phospho data by id_mult found in the protein data
phosdata <- pilot[row.names(pilot) %in% SubtoDEpn$idmult,]

require(limma)
require(statmod)
#Produce the design matrix

fac <- factor(c(1,1,2,2,3,3))##codes the grouping for the ttests
design <- model.matrix(~0 + fac)
dnames <- levels(as.factor(substr(colnames(phosdata), 1, 7))) ##check me out. use 5 digit exp name.
colnames(design) <- dnames

fit <- lmFit(phosdata, design)

#Now to make all pairwise comparisons (group2-1, group3-2, group3-1)
contrast.matrix <- makeContrasts(HL18862-HL18486, HL19160-HL18862, HL19160-HL18486, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

Fvals <- topTableF(fit2, adjust = "BH", n=Inf, sort.by="F")#all F values
sigFvals <- topTableF(fit2, adjust = "BH", n=Inf, sort.by="F", p=.05)#gives 1355 compared to 1549 DE total for separate comparisons

diffphos <- nrow(phosdata[row.names(phosdata) %in% row.names(sigFvals),])
diffphosnorm <- nrow(SubtoDEpn[SubtoDEpn$globalFsigpn == "+",])

##append to SubtoDEpn
SubtoDEpn$diffphos = ifelse(SubtoDEpn$idmult %in% row.names(sigFvals),"+","-")

intersection <- nrow(SubtoDEpn[SubtoDEpn$globalFsigpn == "+" & SubtoDEpn$diffphos == "+",])


#make a double venn
require(VennDiagram)
require(gridExtra)
plot.new()
venn.plot <- draw.pairwise.venn(
  area1 = nrow(SubtoDEpn[SubtoDEpn$globalFsig == "+",]),
  area2 = nrow(SubtoDEpn[SubtoDEpn$globalFsigpn == "+",]),
  cross.area = nrow(SubtoDEpn[SubtoDEpn$globalFsigpn == "+" & SubtoDEpn$globalFsig == "+",]),
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
##better now that the inverse is used!!!

DFs <- list(multExpanded1_withDE, ProtNormalized)
return(DFs)
}
