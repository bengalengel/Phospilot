Enrichment <- function(multExpanded1_withDE){
  #enrichment accepts ME DF and performs categorical enrichment analysis using one sided FE tests of categorical enrichment for reactome and GO terms. It also performs multiple testing correction. Method 'BH'.
  
  #Here I assign multiple annotations/protetin. That is 1 annotation/observed phosphosite measurement (includes multiple phosphorylation 'sites'). Otherwise those proteins who are multiply phosphorylated would have a greater chance of being enriched. The only disadvantage here is that I may be 'diluting' the enrichment effect of a single phosphorylation mark amongst multiple phosphorylated forms of that site. 
  require(GO.db)
  require(reactome.db)
  
  Enrich <- function(x,y){
    #     This function accepts a DE subset and background subset and returns a list of adjusted pvalues for categorical enrichment using a one sided fisher's exact test.
    #     x=background and y=enriched. 
    require(plyr)
    BGtable <- as.matrix(table(x))#the parenthetical should be a factor vector of ids passed to this function
    #remove entries with 0
    BGtable <- BGtable[BGtable!=0,,drop=F]
    DEtable <- as.matrix(table(y))#the parenthetical should be a passed DE factor vector of networKin output
    DEtable <- as.matrix(DEtable[row.names(DEtable) %in% row.names(BGtable),])#removing zeros and all factors not present in BG data. Should not have an affect actually.
    #subset the background table in a similar way to ensure we are making the proper comparisons
    BGtable <- as.matrix(BGtable[row.names(BGtable) %in% row.names(DEtable),])
    NotDE <- BGtable-DEtable
    facttemp <- as.factor(row.names(DEtable))
    pvals <- c()
    frequency <- c()
    ids <- c()
    #for each unique kinase in the DE1 kinases table I need to make a contingency table
    for(i in levels(facttemp)){#this would be the DE dataframe
      #make the first row of the contingency table
      if(DEtable[as.character(i),] & NotDE[as.character(i),] >= 0){#this condition should always be true
        DErow <- c(DEtable[as.character(i),],(sum(DEtable)-DEtable[as.character(i),]))
        NotDErow <- c(NotDE[as.character(i),],(sum(NotDE)-NotDE[as.character(i),]))
        contmatrix <- rbind(DErow,NotDErow)
        #       > contmatrix
        #                     GO:0000002  notGO:0000002
        #       DErow             6      27817
        #       NotDErow         15      51891
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
    DF <- DF[DF$adjpvalue<=.05,]
    return(DF)
  }
  
  ##########################

  #First is GO
  #######################
  #get background
  backgroundGO <- multExpanded1_withDE[multExpanded1_withDE$SubtoDE == "+",]
  #remove duplicate entries of the same peptide (multiplicities). Subset by unique id. id refers to a phosphosite
  # backgroundGO <- backgroundGO[!duplicated(backgroundGO$id),]#the duplicated returns a logical identifying the first instance of duplication from an ordererd query. Therefore it is safe to subset using 'duplicated'
  backgroundGO <- backgroundGO$GOIDs
  BGGO <- sapply(backgroundGO, function(x) strsplit(x, ";"))
  BGGO <- unlist(BGGO)
  index <- which(BGGO=="")#some peptides were not assigned an GOid
  BGGO[index] <- NA
  BGGO <- na.omit(BGGO)#5299 unique
  
  enrichedGO <- multExpanded1_withDE[multExpanded1_withDE$globalFsig == "+",]
  # enrichedGO <- enrichedGO[!duplicated(enrichedGO$id),]
  enrichedGO <- enrichedGO$GOIDs
  DEGO <- sapply(enrichedGO, function(y) strsplit(y, ";"))
  DEGO <- unlist(DEGO)
  index <- which(DEGO=="")
  DEGO[index] <- NA
  DEGO <- na.omit(DEGO)#3474 unique
  
  #for the variance components.
  backgroundGO <- multExpanded1_withDE[multExpanded1_withDE$SubtoVarcomp == "+",]
  #remove duplicate entries of the same peptide (multiplicities). Subset by unique id. id refers to a phosphosite
  # backgroundGO <- backgroundGO[!duplicated(backgroundGO$id),]#the duplicated returns a logical identifying the first instance of duplication from an ordererd query. Therefore it is safe to subset using 'duplicated'
  backgroundGO <- backgroundGO$GOIDs
  BGGOvarcomp <- sapply(backgroundGO, function(x) strsplit(x, ";"))
  BGGOvarcomp <- unlist(BGGOvarcomp)
  index <- which(BGGOvarcomp=="")#some peptides were not assigned an GOid
  BGGOvarcomp[index] <- NA
  BGGOvarcomp <- na.omit(BGGOvarcomp)#5299 unique
  
  #high ind/highbio 
  enrichedGO <- multExpanded1_withDE[multExpanded1_withDE$HighIndVar == "+" & multExpanded1_withDE$HighBioVar == "+",]
  # enrichedGO <- enrichedGO[!duplicated(enrichedGO$id),]
  enrichedGO <- enrichedGO$GOIDs
  varcompHIHB <- sapply(enrichedGO, function(y) strsplit(y, ";"))
  varcompHIHB <- unlist(varcompHIHB)
  index <- which(varcompHIHB=="")
  varcompHIHB[index] <- NA
  varcompHIHB <- na.omit(varcompHIHB)
  
  #low ind/low bio 
  enrichedGO <- multExpanded1_withDE[multExpanded1_withDE$LowIndVar == "+" & multExpanded1_withDE$LowBioVar == "+",]
  # enrichedGO <- enrichedGO[!duplicated(enrichedGO$id),]
  enrichedGO <- enrichedGO$GOIDs
  varcompLILB <- sapply(enrichedGO, function(y) strsplit(y, ";"))
  varcompLILB <- unlist(varcompLILB)
  index <- which(varcompLILB == "")
  varcompLILB[index] <- NA
  varcompLILB <- na.omit(varcompLILB)
  
  #low ind/high bio 
  enrichedGO <- multExpanded1_withDE[multExpanded1_withDE$LowIndVar == "+" & multExpanded1_withDE$HighBioVar == "+",]
  # enrichedGO <- enrichedGO[!duplicated(enrichedGO$id),]
  enrichedGO <- enrichedGO$GOIDs
  varcompLIHB <- sapply(enrichedGO, function(y) strsplit(y, ";"))
  varcompLIHB <- unlist(varcompLIHB)
  index <- which(varcompLIHB == "")
  varcompLIHB[index] <- NA
  varcompLIHB <- na.omit(varcompLIHB)
  
  #high ind/low bio 
  enrichedGO <- multExpanded1_withDE[multExpanded1_withDE$HighIndVar == "+" & multExpanded1_withDE$LowBioVar == "+",]
  # enrichedGO <- enrichedGO[!duplicated(enrichedGO$id),]
  enrichedGO <- enrichedGO$GOIDs
  varcompHILB <- sapply(enrichedGO, function(y) strsplit(y, ";"))
  varcompHILB <- unlist(varcompHILB)
  index <- which(varcompHILB == "")
  varcompHILB[index] <- NA
  varcompHILB <- na.omit(varcompHILB)
  
  
  
  
  #now for the protein normalized DE
  backgroundGO <- multExpanded1_withDE[multExpanded1_withDE$SubtoDEpn =="+",]
  #remove duplicate entries of the same peptide (multiplicities). Subset by unique id. id refers to a phosphosite
  # backgroundGO <- backgroundGO[!duplicated(backgroundGO$id),]#the duplicated returns a logical identifying the first instance of duplication from an ordererd query. Therefore it is safe to subset using 'duplicated'
  backgroundGO <- backgroundGO$ppGOIDs
  BGGOpn <- sapply(backgroundGO, function(x) strsplit(x, ";"))
  BGGOpn <- unlist(BGGOpn)
  index <- which(BGGOpn=="")#some peptides were not assigned an GOid
  BGGOpn[index] <- NA
  BGGOpn <- na.omit(BGGOpn)#to be passed
  
  enrichedGO <- multExpanded1_withDE[multExpanded1_withDE$globalFsigpn == "+",]
  # enrichedGO <- enrichedGO[!duplicated(enrichedGO$id),]
  enrichedGO <- enrichedGO$GOIDs
  DEGOpn <- sapply(enrichedGO, function(y) strsplit(y, ";"))
  DEGOpn <- unlist(DEGOpn)
  index <- which(DEGOpn=="")
  DEGOpn[index] <- NA
  DEGOpn <- na.omit(DEGOpn)

  #for the variance components.
  backgroundGO <- multExpanded1_withDE[multExpanded1_withDE$ppSubtoVarcomp == "+",]
  #remove duplicate entries of the same peptide (multiplicities). Subset by unique id. id refers to a phosphosite
  # backgroundGO <- backgroundGO[!duplicated(backgroundGO$id),]#the duplicated returns a logical identifying the first instance of duplication from an ordererd query. Therefore it is safe to subset using 'duplicated'
  backgroundGO <- backgroundGO$GOIDs
  pnBGGOvarcomp <- sapply(backgroundGO, function(x) strsplit(x, ";"))
  pnBGGOvarcomp <- unlist(pnBGGOvarcomp)
  index <- which(pnBGGOvarcomp=="")#some peptides were not assigned an GOid
  pnBGGOvarcomp[index] <- NA
  pnBGGOvarcomp <- na.omit(pnBGGOvarcomp)#5299 unique
  
  #high ind/highbio 
  enrichedGO <- multExpanded1_withDE[multExpanded1_withDE$pnHighIndVar == "+" & multExpanded1_withDE$pnHighBioVar == "+",]
  # enrichedGO <- enrichedGO[!duplicated(enrichedGO$id),]
  enrichedGO <- enrichedGO$GOIDs
  pnvarcompHIHB <- sapply(enrichedGO, function(y) strsplit(y, ";"))
  pnvarcompHIHB <- unlist(pnvarcompHIHB)
  index <- which(pnvarcompHIHB=="")
  pnvarcompHIHB[index] <- NA
  pnvarcompHIHB <- na.omit(pnvarcompHIHB)
  
  #low ind/low bio 
  enrichedGO <- multExpanded1_withDE[multExpanded1_withDE$pnLowIndVar == "+" & multExpanded1_withDE$pnLowBioVar == "+",]
  # enrichedGO <- enrichedGO[!duplicated(enrichedGO$id),]
  enrichedGO <- enrichedGO$GOIDs
  pnvarcompLILB <- sapply(enrichedGO, function(y) strsplit(y, ";"))
  pnvarcompLILB <- unlist(pnvarcompLILB)
  index <- which(pnvarcompLILB == "")
  pnvarcompLILB[index] <- NA
  pnvarcompLILB <- na.omit(pnvarcompLILB)
  
  #low ind/high bio 
  enrichedGO <- multExpanded1_withDE[multExpanded1_withDE$pnLowIndVar == "+" & multExpanded1_withDE$pnHighBioVar == "+",]
  # enrichedGO <- enrichedGO[!duplicated(enrichedGO$id),]
  enrichedGO <- enrichedGO$GOIDs
  pnvarcompLIHB <- sapply(enrichedGO, function(y) strsplit(y, ";"))
  pnvarcompLIHB <- unlist(pnvarcompLIHB)
  index <- which(pnvarcompLIHB == "")
  pnvarcompLIHB[index] <- NA
  pnvarcompLIHB <- na.omit(pnvarcompLIHB)
  
  #high ind/low bio 
  enrichedGO <- multExpanded1_withDE[multExpanded1_withDE$pnHighIndVar == "+" & multExpanded1_withDE$pnLowBioVar == "+",]
  # enrichedGO <- enrichedGO[!duplicated(enrichedGO$id),]
  enrichedGO <- enrichedGO$GOIDs
  pnvarcompHILB <- sapply(enrichedGO, function(y) strsplit(y, ";"))
  pnvarcompHILB <- unlist(pnvarcompHILB)
  index <- which(pnvarcompHILB == "")
  pnvarcompHILB[index] <- NA
  pnvarcompHILB <- na.omit(pnvarcompHILB)
  
  
  
  #run the enrich function
  ############################################################
  GODF <- Enrich(BGGO,DEGO)#from above
  #assign the annotation using GO.db
  GOannotation <- select(GO.db, keys = as.character(GODF$ID), keytype = "GOID", columns = c("TERM", "ONTOLOGY", "DEFINITION"))
  GOenrichment <- cbind(GODF,GOannotation)
  
  GODF <- Enrich(BGGOvarcomp,varcompHIHB)#HIHB
  #assign the annotation using GO.db
  GOannotation <- select(GO.db, keys = as.character(GODF$ID), keytype = "GOID", columns = c("TERM", "ONTOLOGY", "DEFINITION"))
  GOenrichHIHB <- cbind(GODF,GOannotation)
  
  GODF <- Enrich(BGGOvarcomp,varcompHILB)# NOTHING ENRICHED
  #assign the annotation using GO.db
  GOannotation <- select(GO.db, keys = as.character(GODF$ID), keytype = "GOID", columns = c("TERM", "ONTOLOGY", "DEFINITION"))
  GOenrichHILB <- cbind(GODF,GOannotation)
  
  GODF <- Enrich(BGGOvarcomp,varcompLIHB)#from above
  #assign the annotation using GO.db
  GOannotation <- select(GO.db, keys = as.character(GODF$ID), keytype = "GOID", columns = c("TERM", "ONTOLOGY", "DEFINITION"))
  GOenrichLIHB <- cbind(GODF,GOannotation)
  
  GODF <- Enrich(BGGOvarcomp,varcompLILB)#from above
  #assign the annotation using GO.db
  GOannotation <- select(GO.db, keys = as.character(GODF$ID), keytype = "GOID", columns = c("TERM", "ONTOLOGY", "DEFINITION"))
  GOenrichLILB <- cbind(GODF,GOannotation)
  
  
  
  #run the function for protein normalized data
  GODFpn <- Enrich(BGGOpn,DEGOpn)#from above
  #assign the annotation
  GOannotationpn <- select(GO.db, keys = as.character(GODFpn$ID), keytype = "GOID", columns = c("TERM", "ONTOLOGY", "DEFINITION"))
  GOenrichmentpn <- cbind(GODFpn,GOannotationpn)
  
  GODF <- Enrich(pnBGGOvarcomp,pnvarcompHIHB)#HIHB
  #assign the annotation using GO.db
  GOannotation <- select(GO.db, keys = as.character(GODF$ID), keytype = "GOID", columns = c("TERM", "ONTOLOGY", "DEFINITION"))
  pnGOenrichHIHB <- cbind(GODF,GOannotation)
  
  GODF <- Enrich(BGGOvarcomp,varcompHILB)# NOTHING ENRICHED
  #assign the annotation using GO.db
  GOannotation <- select(GO.db, keys = as.character(GODF$ID), keytype = "GOID", columns = c("TERM", "ONTOLOGY", "DEFINITION"))
  pnGOenrichHILB <- cbind(GODF,GOannotation)
  
  GODF <- Enrich(BGGOvarcomp,varcompLIHB)#from above
  #assign the annotation using GO.db
  GOannotation <- select(GO.db, keys = as.character(GODF$ID), keytype = "GOID", columns = c("TERM", "ONTOLOGY", "DEFINITION"))
  pnGOenrichLIHB <- cbind(GODF,GOannotation)
  
  GODF <- Enrich(BGGOvarcomp,varcompLILB)#from above
  #assign the annotation using GO.db
  GOannotation <- select(GO.db, keys = as.character(GODF$ID), keytype = "GOID", columns = c("TERM", "ONTOLOGY", "DEFINITION"))
  pnGOenrichLILB <- cbind(GODF,GOannotation)  
  
  
  
  
  
  
  
  
  #very big difference between confounded and non-confounded data that may be due to protein level confounding and use of 'leading protein' in phospho data.
  x <- GOenrichment$GOID#101
  y <- GOannotationpn$GOID#76
  
  length(intersect(y,x))#11
  length(setdiff(y,x))#65
  
  
  ################################################
  
  #same steps for reactome
  ###################
  #get background and correct for multiple observations of the same peptides (multiplicities)
  backgroundReact <- multExpanded1_withDE[multExpanded1_withDE$SubtoDE=="+",]
  backgroundReact <- backgroundReact$ReactIDs
  BGRO <- sapply(backgroundReact, function(x) strsplit(x, ";"))
  BGRO <- unlist(BGRO)
  index <- which(BGRO == "")
  BGRO[index] <- NA
  BGRO <- na.omit(BGRO)
  
  enrichedReact <- multExpanded1_withDE[multExpanded1_withDE$globalFsig == "+",]
  enrichedReact <- enrichedReact$ReactIDs
  DERO <- sapply(enrichedReact, function(y) strsplit(y, ";"))
  DERO <- unlist(DERO)
  index <- which(DERO == "")
  DERO[index] <- NA
  DERO <- na.omit(DERO)
  
  #now for the protein normalized DE
  backgroundReact <- multExpanded1_withDE[multExpanded1_withDE$SubtoDEpn == "+",]
  backgroundReact <- backgroundReact$ppReactIDs
  pnBGRO <- sapply(backgroundReact, function(x) strsplit(x, ";"))
  pnBGRO <- unlist(pnBGRO)
  index <- which(pnBGRO == "")
  pnBGRO[index] <- NA
  pnBGRO <- na.omit(pnBGRO)
  
  enrichedReact <- multExpanded1_withDE[multExpanded1_withDE$globalFsigpn == "+",]
  enrichedReact <- enrichedReact$ppReactIDs
  pnDERO <- sapply(enrichedReact, function(y) strsplit(y, ";"))
  pnDERO <- unlist(pnDERO)
  index <- which(pnDERO == "")
  pnDERO[index] <- NA
  pnDERO <- na.omit(pnDERO)
  
  #run the function
  RODF <- Enrich(BGRO,DERO)
  
  #assign the annotation
  ROannotation <- select(reactome.db, keys = as.character(RODF$ID), keytype = "REACTOMEID", columns = "PATHNAME")
  ROenrichment <- cbind(RODF,ROannotation)
  
  #protein normalized data
  RODFpn <- Enrich(pnBGRO,pnDERO)
  
  #assign the annotation
  ROannotationpn <- select(reactome.db, keys = as.character(RODFpn$ID), keytype = "REACTOMEID", columns = "PATHNAME")
  ROenrichmentpn <- cbind(RODFpn,ROannotationpn)  
  
  #return a list of enrichment DFs
  DFs <- list(GOenrichment, GOenrichmentpn, ROenrichment, ROenrichmentpn)
  return(DFs)
}

  
  
  
  
  