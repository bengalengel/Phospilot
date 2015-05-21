Enrichment <- function(multExpanded1_withDE){
  #enrichment accepts ME DF and performs categorical enrichment analysis using one sided FE tests of categorical enrichment for reactome and GO terms. It also performs multiple testing correction. Method 'BH'.
  
  #Here I assign multiple annotations/protetin. That is 1 annotation/observed phosphosite measurement (includes multiple phosphorylation 'sites'). Otherwise those proteins who are multiply phosphorylated would have a greater chance of being enriched. The only disadvantage here is that I may be 'diluting' the enrichment effect of a single phosphorylation mark amongst multiple phosphorylated forms of that site. 
  require(GO.db)
  require(reactome.db)
  
  Enrich <- function(x,y,ontology = c("GO","Reactome")){
    #     This function accepts character vectors of a selected subset and background and returns a DF of adjusted pvalues for categorical enrichment using a one sided fisher's exact test.
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
    for(i in levels(facttemp)){
      #make contingency table and perform FE test one sided for enrichment
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
    DF <- DF[DF$adjpvalue<=.05,]
    if(nrow(DF)>=1){
      if(ontology == "GO" ){
        GOannotation <- select(GO.db, keys = as.character(DF$ID), keytype = "GOID", columns = c("TERM", "ONTOLOGY", "DEFINITION"))
        DF <- cbind(DF,GOannotation)
      }
      if(ontology == "Reactome"){
        ROannotation <- select(reactome.db, keys = as.character(DF$ID), keytype = "REACTOMEID", columns = "PATHNAME")
        DF <- cbind(DF,ROannotation)
      }
    }else{
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
  
  
  #First is GO
  #######################
  #get background
  backgroundGO <- multExpanded1_withDE[multExpanded1_withDE$SubtoDE == "+",]
  backgroundGO <- backgroundGO$GOIDs
  BGGO <- SplitNClean(backgroundGO)
  
  enrichedGO <- multExpanded1_withDE[multExpanded1_withDE$globalFsig == "+",]
  enrichedGO <- enrichedGO$GOIDs
  DEGO <- SplitNClean(enrichedGO)
  
  #for the variance components.
  backgroundGO <- multExpanded1_withDE[multExpanded1_withDE$SubtoVarcomp == "+",]
  backgroundGO <- backgroundGO$GOIDs
  BGGOvarcomp <- SplitNClean(backgroundGO)
  
  #high ind/highbio 
  enrichedGO <- multExpanded1_withDE[multExpanded1_withDE$HighIndVar == "+" & multExpanded1_withDE$HighBioVar == "+",]
  enrichedGO <- enrichedGO$GOIDs
  varcompHIHB <- SplitNClean(enrichedGO)
  
  #low ind/low bio 
  enrichedGO <- multExpanded1_withDE[multExpanded1_withDE$LowIndVar == "+" & multExpanded1_withDE$LowBioVar == "+",]
  enrichedGO <- enrichedGO$GOIDs
  varcompLILB <- SplitNClean(enrichedGO)
  
  #low ind/high bio 
  enrichedGO <- multExpanded1_withDE[multExpanded1_withDE$LowIndVar == "+" & multExpanded1_withDE$HighBioVar == "+",]
  enrichedGO <- enrichedGO$GOIDs
  varcompLIHB <- SplitNClean(enrichedGO)
  
  #high ind/low bio 
  enrichedGO <- multExpanded1_withDE[multExpanded1_withDE$HighIndVar == "+" & multExpanded1_withDE$LowBioVar == "+",]
  enrichedGO <- enrichedGO$GOIDs
  varcompHILB <- SplitNClean(enrichedGO)
  
  #now for the protein normalized DE
  backgroundGO <- multExpanded1_withDE[multExpanded1_withDE$SubtoDEpn =="+",]
  backgroundGO <- backgroundGO$ppGOIDs
  BGGOpn <- SplitNClean(backgroundGO)
    
    
  enrichedGO <- multExpanded1_withDE[multExpanded1_withDE$globalFsigpn == "+",]
  # enrichedGO <- enrichedGO[!duplicated(enrichedGO$id),]
  enrichedGO <- enrichedGO$ppGOIDs
  DEGOpn <- SplitNClean(enrichedGO)

  #for the variance components.
  backgroundGO <- multExpanded1_withDE[multExpanded1_withDE$ppSubtoVarcomp == "+",]
  backgroundGO <- backgroundGO$ppGOIDs
  pnBGGOvarcomp <- SplitNClean(backgroundGO)
  
  #high ind/highbio 
  enrichedGO <- multExpanded1_withDE[multExpanded1_withDE$pnHighIndVar == "+" & multExpanded1_withDE$pnHighBioVar == "+",]
  enrichedGO <- enrichedGO$ppGOIDs
  pnvarcompHIHB <- SplitNClean(enrichedGO)
  
  #low ind/low bio 
  enrichedGO <- multExpanded1_withDE[multExpanded1_withDE$pnLowIndVar == "+" & multExpanded1_withDE$pnLowBioVar == "+",]
  enrichedGO <- enrichedGO$ppGOIDs
  pnvarcompLILB <- SplitNClean(enrichedGO)
  
  #low ind/high bio 
  enrichedGO <- multExpanded1_withDE[multExpanded1_withDE$pnLowIndVar == "+" & multExpanded1_withDE$pnHighBioVar == "+",]
  enrichedGO <- enrichedGO$ppGOIDs
  pnvarcompLIHB <- SplitNClean(enrichedGO)
  
  #high ind/low bio 
  enrichedGO <- multExpanded1_withDE[multExpanded1_withDE$pnHighIndVar == "+" & multExpanded1_withDE$pnLowBioVar == "+",]
  enrichedGO <- enrichedGO$ppGOIDs
  pnvarcompHILB <- SplitNClean(enrichedGO)
  
  
  #run the enrich function
  ############################################################
  GODF <- Enrich(BGGO,DEGO)
  
  #assign the annotation. If no annotation is available, return a character vector "No Enrichment"
  pos_err <- tryCatch(select(GO.db, keys = as.character(GODF$ID), keytype = "GOID", columns = c("TERM", "ONTOLOGY", "DEFINITION")), error=function(e) e)
  if(!inherits(pos_err, "error")){
  GOannotation <- select(GO.db, keys = as.character(GODF$ID), keytype = "GOID", columns = c("TERM", "ONTOLOGY", "DEFINITION"))
  GOenrichment <- cbind(GODF,GOannotation)}else{
    GOenrichment <- "NO Enrichment"}
  
  GODF <- Enrich(BGGOvarcomp,varcompHIHB)#HIHB
  #assign the annotation using GO.db
  GOannotation <- select(GO.db, keys = as.character(GODF$ID), keytype = "GOID", columns = c("TERM", "ONTOLOGY", "DEFINITION"))
  GOenrichHIHB <- cbind(GODF,GOannotation)
  
#   GODF <- Enrich(BGGOvarcomp,varcompHILB)# NOTHING ENRICHED
#   #assign the annotation using GO.db
#   GOannotation <- select(GO.db, keys = as.character(GODF$ID), keytype = "GOID", columns = c("TERM", "ONTOLOGY", "DEFINITION"))
#   GOenrichHILB <- cbind(GODF,GOannotation)
#   
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
  
#   GODF <- Enrich(BGGOvarcomp,varcompHILB)# NOTHING ENRICHED
  #assign the annotation using GO.db
#   GOannotation <- select(GO.db, keys = as.character(GODF$ID), keytype = "GOID", columns = c("TERM", "ONTOLOGY", "DEFINITION"))
#   pnGOenrichHILB <- cbind(GODF,GOannotation)
  
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
  backgroundReact <- multExpanded1_withDE[multExpanded1_withDE$SubtoDE=="+",]
  backgroundReact <- backgroundReact$ReactIDs
  BGRO <- SplitNClean(backgroundReact)
  
  enrichedReact <- multExpanded1_withDE[multExpanded1_withDE$globalFsig == "+",]
  enrichedReact <- enrichedReact$ReactIDs
  DERO <- SplitNClean(enrichedReact)
  
  #for the variance components.
  backgroundReact <- multExpanded1_withDE[multExpanded1_withDE$SubtoVarcomp == "+",]
  backgroundReact <- backgroundReact$ReactIDs
  BGROvarcomp <- SplitNClean(backgroundReact)
  
  #high ind/highbio 
  enrichedReact <- multExpanded1_withDE[multExpanded1_withDE$HighIndVar == "+" & multExpanded1_withDE$HighBioVar == "+",]
  enrichedReact <- enrichedReact$ReactIDs
  ROvarcompHIHB <- SplitNClean(enrichedReact)
  
  #low ind/low bio 
  enrichedReact <- multExpanded1_withDE[multExpanded1_withDE$LowIndVar == "+" & multExpanded1_withDE$LowBioVar == "+",]
  enrichedReact <- enrichedReact$ReactIDs
  ROvarcompLILB <- SplitNClean(enrichedReact)
  
  #low ind/high bio 
  enrichedReact <- multExpanded1_withDE[multExpanded1_withDE$LowIndVar == "+" & multExpanded1_withDE$HighBioVar == "+",]
  enrichedReact <- enrichedReact$ReactIDs
  ROvarcompLIHB <- SplitNClean(enrichedReact)
  
  #high ind/low bio 
  enrichedReact <- multExpanded1_withDE[multExpanded1_withDE$HighIndVar == "+" & multExpanded1_withDE$LowBioVar == "+",]
  enrichedReact <- enrichedReact$ReactIDs
  ROvarcompHILB <- SplitNClean(enrichedReact)
  
  
  #now for the protein normalized DE
  backgroundReact <- multExpanded1_withDE[multExpanded1_withDE$SubtoDEpn == "+",]
  backgroundReact <- backgroundReact$ppReactIDs
  pnBGRO <- SplitNClean(backgroundReact)
  
  enrichedReact <- multExpanded1_withDE[multExpanded1_withDE$globalFsigpn == "+",]
  enrichedReact <- enrichedReact$ppReactIDs
  pnDERO <- SplitNClean(enrichedReact)
  
  #for the variance components.
  backgroundReact <- multExpanded1_withDE[multExpanded1_withDE$ppSubtoVarcomp == "+",]
  backgroundReact <- backgroundReact$ppReactIDs
  pnBGROvarcomp <- SplitNClean(backgroundReact)
  
  #high ind/highbio 
  enrichedReact <- multExpanded1_withDE[multExpanded1_withDE$pnHighIndVar == "+" & multExpanded1_withDE$pnHighBioVar == "+",]
  enrichedReact <- enrichedReact$ppReactIDs
  pnROvarcompHIHB <- SplitNClean(enrichedReact)
  
  #low ind/low bio 
  enrichedReact <- multExpanded1_withDE[multExpanded1_withDE$pnLowIndVar == "+" & multExpanded1_withDE$pnLowBioVar == "+",]
  enrichedReact <- enrichedReact$ppReactIDs
  pnROvarcompLILB <- SplitNClean(enrichedReact)
  
  #low ind/high bio 
  enrichedReact <- multExpanded1_withDE[multExpanded1_withDE$pnLowIndVar == "+" & multExpanded1_withDE$pnHighBioVar == "+",]
  enrichedReact <- enrichedReact$ppReactIDs
  pnROvarcompLIHB <- SplitNClean(enrichedReact)
  
  #high ind/low bio 
  enrichedReact <- multExpanded1_withDE[multExpanded1_withDE$pnHighIndVar == "+" & multExpanded1_withDE$pnLowBioVar == "+",]
  enrichedReact <- enrichedReact$ppReactIDs
  pnROvarcompHILB <- SplitNClean(enrichedReact)
  
  
  
  #run the function
  RODF <- Enrich(BGRO,DERO)
  #assign the annotation
  ROannotation <- select(reactome.db, keys = as.character(RODF$ID), keytype = "REACTOMEID", columns = "PATHNAME")
  ROenrichment <- cbind(RODF,ROannotation)
  
  RODF <- Enrich(BGROvarcomp,ROvarcompHIHB)#HIHB
  #assign the annotation using GO.db
  ROannotation <- select(reactome.db, keys = as.character(RODF$ID), keytype = "REACTOMEID", columns = "PATHNAME")
  ROenrichHIHB <- cbind(RODF,ROannotation)
  
  RODF <- Enrich(BGROvarcomp,ROvarcompHILB)# 
  #assign the annotation using GO.db
  ROannotation <- select(reactome.db, keys = as.character(RODF$ID), keytype = "REACTOMEID", columns = "PATHNAME")
  ROenrichHILB <- cbind(RODF,ROannotation)
  
  RODF <- Enrich(BGROvarcomp,ROvarcompLIHB)#from above
  #assign the annotation using GO.db
  ROannotation <- select(reactome.db, keys = as.character(RODF$ID), keytype = "REACTOMEID", columns = "PATHNAME")
  ROenrichLIHB <- cbind(RODF,ROannotation)
  
#   RODF <- Enrich(BGROvarcomp,ROvarcompLILB)#NOTHING ENRICHED
  #assign the annotation using GO.db
#   ROannotation <- select(reactome.db, keys = as.character(RODF$ID), keytype = "REACTOMEID", columns = "PATHNAME")
#   ROenrichLILB <- cbind(RODF,ROannotation)
  
  
  #protein normalized data
  RODFpn <- Enrich(pnBGRO,pnDERO)
  #assign the annotation
  ROannotationpn <- select(reactome.db, keys = as.character(RODFpn$ID), keytype = "REACTOMEID", columns = "PATHNAME")
  ROenrichmentpn <- cbind(RODFpn,ROannotationpn)  
  
  RODF <- Enrich(pnBGROvarcomp,pnROvarcompHIHB)#HIHB
  #assign the annotation using GO.db
  ROannotation <- select(reactome.db, keys = as.character(RODF$ID), keytype = "REACTOMEID", columns = "PATHNAME")
  pnROenrichHIHB <- cbind(RODF,ROannotation)
  
  RODF <- Enrich(pnBGROvarcomp,pnROvarcompHILB)#
  #assign the annotation using GO.db
  ROannotation <- select(reactome.db, keys = as.character(RODF$ID), keytype = "REACTOMEID", columns = "PATHNAME")
  pnROenrichHILB <- cbind(RODF,ROannotation)
  
  RODF <- Enrich(pnBGROvarcomp,pnROvarcompLIHB)#from above
  #assign the annotation using GO.db
  ROannotation <- select(reactome.db, keys = as.character(RODF$ID), keytype = "REACTOMEID", columns = "PATHNAME")
  pnROenrichLIHB <- cbind(RODF,ROannotation)
  
#   RODF <- Enrich(pnBGROvarcomp,pnROvarcompLILB)#nothing enriched
#   #assign the annotation using GO.db
#   ROannotation <- select(reactome.db, keys = as.character(RODF$ID), keytype = "REACTOMEID", columns = "PATHNAME")
#   pnROenrichLILB <- cbind(pnRODF,pnROannotation)
  
  
  #return a list of enrichment DFs
  DFs <- list(GOenrichment = GOenrichment, GOenrichHIHB = GOenrichHIHB, GOenrichLIHB = GOenrichLIHB, GOenrichLILB = GOenrichLILB, 
              GOenrichmentpn = GOenrichmentpn, pnGOenrichHIHB = pnGOenrichHIHB, pnGOenrichLIHB = pnGOenrichLIHB, pnGOenrichLILB = pnGOenrichLILB,
              ROenrichment = ROenrichment, ROenrichHIHB = ROenrichHIHB, ROenrichHILB = ROenrichHILB, ROenrichLIHB = ROenrichLIHB,
              ROenrichmentpn = ROenrichmentpn, pnROenrichHIHB = pnROenrichHIHB, pnROenrichHILB = pnROenrichHILB, pnROenrichLIHB = pnROenrichLIHB)
  return(DFs)
}

  
  
  
  
  