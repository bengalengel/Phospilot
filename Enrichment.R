Enrichment <- function(multExpanded1_withDE){
  #enrichment accepts ME DF and performs categorical enrichment analysis using one sided FE tests of categorical enrichment for reactome and GO terms. It also performs multiple testing correction. Method 'BH'.
  
  #Here I assign multiple annotations/protetin. That is 1 annotation/observed phosphosite measurement (includes multiple phosphorylation 'sites'). Otherwise those proteins who are multiply phosphorylated would have a greater chance of being enriched. The only disadvantage here is that I may be 'diluting' the enrichment effect of a single phosphorylation mark amongst multiple phosphorylated forms of that site. 
  
  Enrich <- function(x,y,ontology = c("GO","Reactome")){
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
    DF <- DF[DF$adjpvalue<=.05,]
    if(nrow(DF)>=1){
      if(ontology == "GO" ){
        require(GO.db)
        GOannotation <- select(GO.db, keys = as.character(DF$ID), keytype = "GOID", columns = c("TERM", "ONTOLOGY", "DEFINITION"))
        DF <- cbind(DF,GOannotation)
      }
      if(ontology == "Reactome"){
        require(reactome.db)
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
  
  
  #run the enrich function for GO data. 10 enrichments total
  ############################################################
  #omnibus F enrichment
  GOenrichment <- Enrich(BGGO, DEGO, ontology = "GO")
  
  #varcomp enrichments
  GOenrichHIHB <- Enrich(BGGOvarcomp, varcompHIHB, ontology = "GO")
  GOenrichHILB <- Enrich(BGGOvarcomp, varcompHILB, ontology = "GO")# NOTHING ENRICHED
  GOenrichLIHB <- Enrich(BGGOvarcomp, varcompLIHB, ontology = "GO")
  GOenrichLILB <- Enrich(BGGOvarcomp, varcompLILB, ontology = "GO")
  
  #protein normalized data
  
  #omnibus F enrichment
  GOenrichmentpn <- Enrich(BGGOpn, DEGOpn, ontology = "GO")
  
  #varcomp enrichments
  pnGOenrichHIHB <- Enrich(pnBGGOvarcomp, pnvarcompHIHB, ontology = "GO")
  pnGOenrichHILB <- Enrich(pnBGGOvarcomp, pnvarcompHILB, ontology = "GO")# NOTHING ENRICHED
  pnGOenrichLIHB <- Enrich(pnBGGOvarcomp, pnvarcompLIHB, ontology = "GO")
  pnGOenrichLILB <- Enrich(pnBGGOvarcomp, pnvarcompLILB, ontology = "GO")


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
  
  #run the function for Reactome. 10 Enrichments total
  ##################################
  
  #omnibus F
  ROenrichment <- Enrich(BGRO, DERO, ontology = "Reactome")
  #varcomps
  ROenrichHIHB <- Enrich(BGROvarcomp, ROvarcompHIHB, ontology = "Reactome")
  ROenrichHILB <- Enrich(BGROvarcomp, ROvarcompHILB, ontology = "Reactome")
  ROenrichLIHB <- Enrich(BGROvarcomp, ROvarcompLIHB, ontology = "Reactome")#from above
  ROenrichLILB <- Enrich(BGROvarcomp, ROvarcompLILB, ontology = "Reactome")#NOTHING ENRICHED
  
  #protein normalized data
  #omnibus F
  ROenrichmentpn <- Enrich(pnBGRO, pnDERO, ontology = "Reactome")
  #varcomps
  pnROenrichHIHB <- Enrich(pnBGROvarcomp, pnROvarcompHIHB, ontology = "Reactome")#HIHB
  pnROenrichHILB <- Enrich(pnBGROvarcomp, pnROvarcompHILB, ontology = "Reactome")#
  pnROenrichLIHB <- Enrich(pnBGROvarcomp, pnROvarcompLIHB, ontology = "Reactome")#from above
  pnROenrichLILB <- Enrich(pnBGROvarcomp, pnROvarcompLILB, ontology = "Reactome")#nothing enriched
  
  
  #return a list of enrichment DFs
  ###############################
  #final all dataframes in environment with names 'enrich'
  EnrichmentDFs <- setNames(lapply(ls(pattern="enrich"), function(x) {if(class(get(x)) == "data.frame") get(x)}),ls(pattern="enrich")) #returns some NULL list elements
  EnrichmentDFs <- EnrichmentDFs[!sapply(EnrichmentDFs,is.null)]#note there is an 'is.null' function! and needed to use sappply for ligical indexing
  return(EnrichmentDFs)
  
  #note if there is no enrichment nothing will be returned!
}

  
  
  
  
  