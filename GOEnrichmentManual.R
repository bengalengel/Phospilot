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
  #Confounded
  backgroundGO <- multExpanded1_withDE[multExpanded1_withDE$SubtoDEConfounded == "+",]
  backgroundGO <- backgroundGO$GOIDs
  BGGO <- SplitNClean(backgroundGO)
  
  enrichedGO <- multExpanded1_withDE[multExpanded1_withDE$globalFsigConfounded == "+",]
  enrichedGO <- enrichedGO$GOIDs
  DEGO <- SplitNClean(enrichedGO)
  
    
#   #PhosPrep on hold until annotation function is updated
#   backgroundGO <- multExpanded1_withDE[multExpanded1_withDE$SubtoDEPhosProt =="+",]
#   backgroundGO <- backgroundGO$ppGOIDs
#   BGGOPhosPrep <- SplitNClean(backgroundGO)
#   
#   
#   enrichedGO <- multExpanded1_withDE[multExpanded1_withDE$globalFsigPhosProt == "+",]
#   # enrichedGO <- enrichedGO[!duplicated(enrichedGO$id),]
#   enrichedGO <- enrichedGO$ppGOIDs
#   DEGOPhosPrep <- SplitNClean(enrichedGO)
  
  #GelPrep
  backgroundGO <- multExpanded1_withDE[multExpanded1_withDE$SubtoDEGelProt =="+",]
  backgroundGO <- backgroundGO$ppGOIDs
  BGGOGelPrep <- SplitNClean(backgroundGO)
  
  
  enrichedGO <- multExpanded1_withDE[multExpanded1_withDE$globalFsigGelProt == "+",]
  # enrichedGO <- enrichedGO[!duplicated(enrichedGO$id),]
  enrichedGO <- enrichedGO$ppGOIDs
  DEGOGelPrep  <- SplitNClean(enrichedGO)
  
  #run the enrich function for GO data. 10 enrichments total
  ############################################################
  #omnibus F enrichments
  GOenrichment <- Enrich(BGGO, DEGO, ontology = "GO")
  
  #GelPrep protein as covariate data
  #omnibus F enrichment
  GOenrichmentGelPrep <- Enrich(BGGOGelPrep, DEGOGelPrep, ontology = "GO")
  
  
  
  #same steps for reactome
  ###################
  backgroundReact <- multExpanded1_withDE[multExpanded1_withDE$SubtoDEConfounded == "+",]
  backgroundReact <- backgroundReact$ReactIDs
  BGRO <- SplitNClean(backgroundReact)
  
  enrichedReact <- multExpanded1_withDE[multExpanded1_withDE$globalFsigConfounded == "+",]
  enrichedReact <- enrichedReact$ReactIDs
  DERO <- SplitNClean(enrichedReact)
  
  
  #now for the protein normalized DE
  backgroundReact <- multExpanded1_withDE[multExpanded1_withDE$SubtoDEGelProt == "+",]
  backgroundReact <- backgroundReact$ppReactIDs
  BGROGelPrep <- SplitNClean(backgroundReact)
  
  enrichedReact <- multExpanded1_withDE[multExpanded1_withDE$globalFsigGelProt == "+",]
  enrichedReact <- enrichedReact$ppReactIDs
  DEROGelPrep <- SplitNClean(enrichedReact)
 
  
  #run the function for Reactome. 10 Enrichments total
  ##################################
  
  #omnibus F
  ROenrichment <- Enrich(BGRO, DERO, ontology = "Reactome")
  #protein normalized data
  #omnibus F
  ROenrichmentGelPrep <- Enrich(BGROGelPrep, DEROGelPrep, ontology = "Reactome")

  
  #return a list of enrichment DFs
  ###############################
  #final all dataframes in environment with names 'enrich'
  EnrichmentDFs <- setNames(lapply(ls(pattern="enrich"), function(x) {if(class(get(x)) == "data.frame") get(x)}),ls(pattern="enrich")) #returns some NULL list elements
  EnrichmentDFs <- EnrichmentDFs[!sapply(EnrichmentDFs,is.null)]#note there is an 'is.null' function! and needed to use sappply for ligical indexing
  return(EnrichmentDFs)
  
  #note if there is no enrichment nothing will be returned!
}





