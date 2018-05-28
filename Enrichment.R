Enrichment <- function(multExpanded1_withDE){
  #enrichment accepts ME DF and performs categorical enrichment analysis using one sided FE tests of categorical enrichment for reactome, GO, and PSP terms. It also performs multiple testing correction. Method 'BH'.
  
  Enrich <- function(x,y,ontology = c("GO","Reactome", "PSP")){
    # This function accepts character vectors of a selected subset and background and returns a DF of adjusted pvalues for categorical enrichment using a one sided fisher's exact test.
    # x=background and y=enriched. 
    require(plyr)
    BGtable <- as.matrix(table(x))
    #remove entries with 0
    BGtable <- BGtable[BGtable!=0,,drop=F]
    DEtable <- as.matrix(table(y))
    DEtable <- as.matrix(DEtable[row.names(DEtable) %in% row.names(BGtable),])#removing zeros and all factors not present in BG data
    
    #subset the background table
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
  
  
  #GO and reactome enrichments 
  #######################
  
  #confounded
  Confounded.BGGO <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$ConfoundedSubtoDE == "+", "confoundedGOID"]))
  Confounded.DEGO <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$ConfoundedglobalFsig == "+", "confoundedGOID"]))  
  Confounded.BGRO <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$ConfoundedSubtoDE == "+", "confoundedReactIDs"]))
  Confounded.DERO <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$ConfoundedglobalFsig == "+", "confoundedReactIDs"]))  
  
  #phosprep protein covariate
  PhosPrep.BGGO <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$PhosPrepCovSubtoDE == "+", "PhosPrepGOID"]))
  PhosPrep.DEGO <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$PhosPrepCovglobalFsig == "+", "PhosPrepGOID"]))  
  PhosPrep.BGRO <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$PhosPrepCovSubtoDE == "+", "PhosPrepReactIDs"]))
  PhosPrep.DERO <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$PhosPrepCovglobalFsig == "+", "PhosPrepReactIDs"]))  
  
  #gelprep protein covariate
  GelPrep.BGGO <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepCovSubtoDE == "+", "GelPrepGOID"]))
  GelPrep.DEGO <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepCovglobalFsig == "+", "GelPrepGOID"]))  
  GelPrep.BGRO <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepCovSubtoDE == "+", "GelPrepReactIDs"]))
  GelPrep.DERO <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepCovglobalFsig == "+", "GelPrepReactIDs"]))  
  GelPrep.BGPSP <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepCovSubtoDE == "+", "PSPKINASES"]))
  GelPrep.DEPSP <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepCovglobalFsig == "+", "PSPKINASES"]))  
  
  
  
  
  #run the enrich function for GO data. 10 enrichments total
  ############################################################
  #omnibus F enrichment
  Enrich.GO.confounded <- Enrich(Confounded.BGGO, Confounded.DEGO, ontology = "GO")
  Enrich.GO.PhosPrep <- Enrich(PhosPrep.BGGO, PhosPrep.DEGO, ontology = "GO")
  Enrich.GO.GelPrep <- Enrich(GelPrep.BGGO, GelPrep.DEGO, ontology = "GO")
  Enrich.RO.confounded <- Enrich(Confounded.BGRO, Confounded.DERO, ontology = "Reactome")
  Enrich.RO.PhosPrep <- Enrich(PhosPrep.BGRO, PhosPrep.DERO, ontology = "Reactome")
  Enrich.RO.GelPrep <- Enrich(GelPrep.BGRO, GelPrep.DERO, ontology = "Reactome")
  Enrich.PSP.GelPrep <- Enrich(GelPrep.BGPSP, GelPrep.DEPSP, ontology = "PSP")
  
#   NO SIG ENRICHMENT FOR PROT NORMALIZED VALUES. PHOSPREP HAS A BIT BUT THESE ARE LIKELY SPURIOUS
  
  #bonus for gelprep only
  ################################


#phosprep protein covariate
GelPrep.BGGO.cont1 <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepCovSubtoDE == "+", "GelPrepGOID"]))
GelPrep.DEGO.cont1 <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepCovDEcont1 == "+", "GelPrepGOID"]))  
GelPrep.BGRO.cont1 <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepCovSubtoDE == "+", "GelPrepReactIDs"]))
GelPrep.DERO.cont1 <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepCovDEcont1 == "+", "GelPrepReactIDs"]))  
GelPrep.BGPSP.cont1 <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepCovSubtoDE == "+", "PSPKINASES"]))
GelPrep.DEPSP.cont1 <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepCovDEcont1 == "+", "PSPKINASES"]))  


GelPrep.BGGO.cont2 <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepCovSubtoDE == "+", "GelPrepGOID"]))
GelPrep.DEGO.cont2 <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepCovDEcont2 == "+", "GelPrepGOID"]))  
GelPrep.BGRO.cont2 <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepCovSubtoDE == "+", "GelPrepReactIDs"]))
GelPrep.DERO.cont2 <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepCovDEcont2 == "+", "GelPrepReactIDs"])) 
GelPrep.BGPSP.cont2 <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepCovSubtoDE == "+", "PSPKINASES"]))
GelPrep.DEPSP.cont2 <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepCovDEcont2 == "+", "PSPKINASES"])) 



GelPrep.BGGO.cont3 <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepCovSubtoDE == "+", "GelPrepGOID"]))
GelPrep.DEGO.cont3 <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepCovDEcont3 == "+", "GelPrepGOID"]))  
GelPrep.BGRO.cont3 <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepCovSubtoDE == "+", "GelPrepReactIDs"]))
GelPrep.DERO.cont3 <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepCovDEcont3 == "+", "GelPrepReactIDs"]))
GelPrep.BGPSP.cont3 <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepCovSubtoDE == "+", "PSPKINASES"]))
GelPrep.DEPSP.cont3 <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepCovDEcont3 == "+", "PSPKINASES"]))  




#GelPrep normalized contrast enrichments
Enrich.GO.1 <- Enrich(GelPrep.BGGO.cont1, GelPrep.DEGO.cont1, ontology = "GO")
Enrich.GO.2 <- Enrich(GelPrep.BGGO.cont2, GelPrep.DEGO.cont2, ontology = "GO")
Enrich.GO.3 <- Enrich(GelPrep.BGGO.cont3, GelPrep.DEGO.cont3, ontology = "GO")
Enrich.RO.1 <- Enrich(GelPrep.BGRO.cont1, GelPrep.DERO.cont1, ontology = "Reactome")
Enrich.RO.2 <- Enrich(GelPrep.BGRO.cont2, GelPrep.DERO.cont2, ontology = "Reactome")
Enrich.RO.3 <- Enrich(GelPrep.BGRO.cont3, GelPrep.DERO.cont3, ontology = "Reactome")
Enrich.PSP.1 <- Enrich(GelPrep.BGPSP.cont1, GelPrep.DEPSP.cont1, ontology = "PSP")
Enrich.PSP.2 <- Enrich(GelPrep.BGPSP.cont2, GelPrep.DEPSP.cont2, ontology = "PSP")
Enrich.PSP.3 <- Enrich(GelPrep.BGPSP.cont3, GelPrep.DEPSP.cont3, ontology = "PSP")

  
  #return a list of enrichment DFs
  ###############################
  #final all dataframes in environment with names 'enrich'
  EnrichmentDFs <- setNames(lapply(ls(pattern="Enrich."), function(x) {if(class(get(x)) == "data.frame") get(x)}),ls(pattern="Enrich.")) #returns some NULL list elements
  EnrichmentDFs <- EnrichmentDFs[!sapply(EnrichmentDFs,is.null)]#note there is an 'is.null' function! and needed to use sappply for ligical indexing
  return(EnrichmentDFs)
  
  #note if there is no enrichment nothing will be returned!
}

  
  
  
  
  