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
  
  #phosprep protein covariate
  GelPrep.BGGO <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepNormSubtoDE == "+", "GelPrepGOID"]))
  GelPrep.DEGO <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepNormglobalFsig == "+", "GelPrepGOID"]))  
  GelPrep.BGRO <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepNormSubtoDE == "+", "GelPrepReactIDs"]))
  GelPrep.DERO <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepNormglobalFsig == "+", "GelPrepReactIDs"]))  
  
  
  #run the enrich function for GO data. 10 enrichments total
  ############################################################
  #omnibus F enrichment
  Enrich.GO.confounded <- Enrich(Confounded.BGGO, Confounded.DEGO, ontology = "GO")
  Enrich.GO.PhosPrep <- Enrich(PhosPrep.BGGO, PhosPrep.DEGO, ontology = "GO")
  Enrich.GO.GelPrep <- Enrich(GelPrep.BGGO, GelPrep.DEGO, ontology = "GO")
  Enrich.RO.confounded <- Enrich(Confounded.BGRO, Confounded.DERO, ontology = "Reactome")
  Enrich.RO.PhosPrep <- Enrich(PhosPrep.BGRO, PhosPrep.DERO, ontology = "Reactome")
  Enrich.RO.GelPrep <- Enrich(GelPrep.BGRO, GelPrep.DERO, ontology = "Reactome")
  
#   NO SIG ENRICHMENT FOR PROT NORMALIZED VALUES. PHOSPREP HAS A BIT BUT THESE ARE LIKELY SPURIOUS
  
  #bonus for gelprep only
  ################################


#phosprep protein covariate
GelPrep.BGGO.cont1 <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepNormSubtoDE == "+", "GelPrepGOID"]))
GelPrep.DEGO.cont1 <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepNormDEcont1 == "+", "GelPrepGOID"]))  
GelPrep.BGRO.cont1 <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepNormSubtoDE == "+", "GelPrepReactIDs"]))
GelPrep.DERO.cont1 <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepNormDEcont1 == "+", "GelPrepReactIDs"]))  

GelPrep.BGGO.cont2 <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepNormSubtoDE == "+", "GelPrepGOID"]))
GelPrep.DEGO.cont2 <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepNormDEcont2 == "+", "GelPrepGOID"]))  
GelPrep.BGRO.cont2 <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepNormSubtoDE == "+", "GelPrepReactIDs"]))
GelPrep.DERO.cont2 <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepNormDEcont2 == "+", "GelPrepReactIDs"]))  

GelPrep.BGGO.cont3 <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepNormSubtoDE == "+", "GelPrepGOID"]))
GelPrep.DEGO.cont3 <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepNormDEcont3 == "+", "GelPrepGOID"]))  
GelPrep.BGRO.cont3 <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepNormSubtoDE == "+", "GelPrepReactIDs"]))
GelPrep.DERO.cont3 <- SplitNClean(as.character(multExpanded1_withDE[multExpanded1_withDE$GelPrepNormDEcont3 == "+", "GelPrepReactIDs"]))  

#GelPrep normalized contrast enrichments
Enrich.GO.1 <- Enrich(GelPrep.BGGO.cont1, GelPrep.DEGO.cont1, ontology = "GO")
Enrich.GO.2 <- Enrich(GelPrep.BGGO.cont2, GelPrep.DEGO.cont2, ontology = "GO")
Enrich.GO.3 <- Enrich(GelPrep.BGGO.cont3, GelPrep.DEGO.cont3, ontology = "GO")
Enrich.RO.1 <- Enrich(GelPrep.BGRO.cont1, GelPrep.DERO.cont1, ontology = "Reactome")
Enrich.RO.2 <- Enrich(GelPrep.BGRO.cont2, GelPrep.DERO.cont2, ontology = "Reactome")
Enrich.RO.3 <- Enrich(GelPrep.BGRO.cont3, GelPrep.DERO.cont3, ontology = "Reactome")


  
  ####test enrichment using spearman rank correlation vs adjusted p values.
#   
# #   lets look at the top GO enrich term (GO:0046872) for gel prep. This give a frequency of 438/494 (
#   #for each phosphopeptide assign 0/1 depending on abs presence of this GO id
#   ont <- sapply(as.character(multExpanded1_withDE$GelPrepGOID), function(x){
#     tmp <- unlist(strsplit(x, ";"))
#     ifelse(any(tmp %in% "GO:0046872"), 1, 0)
#   })
#   
#   p.vals <- as.numeric(multExpanded1_withDE$GelPrepNormFAdjPval)
#   p.vals <- as.numeric(multExpanded1_withDE$GelPrepNormFPval)
# 
# 
# # Why doesn't the frequency output in the enrichment table match the number of ids isn the subtodiffphos subset for GO:0046872? it is identified 494 times...
# #   Ans: This is due to the GO term being duplicated within a single protein group.
# dups <- sapply(as.character(multExpanded1_withDE$GelPrepGOID), function(x){
#   tmp <- unlist(strsplit(x, ";"))
#   any(duplicated(tmp))
# })
# > table(dups)
# dups
# FALSE  TRUE 
# 12800  4974 
# 
# # same for reactome
# dups2 <- sapply(as.character(multExpanded1_withDE$GelPrepReactIDs), function(x){
#   tmp <- unlist(strsplit(x, ";"))
#   any(duplicated(tmp))
# })
# table(dups2)
# dups2
# FALSE  TRUE 
# 17719    55 
# 
# dups3 <- sapply(ReactIDs, function(x){
#   if(is.character(x)){
#   tmp <- unlist(strsplit(x, ";"))
#   any(duplicated(tmp))
#   }else{
#     NA}
# })
# table(dups3)
# 
# #so for a given phosphosite multiple annotation terms are added by eg the same isoforms within that protein group. 
# 
# # It would not bias enrichment if redundant terms from the same protein group are removed. In fact, including redundant terms contributed by isoforms only clouds the issue and may negatively impact p.values form exact test. Only the unique terms are kept PER SITE. 
# # However, each phosphosite contributes to the enrichment test (in the case of fisher's test) to not bias result in favor of mult phosphorylated proteins. 
# 
# # How does this affect the pvalue binary vector of annotation membership spearman correlation business?...Is there a multiple phos problem?
# # Not if each site's p value is included in the analysis? If the pvalues are random then the results will be negative. if the pvalues are correlated then it is biology driving that enrichment and that's what we want to identify.
# 
# # What about duplicated p values? p.adjust makes the issue worse. what is limma doing to create duplicate pvalues? 
# 
# sum(duplicated(multExpanded1_withDE[multExpanded1_withDE$ConfoundedFPval != "-", "ConfoundedFPval"] ))
# 201
# sum(duplicated(multExpanded1_withDE[multExpanded1_withDE$ConfoundedFAdjPval != "-", "ConfoundedFAdjPval"] ))
# 1342
# 
# # Is this related to multiplicity?
# 
# 
# 
# combined <- cbind(p.vals, ont)
#   combined <- combined[!is.na(combined[,1]),]
#   p.vals <- combined[,1]
#   ont <- combined[,2]
# 
#   plot(log10(p.vals),ont)
#   cor(p.vals, ont, method = "spearman")
#   test <- cor.test(p.vals, ont, method = "spearman", alternative = "two.sided")$p.value
# 
#   
#   
#   ## Sample code for computing correlation p-value.
#   n <- 25
#   x <- rnorm(n)
#   y <- -x + rnorm(n)
#   R <- cor(x,y)
#   p.value <- 2*pt(-abs(R * sqrt(n-2) / sqrt(1-R*R)),df=n-2)
#   cor.test(x,y,)$p.value
#   
  
  
  #return a list of enrichment DFs
  ###############################
  #final all dataframes in environment with names 'enrich'
  EnrichmentDFs <- setNames(lapply(ls(pattern="Enrich."), function(x) {if(class(get(x)) == "data.frame") get(x)}),ls(pattern="Enrich.")) #returns some NULL list elements
  EnrichmentDFs <- EnrichmentDFs[!sapply(EnrichmentDFs,is.null)]#note there is an 'is.null' function! and needed to use sappply for ligical indexing
  return(EnrichmentDFs)
  
  #note if there is no enrichment nothing will be returned!
}

  
  
  
  
  