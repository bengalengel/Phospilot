#' Estimate per peptide variance components 
#'
#' Per gene differential phosporylation analysis controlling for protein levels 
#' quantified using PhosPre. Each model is a mixed effect model with a fixed effect
#' of individual cell lines and a random effect of cell culture; each individual's
#' corresponding protein level for the peptide is used as a covariate.
#' 
#' @param ratios a peptide by sample matrix of H/L ratio
#' @param noMissing remove peptides with samples that contain any missing values.
#'          TRUE or FALSE.
#' @export
#' 
#' @examples
#' NestedVar()
#' dir <- "~/Dropbox/Brett"
#' rdadir <- file.path(dir, "rdas")
#' load(file.path(dir,"melted.RData"))
#' ratios <- melted

NestedVar <- function(ratios, noMissing = TRUE, includeProteinCovariate = FALSE) {
    
    ## Set up
    require(plyr)
    require(dplyr)
    require(reshape2)
    require(MCMCglmm)
    require(qdapRegex)
    
    ##------ Prepare data ------#
    ratios <- as.matrix(ratios)
    
    # Remove rows with any missing values.
    if(noMissing) {
      ratios <- na.omit(ratios)
    }
    
    # Compute a variable indicating which columns
    # are not phospho samples (i.e, protein)
    variable_not_phospho <- grep("_", colnames(ratios), invert = TRUE)
    
    ## If not including protein as a covariate,
    ## check if protein values are also in the data.frame
    if(!includeProteinCovariate) {
      if(length( variable_not_phospho ) > 0) {
        ratios <- ratios[ , -variable_not_phospho]
      } 
    }
    
    ## If including protein as a covariate,
    ## make sure protein values are in the data.frame
    if(includeProteinCovariate) {
      if( length( variable_not_phospho ) != 3 ) {
        stop("Missing covariate vectors", call. = TRUE)
      } 
    }
    
    # Transform the data from one peptide per row
    # to one sample per row
    melted <- melt(ratios, measure.vars = names(ratios))
    
    # Make meta data matrix
    #     sampleLabels <- strsplit(as.character(unique(melted$Var2)), split = "_")
    #     metaData <- data.frame(individual = as.factor(sapply(sampleLabels, "[[", 1)),
    #                            biorep = as.factor(sapply(sampleLabels, "[[", 2)),
    #                            techrep = as.factor(sapply(sampleLabels, "[[", 3)) )
    #     metaData$label <- with(metaData, paste(individual, biorep, techrep, sep = "_"))
    
    
    #identify individual name and add it to the table
    matches <- gregexpr("[0-9]{5}", melted$Var2, perl=T)
    individual <- regmatches(melted$Var2,matches)
    individual <- as.character(individual)
    individual <- as.factor(individual)
    melted$individual <- individual
    
    #identify the biological replicate
    biorep <- rm_between(melted$Var2, "_", "_", extract=TRUE)
    biorep <- as.character(biorep)
    biorep <- as.factor(biorep)
    melted$biorep <- biorep
    
    #identify the technical replicate
    matches <- gregexpr("[0-9]$", melted$Var2, perl=T)
    techrep <- regmatches(melted$Var2,matches)
    techrep <- as.character(techrep)
    techrep <- as.factor(techrep)
    melted$techrep <- techrep
    
    # Append meta data matrix
    #     melted <- cbind(melted, metaData)
    
    ## Create a unique identifier for biological replicates
    melted$biorep_unique <- as.factor(paste(melted$individual, melted$biorep,sep="_"))
    
    
    ##------ MCMCglmm for variance estimation ------#
    mcmcVarcomp <- lapply( levels(melted$Var1), function(id) {
      
      test <- melted[melted$Var1 %in% id,]
      #      test1 <- test[ ,3:7]
      
      if (includeProteinCovariate == FALSE) { 
        stopifnot()  
        fit_try <- tryCatch( MCMCglmm(value ~ 1, 
                                      random = ~ individual + individual:biorep_unique,
                                      data = test, verbose = FALSE),
                             condition = function(c) c)
      } 
      
      if (includeProteinCovariate == TRUE) {
        protein_values <- test[ grep("_", test$Var2, invert = TRUE), ]
        test$protein <- test$value[ match(test$individual, protein_values$individual) ]
        
        # Add uncertains to the protein vector 
        test$protein <- as.numeric( test$protein + runif(NROW(test), 0, 1e-10) )
        
        fit_try <- tryCatch( MCMCglmm(value ~ protein, 
                                      random = ~ individual + individual:biorep_unique,
                                      data = test, verbose = FALSE),
                             condition = function(c) c)
      }
      
      if(inherits(fit_try, "condition")){
        var_foo <- rep(NA, 3) 
        return(var_foo)
      }
      if(!inherits(fit_try, "condition")){
        mcmc_varest <- c(summary(fit_try)$Gcovariances[,1], 
                         summary(fit_try)$Rcovariances[,1])
        mcmc_varest            
      }
      
    })
    mcmcVarcomp <- do.call(rbind, mcmcVarcomp)
    rownames(mcmcVarcomp) <- levels(melted$Var1)
    colnames(mcmcVarcomp) <- c("individual","biorep","residual")
    mcmcVarcomp
}


  
  
  
  