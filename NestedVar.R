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

NestedVar <- function(ratios, noMissing = TRUE){
    
    ## Set up
    require(dplyr)
    require(reshape2)
    require(MCMCglmm)
    
    ##------ Prepare data ------#
    
    ratios <- as.matrix(ratios)
    
    # Remove rows with any missing values.
    if(noMissing) {
        noMissing <- na.omit(noMissing)
    }
    
    # Transform the data from one peptide per row
    # to one sample per row
    melted <- melt(ratios, measure.vars = names(ratios))
    
    # Make meta data matrix
    sampleLabels <- strsplit( colnames(melted$Var2), split = "_")
    metaData <- data.frame(individual = as.factor(sapply(sampleLabels, "[[", 1)),
                           biorep = as.factor(sapply(sampleLabels, "[[", 2)),
                           techrep = as.factor(sapply(sampleLabels, "[[", 3)) )
    metaData$label <- with(metaData, paste(individual, biorep, techrep, sep = "_"))
    
    # Append meta data matrix
    melted <- cbind(melted, metaData)
        
    ## Create a unique identifier for biological replicates
    melted$biorep <- as.factor(paste(melted$individual, melted$biorep,sep="_"))

        
    ##------ MCMCglmm for variance estimation ------#
    mcmcVarcomp <- lapply( levels(melted$Var1), function(id) {
        test <- melted[melted$Var1 %in% id,]
        test1 <- test[,3:6]
        fit_try <- tryCatch( MCMCglmm(value ~ 1, 
                                      random = ~ individual + individual:biorep,
                                      data = test1, verbose = FALSE),
                             condition = function(c) c)
        
        mcmc_varest <- c(c(summary(fit_try)$Gcovariances[,1], 
                           summary(fit_try)$Rcovariances[,1]) )
        
        if(inherits(fit_try, "condition")){
            varFoo <- rep(NA, 3)
            return(var_foo)
        }
        if(!inherits(fit_try, "condition")){
            var_foo <- mcmc_varest
            var_foo
        }
    })
    mcmcVarcomp <- do.call(rbind, mcmcVarcomp)
    rownames(mcmcVarcomp) <- levels(melted$Var1)
    colnames(mcmcVarcomp) <- c("individual","biorep","residual")
    

    ##------ Results ------#
    
    # Plot the variance component distributions
    colnames(mcmcVarcomp) <- c("individual","biorep","residual")
    boxplot(log10(mcmcVarcomp), ylab = "log10 variance component")
    summary(Varcomp)
    head(Varcomp)
    
    # Histograms of log10 variance 
    for (i in 1:ncol(mcmcVarcomp) ) {
        plot( density(log10(mcmcVarcomp[ ,i]) ), xlab = "log10 variance", 
              main = paste(colnames(mcmcVarcomp)[i], "variance") )
    }

    # Scatter plots of log10 variance components
    plot(log10(mcmcVarcomp[,1]),log10(mcmcVarcomp[,3]), 
         main = "log10 variance", xlab = colnames(mcmcVarcomp)[1], ylab = colnames(mcmcVarcomp)[3])
    plot(log10(mcmcVarcomp[,1]),log10(mcmcVarcomp[,2]), 
         main = "log10 variance", xlab = colnames(mcmcVarcomp)[1], ylab = colnames(mcmcVarcomp)[2])
    plot(log10(mcmcVarcomp[,2]),log10(mcmcVarcomp[,3]), 
         main = "log10 variance", xlab = colnames(mcmcVarcomp)[2], ylab = colnames(mcmcVarcomp)[3])
  
    
    # Here we'd like to identify phosphpeptides with little or no variability 
    # at the individual level and at the biological replicate level. To do so, 
    # we standardized the values of the variance components for each phosphopeptides 
    # with respect to its sum of variance components. The standardized variance 
    # components are the proportion of the total variation in each phosphopeptides
    # attributed to individuals, biological replicates, and technical replicates. 

    # Boxplots of the standardized VCs confirm our observations from the raw VC values. 
    # Proportion of variability attributed to biological replicates is the smallest, 
    # followed by technical replicates, with individaul samples contributing the largest 
    # portion of variabilty in expression levels. 
    par(mfrow = c(1,1))
    varprop <- mcmcVarcomp/rowSums(mcmcVarcomp)
    labs = c("individual","biorep","tech")
    boxplot((varprop), axes = F)
    axis(1, at = c(1, 2, 3), labels = labs, col = "white"); axis(2)


    # Heatmap representation of the standardized VCs. 
    require(gplots)
    require(RColorBrewer)
    colnames(varprop) = c("individual","bio","tech")
    heatmap.2(as.matrix(varprop),
              col=brewer.pal(9,"YlGnBu"),
              Colv=F,
              labRow="",
              trace="none",
              srtCol=45,  ,adjCol = c(1,1),
              margins = c(6,5),
              cexCol=1.5,
              key.xlab = "Standardized VC", key.ylab=NULL, key.title = "",
              )

    return(mcmcVarcomp)
}
  
  
  
  
  