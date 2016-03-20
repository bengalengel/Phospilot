#' Differential phosphorylation analysis
#'
#' Per gene differential phosporylation analysis controlling for protein levels 
#' quantified using PhosPre. Each model is a mixed effect model with a fixed effect
#' of individual cell lines and a random effect of cell culture; each individual's
#' corresponding protein level for the peptide is used as a covariate.
#' 
#' @param phosdata Phosphorylated peptides quantifications. Each row uniquely identifies
#'	           a peptide with an individual, cell culture, and technical replicate.
#' @param PhosPrep Protein quantification from PhosPrep workup. Each row uniquely 
#'		   identifies a protein with an individual, cell culture, and
#'		   technical replicate.
#' @param GelPrep Protein quantification from GelPrep workup. Each row uniquely 
#'		  identifies a protein with an individual, cell culture, and
#'		  technical replicate.
#'  
#' @param multExpanded1 
#' @export
#' 
#' @examples
#' 
#' setwd("~/Dropbox/Brett/diffPhose-mean-variance-20160313")
#' 
#' # pQTL PCA regressed data.frame from the gel based workup
#' GelPrep <- readRDS("GelPrep.rds")
#' 
#' # normalized batch effect corrected phosphopeptide H/L ratios
#' phosdata <- readRDS("phosdata.rds")
#'
#' # protein estimates derived from the phospho workup normalized
#' # and BE correct 
#' PhosPrep <- readRDS("phosprepcombatbio.rds")
#' 
#' # Contains the phosphopeptide intensities
#' MultExpanded1 <- readRDS("MultExpanded1.rds")
#'
#' DiffPhos()


DiffPhos <- function(phosdata, PhosPrep, GelPrep, multExpanded1,
                     weight_by_intensities = TRUE){
  # This function accepts phospho and protein matrices runs diffphos 
  # analysis using limma. Columns are appended
  # to the multexpanded1 file according to the presence of DE or not. 
  # As of now combat correction is being performed upstream and protein levels are being fitted as a covariate as opposed to normalizing phosdata and fitting withiout the covariate. It is unclear how normalization would affect fitting biorep as defacto RE.
  
  require(limma)
  require(sva)
  require(statmod)
  require(qdapRegex)
  require(plyr)
  require(reshape2)
  require(matrixStats)
  
  ### If intending to weight peptides H/L ratios with intensity measurements,
  ### Then, use the following function:
  #' @param cov_matrix Matrix of the same dimension as  counts. Contain measurements that are 
  #'                   used to compute weights for the dependent measurment counts. 
        voom_general  <- function (counts, cov_matrix,
                                   span = 0.5, plot = FALSE) 
        {
                design <- matrix(1, ncol(counts), 1)
                rownames(design) <- colnames(counts)
                colnames(design) <- "GrandMean"
                fit <- lmFit(counts, design)
                sx <- log10(rowMeans(cov_matrix, na.rm = TRUE) + 1)
                sy <- sqrt(fit$sigma)
                df <- cbind(sx,sy)
                df <- df[rowSums(is.na(df)) == 0, ]
                l <- lowess(df, f = .5)
                if (plot) {
                        plot(df, xlab = "log10 meausrement",
                             ylab = "Sqrt( standard deviation )",
                             pch = 16, cex = 0.25)
                        title("voom: Mean-variance trend")
                        lines(l, col = "red")
                }
                f <- approxfun(l, rule = 2)
                #fitted.values <- fit$coef %*% t(fit$design)
                
                w <- 1/f(as.matrix(cov_matrix))^4
                dim(w) <- dim(cov_matrix)
                return(w)
        }
        

        
  ### LINEAR MIXED MODEL USING CULTURE AS RANDOM EFFECT ####################  

  ### PREPARE DATA
    
  # Make meta data matrix
  sampleLabels <- strsplit( colnames(phosdata), split = "_")
  metaData <- data.frame(individual = as.factor(sapply(sampleLabels, "[[", 1)),
                          biorep = as.factor(sapply(sampleLabels, "[[", 2)),
                          techrep = as.factor(sapply(sampleLabels, "[[", 3)) )
    
  # Make the design matrix for the individual effect
  designMatrix <- model.matrix(~ 0 + individual, data = metaData)
  
  # Make the design matrix for the random effect of 
  # culture replicate
  block <- as.factor(metaData$biorep)
  
  # Compute correlations between culture replciates
  dupcor <- duplicateCorrelation(phosdata, designMatrix, block = block)
  allCorrelations <- tanh(dupcor$atanh.correlations)
  boxplot(allCorrelations, main = "Per peptide biorep correlation")
  
  # Fit a linear mixed model using biorep as a 
  # random effect where the magnitude of the random effect
  # is estimated to be the same across peptides
  if (weight_by_intensities == FALSE) {
  fit <- lmFit(phosdata, 
               designMatrix, 
               block = block, 
               correlation = dupcor$consensus)
  } 
  
  if (weight_by_intensities == TRUE) {
  intensity_matrix <- data.frame(
      id = MultExpanded1$idmult,
      MultExpanded1[ ,grep("Intensity.1", colnames(MultExpanded1))]) 
  
  ii_phosdata <- match(rownames(phosdata), intensity_matrix$id)
  intensity_phosdata <- intensity_matrix[ii_phosdata, ]
  weights_phosdata <- voom_general(counts = phosdata,
                                   cov_matrix = intensity_phosdata[,c(2:13)],
                                   plot = TRUE)
  fit <- lmFit(phosdata, 
               designMatrix, 
               block = block, 
               correlation = dupcor$consensus,
               weights = weights_phosdata)
  } 
  
  
  # Construct the contrast matrix
  contrastMatrix <- makeContrasts(individualHL18862 - individualHL18486, individualHL19160 - individualHL18862, 
                                   individualHL19160 - individualHL18486, levels = designMatrix)

  fit2 <- contrasts.fit(fit, contrastMatrix)
  
  # eBayes
  ConfoundedFit <- eBayes(fit2)
  

  
    
  ### PHOSPREP AS COVARIATE ################################################
  ### DIFFERENTIAL PHOSPHORYLATION ANALYSIS
  ### using culture replicate as a random effect

  ### PREPARE DATA 
  PhosProt <- merge(phosdata, PhosPrep, by = "row.names", 
                    suffixes = c("_peptide", "_PhosPrep") ) #1308 observations
  rownames(PhosProt) <- PhosProt$Row.names
  PhosProt <- PhosProt[ , -1]
  PhosProt <- as.matrix(PhosProt)
  
  # Make meta data matrix
  sampleLabels <- strsplit( colnames(PhosProt), split = "_", fixed = FALSE)
  metaData <- data.frame(individual = as.factor(sapply(sampleLabels, "[[", 1)),
                         biorep = as.factor(sapply(sampleLabels, "[[", 2)),
                         techrep = as.factor(sapply(sampleLabels, "[[", 3)),
                         dataType = as.factor(sapply(sampleLabels, "[[", 4)) )
  metaData$label <- with(metaData, paste(individual, biorep, techrep, sep = "_"))

  # Imput missing values in the PhosPrep part of the data
  metaData_phos <- metaData[ metaData$dataType == "PhosPrep", ]
  PhosProt_phos <- PhosProt[ , metaData$dataType == "PhosPrep"]
  PhosProt_phos_imput <- 
        lapply( 1:nrow(PhosProt_phos), function(per_peptide) {
                phosMat <- PhosProt_phos[per_peptide, ]
                
                ind_vals <- lapply(unique(metaData_phos$individual), function(per_individual) {
                
                                    bio_vals <- lapply(unique(metaData_phos$biorep), function(per_biorep) {
                                                        iiReplace <- metaData_phos$individual == per_individual &
                                                                        metaData_phos$biorep == per_biorep 
                                                        vals <- phosMat[iiReplace]
                                                        if( sum(is.na(vals) > 0) ) {
                                                            vals[ is.na(vals)] <- vals[ !is.na(vals) ]
                                                            vals
                                                        } else {
                                                            vals    
                                                        }
                                                        }) 
                                    bio_vals <- do.call(c, bio_vals)
                                    bio_vals
                    }) 
                ind_vals <- do.call(c, ind_vals)
                ind_vals
    })
  PhosProt_phos_imput <- do.call(rbind, PhosProt_phos_imput)
  
  ### FIT A LINEAR MODEL FOR ONE PEPTIDE AT A TIME
  ### Use the consensus correlation computed earlier in
  ### the linear mixed model without protein level as
  ### a covariate

  ####################Needed Functions - Temporary!!! ######################
  #' Fitting mixed models of one random effect
  #' 
  #' Adapted from statmod package's mixedModel2Fit function. This function
  #' includes the additional feature of fitting a different design matrix
  #' for each gene (row) of the input data matrix.
  #'
  #' 
  #' @param design
  #' @param block
  #' @param yy
  #'
  #' @keywords Humanzee
  #' 
  #' @export
  #' 
  #' @examples
  #' mixedModel2Fit_multiple_design()
  #' data_matrix <- pheno_data
  #' design <- model.matrix(~ 1+ as.factor(individual) + protein, 
  #'                        data = pheno_data)
  #' block <- as.factor(pheno_data$biorep)
  #' yy <- pheno_data$phos
  
  mixedModel2Fit_multiple_design <- 
    function(design, block, yy) {
      o <- is.finite(yy)
      A <- factor(block[o])
      nobs <- sum(o)
      nblocks <- length(levels(A))
      nbeta <- NCOL(design)
      nafun <- function(e) NA
      
      #        if (nobs > (nbeta + 2) && nblocks > 1 && nblocks < nobs - 1) {
      yy <- yy[o]
      X <- design[o, , drop = FALSE]
      Z <- model.matrix(~0 + A)
      s <- tryCatch(statmod::mixedModel2Fit(yy, X, 
                                            Z, only.varcomp = TRUE, maxit = 20)$varcomp, 
                    error = nafun)
      if (!is.na(s[1])) 
        rho <- s[2]/sum(s)
      return(rho)
      #        }
    }
  
  
  #' limma gls.series for varying covariates across genes
  #' 
  #' This function tests for divergence between two species
  #' in one molecular phenotype at a time. 
  #' 
  #' @param M
  #'
  #' @keywords Humanzee
  #' 
  #' @export
  #'
  #' @examples 
  #' M = phos_data
  #' cov_matrix = protein_pheno_data
  #' individual <- as.numeric(str_extract(colnames(phos_data), "[0-9]+"))
  #' design <- model.matrix(~ 1 + as.factor(individual))
  #' block <- block
  #' correlation <- mrho
  #' ndups = 1; weights = NULL; spacing = 1
  
  gls.series_multiple_designs <- function (M, 
                                           cov_matrix = NULL, design = NULL, ndups = 1, 
                                           spacing = 1, block = NULL, 
                                           correlation = NULL, weights = NULL, ...) {
    
    M <- as.matrix(M)
    
    # narrays: number of samples
    narrays <- ncol(M)
    
    if (is.null(design)) 
      design <- matrix(1, narrays, 1)
    design <- as.matrix(design)
    if (nrow(design) != narrays) 
      stop("Number of rows of design matrix does not match number of arrays")
    if (is.null(correlation)) 
      correlation <- duplicateCorrelation(M, design = design, 
                                          ndups = ndups, spacing = spacing, block = block, 
                                          weights = weights, ...)$consensus.correlation
    if (!is.null(weights)) {
      weights[is.na(weights)] <- 0
      weights <- asMatrixWeights(weights, dim(M))
      M[weights < 1e-15] <- NA
      weights[weights < 1e-15] <- NA
    }
    
    if (!is.null(cov_matrix)) {
      nbeta <- ncol(design) + 1        
      coef.names <- c(colnames(design), "cov")
    } else {
      nbeta <- ncol(design)
      coef.names <- colnames(design)
    }
    
    if (is.null(block)) {
      if (ndups < 2) {
        warning("No duplicates: correlation between duplicates set to zero")
        ndups <- 1
        correlation <- 0
      }
      if (is.null(spacing)) 
        spacing <- 1
      cormatrix <- diag(rep(correlation, len = narrays), nrow = narrays, 
                        ncol = narrays) %x% array(1, c(ndups, ndups))
      M <- unwrapdups(M, ndups = ndups, spacing = spacing)
      if (!is.null(weights)) 
        weights <- unwrapdups(weights, ndups = ndups, spacing = spacing)
      design <- design %x% rep(1, ndups)
      colnames(design) <- coef.names 
    } else {
      if (ndups > 1) {
        stop("Cannot specify ndups>2 and non-null block argument")
      }
      else {
        ndups <- spacing <- 1
      }
      block <- as.vector(block)
      if (length(block) != narrays) 
        stop("Length of block does not match number of arrays")
      ub <- unique(block)
      nblocks <- length(ub)
      Z <- matrix(block, narrays, nblocks) == matrix(ub, narrays, 
                                                     nblocks, byrow = TRUE)
      cormatrix <- Z %*% (correlation * t(Z))
    }
    diag(cormatrix) <- 1
    ngenes <- nrow(M)
    stdev.unscaled <- matrix(NA, ngenes, nbeta, 
                             dimnames = list(rownames(M), 
                                             coef.names))
    # NoProbeWts <- all(is.finite(M)) && (is.null(weights) || !is.null(attr(weights, 
    #                                                                      "arrayweights")))
    #     if (NoProbeWts) {
    #         V <- cormatrix
    #         if (!is.null(weights)) {
    #             wrs <- 1/sqrt(weights[1, ])
    #             V <- wrs * t(wrs * t(V))
    #         }
    #         cholV <- chol(V)
    #         y <- backsolve(cholV, t(M), transpose = TRUE)
    #         dimnames(y) <- rev(dimnames(M))
    #         X <- backsolve(cholV, design, transpose = TRUE)
    #         dimnames(X) <- dimnames(design)
    #         fit <- lm.fit(X, y)
    #         if (fit$df.residual > 0) {
    #             if (is.matrix(fit$effects)) 
    #                 fit$sigma <- sqrt(colMeans(fit$effects[-(1:fit$rank), 
    #                                                        , drop = FALSE]^2))
    #             else fit$sigma <- sqrt(mean(fit$effects[-(1:fit$rank)]^2))
    #         }
    #         else fit$sigma <- rep(NA, ngenes)
    #         fit$fitted.values <- fit$residuals <- fit$effects <- NULL
    #         fit$coefficients <- t(fit$coefficients)
    #         fit$cov.coefficients <- chol2inv(fit$qr$qr, size = fit$qr$rank)
    #         est <- fit$qr$pivot[1:fit$qr$rank]
    #         dimnames(fit$cov.coefficients) <- list(coef.names[est], 
    #                                                coef.names[est])
    #         stdev.unscaled[, est] <- matrix(sqrt(diag(fit$cov.coefficients)), 
    #                                         ngenes, fit$qr$rank, byrow = TRUE)
    #         fit$stdev.unscaled <- stdev.unscaled
    #         fit$df.residual <- rep.int(fit$df.residual, ngenes)
    #         dimnames(fit$stdev.unscaled) <- dimnames(fit$stdev.unscaled) <- dimnames(fit$coefficients)
    #         fit$pivot <- fit$qr$pivot
    #         fit$ndups <- ndups
    #         fit$spacing <- spacing
    #         fit$block <- block
    #         fit$correlation <- correlation
    #         return(fit)
    #     }
    
    beta <- stdev.unscaled
    sigma <- rep(NA, ngenes)
    df.residual <- rep(0, ngenes)
    for (i in 1:ngenes) {
      design_gene <- cbind(design, unlist(cov_matrix[i,]))
      y <- drop(M[i, ])
      o <- is.finite(y)
      y <- y[o]
      n <- length(y)
      if (n > 0) {
        X <- design_gene[o, , drop = FALSE]
        V <- cormatrix[o, o]
        if (!is.null(weights)) {
          wrs <- 1/sqrt(drop(weights[i, o]))
          V <- wrs * t(wrs * t(V))
        }
        cholV <- chol(V)
        y <- backsolve(cholV, y, transpose = TRUE)
        if (all(X == 0)) {
          df.residual[i] <- n
          sigma[i] <- sqrt(array(1/n, c(1, n)) %*% y^2)
        }
        else {
          X <- backsolve(cholV, X, transpose = TRUE)
          out <- lm.fit(X, y)
          est <- !is.na(out$coefficients)
          beta[i, ] <- out$coefficients
          stdev.unscaled[i, est] <- sqrt(diag(chol2inv(out$qr$qr, 
                                                       size = out$rank)))
          df.residual[i] <- out$df.residual
          if (df.residual[i] > 0) 
            sigma[i] <- sqrt(array(1/out$df.residual, c(1, 
                                                        n)) %*% out$residuals^2)
        }
      }
    }
    cholV <- chol(cormatrix)
    QR <- qr(backsolve(cholV, design_gene, transpose = TRUE))
    cov.coef <- chol2inv(QR$qr, size = QR$rank)
    est <- QR$pivot[1:QR$rank]
    dimnames(cov.coef) <- list(coef.names[est], coef.names[est])
    
    list(coefficients = beta, stdev.unscaled = stdev.unscaled, 
         sigma = sigma, df.residual = df.residual, ndups = ndups, 
         spacing = spacing, block = block, correlation = correlation, 
         cov.coefficients = cov.coef, pivot = QR$pivot, rank = QR$rank,
         design = design)
  }
  
  
  #######################################################
  
  
  
  # Compute correlation due to the random effect of
  # biological replicate
  corrsPhosProt <- sapply(1:nrow(PhosProt), function(per_peptide) {
      
      # Make a data matrix for peptide i information
      per_phenoData <- cbind(metaData[ metaData$dataType == "peptide", ],
                             peptide = PhosProt[per_peptide, metaData$dataType == "peptide"],
                             phosprep = PhosProt_phos_imput[per_peptide, ] )
      designMatrix <- model.matrix(~ 0 + individual + phosprep, 
                                    data = per_phenoData)
      block <- per_phenoData$biorep
      yy <- per_phenoData$peptide
      
      mixedModel2Fit_multiple_design(design = designMatrix, block = block, yy = yy)
  })
  corrsPhosProt_mean <- mean(corrsPhosProt)
  
  # Per gene limma mixed model 
  ii_peptide <- metaData$dataType == "peptide"
  metaData_peptide <- metaData[ ii_peptide, ]
  designMatrix <- model.matrix(~ 0 + individual, data = metaData_peptide )
  
  if (weight_by_intensities == FALSE) {
  glsRes <- gls.series_multiple_designs(M = PhosProt[ ,ii_peptide ], 
                                        cov_matrix = PhosProt_phos_imput, 
                                        design = designMatrix, 
                                        ndups = 1, 
                                        spacing = 1, 
                                        block = metaData[ ii_peptide, ]$biorep, 
                                        correlation = corrsPhosProt_mean)
  }
  if (weight_by_intensities == TRUE) {
  intensity_matrix <- data.frame(
          id = MultExpanded1$idmult,
          MultExpanded1[ ,grep("Intensity.1", colnames(MultExpanded1))]) 
  
  ii_phosdata <- match(rownames(PhosProt[,ii_peptide]), intensity_matrix$id)
  intensity_phosdata <- intensity_matrix[ii_phosdata, ]
  weights_phosdata <- voom_general(counts = PhosProt[,ii_peptide],
                                   cov_matrix = intensity_phosdata[,c(2:13)],
                                   plot = TRUE)
          
  glsRes <- gls.series_multiple_designs(M = PhosProt[ ,ii_peptide ], 
                                        cov_matrix = PhosProt_phos_imput, 
                                        design = designMatrix, 
                                        ndups = 1, 
                                        spacing = 1, 
                                        block = metaData[ ii_peptide, ]$biorep, 
                                        correlation = corrsPhosProt_mean,
                                        weights = weights_phosdata)
  }
  
  # Contrast test
  designMatrixCov <- model.matrix(~ 0 + individual + pseudoCov, 
                                  data = data.frame(metaData_peptide, pseudoCov = runif(12) ) )
  contrastMatrix <- makeContrasts(individualHL18862 - individualHL18486, individualHL19160 - individualHL18862, 
                                   individualHL19160 - individualHL18486, 
                                   levels = designMatrixCov)
  rownames(contrastMatrix) <- colnames(glsRes$coefficients)
  
  fit2 <- contrasts.fit(glsRes, contrastMatrix)
  PhosPrepCovFit <- eBayes(fit2)
  

  
  ### GELPREP AS COVARIATE ####################################
  ### DIFFERENTIAL PHOSPHORYLATION ANALYSIS
  ### using culture replicate as a random effect
  
  ### PREPARE DATA 
  colnames(GelPrep) <- c("HL18862", "HL18486", "HL19160")
  PhosProtGel <- merge(phosdata, GelPrep, by = "row.names", 
                    suffixes = c("_peptide", "_GelPrep") ) #3257 observations
  rownames(PhosProtGel) <- PhosProtGel$Row.names
  PhosProtGel <- PhosProtGel[ , -1]
  PhosProtGel <- as.matrix(PhosProtGel)
  

  # Make meta data matrix
  sampleLabels <- strsplit( colnames(PhosProtGel)[1:12], split = "_", fixed = FALSE)
  metaData <- data.frame(individual = as.factor(sapply(sampleLabels, "[[", 1)),
                         biorep = as.factor(sapply(sampleLabels, "[[", 2)),
                         techrep = as.factor(sapply(sampleLabels, "[[", 3)) )
  metaData <- rbind(metaData,
                    data.frame(individual = colnames(PhosProtGel)[13:15],
                               biorep = "NA",
                               techrep = "NA") )
  metaData$dataType <- c(rep("peptide", 12), rep("GelPrep", 3) )
  metaData$label <- with(metaData, paste(individual, biorep, techrep, sep = "_"))
  

  ### FIT A LINEAR MODEL FOR ONE PEPTIDE AT A TIME

  # First, add random perturbation to the covariate matrix
  ii_phosprotgel <- metaData$dataType == "GelPrep"
  randomPhosProtGel <- PhosProtGel[ , ii_phosprotgel] + 
                        matrix( runif(NROW(PhosProtGel[ , ii_phosprotgel])*NCOL(PhosProtGel[ , ii_phosprotgel]), 
                                0, 10^(-6)), 
                                nrow(PhosProtGel[ , ii_phosprotgel]) )
  ii_peptide <- metaData$dataType == "peptide"
  metaData_peptide <- metaData[ ii_peptide, ]
  order_phosgel <- match(metaData_peptide$individual, colnames(randomPhosProtGel) )
  randomPhosProtGel <- randomPhosProtGel[, order_phosgel]
  
  # Compute correlation due to the random effect of
  # biological replicate
  corrsPhosProtGel <- sapply(1:nrow(PhosProtGel), function(per_peptide) {
      
      # Make a data matrix for peptide i information
      per_phenoData <- cbind(metaData[ metaData$dataType == "peptide", ],
                             peptide = PhosProtGel[per_peptide, metaData$dataType == "peptide"],
                             phosgel = randomPhosProtGel[per_peptide, ])

      designMatrix <- model.matrix(~ 0 + individual + phosgel, 
                                   data = per_phenoData)
      block <- per_phenoData$biorep
      yy <- per_phenoData$peptide
      
      mixedModel2Fit_multiple_design(design = designMatrix, block = block, yy = yy)
      })
  corrsPhosProtGel_mean <- mean(corrsPhosProtGel)
  
  # Per gene limma mixed model 
  designMatrix <- model.matrix(~ 0 + individual, data = metaData_peptide )
  if (weight_by_intensities == FALSE) {
  glsRes <- gls.series_multiple_designs(M = PhosProtGel[ ,ii_peptide ], 
                                        cov_matrix = randomPhosProtGel, 
                                        design = designMatrix, 
                                        ndups = 1, 
                                        spacing = 1, 
                                        block = metaData[ ii_peptide, ]$biorep, 
                                        correlation = corrsPhosProtGel_mean)
  }
  if (weight_by_intensities == TRUE) {
  intensity_matrix <- data.frame(
          id = MultExpanded1$idmult,
          MultExpanded1[ ,grep("Intensity.1", colnames(MultExpanded1))]) 
  
  ii_phosdata <- match(rownames(PhosProtGel), intensity_matrix$id)
  intensity_phosdata <- intensity_matrix[ii_phosdata, ]
  weights_phosdata <- voom_general(counts = PhosProtGel[ , ii_peptide],
                                   cov_matrix = intensity_phosdata[,c(2:13)],
                                   plot = TRUE)
          
  glsRes <- gls.series_multiple_designs(M = PhosProtGel[ ,ii_peptide ], 
                                        cov_matrix = randomPhosProtGel, 
                                        design = designMatrix, 
                                        ndups = 1, 
                                        spacing = 1, 
                                        block = metaData[ ii_peptide, ]$biorep, 
                                        correlation = corrsPhosProtGel_mean,
                                        weights = weights_phosdata)
  }
  

  # Contrast test
  designMatrixCov <- model.matrix(~ 0 + individual + pseudoCov, 
                                  data = data.frame(metaData_peptide, pseudoCov = runif(12) ) )
  contrastMatrix <- makeContrasts(individualHL18862 - individualHL18486, individualHL19160 - individualHL18862, 
                                  individualHL19160 - individualHL18486, 
                                  levels = designMatrixCov)
  rownames(contrastMatrix) <- colnames(glsRes$coefficients)
  
  fit2 <- contrasts.fit(glsRes, contrastMatrix)
  GelPrepCovFit <- eBayes(fit2)



  ############# Collect the sig hits and annotate ME df --------------

ProcessFit <- function(fit2, header, FitData, multExpanded1){
    # This function accepts the ebays moderated fit data from a given processing choice, 
    # returns charts and annotates the ME dataframe.
    # It requires the ebays modifed contrast fits, a header string (such as "confounded") 
    # and the matrix passed to limma.
    
    # For convert fit2 to a limma object
    fit2 <- new("MArrayLM", fit2)

    # Extract contrasts tests
    contrastFits <- colnames(contrastAsCoef(fit2))
    
    topLists <- lapply(1:3, function(i) {
        topTable(fit2, coef = i, adjust = "BH", n=Inf)#sorts by adjusted p up to the threshold of .
    })
    
    for (i in 1:3) {
        hist(topLists[[i]]$P.Value, nc=40, xlab="P values", main = contrastFits[i])
    }

    for (i in 1:3) {
        plot(topLists[[i]]$logFC, -log10(topLists[[i]]$P.Value), 
             xlab = contrastFits[i], 
             pch = 20, ylab = "-log10(P)",xlim = c(-5, 5))
        #sites with sig difference in comparison 1
        index <- which(topLists[[i]]$adj.P.Val < .05)
        points(topLists[[i]]$logFC[index], -log10(topLists[[i]]$P.Value)[index], 
               col="red3", pch = 20)
    }

    ## Summarize. Results here are shown at a 5% FDR
    results <- decideTests(fit2, adjust.method = "BH", method = "separate")#results is a 'TestResults' matrix
    #separate compares each sample individually and is the default approach
    colnames(results) <- c("18862 - 18486", "19160 - 18862", "19160 - 18486")
    summary(results)
    vennDiagram(results, cex=c(1.2,1,0.7)) #good DE across conditions
    
    #make this venn diagram prettier.
    require(VennDiagram)
    result.matrix <- as.data.frame(results)
    id1 <- result.matrix[abs(result.matrix[[1]]) == 1 & abs(result.matrix[[2]]) == 0 & abs(result.matrix[[3]]) == 0,]
    id2 <- result.matrix[abs(result.matrix[[1]]) == 0 & abs(result.matrix[[2]]) == 1 & abs(result.matrix[[3]]) == 0,]
    id3 <- result.matrix[abs(result.matrix[[1]]) == 0 & abs(result.matrix[[2]]) == 0 & abs(result.matrix[[3]]) == 1,]
    id12 <- result.matrix[abs(result.matrix[[1]]) == 1 & abs(result.matrix[[2]]) == 1 & abs(result.matrix[[3]]) == 0,]
    id13 <- result.matrix[abs(result.matrix[[1]]) == 1 & abs(result.matrix[[2]]) == 0 & abs(result.matrix[[3]]) == 1,]
    id23 <- result.matrix[abs(result.matrix[[1]]) == 0 & abs(result.matrix[[2]]) == 1 & abs(result.matrix[[3]]) == 1,]
    id123 <- result.matrix[abs(result.matrix[[1]]) == 1 & abs(result.matrix[[2]]) == 1 & abs(result.matrix[[3]]) == 1,]
    
    pdf(paste(header, "TripleVenn.pdf", sep =""), 9, 8)
    plot.new()
    venn.plot <- draw.triple.venn(
      area1 = sum(abs(result.matrix[ , 1])),
      area2 = sum(abs(result.matrix[ , 2])),
      area3 = sum(abs(result.matrix[ , 3])),
      n12 = nrow(id12)+nrow(id123),
      n23 = nrow(id23)+nrow(id123),
      n13 = nrow(id13)+nrow(id123),
      n123 = nrow(id123),
      category = colnames(results),
      fill = c("red", "gray", "blue"),
      lty = "blank",
      cex = 2,
      cat.cex = 2,
      cat.col = rep("black",3), 
      margin = .1,
    )
    dev.off()
    
    vennDiagram(results, cex=c(1.2,1,0.7), include = "up") #good DE across conditions
    vennDiagram(results, cex=c(1.2,1,0.7), include = "down") #good DE across conditions
    vennDiagram(results, cex=c(1.2,1,0.7), include = c("up", "down")) #good DE across conditions


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

    #F stats by contrast type
    
    # #Sorting F Values
    # test <- topTable(fit2, coef = c(1,2,3), adjust = "BH", n=Inf, sort.by="F", p=.05)#equivalent to below
    # test2 <- topTableF(fit2, adjust = "BH", n=Inf, sort.by="F", p=.05)
    
    Fvals <- topTableF(fit2, adjust = "BH", n=Inf, sort.by="F")#all F values
    sigFvals <- topTableF(fit2, adjust = "BH", n=Inf, sort.by="F", p=.05)#gives 1355 compared to 1549 DE total for separate comparisons
    
    #below gives adjusted pvalue
    FDE1 <- Fvals[match(row.names(DE1), row.names(Fvals), nomatch = F), 4:6]
    FDE2 <- Fvals[match(row.names(DE2), row.names(Fvals), nomatch = F), 4:6]
    FDE3 <- Fvals[match(row.names(DE3), row.names(Fvals), nomatch = F), 4:6]
    
    #boxplot(FDE1$F,FDE2$F,FDE3$F)
    boxplot(log10(FDE1$F),log10(FDE2$F),log10(FDE3$F))
    summary(FDE1$F)
    summary(FDE2$F)
    summary(FDE3$F)
    
    plot(density(log10(FDE1$F)),xlim = c(0,3), main = "Sig F Stats Cut by Number of DiffPhos Contrasts")
    lines(density(log10(FDE2$F)), col = 2)
    lines(density(log10(FDE3$F)), col = 3)
    
    #add annotation to multexpanded DF
    multExpanded1$SubtoDE = ifelse(multExpanded1$idmult %in% row.names(FitData),"+","-")
    
    #add F test values to the table
    multExpanded1$globalFsig = ifelse(multExpanded1$idmult %in% row.names(sigFvals),"+","-")
    
    multExpanded1$globalFstat  <-  sapply(as.character(multExpanded1$idmult), function(x){
      ifelse(x %in% row.names(Fvals), Fvals[which(row.names(Fvals) == x), "F"], "-")
    })
    multExpanded1$globalFPval  <-  sapply(as.character(multExpanded1$idmult), function(x){
      ifelse(x %in% row.names(Fvals), Fvals[which(row.names(Fvals) == x), "P.Value"], "-")
    })
    multExpanded1$globalFAdjPval  <- sapply(as.character(multExpanded1$idmult), function(x){
      ifelse(x %in% row.names(Fvals), Fvals[which(row.names(Fvals) == x), "adj.P.Val"], "-")
    })
    
    #add DE to table. 5% FDR
    sig1 <- topTable(fit2, coef = 1, adjust = "BH", n = Inf, sort.by = "P", p = .05)
    sig2 <- topTable(fit2, coef = 2, adjust = "BH", n = Inf, sort.by = "P", p = .05)
    sig3 <- topTable(fit2, coef = 3, adjust = "BH", n = Inf, sort.by = "P", p = .05)
    
    c1up  <- sig1[sig1$logFC > 0,]
    c1down <- sig1[sig1$logFC < 0,]
    c2up <- sig2[sig2$logFC > 0,]
    c2down <- sig2[sig2$logFC < 0,]
    c3up <- sig3[sig3$logFC > 0,]
    c3down <- sig3[sig3$logFC < 0,]
    
    multExpanded1$DEcont1 = ifelse(multExpanded1$idmult %in% row.names(sig1),"+","-")
    multExpanded1$DEcont2 = ifelse(multExpanded1$idmult %in% row.names(sig2),"+","-")
    multExpanded1$DEcont3 = ifelse(multExpanded1$idmult %in% row.names(sig3),"+","-")
    
    #add DE direction to table
    multExpanded1$cont1up = ifelse(multExpanded1$idmult %in% row.names(c1up),"+","-")
    multExpanded1$cont1down = ifelse(multExpanded1$idmult %in% row.names(c1down),"+","-")
    multExpanded1$cont2up = ifelse(multExpanded1$idmult %in% row.names(c2up),"+","-")
    multExpanded1$cont2down = ifelse(multExpanded1$idmult %in% row.names(c2down),"+","-")
    multExpanded1$cont3up = ifelse(multExpanded1$idmult %in% row.names(c3up),"+","-")
    multExpanded1$cont3down = ifelse(multExpanded1$idmult %in% row.names(c3down),"+","-")
    
    #replace the names
    propernames <- c(paste(header, "SubtoDE", sep = ""), paste(header, "globalFsig", sep = ""), paste(header, "Fstat", sep = ""),
                     paste(header, "FPval", sep = ""), paste(header, "FAdjPval", sep = ""), paste(header, "DEcont1", sep = ""),
                     paste(header, "DEcont2", sep = ""), paste(header, "DEcont3", sep = ""), paste(header, "cont1up", sep = ""),
                     paste(header, "cont1down", sep = ""), paste(header, "cont2up", sep = ""), paste(header, "cont2down", sep = ""),
                     paste(header, "cont3up", sep = ""), paste(header, "cont3down", sep = ""))
    names(multExpanded1)[(length(names(multExpanded1))-(length(propernames) - 1)):length(names(multExpanded1))] <- propernames
    return(multExpanded1)
    }

    #process fits
    multExpanded1 <- ProcessFit(fit2 = ConfoundedFit, header = "Confounded", 
                                FitData = phosdata, multExpanded1)
    multExpanded1 <- ProcessFit(fit2 = PhosPrepCovFit, header = "PhosPrepCov", 
                                FitData = PhosProt, multExpanded1)
    multExpanded1 <- ProcessFit(fit2 = GelPrepCovFit, header = "GelPrepCov", 
                                FitData = PhosProtGel, multExpanded1)

    ############# Annotate the ME dataframe --------------
  
  #
  # NOTE THE FOLLOWING WHEN COMPARING ALL THE CONTRASTS AT ONCE PG 62 IN LIMMA USER GUIDE - decideTests method global should be used here.
  
  # method="global" is recommended when a set of closely related contrasts are being tested. This
  # method simply appends all the tests together into one long vector of tests, i.e., it treats all the tests
  # as equivalent regardless of which probe or contrast they relate to. An advantage is that the raw
  # p-value cutoff is consistent across all contrasts. For this reason, method="global" is recommended if
  # you want to compare the number of DE genes found for different contrasts, for example interpreting
  # the number of DE genes as representing the strength of the contrast. However users need to be aware
  # that the number of DE genes for any particular contrasts will depend on which other contrasts are
  # tested at the same time. Hence one should include only those contrasts which are closely related to
  # the question at hand. Unnecessary contrasts should be excluded as these would affect the results for
  # the contrasts of interest. Another more theoretical issue is that there is no theorem which proves that
  # adjust.method="BH" in combination with method="global" will correctly control the false discovery
  # rate for combinations of negatively correlated contrasts, however simulations, experience and some
  # theory suggest that the method is safe in practice.
  
 
  
  #write the multExpanded table with DE information
  write.table(multExpanded1,"multExpanded1_withDE1.csv", sep=',',col.names=T, row.names=F)
  saveRDS(multExpanded1, file = "multExpanded1_withDE.rds")
  return(multExpanded1)
}
