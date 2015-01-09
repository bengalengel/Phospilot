breakdown <- function(protein, phospho, multExpanded, cls=TRUE){
  require(plyr)
  require(VennDiagram)
  require(gridExtra)
  
  ##This function prints multiple tables and venn diagrams for phospho and protein files. 
  ##It accepts a protein file, a phospho file, and an  expanded phospho file (each row is an experimental observation), 
  ##and a class designation for the expanded phospho file. The default is class 1 (cls=TRUE).
  
  # General functions
  nmeasure <- function(x) sum(!(is.na(x)))#function to count valid values
  totalgt0 <- function(x) sum(x > 0, na.rm = TRUE)#function that counts total greater than 0
  #create vectors of cumulative percentages for overlay graphic called within 'if' statements below
  cumulative <- function(x) {
    s <- as.numeric(vector(length = nrow(x)))
    s[1] <- x[1]
    for (i in seq_along(x)) {
      if(i>1){
        s[i] <- s[i-1] + x[i]
      }
    }
    return(s*100)
  }
  
  # Functions and definitions for protein work
  expCol <- grep("H.L.(.*)", colnames(protein1))
  protein2 <- protein[rowSums(is.na(protein[,expCol]))!=length(expCol),]
  # change protein2 names
  names(protein2)[expCol] <- sub(names(protein2)[expCol], pattern = "Ratio.H.L.normalized.", replacement = "HL")
  
  # barplot of number of unique protein groups that are quantified in at least one sample
  uniqueQuantEvents <- colSums(!is.na(protein2[expCol]))
  par(mar=c(6,4,4,4))
  barplot(uniqueQuantEvents, las=2, cex.names = .8, ylim = c(0,max(uniqueQuantEvents)+500), 
          main = "proteins per replicate")##number of quant events
  
  # barplot of number of overlapping protein groups common each sample
  ExpOverlap <- table(rowSums(!is.na(protein2[,expCol])))
  ExpOverlap <- rev(ExpOverlap)
  barplot(ExpOverlap, las = 1, cex.names = .8, ylab = "Proteins", main = 
            "Overlap between N samples")
  
  ##barplot with summary line overlay
  
  #vector of percentages
  percentExp <- ExpOverlap/sum(ExpOverlap)
  percentTotalExp <- cumulative(percentExp)##as percentage
  ##Overlaid graphic of sample overlap and cumulative percentage using base graphics. Must be aligned later
  BarCumOverlay <- function(overlap,cumPercent){
    bp <- barplot(overlap)
    bp <- barplot(overlap,las=1, cex.names = .8, ann=FALSE, xlim = c(0,max(bp)+1), ylim = c(0,max(overlap)+500), ylab = "# of overlapping proteins", xlab = "overlap between N samples")##note needs a matrix as input and other variables used
    par(new=TRUE)
    par(mar=c(5,4,4,4))
    plot(bp,cumPercent,axes="FALSE", ann=FALSE, xlim = c(0,max(bp)+1), ylim = c(0,100), col = "red", type = "b", pch=19)##note the same coordinate ranges 'xlim' so that the points are in the center of the barchart; type b is points connected by lines.
    mtext("% of total proteins",side=4,line=2)
    axis(4,at=seq(0,100,10), las=1)
    box()
  }
  BarCumOverlay(ExpOverlap,percentTotalExp)
  
  
  
  #   Functions and definitions for phosphosite and phosphomeasurement calculations
  expCol <- grep("HL(.*)", colnames(multExpanded))
  #newnames <- colnames(multExpanded)[grep("_", colnames(multExpanded))]#assumes experimental obs have an underscore
  newnames <- colnames(multExpanded)[expCol]
  idBreakdown <- ddply(multExpanded,.(id), colwise(nmeasure,newnames))##breaks down by id number of mults observed per sample
  #the subset of each vector > 0 gives the unique ids
  uniqueids <- colwise(totalgt0,newnames)(idBreakdown)##total number of unique ids
  # number of overlapping experimental observations common to all experiments
  ExpOverlap <- table(rowSums(!is.na(multExpanded[,expCol])))##removes rows containing all 'NA's using the sums of the logical per row                        
  ExpOverlap <- rev(ExpOverlap)
  # number of overlapping phosphosites common to all experiments
  idOverlap <- table(rowSums(idBreakdown[,newnames] > 0))
  idOverlap <- rev(idOverlap)
  #vector of percentages
  percentExp <- ExpOverlap/sum(ExpOverlap)
  percentId <- idOverlap/sum(idOverlap)
  percentTotalExp <- cumulative(percentExp)##as percentage
  percentTotalId <- cumulative(percentId)##as percentage
  
  
  if(cls){
    #   Overlap for all class one phospho MEASUREMENTS per sample
    uniqueQuantEvents <- colSums(!is.na(multExpanded[expCol]))
    barplot(uniqueQuantEvents, las=1, cex.names = .80, ylab = "class 1 phosphomeasurements",
            main = "class 1 phosphomeasurements per sample")##number of quant events
    #     Overlap for all class one phosphosites per sample
    barplot(as.matrix(uniqueids),las=1, cex.names = .80, ylab = "Class 1 phosphosites", main = 
              "Class 1 phosphosites per sample")##note needs a matrix as input and other variables used
    # barplot of number of overlapping experimental observations and sites common to all experiments
    barplot(ExpOverlap, las = 1, cex.names = .80, ylab = "Class 1 phosphomeasurements", main = 
              "Overlap between N samples")
    barplot(idOverlap, las = 1, cex.names = .80, ylab = "Class 1 phosphosites", main = 
              "Overlap between N samples")
    ##Overlaid graphic of sample overlap and cumulative percentage using base graphics. Must be aligned later
    BarCumOverlay <- function(overlap,cumPercent, site="yes"){
      bp <- barplot(overlap)
      if(site=="yes"){
        bp <- barplot(overlap,las=1, cex.axis = 1, cex.names = 1, cex.lab = 1, ann=FALSE, xlim = c(0,max(bp)+1), ylim = c(0,max(overlap)+500), ylab = "# of overlapping class 1 phosphosites", xlab = "overlap between N samples")##note needs a matrix as input and other variables used
        par(new=TRUE)
        par(mar=c(5,4,4,4))
        plot(bp,cumPercent,axes="FALSE", ann=FALSE, xlim = c(0,max(bp)+1), ylim = c(0,100), col = "red", type = "b", pch=19)##note the same coordinate ranges 'xlim' so that the points are in the center of the barchart; type b is points connected by lines.
        mtext("% of total phosphosites",side=4,line=2, cex.axis = 1)
        axis(4,at=seq(0,100,10), las=1, cex.axis = 1)
        box()
      }
      if(site=="no"){
        bp <- barplot(overlap,las=1, cex.axis = 1, cex.names = 1, cex.lab = 1, ann=FALSE, xlim = c(0,max(bp)+1), ylim = 
                        c(0,max(overlap)+500), ylab = "# of overlapping class 1 phosphomeasurements", 
                      xlab = "overlap between N samples")##note needs a matrix as input and other variables used
        par(new=TRUE)
        par(mar=c(5,4,4,4))
        plot(bp,cumPercent,axes="FALSE", ann=FALSE, xlim = c(0,max(bp)+1), ylim = c(0,100), col = "red", type = "b", pch=19)##note the same coordinate ranges 'xlim' so that the points are in the center of the barchart; type b is points connected by lines.
        mtext("% of total phosphomeasurements",side=4,line=2, cex.axis = 1)
        axis(4,at=seq(0,100,10), las=1, cex.axis = 1)
        box()
      }
    }
    BarCumOverlay(ExpOverlap,percentTotalExp,site="no")
    BarCumOverlay(idOverlap,percentTotalId)
    
    
    ##number of modifications per protein. Here I can use the leading razor protein associated with each site.
    #     or I can use the protein groups file.
    #Must go from the sites file so that each site is used only once as opposed to each group used only once with the same site assigned multiple times
    
    hist(table(phospho$Protein), breaks = max(table(phospho$Protein)), xlim = c(0,20), xlab = "Number of phosphosites", ylab = "Number of Proteins", main = "number of C1 sites per protein")
    
    ##pie chart of AAs with percentages of total IDd CLASS 1 (but not necessarily quantified!)
    mytable <- table(phospho$Amino.acid)
    lbls <- paste(names(mytable), mytable, sep=" ")##pastes the labels and numbers from the table
    pct <- round(mytable/(sum(mytable)),3)*100 ##calculates percentages
    pct <- paste0(pct,"%") ##adds % sign
    pct <- paste("(",pct,")",sep="") ##adds parentheses
    lbls <- paste(lbls, pct,sep=" ") ##combines
    pie(mytable, labels = lbls,
        main="Class 1 Amino Acid breakdown")
    
    ##pie chart of multiplicity with percentages
    mytable <- table(multExpanded$multiplicity)
    lbls <- paste(names(mytable), mytable, sep=" ")##pastes the labels and numbers from the table
    pct <- round(mytable/(sum(mytable)),3)*100 ##calculates percentages
    pct <- paste0(pct,"%") ##adds % sign
    pct <- paste("(",pct,")",sep="") ##adds parentheses
    lbls <- paste(lbls, pct,sep=" ") ##combines
    pie(mytable, labels = lbls,
        main="Multiplicity states of class 1 quantifications")
    
    # Class 1 phospho breakdowns by multiplicity;can be used later to subset analyses by confidence of quantification.
    z <- table(multExpanded$id,multExpanded$multiplicity) ##two dimensional table
    #colSums(z)
    ##should be a better way to do this but I will just subset
    colnames(z) <- c("m1","m2","m3")
    dfz <- as.data.frame.matrix(z)##converts table into a dataframe for 
    
    id1 <- dfz[dfz$m1==1 & dfz$m2==0 & dfz$m3==0,]
    nrow(id1)##sum also works
    
    id2 <- dfz[dfz$m1==0 & dfz$m2==1 & dfz$m3==0,]##which also works here
    nrow(id2)##sum also works
    
    id3 <- dfz[dfz$m1==0 & dfz$m2==0 & dfz$m3==1,]##which also works here
    nrow(id3)##sum also works
    
    id12 <- dfz[dfz$m1==1 & dfz$m2==1 & dfz$m3==0,]
    nrow(id12)##sum also works
    
    id13 <- dfz[dfz$m1==1 & dfz$m2==0 & dfz$m3==1,]##which also works here
    nrow(id13)##sum also works
    
    id23 <- dfz[dfz$m1==0 & dfz$m2==1 & dfz$m3==1,]##which also works here
    nrow(id23)##sum also works
    
    id123 <- dfz[dfz$m1==1 & dfz$m2==1 & dfz$m3==1,]##which also works here
    nrow(id123)##sum also works
    
    
    ##output as table arranged in decending order
    combos <- c(nrow(id1), nrow(id2), nrow(id3), nrow(id12), nrow(id13), nrow(id23), nrow(id123)) ##vector of counts
    combos <- as.data.frame(combos)
    rownames(combos) <- c(1,2,3,12,13,23,123)
    combos <- cbind(combos,prop.table(combos))
    colnames(combos) <- c("count","%")
    combos=cbind(multiplicity=row.names(combos), combos)
    combos <- arrange(combos,desc(`%`))##note the backticks
    
    # #make the venn diagram
    plot.new()
    venn.plot <- draw.triple.venn(
      area1 = sum(dfz$m1),
      area2 = sum(dfz$m2),
      area3 = sum(dfz$m3),
      n12 = nrow(id12)+nrow(id123),
      n23 = nrow(id23)+nrow(id123),
      n13 = nrow(id13)+nrow(id123),
      n123 = nrow(id123),
      category = c("Singly", "Doubly", "Triply"),
      fill = c("orange", "green", "blue"),
      lty = "blank",
      cex = 2,
      cat.cex = 2,
      cat.col = c("orange", "green", "blue"), 
      margin = .1,
      main="test"
    )
    plot.new()
    grid.arrange(gTree(children=venn.plot), main="Class 1 phosphopeptide multiplicity")
  }
  
  
  
  if(!cls){
    #   Overlap for all phospho MEASUREMENTS per sample
    uniqueQuantEvents <- colSums(!is.na(multExpanded[expCol]))
    barplot(uniqueQuantEvents, las=1, cex.names = .80, ylab = "phosphomeasurements",
            main = "phosphomeasurements per sample")##number of quant events
    #     Overlap for all class one phosphosites per sample
    barplot(as.matrix(uniqueids),las=1, cex.names = .80, ylab = "Phosphosites", main = 
              "Phosphosites per sample")##note needs a matrix as input and other variables used
    # barplot of number of overlapping experimental observations and sites common to all experiments
    barplot(ExpOverlap, las = 1, cex.names = .80, ylab = "Phosphomeasurements", main = 
              "Overlap between N samples")
    barplot(idOverlap, las = 1, cex.names = .80, ylab = "Phosphosites", main = 
              "Overlap between N samples")
    
    ##Overlaid graphic of sample overlap and cumulative percentage using base graphics. Must be aligned later
    BarCumOverlay <- function(overlap,cumPercent, site="yes"){
      bp <- barplot(overlap)
      if(site=="yes"){
        bp <- barplot(overlap,las=1, cex.axis = 1, cex.names = 1, cex.lab = 1, ann=FALSE, xlim = c(0,max(bp)+1), ylim = c(0,max(overlap)+500), ylab = "# of overlapping phosphosites", xlab = "overlap between N samples")##note needs a matrix as input and other variables used
        par(new=TRUE)
        par(mar=c(5,4,4,4))
        plot(bp,cumPercent,axes="FALSE", ann=FALSE, xlim = c(0,max(bp)+1), ylim = c(0,100), col = "red", type = "b", pch=19)##note the same coordinate ranges 'xlim' so that the points are in the center of the barchart; type b is points connected by lines.
        mtext("% of total phosphosites",side=4,line=2, cex.axis = 1)
        axis(4,at=seq(0,100,10), las=1, cex.axis = 1)
        box()
      }
      if(site=="no"){
        bp <- barplot(overlap,las=1, cex.axis = 1, cex.names = 1, cex.lab = 1, ann=FALSE, xlim = c(0,max(bp)+1), ylim = 
                        c(0,max(overlap)+500), ylab = "# of overlapping phosphomeasurements", 
                      xlab = "overlap between N samples")##note needs a matrix as input and other variables used
        par(new=TRUE)
        par(mar=c(5,4,4,4))
        plot(bp,cumPercent,axes="FALSE", ann=FALSE, xlim = c(0,max(bp)+1), ylim = c(0,100), col = "red", type = "b", pch=19)##note the same coordinate ranges 'xlim' so that the points are in the center of the barchart; type b is points connected by lines.
        mtext("% of total phosphomeasurements",side=4,line=2, cex.axis = 1)
        axis(4,at=seq(0,100,10), las=1, cex.axis = 1)
        box()
      }
    }
    BarCumOverlay(ExpOverlap,percentTotalExp,site="no")
    BarCumOverlay(idOverlap,percentTotalId)
    
    ##number of modifications per protein. Here I can use the leading razor protein associated with each site.
    #     or I can use the protein groups file.
    #Must go from the sites file so that each site is used only once as opposed to each group used only once with the same site assigned multiple times
    
    hist(table(phospho$Protein), breaks = max(table(phospho$Protein)), xlim = c(0,20), xlab = "Number of phosphosites", ylab = "Number of Proteins", main = "number of sites per protein")
    
    ##pie chart of AAs with percentages of total IDd phosphosites (but not necessarily quantified!)
    mytable <- table(phospho$Amino.acid)
    lbls <- paste(names(mytable), mytable, sep=" ")##pastes the labels and numbers from the table
    pct <- round(mytable/(sum(mytable)),3)*100 ##calculates percentages
    pct <- paste0(pct,"%") ##adds % sign
    pct <- paste("(",pct,")",sep="") ##adds parentheses
    lbls <- paste(lbls, pct,sep=" ") ##combines
    pie(mytable, labels = lbls,
        main="Amino Acid breakdown")
    
    ##pie chart of multiplicity with percentages
    mytable <- table(multExpanded$multiplicity)
    lbls <- paste(names(mytable), mytable, sep=" ")##pastes the labels and numbers from the table
    pct <- round(mytable/(sum(mytable)),3)*100 ##calculates percentages
    pct <- paste0(pct,"%") ##adds % sign
    pct <- paste("(",pct,")",sep="") ##adds parentheses
    lbls <- paste(lbls, pct,sep=" ") ##combines
    pie(mytable, labels = lbls,
        main="Multiplicity states of phosphopeptide quantifications")
    
    # Class 1 phospho breakdowns by multiplicity;can be used later to subset analyses by confidence of quantification.
    z <- table(multExpanded$id,multExpanded$multiplicity) ##two dimensional table
    #colSums(z)
    ##should be a better way to do this but I will just subset
    colnames(z) <- c("m1","m2","m3")
    dfz <- as.data.frame.matrix(z)##converts table into a dataframe for 
    
    id1 <- dfz[dfz$m1==1 & dfz$m2==0 & dfz$m3==0,]
    nrow(id1)##sum also works
    
    id2 <- dfz[dfz$m1==0 & dfz$m2==1 & dfz$m3==0,]##which also works here
    nrow(id2)##sum also works
    
    id3 <- dfz[dfz$m1==0 & dfz$m2==0 & dfz$m3==1,]##which also works here
    nrow(id3)##sum also works
    
    id12 <- dfz[dfz$m1==1 & dfz$m2==1 & dfz$m3==0,]
    nrow(id12)##sum also works
    
    id13 <- dfz[dfz$m1==1 & dfz$m2==0 & dfz$m3==1,]##which also works here
    nrow(id13)##sum also works
    
    id23 <- dfz[dfz$m1==0 & dfz$m2==1 & dfz$m3==1,]##which also works here
    nrow(id23)##sum also works
    
    id123 <- dfz[dfz$m1==1 & dfz$m2==1 & dfz$m3==1,]##which also works here
    nrow(id123)##sum also works
    
    
    ##output as table arranged in decending order
    combos <- c(nrow(id1), nrow(id2), nrow(id3), nrow(id12), nrow(id13), nrow(id23), nrow(id123)) ##vector of counts
    combos <- as.data.frame(combos)
    rownames(combos) <- c(1,2,3,12,13,23,123)
    combos <- cbind(combos,prop.table(combos))
    colnames(combos) <- c("count","%")
    combos=cbind(multiplicity=row.names(combos), combos)
    combos <- arrange(combos,desc(`%`))##note the backticks
    
    #make the venn diagram
    plot.new()
    venn.plot <- draw.triple.venn(
      area1 = sum(dfz$m1),
      area2 = sum(dfz$m2),
      area3 = sum(dfz$m3),
      n12 = nrow(id12)+nrow(id123),
      n23 = nrow(id23)+nrow(id123),
      n13 = nrow(id13)+nrow(id123),
      n123 = nrow(id123),
      category = c("Singly", "Doubly", "Triply"),
      fill = c("orange", "green", "blue"),
      lty = "blank",
      cex = 2,
      cat.cex = 2,
      cat.col = c("orange", "green", "blue"), 
      margin = .1,
      main="test"
    )
    plot.new()
    grid.arrange(gTree(children=venn.plot), main="Phosphopeptide multiplicity")
  }
}
  