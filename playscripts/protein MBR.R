# First perform all processing steps using plyr and related tools.
# load required libraries
library(reshape2)
library(stringr)
library(plyr)

# Protein - all non reverse and contaminant proteins idd
# Protein1 - protein - "only Idd by site"
# Protein2 - protein1 with quantification in at least one sample
# Protein3 - protein1 with quant in every biological rep
# Protein4 - protein1 with quant in every sample

# Read in phospho table. Note the quote option is key.
protein <- read.table("./MQ output/7_29_output/ProteinGroups.txt", sep = "\t", header=T, fill = T, quote = "")

# subset those hits that are not contaminants and not reverse.
protein <- protein[(protein$Potential.contaminant != "+" & protein$Reverse != "+"),]## here I am specifying the rows using column variables and logical operators while not specifying any particular columns

# "only identified by site" hits are removed because they tend to have lower PEPs (wouldn't pass the FDR TH anyway) and can't be
# quantified since they are not idd by non-modified peptides. Note there are some high probability proteins here given some proteins
# are idd by 20+ phosphopeptides. eg is A6NKT7 (PEP = 2.23E-70)
protein1 <- protein[(protein$Only.identified.by.site != "+"),]



# Produce a trimmed matrix like before for quantification. Note unsure what's going on with sequence coverage
# and unique + razore variations. 


vars <- c("id", "Protein.IDs", "Majority.protein.IDs",  "Protein.names", "Gene.names", "Number.of.proteins", "Peptides", "Razor...unique.peptides", "Unique.peptides", "Sequence.coverage....", "Mol..weight..kDa.", "Sequence.length", "PEP", "Peptide.IDs", "Mod..peptide.IDs", "Phospho..STY..site.IDs")

other_data <- protein1[,vars]

##dataframe that collects only the relevent expression columns. NOTE THE NEED TO USE REP!!!!!
##The sample number precedes 'Rep' (technical replicate) and the triple underscore denotes the multiplicity 
expression <- protein1[,grep("Ratio.H.L.normalized(.*)_[12]", colnames(protein1))]

# Replace the column names
names(expression) <- sub(names(expression), pattern ="_", replacement = "Rep")

##combine the two
protein1 <- cbind(expression,other_data)

## remove rows with only NAs in expression columns

expCol <- grep("Ratio.H.L.(.*)", colnames(protein1))


protein2 <- protein1[rowSums(is.na(protein1[,expCol]))!=length(expCol),]##removes rows containing all 'NA's using the sums of the logical per row  


# remove sites if not observed in each sample (not sure how to automate this for larger datasets)
protein3 <- protein2[rowSums(is.na(protein2[ , 1:2])) < 2 & rowSums(is.na(protein2[ , 3:4])) < 2 & rowSums(is.na(protein2[ , 5:6])) < 2, ]

# remove sites if not observed in each sample (tech and biological) 
protein4 <- na.omit(protein3)

##############total sites and protein groups cut by samples and replicates #########################################################

##number of unique protein groups
proteingroups <- nrow(protein)

# number of unique protein groups with non-modified peptides used for id
proteingroups1 <- nrow(protein1)

# number of unique protein groups with quantification in at least one sample
proteingroups2 <- nrow(protein2)

# number of unique protein groups with quantification in each bio rep
proteingroups3 <- nrow(protein3)

# number of unique protein groups with quantification in each sample (tech and bio)
proteingroups4 <- nrow(protein4)

# percent of unique protein groups with phosphosite identification

percent_phos <- count(protein$Phospho..STY..site.IDs != "")[2,2]/proteingroups #total groups

percent_phos1 <- count(protein1$Phospho..STY..site.IDs != "")[2,2]/proteingroups1 #without idd by site

percent_phos2 <- count(protein2$Phospho..STY..site.IDs != "")[2,2]/proteingroups2 #quant in at least one sample

percent_phos3 <- count(protein3$Phospho..STY..site.IDs != "")[2,2]/proteingroups3 #quant in all 3 bio reps

percent_phos4 <- count(protein4$Phospho..STY..site.IDs != "")[2,2]/proteingroups4 #quant in all samples

# protein4$Phospho..STY..site.IDs[protein4$Phospho..STY..site.IDs == ""]=NA may be helpful in the future should probably melt then cast here

# above outputted as a table
out=data.frame(proteingroups,proteingroups1,proteingroups2,proteingroups3,proteingroups4, percent_phos, percent_phos1,percent_phos2,percent_phos3,percent_phos4)

#write.table(t(colnames(out)),"fitsheader.csv",sep=',',col.names=F,row.names=F,append=T)
write.table(out,"proteinstats.csv",sep=',',col.names=T,row.names=F)


# change protein2 names
names(protein2)[expCol] <- sub(names(protein2)[expCol], pattern = "Ratio.H.L.normalized", replacement = "HL")

# barplot of number of unique protein groups that are quantified in at least one sample
uniqueQuantEvents <- colSums(!is.na(protein2[expCol]))
par(mar=c(6,4,4,4))
barplot(uniqueQuantEvents, las=2, cex.names = .8, ylim = c(0,max(uniqueQuantEvents)+500), main = "proteins per replicate")##number of quant events


# barplot of number of overlapping protein groups common to 6,5,4,3,2,1 etc replicate
ExpOverlap <- table(rowSums(!is.na(protein2[,expCol])))##removes rows containing all 'NA's using the sums of the logical per row                        
ExpOverlap <- rev(ExpOverlap)
barplot(ExpOverlap, las = 1)

##barplot with summary line overlay

#vector of percentages
percentExp <- ExpOverlap/sum(ExpOverlap)

#vector of cumulative percentages
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

percentTotalExp <- cumulative(percentExp)##as percentage


##Overlaid graphic of sample overlap and cumulative percentage using base graphics. Must be aligned later
BarCumOverlay <- function(overlap,cumPercent){
  bp <- barplot(overlap)
  bp <- barplot(overlap,las=1, cex.names = 1, ann=FALSE, xlim = c(0,max(bp)+1), ylim = c(0,max(overlap)+500), ylab = "# of overlapping proteins", xlab = "overlap between N samples")##note needs a matrix as input and other variables used
  par(new=TRUE)
  par(mar=c(5,4,4,4))
  plot(bp,cumPercent,axes="FALSE", ann=FALSE, xlim = c(0,max(bp)+1), ylim = c(0,100), col = "red", type = "b", pch=19)##note the same coordinate ranges 'xlim' so that the points are in the center of the barchart; type b is points connected by lines.
  mtext("% of total proteins",side=4,line=2)
  axis(4,at=seq(0,100,10), las=1)
  box()
}
BarCumOverlay(ExpOverlap,percentTotalExp)


#################################################################################################################
# Data Analysis. Histograms, Dendograms/Clustering/Heatmaps, PCA, DE, QQplots for 'protein3' dataset


# Filter so that there is at least one valid value in each sample
data <- protein2[,expCol]

#replace row names with protein group id
row.names(data) <- protein2$id

#log transform the data
data <- log2(data)

#normalize by col median using columnwise (when to normalize??)

median.subtract <- function(x){ x - median(x, na.rm = TRUE)}##create a wrapper for median subtraction


datanorm <- colwise(median.subtract)(data) #create median subtracted data but loose the row names here...

row.names(datanorm) <- row.names(data)##add back the row names


# remove sites if not observed in each sample (not sure how to automate this for larger datasets)
datanorm <- datanorm[rowSums(is.na(datanorm[ , 1:2])) < 2 & rowSums(is.na(datanorm[ , 3:4])) < 2 & rowSums(is.na(datanorm[ , 5:6])) < 2, ]                    

# note the dimensions here are 6071 where n>=3 for each phosphosite




# Histograms of datanorm
par(mfrow = c(3,2))
proth1 <- hist(datanorm$HL.16770Rep1, breaks = 60)
proth2 <- hist(datanorm$HL.16770Rep2, breaks = 60)
proth3 <- hist(datanorm$HL.16778Rep1, breaks = 60)
proth4 <- hist(datanorm$HL.16778Rep2, breaks = 60)
proth5 <- hist(datanorm$HL.16788Rep1, breaks = 60)
proth6 <- hist(datanorm$HL.16788Rep2, breaks = 60)


################################################ Dendograms and Clustering  ##########################################

#Zscore scale the samples. Note that NAs have an effect on the standardization and the distance metrics so they are removed
# Z score subtracts the mean of each vector from each element, then divides by the sd of the vector
# Z score makes each vector mean = 0 and sd = 1
datanorm <- na.omit(datanorm)#now I have limited to observed in all cases.NOTE WAS NORMALIZED USING MORE DATA POINTS



dataZ <- scale(datanorm)##Z-scored column wise

# now all data excepting complete cases (note that the sample dendograms look the same)
hist(dataZ[,6], breaks = 100)

# dendogram using euclidian distance (default) and ward or complete agglomeration
dend.ward<- as.dendrogram(hclust(dist(t(dataZ)),method="ward"))
dend.complete<- as.dendrogram(hclust(dist(t(dataZ))))

ward.o<- order.dendrogram(dend.ward)
complete.o<- order.dendrogram(dend.complete)

par(mar = c(6,4,4,4))
plot(dend.complete,ylab="height", main = "Euclidian/Complete")
plot(dend.ward, leaflab = "perpendicular", ylab = "height", main = "Euclidian/Ward")

dev.off()##turns off plot

# Cluster using euclidian distance and ward linkage for both sites(rows) and samples (columns)
# Note that both dendograms are created independently and row Z scores are presented in the heatmap

# row scaled
r <- t(scale(t(datanorm)))#transpose to zscale the rows then transpose back to original format

# sample scaled
c <- scale(datanorm)


# install heatmap.2 package
install.packages("gplots")
library(gplots)

# Create dendrogram using the data without NAs
feature.dend<- as.dendrogram(hclust(dist(r),method="ward"))
sample.dend<- as.dendrogram(hclust(dist(t(c)),method="ward"))##note that dist caclculates distance between rows by default


##produce the heatmap. Note that the help page has a nice section on identifying subregions by color. Although I will likely have to cut the dendogram to id clusters of interest

heatmap.2(
  r,#row Z scores
  Colv=sample.dend,
  Rowv=feature.dend,
  col=bluered(25),
  scale="none",
  trace="none",
  key.xlab = "Row Z scores", key.ylab=NULL, key.title = "",
  srtCol=45,  ,adjCol = c(1,1),
  margins = c(6,5),
  cexCol=1,
  labRow = NA#remove row labels
)


dev.off()


###########################################################################################################################

##PCA of data from quick R (I will need to use the column standardized data next)

# fit <- princomp(datanorm, cor=TRUE)
# summary(fit) # print variance accounted for
# loadings(fit) # pc loadings
# plot(fit,type="lines") # scree plot
# fit$scores # the principal components
# biplot(fit) 
# plot(fit)
# 
# 
# # PCA two using prcomp
# pca.res <- prcomp(datanorm, retx=TRUE)
# pca.res
# summary(pca.res)
# 
# 
# # Get principal component vectors using prcomp instead of princomp
# pc <- prcomp(dataZ)
# 
# # First for principal components
# comp <- data.frame(pc$x[,1:4])
# # Plot
# plot(comp, pch=16, col=rgb(0,0,0,0.5))







# Rafa PCA plots!
x <- t(datanorm)#samples are the rows of the column matrix
pc <- prcomp(x)#scale = T, center = T) as of now I am not scaling

names(pc)

cols <- as.factor(substr(colnames(datanorm), 4, 8))##check me out. use 5 digit exp name.
plot(pc$x[, 1], pc$x[, 2], col=as.numeric(cols), main = "PCA", xlab = "PC1", ylab = "PC2")
legend("bottomright", levels(cols), col = seq(along=levels(cols)), pch = 1)


summary(pc)

#SVD for calculating variance explained; see Rafa's notes for an explaination
cx <- sweep(x, 2, colMeans(x), "-")
sv <- svd(cx)
names(sv)
plot(sv$u[, 1], sv$u[, 2], col = as.numeric(cols), main = "SVD", xlab = "U1", ylab = "U2")


plot(sv$d^2/sum(sv$d^2), xlim = c(1, 6), type = "b", pch = 16, xlab = "principal components", 
     ylab = "variance explained")


# scaled? analysis...
# x <- t(dataZ)
# pc <- prcomp(x)
# 
# names(pc)
# 
# cols <- as.factor(substr(colnames(datanorm), 3, 7))##check me out. use 5 digit exp name.
# plot(pc$x[, 1], pc$x[, 2], col=as.numeric(cols), main = "PCA", xlab = "PC1", ylab = "PC2")
# legend("bottomleft", levels(cols), col = seq(along=levels(cols)), pch = 1)






# MDS with no Zscaling (MDS using euclidian distance is equivalent to PCA first two components)
d <- dist(t(datanorm))##note that their is no Z scaling here ?!
mds <- cmdscale(d)
cols <- as.factor(substr(colnames(datanorm), 3, 7))##check me out. use 5 digit exp name..
plot(mds, col=as.numeric(cols), lwd = 1.5)
legend("topleft", levels(cols), col = seq(along=levels(cols)), pch = 1)

# MDS with Zscaling I think...
d <- dist(t(dataZ))
mds <- cmdscale(d)
cols <- as.factor(substr(colnames(datanorm), 3, 7))##check me out. got this muthatrucka. use 5 digit exp name..
plot(mds, col=as.numeric(cols), lwd = 1.5)
legend("topleft", levels(cols), col = seq(along=levels(cols)), pch = 1)









#################################################################################################################




