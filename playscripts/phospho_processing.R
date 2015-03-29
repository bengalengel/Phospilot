# This program will produce a series of charts relating to the phospho(STY) output table from maxquant. 
# It is IMPORTANT that the experiments be processed in MQ as "XXXXXRepY", where the first set of numbers
# is the cell line and the last number is the replicate number. Many of the regular expressions that I use require this
# type of annotation.


# First perform all processing steps using plyr and related tools.
# load required libraries
library(reshape2)
library(stringr)
library(plyr)


# Read in phospho table. Note the quote option is key.
phospho <- read.table("./MQ output/7_17_output/Phospho (STY)Sites.txt", sep = "\t", header=T, fill = T, quote = "")

# subset those hits that are not contaminants and not reverse hits.

phospho <- phospho[(phospho$Contaminant != "+" & phospho$Reverse != "+"),]## here I am specifying the rows using column variables and logical operators while not specifying any particular columns
# I can come back here to get phospho ID information if needed


# subset to class 1
phospho1 <- phospho[(phospho$Localization.prob >= .75),]##Why just one localization probability?...
  # Here the BEST localization probablility across all samples is used to justify identification and use in each. What doesn't make sense is 
# multiplicity, where a single probability per sample is assigned for each multiplicity state

#####now must condense DF and 'melt' the dataframe so that each ratio for each multiplicity state has its own observation

##first condense dataframe to include variables of interest. colect non-expression variables first
other_data <- phospho1[,c("id","Amino.acid","Charge","Reverse","Contaminant","Proteins","Positions.within.proteins","Leading.proteins",
                          "Sequence.window","Modified.sequence","Localization.prob","PEP", "Score", "Delta.score", "Score.for.localization", 
                          "m.z", "Mass.error..ppm.", "Intensity", "Intensity.L", "Intensity.H", "Position", "Number.of.Phospho..STY.", 
                          "Protein.group.IDs", "Protein") ]

##dataframe that collects only the relevent expression columns. NOTE THE NEED TO USE REP!!!!!
##The sample number precedes 'Rep' (technical replicate) and the triple underscore denotes the multiplicity 
expression <- phospho1[,grep("Ratio.H.L.normalized(.*)Rep.___", colnames(phospho1))]

##combine the two
phospho2 <- cbind(expression,other_data)

##now to ensure each multiplicity has its own row and variables are condensed. I am converting from wide to long format to ensure that each
# observation is uniquely represented in the data

# Use reshape2 'melt' function 
##lowercase names if I want
# names(phospho2) <- tolower(names(phospho2))

##melt all expression variables into one column

#names(expression)
#melted <- melt(phospho2, measure.vars = names(expression), variable.name = "sample_rep_mult", value.name = "H/L Ratio")

melted <- melt(phospho2, measure.vars = names(expression))

# Here I split (that is add a variable identifier) the melted 'variable' column so that the sample, replicate, and multiplicity are now explicit
# for each "variable"
melted <- cbind(melted, colsplit(melted$variable, "Rep", c("sample", "replicate_mult"))) ##first split
melted <- cbind(melted, colsplit(melted$replicate_mult, "___", c("replicate","multiplicity"))) ##second split
melted$sample <- gsub(melted$sample, pattern = "Ratio.H.L.normalized.", replacement = "") ##remove redundant information next 3 lines
drop <- "replicate_mult" 
melted <- melted[,!(names(melted) %in% drop)]


##cast data so that each unique 'sample/replicate' combination has its own column populated by the measurement 'currently the value column'.

# testnew <- dcast(test, ... ~ variable)##will recast the entire data frame
# testnew1 <- dcast(test, ... ~ sample, value.var="value") ##casts by sample

casted <- dcast(melted, ... ~ sample + replicate, value.var="value") ##close but creates extra rows


##produce the multiplicity explicit table **********************************

##gives index of experiment and replicate
data <- grep("_", colnames(casted))

##gives string of experiment and replicate
data2 <- colnames(casted)[grep("_", colnames(casted))]

#produces a new string with proper alpha leading R variable names
newnames <- paste0("HL",data2)

# rename the experiment variables within the dataframe
colnames(casted)[data] <- newnames

## columnwise application of mean to condense the dataframe. PRODUCES CORRECT NUMBERS!
out <- ddply(casted, .(id, multiplicity), colwise(mean,newnames,na.rm=T))

#merge with identifying information by id to produce the multiplicity expanded table (each obs has a row)
multExpanded <- merge(other_data, out, by="id")

## remove rows with only NAs in expression columns

expCol <- grep("HL(.*)", colnames(multExpanded))

multExpanded <- multExpanded[rowSums(is.na(multExpanded[,expCol]))!=length(expCol),]##removes rows containing all 'NA's using the sums of the logical per row                        

##############total sites and protein groups cut by samples and replicates #########################################################

##number of unique sites
sites <- nrow(phospho)

##number of unique class 1 sites
class1 <- nrow(phospho1)

##number of unique protein groups associated with class 1 sites
pgroups <- nrow(table(phospho1$Proteins))

leadingp <- nrow(table(phospho1$Leading.proteins))##These are the proteins from each group with the most Ids (first on the list anyway)

protein <- nrow(table(phospho1$Protein))##This is the leading razor protein (the one with the most ids amongst all the associated proteins)

##number of unique class 1 sites with quantification in at least one sample
##note that the multiplicity references the phosphorylation state of the peptide (singly,doubly,triply) when this site was quantified.
quantClass1 <- nrow(table(multExpanded$id))

##protein groups for quantified class 1
pgroupsClass1 <- nrow(table(multExpanded$Proteins))

##number of leading proteins class 1
leadProteinClass1 <- nrow(table(multExpanded$Leading.proteins))

# number of sites per class 1 replicate 
uniqueQuantEvents <- colSums(!is.na(multExpanded[expCol]))
barplot(uniqueQuantEvents, las=1, cex.names = .80)##number of quant events

# number of unique ids per experiment (CERTAINLY A MUCH BETTER WAY TO DO THIS!)

nmeasure <- function(x) sum(!(is.na(x)))#function to count valid values

#nmeasure(multExpanded$HL16770_1)##works
#nmeasure(multExpanded[,expCol]) ##sums over all the columns! must use colwise (plyr) see below

#t <- colwise(nmeasure,newnames)(multExpanded)#same as unique events above

idBreakdown <- ddply(multExpanded,.(id), colwise(nmeasure,newnames))##breaks down by id number of mults observed per sample

#now I just need the subset of each vector which is greater than 0
totalgt0 <- function(x) sum(x > 0, na.rm = TRUE)#function that counts total greater than 0

uniqueids <- colwise(totalgt0,newnames)(idBreakdown)##total number of unique ids

barplot(as.matrix(uniqueids),las=1, cex.names = .80)##note needs a matrix as input and other variables used


# barplot of number of overlapping experimental observations (not sites) common to 6,5,4,3,2,1 etc replicate
ExpOverlap <- table(rowSums(!is.na(multExpanded[,expCol])))##removes rows containing all 'NA's using the sums of the logical per row                        
ExpOverlap <- rev(ExpOverlap)
barplot(ExpOverlap, las = 1)

# barplot of number of overlapping sites common to replicates
idOverlap <- table(rowSums(idBreakdown[,newnames] > 0))
idOverlap <- rev(idOverlap)
barplot(idOverlap, las = 1)

##barplot with summary line overlay

#vector of percentages
percentExp <- ExpOverlap/sum(ExpOverlap)
percentId <- idOverlap/sum(idOverlap)

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
percentTotalId <- cumulative(percentId)##as percentage


##Overlaid graphic of sample overlap and cumulative percentage using base graphics. Must be aligned later
BarCumOverlay <- function(overlap,cumPercent){
bp <- barplot(overlap)
bp <- barplot(overlap,las=1, cex.names = 1, ann=FALSE, xlim = c(0,max(bp)+1), ylim = c(0,max(overlap)+500), ylab = "# of overlapping phosphosites", xlab = "overlap between N samples")##note needs a matrix as input and other variables used
par(new=TRUE)
par(mar=c(5,4,4,4))
plot(bp,cumPercent,axes="FALSE", ann=FALSE, xlim = c(0,max(bp)+1), ylim = c(0,100), col = "red", type = "b", pch=19)##note the same coordinate ranges 'xlim' so that the points are in the center of the barchart; type b is points connected by lines.
mtext("% of total phosphosites",side=4,line=2)
axis(4,at=seq(0,100,10), las=1)
box()
}
BarCumOverlay(ExpOverlap,percentTotalExp)
BarCumOverlay(idOverlap,percentTotalId)




##number of modifications per protein. Here I can use the leading razor protein associated with each site or I can use the protein groups file.
#Must go from the sites file so that each site is used only once as opposed to each group used only once with the same site assigned multiple times

hist(table(phospho1$Protein), breaks = max(table(phospho1$Protein)), xlim = c(0,20), xlab = "Number of phosphosites", ylab = "Number of Proteins", main = "number of C1 sites per protein")





#################################################################################################################################

###################################################PIES in separate script later

##pie chart of AAs wither percentages of total IDd (but not necessarily quantified) FROM PHOSPHO TABLE!
mytable <- table(phospho$Amino.acid)
lbls <- paste(names(mytable), mytable, sep=" ")##pastes the labels and numbers from the table
pct <- round(mytable/(sum(mytable)),3)*100 ##calculates percentages
pct <- paste0(pct,"%") ##adds % sign
pct <- paste("(",pct,")",sep="") ##adds parentheses
lbls <- paste(lbls, pct,sep=" ") ##combines
pie(mytable, labels = lbls,
    main="Amino Acid breakdown")


##pie chart of AAs wither percentages of total IDd CLASS 1 (but not necessarily quantified) FROM PHOSPHO1 TABLE!
mytable <- table(phospho1$Amino.acid)
lbls <- paste(names(mytable), mytable, sep=" ")##pastes the labels and numbers from the table
pct <- round(mytable/(sum(mytable)),3)*100 ##calculates percentages
pct <- paste0(pct,"%") ##adds % sign
pct <- paste("(",pct,")",sep="") ##adds parentheses
lbls <- paste(lbls, pct,sep=" ") ##combines
pie(mytable, labels = lbls,
    main="Class 1 Amino Acid breakdown")

##pie chart of AAs wither percentages of total IDd and quantified CLASS 1  FROM PHOSPHO1 TABLE!
mytable <- table(phospho1$Amino.acid)
lbls <- paste(names(mytable), mytable, sep=" ")##pastes the labels and numbers from the table
pct <- round(mytable/(sum(mytable)),3)*100 ##calculates percentages
pct <- paste0(pct,"%") ##adds % sign
pct <- paste("(",pct,")",sep="") ##adds parentheses
lbls <- paste(lbls, pct,sep=" ") ##combines
pie(mytable, labels = lbls,
    main="Class 1 Quantified Amino Acid breakdown")


##pie chart of multiplicity with percentages
mytable <- table(multExpanded$multiplicity)
lbls <- paste(names(mytable), mytable, sep=" ")##pastes the labels and numbers from the table
pct <- round(mytable/(sum(mytable)),3)*100 ##calculates percentages
pct <- paste0(pct,"%") ##adds % sign
pct <- paste("(",pct,")",sep="") ##adds parentheses
lbls <- paste(lbls, pct,sep=" ") ##combines
pie(mytable, labels = lbls,
main="Peptide multiplicity states of class 1 quantifications")




############################################################################################################

# Class 1 phospho breakdowns by multiplicity;can be used later to subset analyses by confidence of quantification!!!!
z <- table(multExpanded$id,multExpanded$multiplicity) ##two dimensional table
colSums(z)
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

##make a venn diagram
install.packages("VennDiagram")
library(VennDiagram)


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
  cat.col = c("orange", "green", "blue"), margin = .1
);


#################################################################################################################
# Data Analysis. Histograms, Dendograms/Clustering/Heatmaps, PCA, DE, QQplots

##log2 the expression values (updated in main table)
multExpanded[newnames] <- log2(multExpanded[,newnames])##log2 transform

##a quick summary
summary(multExpanded[,newnames])

# Filter so that there is at least one valid value in each sample
data <- multExpanded[,newnames]

#replace row names with phospho id
row.names(data) <- multExpanded$id

# remove sites if not observed in each sample (not sure how to automate this for larger datasets)
data <- data[rowSums(is.na(data[ , 1:2])) < 2 & rowSums(is.na(data[ , 3:4])) < 2 & rowSums(is.na(data[ , 5:6])) < 2, ]                    

# note the dimensions here are 6071 where n>=3 for each phosphosite

#normalize by col median using columnwise (when to normalize??)

median.subtract <- function(x){ x - median(x, na.rm = TRUE)}##create a wrapper for median subtraction
datanorm <- colwise(median.subtract, newnames)(data) #create median subtracted data but loose the row names here...

row.names(datanorm) <- row.names(data)##add back the row names

# Histograms of data

##################################################### LIMMA for DE #####################################################
#Biological replication is needed for a valid comparison 

library(limma)


install.packages("statmod")
library(statmod)

# Calculate the correlation between technical replicates?...
# biolrep <- c(1, 1, 2, 2, 3, 3) 
# corfit <- duplicateCorrelation(datanorm, ndups = 1, block = biolrep)


# Produce dataframe from sample means ignoring missing data

HL16770 <- rowMeans(datanorm[,1:2], na.rm = T)
HL16778 <- rowMeans(datanorm[,3:4], na.rm = T)
HL16788 <- rowMeans(datanorm[,5:6], na.rm = T)

pilot <- cbind(HL16770, HL16778, HL16788)



#Produce the design matrix

fac <- factor(rep(1:3, each = 1))##codes the grouping for the ttests
design <- model.matrix(~0 + fac)
dnames <- levels(as.factor(substr(colnames(datanorm), 1, 7))) ##check me out. use 5 digit exp name.
colnames(design) <- dnames

#limma fit 
fit <- lmFit(pilot, design)

#Now to make all pairwise comparisons (from Smyth pg 14)
contrast.matrix <- makeContrasts(HL16778-HL16770, HL16788-HL16778, HL16788-HL16770, levels = design) 
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)





#Look at pairwise DE using toptable and the coef parameter to id which genes you are interested in 
topTable(fit2, coef = 1, adjust = "fdr")















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

cols <- as.factor(substr(colnames(datanorm), 3, 7))##check me out. use 5 digit exp name.
plot(pc$x[, 1], pc$x[, 2], col=as.numeric(cols), main = "PCA", xlab = "PC1", ylab = "PC2")
legend("bottomleft", levels(cols), col = seq(along=levels(cols)), pch = 1)


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

##other shit
#same as above see the quick R page for frequencies and crosstabs
margin.table(z,2) ##marginal frequencies across the 2nd (column) dimension
prop.table(z,2)  ##proportion across the 2nd dimension (the first dimension is quite interesting)

mytable <- xtabs(~A+B+c, data=mydata)
q <- xtabs(~id+multiplicity, data=multExpanded)##this is the same as table
summary(q)##chi-squared test of independence



multbd <- xtabs(id ~ multiplicity, data = multExpanded)

table(multExpanded$id)

quantile(multExpanded$multiplicity)



v <- as.factor(multExpanded$multiplicity)

str(v)
str(multExpanded$multiplicity)






multExpanded <- multExpanded[complete.cases(multExpanded),]#removes rows with ANY NAs



#hit 






























































##gives a dataframe but adds 'NaN' to summarized column and adds NA to other columns as well (NaN vs NA in original data). is way too long

sample <- testnew2[testnew2$id==42,25:31]
sample <- cbind(rep(42),sample)

colnames(sample)[1] <- "id" ##change a particular column name using colnames!!!


ddply(sample, "id", summarize)

ddply(sample, .(id, multiplicity), summarize, mean = mean(`16770_1`,na.rm = T))

##gives index
data <- grep("_", colnames(sample))

##gives string
data2 <- colnames(sample)[grep("_", colnames(sample))]

ddply(sample, .(id, multiplicity), colwise(mean),na.rm=T)##works

ddply(sample, .(id, multiplicity), colwise(mean,.(data2)),na.rm=T)##should work if underscore not there. replace with camelcase? whatever lets paste a letter in front of the character vector and column names?

newnames <- paste0("HL",data2)

##rename to proper alpha first
colnames(sample)[data] <- newnames


## columnwise summary time!...
ddply(sample, .(id, multiplicity), colwise(mean,newnames,na.rm=T))












# repeat each row of other data three times?
df[rep(seq_len(nrow(df)), each=2),]

otherdata3 <- other_data[rep(seq_len(nrow(other_data)), each=3),]

test <- merge(other_data, out, by="id")

colnames(otherdata3)






out2 <- ddply(testnew2, .(id, multiplicity), colwise(mean,newnames,na.rm=T))


length(which(duplicated(out)))
length(which(duplicated(testnew2)))


# now must merge with the rest of the dataframe...

merged <- merge(testnew2,out,by="id",all = TRUE)


try <- aggregate(. ~ id,
          data=merge(testnew2, out, by="id", all=TRUE), # Merged data, including NAs
          na.action=na.pass,              # Aggregate rows with missing values...
          FUN=sum, na.rm=TRUE)            # ...but instruct "sum" to ignore them.





# adjacency[,j] <- ifelse((!is.na(peptides[i]) & !is.na(peptides[i+1]) == "TRUE"),1,0)## winner winner. note && vs &!


a <- testnew2[(!is.na(testnew2$'16770_2') & !is.na(testnew2$'16770_1') == "TRUE"),] ##gives nothing

b <- testnew2[(!is.na(testnew2$multiplicity) & !is.na(testnew2$position) == "TRUE"),] ##gives nothing




##thoughts about a unique ID from which to parse DF
> any(duplicated(test$Modified.sequence))
[1] TRUE
> any(duplicated(phospho2$Modified.sequence))
[1] TRUE




##split this column by multiplicity
test2 <- colsplit(test$sample_rep_mult, pattern = "[0-9]{5}", c("a","b"))
##cast
test2 <- cbind(test, colsplit(test$sample_rep_mult, pattern = "___", c("sample","multiplicity")))

tail(test2$multiplicity)
test3 <- dcast(test2, multiplicity ~ var)
     
# potentially helpful 
# Here I split (that is add a variable identifier) the melted 'variable' column so that the sample, replicate, and multiplicity are now explicit
# for each "variable"
test5 <-   colsplit(names(expression), c("Rep"), c("sample","replicate_mult")) ##first split
test5 <- cbind(test5, colsplit(test5$replicate_mult, "___", c("replicate","multiplicity"))) ##second split
test5$sample <- gsub(test5$sample,pattern = "Ratio.H.L.normalized.", replacement = "") ##remove redundant information next 3 lines
drop <- "replicate_mult" 
test5 <- test5[,!(names(test5) %in% drop)]
test5

##barcharts of class 1 peptides 

# 1) of total class 1 per LC-MS run
# 2) of total class 1 per biolgical sample
# 3) 



