##zia protein workup. goal is to get protein measurements for all 60 human samples. Normalize them against each other (median then quantile normalized) and then subtract (log scale) the relative protein concentration measurements from the relative phospho measurements. Issues will be in isoform quantification and batch correction across the protein measurements.

#these measurements will then be subjected to the same workflow as the confounded phosphomeasurements. 

rm(list=ls(all=TRUE)) #start with empty workspace if desired


# First perform all processing steps using plyr and related tools.
# load required libraries
library(reshape2)
library(stringr)
library(plyr)
require(limma)
require(sva)
require(swamp)
require(statmod)
source("loadMQZ.R")
#source("ExpandPhos.R")
#source("counts.R")
#source("breakdown.R")

# load protein files with particular variables populated using "loadMQ"
protein <- load.MQZ(directory = "D:/November Zia MBR MQ analysis/txt/")#7710 protein groups/81 variables

# load protein files with particular variables populated using "loadMQ" at home
protein <- load.MQZ(directory = "E:/My Documents/Pilot/November Zia MBR MQ analysis/txt/")#7,710 protein groups


# remove contaminants and reverse database hits
protein <- protein[(protein$Potential.contaminant != "+" & protein$Reverse != "+"),]#7354/78


# "only identified by site" hits CAN BE removed because they tend to have lower PEPs (wouldn't pass the FDR TH anyway) and can't be quantified since they are not idd by non-modified peptides. 
# Note there are some high probability proteins here given some proteins are idd by 20+ phosphopeptides.
# eg is A6NKT7 (PEP = 2.23E-70)
protein1 <- protein[(protein$Only.identified.by.site != "+"),]#6705 (really getting down there...)


colnames(protein1)<- gsub(colnames(protein1), pattern = "Ratio.H.L.normalized.", replacement = "HL") ##remove redundant information 

#some strangeness sometimes there are two extra rows!
#protein1 <- protein1[2:length(protein1)]

#remove proteins if not quantified in at least one sample
expCol <- grep("HL(.*)", colnames(protein1))

protein1 <- protein1[rowSums(is.na(protein1[,expCol]))!=length(expCol),]##removes rows containing all 
#'NA's using the sums of the logical per row 
#6421

data <- protein1[,expCol]#60 cell lines

row.names(data) <- protein1$id
data <- 1/data
data <- log2(data)

#inverse because I want light heavy here because the light sample is the standard 


#18871 is all over the place. actually its 18871?
boxplot(data)
which.min(sapply(data,median,na.rm=T))

#I am going to remove it for now
drop <- names(which.min(sapply(data,median,na.rm=T)))
data <- data[,!(names(data) %in% drop)]

#median normalize
names <- colnames(data)
median.subtract <- function(x){ x - median(x, na.rm = TRUE)}##create a wrapper for median subtraction
data <- colwise(median.subtract, names)(data) #create median subtracted data but loose intensity and the row names here...

#add back protien ids
row.names(data) <- protein1$id

#summaries
summary(data)
boxplot(data)#

#subset data to only the samples of interest
data <- data[,c("HL18862","HL18486","HL19160")]

#remove if protein group not found in all samples
data <- na.omit(data)#4270
boxplot(data)#differences in distribution shape for sure with HL18486 and HL19160
par(mfrow = c(1, 1))
for (i in 1:(ncol(data))){
  if(i==1) plot(density(data[, i], na.rm=T), col = i, ylim = c(0,2))
  else lines(density(data[, i], na.rm=T), col = i)
}


#quantile normalize data being compared
quantiled <- normalizeQuantiles(data,ties = T)#ties are all assigned the same value for the common quantile
summary(quantiled)
boxplot(data)
boxplot(quantiled)
# density plots all look the same
plot.new()
par(mfrow = c(1, 1))
for (i in 1:(ncol(quantiled))){
  if(i==1) plot(density(quantiled[, i], na.rm=T), col = i, ylim = c(0,1.9))
  else lines(density(quantiled[, i], na.rm=T), col = i)
}

##now to try some SVA. I need to create the model matrix (adjustment variables and variables of interest) and the null model matrix (only adjustment variables).
# mod <- model.matrix(~levels(as.factor(colnames(quantiled))), data=quantiled)
# colnames(mod) <- levels(as.factor(colnames(quantiled)))
# 
# tmp <- data.frame(x=c(1,1,1))
# mod0 <- model.matrix(~1, data=tmp)#only an intercept is included since we are not adjusting for any other variables...
# 
# #how many 'latent factors' are present in the protein data?
# num.sv(quantiled,mod,method = "be")
# 
# svobj = sva(quantile,mod,mod0)
# 
# 
# 
# fac <- factor(c(1,1,2,2,3,3))##codes the grouping for the ttests
#   design <- model.matrix(~0 + fac)
#   dnames <- levels(as.factor(substr(colnames(pilot), 1, 7))) ##check me out. use 5 digit exp name.
#   colnames(design) <- dnames

##normalized protein ratios from the three samples of interest
#quantiled <- quantiled[,c("HL18862","HL18486","HL19160")]

#and with quants in all three
#datacomp <- na.omit(quantiled)#3,925 now this is 4270? Confirm from home
#boxplot(datacomp)

#add an indicator to protein1 and then subset the gene names or protein names or whatever to get a count of of how many proteins I can use to normalize the phosphosites with...

#subset protein1 by id
tmp <- protein1[protein1$id %in% row.names(quantiled),]

##add the gene names etc of interest

datacomp <- cbind(quantiled,tmp[c("Protein.IDs","Majority.protein.IDs","Protein.names","Gene.names","Sequence.coverage....",
                                 "Number.of.proteins", "Peptides", "Razor...unique.peptides", "Unique.peptides",
                                 "Razor...unique.peptides.18862", "Razor...unique.peptides.18486", 
                                 "Razor...unique.peptides.19160")])


#Protein quants need to be assigned to phosphosites.
#*******************************************************************************************************************

#now how many unique proteins in the phospho data are subjected to DE analysis? 

#LOAD FILE ONLY IF STARTING HERE
#subset of multExpanded that are subjected to DE (load a local file with DE information)
multExpanded1_withDE <- read.csv("multExpanded1_withDE.csv", header=T)

#names(multExpanded1_withDE) <- sub("X.","",names(multExpanded1_withDE))#replace Xs
#names(multExpanded1_withDE) <- sub(".$","",names(multExpanded1_withDE))#replace dot at the end


SubtoDE <- multExpanded1_withDE[as.character(multExpanded1_withDE$SubtoDE) == "+",]

#DEany <- multExpanded1_withDE[multExpanded1_withDE$DEcont1=="+"|multExpanded1_withDE$DEcont2=="+"|multExpanded1_withDE$DEcont3=="+",]

#remove the reverse for now
SubtoDE <- SubtoDE[!grepl(SubtoDE$Protein, pattern = "REV"),]#the logical subsets the DF by row. 4991 observations.



# #remove the reverse for now
# SubtoDE$Protein <- as.character(SubtoDE$Protein)
# SubtoDE <- SubtoDE[SubtoDE$Protein[!grepl("REV",SubtoDE$Protein)],]#4991 observations but messes up id_mult. change to character results in all NAs.

#How many protiens subject to DE?
SubtoDEtable <- as.matrix(table(SubtoDE$Protein))
SubtoDEtable <- SubtoDEtable[SubtoDEtable!=0,,drop=F]#1991 proteins

#Remove isoform designation (this should be noted in the discussion of the paper. Unless there is a specific sequence aligning the peptide to a particular isoform the peptide could reasonably belong to any of the isoforms. Is there a designator in the table indicating if this is so?)

row.names(SubtoDEtable) <- substr(row.names(SubtoDEtable),1,6)#eventually 
SubtoDEproteins <- row.names(SubtoDEtable)
SubtoDEproteins <- unique(SubtoDEproteins)#1968 unique proteins (excluding isoforms) subjected to DE

#*******************************************************************************

#table for comparison with protein. 

# I need protein ID without isoform and ID_multiplicity for phosphopeptide. For matches with Ziaproteins, a new column will be added to that datatable, "norm_id_mult" that maps to the ids this protein group will normalize. For proteins that map to multiple ID_mults they will be added as a semicolon separated list. 

# This Ziaprotein table will be crossed with the original multexpanded table to identify any sites with matching ids. if there is a match, the three measurments of that protein group will be appended to the data table. Three new columns will then be populated containing the "normalized" measurements. 

#How many protiens subject to DE?
SubtoDEtable <- SubtoDE[c("Protein","idmult")]
SubtoDEtable$Protein <- substr(SubtoDEtable$Protein,1,6)#remove isoform designation






#****************************************************************************************************************
#How many proteins that are subjected to DE analysis are also IDd and quantified by proteomic analysis (Zia)? Here I will use the majority protein IDs. I will use the majority protein IDs for each protein group quantification (these proteins have at least half the peptides of the leading protein within the group)
Ziaproteins <- datacomp$Majority.protein.IDs

#unparse and turn into a long string for comparison
#now I need to design a loop to dig into the phospho file 
test <- strsplit(as.character(Ziaproteins), ";")
test <- as.character(unlist(test))
test <- unique(test)
any(duplicated(test))#no duplicates with 11500 proteins (now 12642 proteins with proper workflow)

test <- substr(test,1,6)#now I have some duplicates (3500 duplicates)
test <- unique(test) #Now I have 8885 unique protein with at least 1/2 of the peptides of the majority proteins within 4270 groups.

table(SubtoDEproteins%in%test)
#1206 of 1968 (61.2% are quantified in all three samples in Zia's work)



#how many of the overlapping proteins map to unique protein groups in the Zia dataset? (the same run gets around 40% of the protein groups ) I think this should be all of them.

#from Mann silac ratios 2010 paper:
# "For the phosphopeptides shared between multiple protein identifiers, the
# identifier with the maximum occupancy stoichiometry was used in the
# analysis."


#look for a uniuque identifier within a list of colon separated identifiers (grep?)
#row.names(DEtable)%in%as.character(Ziaproteins)#fix isoforms

# Ziaproteins2 <- substr(Ziaproteins,1,6)
# row.names(DEtable)%in%as.character(Ziaproteins2)#better but removes the other majority proteins
# table(row.names(DEtable)%in%as.character(Ziaproteins2))#

# Ziaproteins <- datacomp[c("Majority.protein.IDs"v)]
# Ziaproteins3 <- gsub("-.", "", Ziaproteins)#removes the isoform indicator
# Ziaproteins$id <- row.names(datacomp)
# Ziaproteins <- cbind(Ziaproteins,datacomp[c("Razor...unique.peptides", "Unique.peptides"             "Razor...unique.peptides.18862", "Razor...unique.peptides.18486", "Razor...unique.peptides.19160")]

Ziaproteins <- datacomp[c("Majority.protein.IDs","Razor...unique.peptides", "Unique.peptides", "Razor...unique.peptides.18862", "Razor...unique.peptides.18486", "Razor...unique.peptides.19160")]
Ziaproteins$id <- row.names(datacomp)
Ziaproteins$Majority.protein.IDs <- gsub("-.", "", Ziaproteins$Majority.protein.IDs)#removes the isoform indicator




# a quick for loop for each level of rownames that returns the number of hits and the protein ids that they match
# facttemp <- as.factor(row.names(SubtoDEtable))


## I need protein ID without isoform and ID_multiplicity from the phosphopeptide table (multexpanded). For matches with Ziaproteins, a new column will be added to that datatable, "norm_id_mult" that maps to the ids this protein group will normalize. For Ziaproteins that map to multiple phosphoobservations they will be added as a semicolon separated list. 

#for every protein linked to an id_mult from the phosphotable, a paired protein group from the ziaproteins table is found (if present) using any of the majority protein ids within that group. If the phospho id maps to multiple protein groups, the one with the most peptides is used.

#are there any problems with this? Nope

proteinindex <- c()
morethan1 <- c()#need to come back to this

for(i in seq_along(SubtoDEtable[,1])){
  tmp <- grep(SubtoDEtable$Protein[i], Ziaproteins$Majority.protein.IDs)
  #more than one value? This can happen with isoforms
  if(length(tmp)>1){
    #compare the razor plus unique count across the two matches and choose the one with the most   matches
    counts <- Ziaproteins$Razor...unique.peptides[tmp]
    proteinindex[i] <- tmp[which.max(counts)]
    morethan1 <- c(morethan1,tmp)
  }
  #if length of tmp >0 add to index
  if(length(tmp) == 1){
    proteinindex[i] <- tmp
  }
}

SubtoDE$ziaindex <- proteinindex

##for each index value add the three values and the majority ids. next for loop performs the normalization.


#declare empty dataframe with proper names

#declare empty data frame and note issues below..
protein_norm <- data.frame()
for(i in seq_along(proteinindex)){
  if(!is.na(proteinindex[i])){
    tmp <- datacomp[c("Majority.protein.IDs","HL18862", "HL18486", "HL19160")][proteinindex[i],]
    protein_norm <- rbind(protein_norm,tmp)
  }
  if(is.na(proteinindex[i])){
    tmp <- rep(NA,4)
    protein_norm <- rbind(protein_norm,tmp)
    if(dim(protein_norm)[1]==1){#reset names etc
      names(protein_norm) <- c("Majority.protein.IDs","HL18862", "HL18486", "HL19160")
      protein_norm$Majority.protein.IDs <- as.factor(protein_norm$Majority.protein.IDs)
      protein_norm$HL18862 <- as.numeric(protein_norm$HL18862)
      protein_norm$HL18486 <- as.numeric(protein_norm$HL18486)
      protein_norm$HL19160 <- as.numeric(protein_norm$HL19160)
    }
  }
}    

#***************************************************************************************************************88
#data frame for normalization

#phospho information. id, id_mult, protein and "pilot dataframe"
normphos <- SubtoDE[,c("id","idmult")]
normphos$PhosphoProtein <- SubtoDE$Protein

#add the pilot data after removing the REV proteins
rowsofbadness <- which(grepl(multExpanded1_withDE$Protein, pattern = "REV") & multExpanded1_withDE$SubtoDE == "+")
idmultsofbadness <- as.character(multExpanded1_withDE$idmult[rowsofbadness])

#subset out normalized batch corrected pilot dataframe
pilot2 <- pilot[!rownames(pilot) %in% idmultsofbadness,]

#combine the two
normphos <- cbind(normphos,pilot2)

#now add the protein data to make a table of normalized/processed data from both molecular phenotypes:
Phos_Protein <- cbind(normphos,protein_norm)

#write out this table
write.csv(Phos_Protein, "Phos_Protein.csv", row.names=F)

##now subset and normalize to make a final table for limma DE and all other downstream analysis
expCol <- grep("HL(.*)", colnames(Phos_Protein))
data <- Phos_Protein[,expCol]
row.names(data) <- Phos_Protein$idmult
data <- na.omit(data)

#perform the normalization
HL18486_1norm <- data$HL18486_1-data$HL18486
HL18486_2norm <- data$HL18486_2-data$HL18486
HL18862_1norm <- data$HL18862_1-data$HL18862
HL18862_2norm <- data$HL18862_2-data$HL18862
HL19160_1norm <- data$HL19160_1-data$HL19160
HL19160_2norm <- data$HL19160_2-data$HL19160

ProtNormalized <- cbind(HL18486_1norm, HL18486_2norm, HL18862_1norm, HL18862_2norm, HL19160_1norm, HL19160_2norm)
row.names(ProtNormalized) <- row.names(data)

boxplot(ProtNormalized)
summary(ProtNormalized)
par(mfrow = c(1, 1))
for (i in 1:(ncol(ProtNormalized))){
  if(i==1) plot(density(ProtNormalized[, i], na.rm=T), col = i, ylim = c(0,1.5))
  else lines(density(ProtNormalized[, i], na.rm=T), col = i)
}

#do these still cluster?..scary YES!!!

dataZ <- scale(ProtNormalized)##Z-scored column wise the complete data matrix

# now all data excepting complete cases (note that the sample dendograms look the same)
#hist(dataZ[,6], breaks = 100)

# dendogram using euclidian distance (default) and ward or complete agglomeration
dend.ward<- as.dendrogram(hclust(dist(t(dataZ)),method="ward"))
dend.complete<- as.dendrogram(hclust(dist(t(dataZ))))

ward.o<- order.dendrogram(dend.ward)
complete.o<- order.dendrogram(dend.complete)

plot(dend.complete,ylab="height", main = "Euclidian/Complete")
plot(dend.ward, leaflab = "perpendicular", ylab = "height", main = "Euclidian/Ward")



# row scaled
r <- t(scale(t(ProtNormalized)))#transpose to zscale the rows then transpose back to original format

# sample scaled
c <- scale(ProtNormalized)


# install heatmap.2 package
# install.packages("gplots")
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
  density.info="none",
  key.xlab = "Row Z scores", key.ylab=NULL, key.title = "",
  srtCol=45,  ,adjCol = c(1,1),
  margins = c(6,5),
  cexCol=1,
  labRow = NA#remove row labels
)


# plot.new()

#PCA analysis 
x <- t(ProtNormalized)#samples are the rows of the column matrix
pc <- prcomp(x)#scale = T, center = T) as of now I am not scaling

cols <- as.factor(substr(colnames(ProtNormalized), 3, 7))##check me out. use 5 digit exp name.
plot(pc$x[, 1], pc$x[, 2], col=as.numeric(cols), main = "PCA", xlab = "PC1", ylab = "PC2")
legend("bottomleft", levels(cols), col = seq(along=levels(cols)), pch = 1)


summary(pc)

#SVD for calculating variance explained; see Rafa's notes for an explaination
cx <- sweep(x, 2, colMeans(x), "-")
sv <- svd(cx)
names(sv)
plot(sv$u[, 1], sv$u[, 2], col = as.numeric(cols), main = "SVD", xlab = "U1", ylab = "U2")


plot(sv$d^2/sum(sv$d^2), xlim = c(1, 12), type = "b", pch = 16, xlab = "principal components", 
     ylab = "variance explained")


#***Differential Expression******************************************************

fac <- factor(c(1,1,2,2,3,3))##codes the grouping for the ttests
design <- model.matrix(~0 + fac)
dnames <- levels(as.factor(substr(colnames(ProtNormalized), 1, 7))) ##check me out. use 5 digit exp name.
colnames(design) <- dnames

fit <- lmFit(ProtNormalized, design)
# Now to make all pairwise comparisons (group2-1, group3-2, group3-1)
contrast.matrix <- makeContrasts(HL18862-HL18486, HL19160-HL18862, HL19160-HL18486, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

sig1 <- topTable(fit2, coef = 1, adjust = "BH", n=Inf, sort="p", p=.05)#sorts by adjusted p up to the threshold of .05, which is the default FDR chosen for differential expression ("results" function). This actually seems a conservative way to sort.
sig2 <- topTable(fit2, coef = 2, adjust = "BH", n=Inf, sort="p", p=.05)
sig3 <- topTable(fit2, coef = 3, adjust = "BH", n=Inf, sort="p", p=.05)

# sig1 - 18862-18486
# sig2 - 19160-18862
# sig3 - 19160-18486

c1up  <- sig1[sig1$logFC > 0,]
c1down <- sig1[sig1$logFC < 0,]
c2up <- sig2[sig2$logFC > 0,]
c2down <- sig2[sig2$logFC < 0,]
c3up <- sig3[sig3$logFC > 0,]
c3down <- sig3[sig3$logFC < 0,]


tt1 <- topTable(fit2, coef = 1, adjust = "BH", n=Inf)#sorts by adjusted p up to the threshold of .
tt2 <- topTable(fit2, coef = 2, adjust = "BH", n=Inf)#sorts by adjusted p up to the threshold of .
tt3 <- topTable(fit2, coef = 3, adjust = "BH", n=Inf)#sorts by adjusted p up to the threshold of .

hist(tt1$P.Value, nc=40, xlab="P values", main = colnames(contrast.matrix)[1])
hist(tt2$P.Value, nc=40, xlab="P values", main = colnames(contrast.matrix)[2])
hist(tt3$P.Value, nc=40, xlab="P values", main = colnames(contrast.matrix)[3])

plot(tt1$logFC,-log10(tt1$P.Value), xlab = colnames(contrast.matrix)[1], pch = 20, ylab = "-log10(P)",xlim = c(-6, 6))
#sites with sig difference in comparison 1
names <- row.names(sig1)
names2 <- row.names(tt1)
index <- which(names2 %in% names)
points(tt1$logFC[index],-log10(tt1$P.Value)[index], col="red3", pch = 20)


plot(tt2$logFC,-log10(tt2$P.Value), xlab = colnames(contrast.matrix)[2], pch = 20, ylab = "-log10(P)",xlim = c(-6, 6))
#sites with sig difference in comparison 1
names <- row.names(sig2)
names2 <- row.names(tt2)
index <- which(names2 %in% names)
points(tt2$logFC[index],-log10(tt2$P.Value)[index], col="red3", pch = 20)


plot(tt3$logFC,-log10(tt3$P.Value), xlab = colnames(contrast.matrix)[3], pch = 20, ylab = "-log10(P)",xlim = c(-6, 6))
#sites with sig difference in comparison 1
names <- row.names(sig3)
names2 <- row.names(tt3)
index <- which(names2 %in% names)
points(tt3$logFC[index],-log10(tt3$P.Value)[index], col="red3", pch = 20)




results <- decideTests(fit2, adjust.method = "BH", method = "separate")#results is a 'TestResults' matrix
#separate compares each sample individually and is the default approach
summary(results)


vennDiagram(results, cex=c(1.2,1,0.7)) #good DE across conditions
vennDiagram(results, cex=c(1.2,1,0.7), include = "up")
vennDiagram(results, cex=c(1.2,1,0.7), include = "down") 
vennDiagram(results, cex=c(1.2,1,0.7), include = c("up", "down")) #reproduces 'signature'


#F statistics
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



Fvals <- topTableF(fit2, adjust = "BH", n=Inf, sort.by="F")#all F values
sigFvals <- topTableF(fit2, adjust = "BH", n=Inf, sort.by="F", p=.05)#gives 1355 compared to 1549 DE total for separate comparisons

#subsets by contrast specific DE
FDE1 <- Fvals[which(row.names(DE1)%in%row.names(Fvals)),5:6]
FDE2 <- Fvals[which(row.names(DE2)%in%row.names(Fvals)),5:6]
FDE3 <- Fvals[which(row.names(DE3)%in%row.names(Fvals)),5:6]

#below gives adjusted pvalue
FDE1 <- Fvals[match(row.names(DE1), row.names(Fvals), nomatch = F),5:7]
FDE2 <- Fvals[match(row.names(DE2), row.names(Fvals), nomatch = F),5:7]
FDE3 <- Fvals[match(row.names(DE3), row.names(Fvals), nomatch = F),5:7]

#boxplot(FDE1$F,FDE2$F,FDE3$F)
boxplot(log10(FDE1$F),log10(FDE2$F),log10(FDE3$F))
summary(FDE1$F)
summary(FDE2$F)
summary(FDE3$F)

plot(density(log10(FDE1$F)),xlim = c(0,3))
lines(density(log10(FDE2$F)), col = 2)
lines(density(log10(FDE3$F)), col = 3)





##overlap with phospho dataset
#add annotation to multExpanded withDE and save the new file


#add annotation to multexpanded DF
multExpanded1_withDE$SubtoDEpn = ifelse(multExpanded1_withDE$idmult %in% row.names(ProtNormalized),"+","-")

#add F test values to the table
multExpanded1_withDE$globalFsigpn = ifelse(multExpanded1_withDE$idmult %in% row.names(sigFvals),"+","-")

#add DE to table
multExpanded1_withDE$DEcont1pn = ifelse(multExpanded1_withDE$idmult %in% row.names(sig1),"+","-")
multExpanded1_withDE$DEcont2pn = ifelse(multExpanded1_withDE$idmult %in% row.names(sig2),"+","-")
multExpanded1_withDE$DEcont3pn = ifelse(multExpanded1_withDE$idmult %in% row.names(sig3),"+","-")

#add DE direction to table
multExpanded1_withDE$cont1uppn = ifelse(multExpanded1_withDE$idmult %in% row.names(c1up),"+","-")
multExpanded1_withDE$cont1downpn = ifelse(multExpanded1_withDE$idmult %in% row.names(c1down),"+","-")
multExpanded1_withDE$cont2uppn = ifelse(multExpanded1_withDE$idmult %in% row.names(c2up),"+","-")
multExpanded1_withDE$cont2downpn = ifelse(multExpanded1_withDE$idmult %in% row.names(c2down),"+","-")
multExpanded1_withDE$cont3uppn = ifelse(multExpanded1_withDE$idmult %in% row.names(c3up),"+","-")
multExpanded1_withDE$cont3downpn = ifelse(multExpanded1_withDE$idmult %in% row.names(c3down),"+","-")

#write the multExpanded table with DE information
write.table(multExpanded1_withDE,"multExpanded1_withDE_protein.csv", sep=',',col.names=T, row.names=F)


##check for overlap using omnibus F

#subset
SubtoDEpn <- multExpanded1_withDE[multExpanded1_withDE$SubtoDEpn == "+",]

#for an apples to apples comparison I will need to subject the same subset to DE in both cases.

#limma on common subset for phosphodata

#subset the phospho data by id_mult found in the protein data
phosdata <- pilot[row.names(pilot) %in% SubtoDEpn$idmult,]

require(limma)
require(statmod)
#Produce the design matrix

fac <- factor(c(1,1,2,2,3,3))##codes the grouping for the ttests
design <- model.matrix(~0 + fac)
dnames <- levels(as.factor(substr(colnames(phosdata), 1, 7))) ##check me out. use 5 digit exp name.
colnames(design) <- dnames

fit <- lmFit(phosdata, design)

#Now to make all pairwise comparisons (group2-1, group3-2, group3-1)
contrast.matrix <- makeContrasts(HL18862-HL18486, HL19160-HL18862, HL19160-HL18486, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

Fvals <- topTableF(fit2, adjust = "BH", n=Inf, sort.by="F")#all F values
sigFvals <- topTableF(fit2, adjust = "BH", n=Inf, sort.by="F", p=.05)#gives 1355 compared to 1549 DE total for separate comparisons

diffphos <- nrow(phosdata[row.names(phosdata) %in% row.names(sigFvals),])
diffphosnorm <- nrow(SubtoDEpn[SubtoDEpn$globalFsigpn == "+",])

##append to SubtoDEpn
SubtoDEpn$diffphos = ifelse(SubtoDEpn$idmult %in% row.names(sigFvals),"+","-")

intersection <- nrow(SubtoDEpn[SubtoDEpn$globalFsigpn == "+" & SubtoDEpn$diffphos == "+",])


#make a double venn
require(VennDiagram)
require(gridExtra)
plot.new()
venn.plot <- draw.pairwise.venn(
  area1 = nrow(SubtoDEpn[SubtoDEpn$globalFsig == "+",]),
  area2 = nrow(SubtoDEpn[SubtoDEpn$globalFsigpn == "+",]),
  cross.area = nrow(SubtoDEpn[SubtoDEpn$globalFsigpn == "+" & SubtoDEpn$globalFsig == "+",]),
  category = c("PhosDE", "normPhosDE"),
  fill = c("green", "blue"),
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.col = c("green", "blue"), 
  margin = .1,
  main="test"
)
plot.new()
grid.arrange(gTree(children=venn.plot), main="Differential Phosphorylation Overlap")


##almost all of the phosDE is picked up. But there is just as many new DE in the non-confounded data!...
##better now that the inverse is used!!!

