ProtAssignment <- function(proteinfull, proteinnorm, multExpanded1_withDE, phosphonorm){
##This program assigns protein groups (from previous proteomic analysis) to phosphosites from the SCX-TiO2 workflow for normalization. 
  
##proteinfull is the 60 estimates from zias work. proteinnorm is the quantile normalized protein data.  phosphonorm is the normalized/batch corrected /confounded phospho data.multExpanded1_withDE is the parent dataframe for the class one sites with DiffPhos annotation from using the confounded data.  

#subset proteinfull by protein group id found in the three normalized samples of interest. 
proteinfull <- proteinfull[proteinfull$id %in% row.names(proteinnorm),]

#revert names of proteinnorm back to 'HL' for continuity with the structure below. (Changed back at the end)
names(proteinnorm) <- gsub(names(proteinnorm), pattern = "LH", replacement = "HL")

##combine the normalized protein information and the annotation data into a common dataframe for matching to the phospho data
datacomp <- cbind(proteinnorm,proteinfull[c("Protein.IDs","Majority.protein.IDs","Protein.names","Gene.names","Sequence.coverage....",
                                  "Number.of.proteins", "Peptides", "Razor...unique.peptides", "Unique.peptides",
                                  "Razor...unique.peptides.18862", "Razor...unique.peptides.18486", 
                                  "Razor...unique.peptides.19160")])

#Descriptive
#how many unique proteins in the phospho data are identified and subjected to DE analysis? Below are counts with and without isoform designation.
######################
SubtoDE <- multExpanded1_withDE[as.character(multExpanded1_withDE$SubtoDE) == "+",]

#remove the reverse (ME should be passed this way) and add 'AllPhos', a data frame for all phosphopeptides.
SubtoDE <- SubtoDE[!grepl(SubtoDE$Protein, pattern = "REV"),]#the logical subsets the DF by row. 4991 observations.
AllPhos <- multExpanded1_withDE[!grepl(multExpanded1_withDE$Protein, pattern = "REV"),] #18238
  
#How many 'protiens' subject to DE? These are majority protein ids from the protein group assigned to the phosphopeptide.
SubtoDEtable <- as.matrix(table(SubtoDE$Protein))
SubtoDEtable <- SubtoDEtable[SubtoDEtable!=0,,drop=F]#1991 proteins

#How many 'proteins' in the entire class1 dataset?
AllPhostable <- as.matrix(table(AllPhos$Protein))
AllPhostable <- AllPhostable[AllPhostable!=0,,drop=F]#4184 proteins

#Remove isoform designation (this should be noted in the discussion of the paper. Unless there is a specific sequence aligning the peptide to a particular isoform the peptide could reasonably belong to any of the isoforms. This workflow (like those that use only the annotated uniprot databases) assumes this is not the case. Only a minor issue using this approach anyway 1991 to 1968 and 4184 to 4127.

row.names(SubtoDEtable) <- substr(row.names(SubtoDEtable),1,6)#eventually 
SubtoDEproteins <- row.names(SubtoDEtable)#1991
SubtoDEproteins <- unique(SubtoDEproteins)#1968 unique proteins (excluding isoforms) subjected to DE

row.names(AllPhostable) <- substr(row.names(AllPhostable),1,6)#eventually 
AllPhosproteins <- row.names(AllPhostable)#4184
AllPhosproteins <- unique(AllPhosproteins)#4127 unique proteins (excluding isoforms) subjected to DE
###############

#Descriptive
#How many proteins subjected to DE analysis are also IDd and quantified by proteomic analysis (Zia)? 
###############################################
#Here I will use the majority protein IDs. I will use the majority protein IDs for each protein group quantification (these proteins have at least half the peptides of the leading protein within the group)
Ziaproteins <- datacomp$Majority.protein.IDs

#unparse and turn into a long string for comparison
#now I need to design a loop to dig into the phospho file 
pcount <- strsplit(as.character(Ziaproteins), ";")
pcount <- as.character(unlist(pcount))
pcount <- unique(pcount)
any(duplicated(pcount))#12642 proteins identified with at least 1/2 of the peptides of the majority protein(s) in Zias work within 4270 protein groups

pcount <- substr(pcount,1,6)#now I have some duplicates due to isoform designation (3500 duplicates)
pcount <- unique(pcount) #Now I have 8885 unique protein with at least 1/2 of the peptides of the majority protein(s) within 4270 groups.

table(SubtoDEproteins%in%pcount)#isoform free
#1206 of 1968  sub to diffphos (61.2% are quantified in all three samples in Zia's work)

table(AllPhosproteins%in%pcount)
#2117 of 4127 (51% of proteins)
#####################

#Assigning protein group quantifications to phosphopeptide quantifications
#####################################################################
Ziaproteins <- datacomp[c("Majority.protein.IDs","Razor...unique.peptides", "Unique.peptides", "Razor...unique.peptides.18862", "Razor...unique.peptides.18486", "Razor...unique.peptides.19160")]
Ziaproteins$id <- row.names(datacomp)#this is a bit redundant
Ziaproteins$Majority.protein.IDs <- gsub("-.", "", Ziaproteins$Majority.protein.IDs)#removes the isoform indicator

#tables for comparison with protein with protein data. 
SubtoDEtable <- SubtoDE[c("Protein","idmult")]#note that reverse entries have already been removed above
SubtoDEtable$Protein <- substr(SubtoDEtable$Protein,1,6)#remove isoform designation
AllPhostable <- AllPhos[c("Protein","idmult")]#note that reverse entries have already been removed above
AllPhostable$Protein <- substr(AllPhostable$Protein,1,6)#remove isoform designation

#For every protein assigned to an id_mult from the phosphotable, a paired protein group from the ziaproteins table is found (if present) using any of the majority protein ids within that group. If the phospho id maps to multiple protein groups, the one with the most peptides is used.


#all phos loops. I can pre-allocate memory in this loop by defining the length of the poroteinindex
nPTable <- length(AllPhostable[,1])
proteinindex <- integer(nPTable)
#morethan1 <- c()#need to come back to this because this isn't work properly
for(i in seq_along(AllPhostable[,1])){#for every protein assigned to a phosphopeptide ask is this protein found in the zia dataset? 
  tmp <- grep(AllPhostable$Protein[i], Ziaproteins$Majority.protein.IDs)#Searches all IDs in comma separated majority protein list.
  #more than one value? This can happen with isoforms
  if(length(tmp)>1){
    #compare the razor plus unique count across the two matches and choose the one with the most   matches
    counts <- Ziaproteins$Razor...unique.peptides[tmp]
    proteinindex[i] <- tmp[which.max(counts)]
    #morethan1 <- c(morethan1,tmp)
  }
  #if length of tmp >0 add to index
  if(length(tmp) == 1){
    proteinindex[i] <- tmp
  }
  else{
    proteinindex[i] <- NA
  }
}

AllPhostable$ziaindex <- proteinindex#add the index to the parent dataframe

##for each proteinindex value add the three protein quantification (H/L) values and the majority ids.
protein_norm2 <- data.frame()
for(i in seq_along(proteinindex)){
  if(!is.na(proteinindex[i])){
    tmp <- datacomp[c("Majority.protein.IDs","HL18862", "HL18486", "HL19160")][proteinindex[i],]#uses the quantile normalized values
    protein_norm2 <- rbind(protein_norm2,tmp)
  }
  if(is.na(proteinindex[i])){
    tmp <- rep(NA,4)
    protein_norm2 <- rbind(protein_norm2,tmp)
    if(dim(protein_norm2)[1]==1){#reset names of dataframe in the event the first loop is an NA
      names(protein_norm2) <- c("Majority.protein.IDs","HL18862", "HL18486", "HL19160")
      protein_norm2$Majority.protein.IDs <- as.factor(protein_norm2$Majority.protein.IDs)
      protein_norm2$HL18862 <- as.numeric(protein_norm2$HL18862)
      protein_norm2$HL18486 <- as.numeric(protein_norm2$HL18486)
      protein_norm2$HL19160 <- as.numeric(protein_norm2$HL19160)
    }
  }
}    
#link the protein quants to the phospho ids to make a dataframe with normalized protein quants appended. Note "REV_" entries are removed.
AllPhos <- cbind(AllPhos,protein_norm2)
colnames(AllPhos)[76] <- "ProtPrep Majority Protein IDs"
#change to proper designation
names(AllPhos)[77:79] <- gsub(names(AllPhos)[77:79], pattern = "HL", replacement = "LH")
#AllPhos is to be returned
#############################

#Normalize the passed phospho dataframe for return and produce EDA plots on normalized dataframe.
################################################################
#subset phosphonorm s.t. all phosphopeptides were mapped to a protein in Zia's data

#first generate the ids for subsetting
normphos <- AllPhos[,c("idmult", "Protein", "Leading.proteins")]
normphos <- cbind(normphos,protein_norm2)
colnames(normphos)[4] <- "ProtPrep Majority Protein IDs"
MappedPhos <- na.omit(normphos)#11641 phosphopeptides are mapped to a protein in zia's dataset

#subset phosphonorm
phosphonorm <- phosphonorm[row.names(phosphonorm) %in% MappedPhos$idmult, ]

#now subset the MappedPhospho
MappedPhos <- MappedPhos[MappedPhos$idmult %in% row.names(phosphonorm), ]

#combine. I should change back to L/H to prevent confusion later on
phosphonorm <- cbind(phosphonorm,MappedPhos)

#normalize and subject to EDA and varcomp. protein will be run as a covariate for DiffPhos.

##now subset and normalize
expCol <- grep("HL(.*)", colnames(phosphonorm))
data <- phosphonorm[,expCol]
row.names(data) <- phosphonorm$idmult

#perform the normalization
HL18486_1_1norm <- data$HL18486_1_1-data$HL18486
HL18486_1_2norm <- data$HL18486_1_2-data$HL18486
HL18486_2_1norm <- data$HL18486_2_1-data$HL18486
HL18486_2_2norm <- data$HL18486_2_2-data$HL18486
HL18862_1_1norm <- data$HL18862_1_1-data$HL18862
HL18862_1_2norm <- data$HL18862_1_2-data$HL18862
HL18862_2_1norm <- data$HL18862_2_1-data$HL18862
HL18862_2_2norm <- data$HL18862_2_2-data$HL18862
HL19160_1_1norm <- data$HL19160_1_1-data$HL19160
HL19160_1_2norm <- data$HL19160_1_2-data$HL19160
HL19160_2_1norm <- data$HL19160_2_1-data$HL19160
HL19160_2_2norm <- data$HL19160_2_2-data$HL19160

ProtNormalized <- cbind(HL18486_1_1norm, HL18486_1_2norm, HL18486_2_1norm, HL18486_2_2norm, HL18862_1_1norm, HL18862_1_2norm, HL18862_2_1norm,
                        HL18862_2_2norm, HL19160_1_1norm, HL19160_1_2norm, HL19160_2_1norm, HL19160_2_2norm)
row.names(ProtNormalized) <- row.names(data)
##############

#summary plots
###############
boxplot(ProtNormalized)
par(mfrow = c(1, 1))
for (i in 1:(ncol(ProtNormalized))){
  if(i==1) plot(density(ProtNormalized[, i], na.rm=T), col = i, ylim = c(0,1.5))
  else lines(density(ProtNormalized[, i], na.rm=T), col = i)
}

#clustering and EDA using complete cases
ProtNormalized2 <- na.omit(ProtNormalized)

dataZ <- scale(ProtNormalized2)##Z-scored column wise the complete data matrix

# dendogram using euclidian distance (default) and ward or complete agglomeration
dend.ward<- as.dendrogram(hclust(dist(t(dataZ)),method="ward"))
dend.complete<- as.dendrogram(hclust(dist(t(dataZ))))

ward.o<- order.dendrogram(dend.ward)
complete.o<- order.dendrogram(dend.complete)

plot(dend.complete,ylab="height", main = "Euclidian/Complete")
plot(dend.ward, leaflab = "perpendicular", ylab = "height", main = "Euclidian/Ward")

# row scaled
r <- t(scale(t(ProtNormalized2)))#transpose to zscale the rows then transpose back to original format

# sample scaled
c <- scale(ProtNormalized2)

require(gplots)

# Create dendrogram using the data without NAs
feature.dend<- as.dendrogram(hclust(dist(r),method="ward"))
sample.dend<- as.dendrogram(hclust(dist(t(c)),method="ward"))##note that dist caclculates distance between rows by default


##produce the heatmap. 
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
  margins = c(7,6),
  cexCol=1,
  labRow = NA#remove row labels
)
# plot.new()

#PCA analysis 
x <- t(ProtNormalized2)#samples are the rows of the column matrix
pc <- prcomp(x)#scale = T, center = T) as of now I am not scaling

cols <- as.factor(substr(colnames(ProtNormalized2), 3, 7))##use 5 digit exp name.
plot(pc$x[, 1], pc$x[, 2], col=as.numeric(cols), main = "PCA", xlab = "PC1", ylab = "PC2")
legend("bottomleft", levels(cols), col = seq(along=levels(cols)), pch = 1)


summary(pc)

#SVD for calculating variance explained; see Rafa's notes for an explaination
cx <- sweep(x, 2, colMeans(x), "-")
sv <- svd(cx)
names(sv)
#plot(sv$u[, 1], sv$u[, 2], col = as.numeric(cols), main = "SVD", xlab = "U1", ylab = "U2")
plot(sv$d^2/sum(sv$d^2), xlim = c(1, 12), type = "b", pch = 16, xlab = "principal components", 
     ylab = "variance explained")


DFs <- list(AllPhos,ProtNormalized)

return(DFs)
}