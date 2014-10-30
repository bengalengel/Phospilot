# Read in phospho table. Note the quote option is key.
phospho <- read.table("./MQ output/7_29_output/Phospho (STY)Sites.txt", sep = "\t", header=T, fill = T, quote = "")

# subset those hits that are not contaminants and not reverse hits.

phospho <- phospho[(phospho$Potential.contaminant != "+" & phospho$Reverse != "+"),]## here I am specifying the rows using column variables and logical operators while not specifying any particular columns
# I can come back here to get phospho ID information if needed


# subset to class 1
phospho1 <- phospho[(phospho$Localization.prob >= .75),]##Why just one localization probability?...


vars <- c("id","Amino.acid","Charge","Reverse","Potential.contaminant","Proteins","Positions.within.proteins","Leading.proteins",
          "Sequence.window","Phospho..STY..Probabilities","Localization.prob","PEP", "Score", "Delta.score", "Score.for.localization", 
          "Mass.error..ppm.", "Intensity", "Intensity.L", "Intensity.H", "Position", "Number.of.Phospho..STY.", 
          "Protein.group.IDs", "Protein")

other_data <- phospho1[,vars]

##dataframe that collects only the relevent expression columns. NOTE THE NEED TO USE REP!!!!!
##The sample number precedes 'Rep' (technical replicate) and the triple underscore denotes the multiplicity 

expression <- phospho1[,grep("Ratio.H.L.normalized(.*)_[12]___", colnames(phospho1))]

expression2 <- phospho1[,grep("Ratio.H.L.[0-9]+_[12]___", colnames(phospho1))]


names(expression) <- sub(names(expression), pattern ="_", replacement = "Rep")


##combine the two
phospho2 <- cbind(expression2, expression,other_data)

boxplot(expression2)


# They are not median centered. I am not sure what normalization means here
nonnorm <- log2(expression2)
norm <- log2(expression)

summary(nonnorm)
summary(norm)

names(nonnorm) <- sub(names(expression), pattern ="_", replacement = "Rep")
names(norm) <- sub(names(expression), pattern ="_", replacement = "Rep")


# melt to see if the normalization was performed across all multiplicities
require(plyr)
require(reshape2)

melted <- melt(nonnorm)

melted <- cbind(melted, colsplit(melted$variable, "Rep", c("sample", "replicate_mult"))) ##first split
melted <- cbind(melted, colsplit(melted$replicate_mult, "___", c("replicate","multiplicity"))) ##second split
melted$sample <- gsub(melted$sample, pattern = "Ratio.H.L.", replacement = "") ##remove redundant information next 3 lines
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


