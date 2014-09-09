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
                          "Protein.group.IDs") ]

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

multExpanded <- multExpanded[rowSums(is.na(multExpanded[,expCol]))!=length(expCol),]##removes all 'NA' rows using the sums of the logical per row                        

##total class 1 sites and protein groups? 
# How are phosphosites attached to a protein?

##number of protein groups associated with class 1 sites



##number of unique class 1 sites
nrow(phospho1)


##number of unique class 1 sites with valid quantification
nrow(table(multExpanded$id))








##pie chart of AAs wither percentages
mytable <- table(multExpanded$Amino.acid)
lbls <- paste(names(mytable), mytable, sep=" ")##pastes the labels and numbers from the table
pct <- round(mytable/(sum(mytable)),3)*100 ##calculates percentages
pct <- paste0(pct,"%") ##adds % sign
pct <- paste("(",pct,")",sep="") ##adds parentheses
lbls <- paste(lbls, pct,sep=" ") ##combines
pie(mytable, labels = lbls,
    main="Amino Acid breakdow")


##pie chart of multiplicity with percentages
mytable <- table(multExpanded$multiplicity)
lbls <- paste(names(mytable), mytable, sep=" ")##pastes the labels and numbers from the table
pct <- round(mytable/(sum(mytable)),3)*100 ##calculates percentages
pct <- paste0(pct,"%") ##adds % sign
pct <- paste("(",pct,")",sep="") ##adds parentheses
lbls <- paste(lbls, pct,sep=" ") ##combines
pie(mytable, labels = lbls,
main="Number of Phosphorylation sites per peptide")





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



