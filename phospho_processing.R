##barcharts of class 1 peptides 

# 1) of total class 1 per LC-MS run
# 2) of total class 1 per biolgical sample
# 3) 




# First reperform all preseus processing steps using plyr and related tools.

# Read in phospho table. Note the quote option is key.
phospho <- read.table("./MQ output/7_17_output/Phospho (STY)Sites.txt", sep = "\t", header=T, fill = T, quote = "")

# subset those hits that are not contaminants and not reverse hits.

phospho <- phospho[(phospho$Contaminant != "+" & phospho$Reverse != "+"),]## here I am specifying the rows using column variables and logical operators while not specifying any particular columns
# I can come back here to get phospho ID information if needed


# subset to class 1
phospho1 <- phospho[(phospho$Localization.prob >= .75),]##Why just one localization probability?...

#####now must condense DF and 'melt' the dataframe so that each ratio for each multiplicity state has its own observation

##first condense dataframe to include variables of interest. colect non-expression variables first
other_data <- phospho1[,c("Amino.acid","Charge","Reverse","Contaminant","Proteins","Positions.within.proteins","Leading.proteins","Sequence.window",
                    "Modified.sequence","Localization.prob","PEP", "Score", "Delta.score", "Score.for.localization", "m.z", "Mass.error..ppm.",
                    "Intensity", "Intensity.L", "Intensity.H", "Position", "Number.of.Phospho..STY.", "Protein.group.IDs") ]

##dataframe that collects only the relevent expression columns. NOTE THE NEED TO USE REP!!!!!
##The sample number precedes 'Rep' (technical replicate) and the triple underscore denotes the multiplicity 
expression <- phospho1[,grep("Ratio.H.L.normalized(.*)Rep.___", colnames(phospho1))]

##combine the two
phospho2 <- cbind(expression,other_data)

##now to ensure each multiplicity has its own row and variables are condensed. I am converting from wide to long format to ensure that each
# observation is uniquely represented in the data

# Use reshape2 'melt' function 
library(reshape2)
library(stringr)
##lowercase names if I want
# names(phospho2) <- tolower(names(phospho2))

##melt all expression variables into one column
names(expression)
#test <- melt(phospho2, measure.vars = names(expression), variable.name = "sample_rep_mult", value.name = "H/L Ratio")

test <- melt(phospho2, measure.vars = names(expression))

# Here I split (that is add a variable identifier) the melted 'variable' column so that the sample, replicate, and multiplicity are now explicit
# for each "variable"
test <- cbind(test, colsplit(test$variable, "Rep", c("sample", "replicate_mult"))) ##first split
test <- cbind(test, colsplit(test$replicate_mult, "___", c("replicate","multiplicity"))) ##second split
test$sample <- gsub(test$sample, pattern = "Ratio.H.L.normalized.", replacement = "") ##remove redundant information next 3 lines
drop <- "replicate_mult" 
test <- test[,!(names(test) %in% drop)]


##cast data so that each unique 'sample/replicate' combination has its own column populated by the measurement 'currently the value column'.

testnew <- dcast(test, ... ~ variable)##will recast the entire data frame

testnew1 <- dcast(test, ... ~ sample, value.var="value") ##casts by sample

testnew2 <- dcast(test, ... ~ sample + replicate, value.var="value") ##close




Name, value.var="Value"





test5 <-   colsplit(names(expression), c("Rep"), c("sample","replicate_mult")) ##first split
test5 <- cbind(test5, colsplit(test5$replicate_mult, "___", c("replicate","multiplicity"))) ##second split
test5$sample <- gsub(test5$sample,pattern = "Ratio.H.L.normalized.", replacement = "") ##remove redundant information next 3 lines
drop <- "replicate_mult" 
test5 <- test5[,!(names(test5) %in% drop)]
test5



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




