##barcharts of class 1 peptides 

# 1) of total class 1 per LC-MS run
# 2) of total class 1 per biolgical sample
# 3) 




# First repreform all preseus processing steps using plyr and related tools.

# Read in phospho table. Note the quote option is key.
phospho <- read.table("./MQ output/7_17_output/Phospho (STY)Sites.txt", sep = "\t", header=T, fill = T, quote = "")

# subset those hits that are not contaminants and not reverse hits.

phospho <- phospho[(phospho$Contaminant != "+" & phospho$Reverse != "+"),]## here I am specifying the rows using column variables and logical operators while not specifying any particular columns
# I can come back here to get phospho ID information if needed


# subset to class 1
phospho1 <- phospho[(phospho$Localization.prob >= .75),]##Why just one localization probability?

#####now must condense DF and 'melt' the dataframe so that each ratio for each multiplicity state has its own observation

##first condense
test <- phospho1[,c("Ratio.H.L.normalized___1","Ratio.H.L.normalized___2") ]

##dataframe that collects only the relevent expression columns. NOTE THE NEED TO USE REP!!!!!
##The sample number precedes 'Rep' (technical replicate) and the triple underscore denotes the multiplicity 
test <- phospho1[,grep("Ratio.H.L.normalized(.*)Rep.___", colnames(phospho1))]


head(test)
str(test)
