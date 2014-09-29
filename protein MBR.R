# First perform all processing steps using plyr and related tools.
# load required libraries
library(reshape2)
library(stringr)
library(plyr)


# Read in phospho table. Note the quote option is key.
protein <- read.table("./MQ output/7_29_output/ProteinGroups.txt", sep = "\t", header=T, fill = T, quote = "")

# subset those hits that are not contaminants and not reverse.
protein <- protein[(protein$Potential.contaminant != "+" & protein$Reverse != "+"),]## here I am specifying the rows using column variables and logical operators while not specifying any particular columns

# "only identified by site" hits are removed because they tend to have lower PEPs (wouldn't pass the FDR TH anyway) and can't be
# quantified since they are not idd by non-modified peptides. Note there are some high probability proteins here given some proteins
# are idd by 20+ phosphopeptides. eg is A6NKT7 (PEP = 2.23E-70)
protein2 <- protein[(protein$Only.identified.by.site != "+"),]



# Produce a trimmed matrix like before for quantification. Note unsure what's going on with sequence coverage
# and unique + razore variations. 


vars <- c("id", "Protein.IDs", "Majority.protein.IDs",  "Protein.names", "Gene.names", "Number.of.proteins", "Peptides", "Razor...unique.peptides", "Unique.peptides", "Sequence.coverage....", "Mol..weight..kDa.", "Sequence.length", "PEP", "Peptide.IDs", "Mod..peptide.IDs", "Evidence.IDs")

other_data <- protein2[,vars]

##dataframe that collects only the relevent expression columns. NOTE THE NEED TO USE REP!!!!!
##The sample number precedes 'Rep' (technical replicate) and the triple underscore denotes the multiplicity 
expression <- protein2[,grep("Ratio.H.L.normalized(.*)_[12]___", colnames(protein2))]

# Replace the column names
names(expression) <- sub(names(expression), pattern ="_", replacement = "Rep")

##combine the two
protein2 <- cbind(expression,other_data)
