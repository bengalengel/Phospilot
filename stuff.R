##phospho loop to assign protein ids, H/L values and ibaq values to phosphosites for normalization from phospho workup

#decision is when there are multiple protein groups that could match the phosphosite. Here I must choose the one with the most unique/razor ids

# also want to remove names from the majority protein ids that couldn't possible contain the phosphopeptide in the future

#use sapply when possible. I will work from me1 dataframe

#collect ids from multexpanded
ids <- multExpanded1$Protein.group.IDs
ids <- as.character(ids)

#for those ids that contain multiple groups. pick the one with the most unique plus razor 
index <- grep(";", ids)
suspect <- ids[index]


x <- suspect[1]
y <- as.numeric(unlist(strsplit(x,split = ";")))
#does either id match a protein group id in the protein1 dataframe?
require(plyr)

hits <- ddply(protein1, "id", function(x) any(y %in% x))
#subset
matches <- protein1[hits$V1, c("id","Razor...unique.peptides")]
#many false matches!


#this works! but ddply doesn't?
listofids <- as.list(protein1$id)
hits2 <- sapply(listofids, function(x) any(y %in% x))
matches <- protein1[hits2, c("id","Razor...unique.peptides")]

#it also works in vector format
vectorofids <- protein1$id
hits3 <- sapply(vectorofids, function(x) any(y %in% x))
matches <- protein1[hits3, c("id","Razor...unique.peptides")] 

#keep the rows of the protein1 data frame that match the id
keeps <- apply(

y





#find the rows of the protein1 DF that contain either of these protein group ids

y <- paste(y, collapse = "|")


matches <- protein1[grep(y, protein1$id), c("id","Razor...unique.peptides")]
                                        
                                        
counts <- matchingGroups$Razor...unique.peptides
      assignedProtein <- matchingGroups[which.max(counts),] #need to fix this to return desired numbers
      