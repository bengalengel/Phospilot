##Loop to assign protein ids, H/L values and ibaq values to phosphosites for normalization using phospho workup

#decision is when there are multiple protein groups that could match the phosphosite. Here I must choose the one with the most unique/razor ids

# also want to remove names from the majority protein ids that couldn't possible contain the phosphopeptide if I am to use these names for enrichment analysis.

#use sapply when possible. I will work from me1 dataframe

#collect ids from multexpanded
ids <- multExpanded1$Protein.group.IDs
ids <- as.character(ids)

#for those ids that contain multiple groups. pick the one with the most unique plus razor 
index <- grep(";", ids)
MultiMatch <- ids[index]

#for each of the suspect ids, match the proteingroup id with the most razor + unique peptides
proteinids <- protein1$id#for the loops
FinalIDMultMatch <- sapply(MultiMatch,function(x){
  y <- as.numeric(unlist(strsplit(x,split = ";")))
  #this works! but ddply doesn't?
  hits <- sapply(proteinids, function(z) any(y %in% z))
  matches <- protein1[hits, c("id","Razor...unique.peptides")] 
  #assign id with greatest amount of razor + unique peptides
  finalid <- matches$id[which.max(matches$Razor...unique.peptides)]
})

#assign these final ids to the main id vector
ids[index] <- FinalIDMultMatch
ids <- as.numeric(as.character(ids))#but why?...length 1 loss with unlist
#one na introduced by coercion

#add ids to multexpanded df
multExpanded1$matchids <- ids

#retrieve the H/L, ibaq and majority protein id information from protein1
ibaqNames <- grep("iBAQ.H.*",names(protein1), value = T)
ratios <- grep("Ratio*",names(protein1), value = T)
info <- protein1[,c("id","Majority.protein.IDs",ratios,ibaqNames)]
names(info)[1] <- "phosprepProteinID"
#merge by 'matchids' and 'id'. I don't want new columns for non-matches
test <- merge(multExpanded1,info, by.x = "matchids", by.y = "phosprepProteinID")#16986 from 17774

test <- merge(multExpanded1,info, by.x = "matchids", by.y = "phosprepProteinID", all.x = T, all.y = F)#16986 from 17774



159 1103



16383 is.na in original ids dataframe. perhaps a tie!?



#this works! but ddply doesn't?
hits <- sapply(proteinids, function(x) any(y %in% x))
matches <- protein[hits, c("id","Razor...unique.peptides")] 

#assign id with greatest amount of razor + unique peptides
finalid <- matches$id[which.max(matches$Razor...unique.peptides)]








