#Here I will extract all sequence windows and the modified residues 7 AAs on either side.Subsetted by modified AA of course

#Subjected to DE
DE <- multExpanded1[multExpanded1$SubtoDE == "+",]#4996

#DE contrasts
DE1 <- multExpanded1[multExpanded1$DEcont1 == "+",]#945
DE1S <- DE1[DE1$Amino.acid == "S",]#905
DE1T <- DE1[DE1$Amino.acid == "T",]

DE2 <- multExpanded1[multExpanded1$DEcont2 == "+",]#945
DE2S <- DE2[DE2$Amino.acid == "S",]#905
DE2T <- DE2[DE2$Amino.acid == "T",]

DE3 <- multExpanded1[multExpanded1$DEcont3 == "+",]#945
DE3S <- DE3[DE3$Amino.acid == "S",]#905
DE3T <- DE3[DE3$Amino.acid == "T",]


#Directional DE contrasts
DE1up <- multExpanded1[multExpanded1$cont1up == "+",]#945
DE1down <- multExpanded1[multExpanded1$cont1down == "+",]#945
DE2up <- multExpanded1[multExpanded1$cont2up == "+",]#945
DE2down <- multExpanded1[multExpanded1$cont2down == "+",]#945
DE3up <- multExpanded1[multExpanded1$cont3up == "+",]#945
DE3down <- multExpanded1[multExpanded1$cont3down == "+",]#945


#for networKIN

DE <- DE[,c("Protein","Position","Amino.acid")]
DE1nk <- DE1[,c("Protein","Position","Amino.acid")]
DE2nk <- DE2[,c("Protein","Position","Amino.acid")]
DE3nk <- DE3[,c("Protein","Position","Amino.acid")]
DE1nkup <- DE1up[,c("Protein","Position","Amino.acid")]
DE1nkdown <- DE1down[,c("Protein","Position","Amino.acid")]
DE2nkup <- DE2up[,c("Protein","Position","Amino.acid")]
DE2nkdown <- DE2down[,c("Protein","Position","Amino.acid")]
DE3nkup <- DE3up[,c("Protein","Position","Amino.acid")]
DE3nkdown <- DE3down[,c("Protein","Position","Amino.acid")]

write.csv(DE,"subtoDE.csv", row.names=F)
write.csv(DE1nk,"DE1nk.csv", row.names=F)
write.csv(DE2nk,"DE2nk.csv", row.names=F)
write.csv(DE3nk,"DE3nk.csv", row.names=F)

write.csv(DE1nkup,"DE1nkup.csv", row.names=F)
write.csv(DE1nkdown,"DE1nkdown.csv", row.names=F)
write.csv(DE2nkup,"DE2nkup.csv", row.names=F)
write.csv(DE2nkdown,"DE2nkdown.csv", row.names=F)
write.csv(DE3nkup,"DE3nkup.csv", row.names=F)
write.csv(DE3nkdown,"DE3nkdown.csv", row.names=F)




#pick one and then run through below. I should make this into a function at some point.
#sequences1 <- DE3S$Sequence.window
#sequences1 <- DE3T$Sequence.window

sequences1 <- DE3down$Sequence.window

##how long is the sequence window
sequences1 <- as.character(sequences1)
nchar(sequences1[1])#31 with 15 AAs on either side of the modified peptide.

##remove alternative protein sequences
sequences2 <- c()
for(i in seq(along=sequences1)){
temp <- sub(";.*","",sequences1[i])
sequences2 <- c(sequences2,temp)}

# I need to extract character 16 with 7 on either side. 
sequences3 <- c()
for(i in seq(along=sequences2)){
  temp <- substr(sequences2[i],9,23)
  sequences3 <- c(sequences3,temp)}

#remove the serine or threonine
sequences4 <- c()
for(i in seq(along=sequences3)){
  temp <- paste(substr(sequences3[i],1,7),substr(sequences3[i],9,nchar(sequences3[i])), sep='')
  sequences4 <- c(sequences4,temp)}

#write output with serines
#DE1serine <- unique(sequences3)
#write.table(DE1serine,"DE1serine.csv",sep=',',col.names=F,row.names=F)

DE1serine <- unique(sequences3)
write.table(DE1serine,"DE3downserine.csv",sep=',',col.names=F,row.names=F)



#write output without serines
DE1serine2 <- unique(sequences4)
write.table(DE1serine2,"DE1serine2.csv",sep=',',col.names=F,row.names=F)

#########OR!

#write output with threonines
DE1thr <- unique(sequences3)
write.table(DE1thr,"DE1thr.csv",sep=',',col.names=F,row.names=F)

#write output without threonines
DE1thr2 <- unique(sequences4)
write.table(DE1thr2,"DE1thr2.csv",sep=',',col.names=F,row.names=F)






#write output with serines
DE2serine <- unique(sequences3)
write.table(DE2serine,"DE2serine.csv",sep=',',col.names=F,row.names=F)

#write output without serines
DE2serine2 <- unique(sequences4)
write.table(DE2serine2,"DE2serine2.csv",sep=',',col.names=F,row.names=F)

#########OR!

#write output with threonines
DE2thr <- unique(sequences3)
write.table(DE2thr,"DE2thr.csv",sep=',',col.names=F,row.names=F)

#write output without threonines
DE2thr2 <- unique(sequences4)
write.table(DE2thr2,"DE2thr2.csv",sep=',',col.names=F,row.names=F)







#write output with serines
DE3serine <- unique(sequences3)
write.table(DE3serine,"DE3serine.csv",sep=',',col.names=F,row.names=F)

#write output without serines
DE3serine2 <- unique(sequences4)
write.table(DE3serine2,"DE3serine2.csv",sep=',',col.names=F,row.names=F)

#########OR!

#write output with threonines
DE3thr <- unique(sequences3)
write.table(DE3thr,"DE3thr.csv",sep=',',col.names=F,row.names=F)

#write output without threonines
DE3thr2 <- unique(sequences4)
write.table(DE3thr2,"DE3thr2.csv",sep=',',col.names=F,row.names=F)


