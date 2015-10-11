
#Table 1 info

# the number of unique proteins is currently being reported as 'proteins' column but I may want to change this to protein ids. 

expcol <- grep("HL.*_", names(multExpanded1))
me1 <- multExpanded1[, c(expcol, 18)] 

#C1 phosphopeptides in all lcls
me.all <- me1[rowSums(is.na(me1[ , 1:4])) < 4 & rowSums(is.na(me1[ , 5:8])) < 4 & rowSums(is.na(me1[ , 9:12])) < 4,]    
#C1 proteins
length(unique(me.all$Proteins))
# 3514

me.all.bio <- me1[rowSums(is.na(me1[ , 1:2])) < 2 & rowSums(is.na(me1[ , 3:4])) < 2 & rowSums(is.na(me1[ , 5:6])) < 2 &
                  rowSums(is.na(me1[ , 7:8])) < 2 & rowSums(is.na(me1[ , 9:10])) < 2 & rowSums(is.na(me1[ , 11:12])) < 2,]

length(unique(me.all.bio$Proteins))
#2073

