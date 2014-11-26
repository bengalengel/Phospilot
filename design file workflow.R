#read design table and fill in rows based on annotation
design <- read.table("./experimentalDesignTemplate.txt", sep = "\t", header=T, fill = T, quote = "")

#identify experiment name
Names <- as.character(design$Name)
#zero length negative look behind (?<!MSB)
matches <- gregexpr("(?<!MSB)[0-9]{5}", Names, perl=T)
experiment <- regmatches(Names,matches)
experiment <- as.character(experiment)
design$Experiment <- experiment


#identify fraction number
#remove last character from names if it is a letter
#the last two numbers of each string is the fraction
matches <- gregexpr("[0-9]{2}$|[0-9]{2}.$", Names, perl=T) #The OR statement is included because of the final letter
fraction <- regmatches(Names,matches)
fraction <- as.character(fraction)
design$Fraction <- fraction
design$Fraction <- gsub(design$Fraction, pattern = "[ABCD]",replacement = "")
                
design$Fraction <- as.integer(design$Fraction)
design$Experiment <- as.integer(design$Experiment)

write.table(design, "design2.txt", sep = "\t", row.names=F)




s <- "nsfghstighsl44.11.36.00-1vsdfgh"
m <- gregexpr("(?<![0-9]{1})[0-9]{2}\\.[0-9]{2}\\.[0-9]{2}\\.[0-9]{2}-[0-9]{1}(?![0-9]{1})", s, perl=TRUE)
regmatches(s, m)
[1] "44.11.36.00-1"



data <- grep("[0-9]{5}", Names, value = T)


##gives index of experiment and replicate
data <- grep("^[0-9]+_[12]_[12]", colnames(casted))

##gives string of experiment and replicate
data2 <- colnames(casted)[grep("^[0-9]+_[12]_[12]", colnames(casted))]


{min,max}

design$Experiment <- 
