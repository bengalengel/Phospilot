#quantile vs loess comparison

DE <- multExpanded1[multExpanded1$SubtoDE == "+" | multExpanded1$SubtoDEloess == "+",]#4996

ANYDE <- multExpanded1[multExpanded1$DEcont1 == "+" | multExpanded1$DEcont1loess == "+" | 
                         multExpanded1$DEcont2 == "+" | multExpanded1$DEcont2loess == "+" | 
                         multExpanded1$DEcont3 == "+" | multExpanded1$DEcont3loess == "+",]#945



#DE contrasts either direction
DE1 <- multExpanded1[multExpanded1$DEcont1 == "+" | multExpanded1$DEcont1loess == "+",]#945

DE2 <- multExpanded1[multExpanded1$DEcont2 == "+" | multExpanded1$DEcont2loess == "+",]#945

DE3 <- multExpanded1[multExpanded1$DEcont3 == "+" | multExpanded1$DEcont3loess == "+",]#945


##now get the id_mults that are unique to each approach at the global and specific contrast level, as well as their percentages.


idsquant <- ANYDE[ANYDE$DEcont1 == "+" | ANYDE$DEcont2 == "+" | ANYDE$DEcont3 == "+", "idmult"]
idsloess <- ANYDE[ANYDE$DEcont1loess == "+" | ANYDE$DEcont2loess == "+" | ANYDE$DEcont3loess == "+", "idmult"]

idsquant <- as.character(idsquant)
idsloess <- as.character(idsloess)

count(idsquant %in% idsloess) #92% coverage of loess data
count(idsloess %in% idsquant) #94% coverage of quantile data





#Directional DE contrasts
DE1up <- multExpanded1[multExpanded1$cont1up == "+",]#945
DE1down <- multExpanded1[multExpanded1$cont1down == "+",]#945
DE2up <- multExpanded1[multExpanded1$cont2up == "+",]#945
DE2down <- multExpanded1[multExpanded1$cont2down == "+",]#945
DE3up <- multExpanded1[multExpanded1$cont3up == "+",]#945
DE3down <- multExpanded1[multExpanded1$cont3down == "+",]#945



