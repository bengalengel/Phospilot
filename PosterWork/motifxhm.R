##melt expand then produce a heatmap from the motif-x output

require(plyr)
require(reshape2)
require(gplots)
require(RColorBrewer)


motifx <- read.csv("motifxmelted.csv")


casted <- dcast(melted, ... ~ sample + bio + tech, value.var="value") ##close but creates extra rows


motifxcast <- dcast(motifx, ... ~ Contrast)
names <- as.character(motifxcast$Motif)
names
motifxcast <- motifxcast[,2:ncol(motifxcast)]
motifxcast <- as.matrix(motifxcast)
row.names(motifxcast) <- names
motifxcast

heatmap.2(motifxcast,
          Colv="NA",
          Rowv="NA",
          col=bluered(25),
          scale="none",
          trace="none",
          density.info="none",
          na.color="darkgrey",
          key.xlab = "-log(P)", key.title = "Motif Enrichment",
          margins = c(7,8),
          cexCol=1,
          cexRow=1
          )



# 
# 
# heatmap.2(
#   r,#row Z scores
#   Colv=sample.dend,
#   Rowv=feature.dend,
#   col=bluered(25),
#   scale="none",
#   trace="none",
#   key.xlab = "Row Z scores", key.ylab=NULL, key.title = "",
#   srtCol=45,  ,adjCol = c(1,1),
#   margins = c(6,5),
#   cexCol=1,
#   labRow = NA#remove row labels
# )
