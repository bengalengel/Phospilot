
# double venn blocked random vs averaged ----------------------------------
#both datasets have been BE corrected using combat
#Comparing global signficant calls across all three contrasts using limma moderated F statistics

#first set size
pilotFsig <- nrow(sigFvalspilot)
#second set size
blocknestedFsig <- nrow(sigFvalsadata)
#intersection
intsect <- length(intersect(row.names(sigFvalsadata), row.names(sigFvalspilot)))

#make a double venn
plot.new()
venn.plot <- draw.pairwise.venn(
  area1 = pilotFsig,
  area2 = blocknestedFsig,
  cross.area = intsect,
  category = c("Averaged Bioreps", "BlockedREmodel"),
  fill = c("green", "blue"),
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.col = c("green", "blue"), 
  margin = .1,
  main="test"
)
grid.newpage();
grid.arrange(gTree(children=venn.plot), main="Design Matrix Impact on Limma Differential Phosphorylation")


# double venn combat vs medianquantilenormalized with batch as a covariate --------------------------

#first set size
batchcovar <- nrow(sigFvalBREBatchCov)
#second set size
combat <- nrow(sigFvalsadata)
#intersection
intsect <- length(intersect(row.names(sigFvalsadata), row.names(sigFvalBREBatchCov)))

#make a double venn
plot.new()
venn.plot <- draw.pairwise.venn(
  area1 = batchcovar,
  area2 = combat,
  cross.area = intsect,
  category = c("BatchCovariate", "Combat"),
  fill = c("green", "blue"),
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.col = c("green", "blue"), 
  margin = .1,
  main="test"
)
grid.newpage();
grid.arrange(gTree(children=venn.plot), main="Batch Correction Apporach Impact on Limma Differential Phosphorylation")


# phosprep protein as cov vs confounded -----------------------------------

#first set size
PhosProtcovar <- nrow(sigFvalsPhosPrepProt)
#second set size
Confounded <- nrow(sigFvalsCombatCommon)
#intersection
intsect <- length(intersect(row.names(sigFvalsPhosPrepProt), row.names(sigFvalsCombatCommon)))

#make a double venn
plot.new()
venn.plot <- draw.pairwise.venn(
  area1 = PhosProtcovar,
  area2 = Confounded,
  cross.area = intsect,
  category = c("PhosPrep Protein Covariate", "Confounded"),
  fill = c("green", "blue"),
  lty = "blank",
  cex = 1.5,
  cat.cex = 1.5,
  cat.col = c("green", "blue"), 
  margin = .2,
  main="test"
)
grid.newpage();
grid.arrange(gTree(children=venn.plot), main="Differential Phosphorylation Overlap \n Confounded vs PhosPrep Protein as a covariate")



# GelPrep protein as cov vs confounded -----------------------------------

#first set size
GelProtcovar <- nrow(sigFvalsGelPrepProt)
#second set size
Confounded <- nrow(sigFvalsCombatCommonGel)
#intersection
intsect <- length(intersect(row.names(sigFvalsGelPrepProt), row.names(sigFvalsCombatCommonGel)))

#make a double venn
plot.new()
venn.plot <- draw.pairwise.venn(
  area1 = GelProtcovar,
  area2 = Confounded,
  cross.area = intsect,
  category = c("GelPrep Protein Covariate", "Confounded"),
  fill = c("green", "blue"),
  lty = "blank",
  cex = 1.5,
  cat.cex = 1.5,
  cat.col = c("green", "blue"), 
  margin = .2,
  main="test"
)
grid.newpage();
grid.arrange(gTree(children=venn.plot), main="Differential Phosphorylation Overlap \n Confounded vs GelPrep Protein as a covariate")




