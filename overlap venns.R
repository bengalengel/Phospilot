
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
#each one uses blocked

