#This script calculates the longest contiguous stretch of disordered residues for each protein within the Iupred list

#function for one protein

test.protein <- Iupred[[5555]]
prediction <- test.protein$DisProb >= .5
run.lengths <-  rle(prediction)

with(run.lengths, {
  max(lengths[values == T])
})


Iupred.disorder.contiguous.max <- sapply(Iupred, function(x) {
  prediction <- x$DisProb >= .5
  run.lengths <-  rle(prediction)
  with(run.lengths, {
    max(lengths[values == T])
  })
})
