#' @export
simulate = function(numSims = 1, aliqot = 20/66, extract = 0.6, numCells = 20, PCReff = 0.83, stutter = TRUE, alpha = c(2,200),  numLoci = 15, print = TRUE){
  x = new(PCRSim)
  x$aliqot = aliqot
  x$extract = extract
  x$numCells = numCells
  x$PCReff = PCReff
  x$stutter = stutter
  x$setAlpha(alpha)
  x$numLoci = numLoci
  x$numSims = numSims


  if(print) x$print()

  results = x$simulate()
  results$stutters[which(results$stutters < 1, arr.ind = TRUE)] = 0

  return(results)
}


