
GeneralMTDT <- function(MultiAssayExperiment,
                      tierList,fixedTier = NULL,
                      ssercutoffList,tierUnitCosts = c(100, 500, 1000), 
                      performanceType = "Sample Error", runtestorruntests = "runtest",
                      classes = NULL, crossValParams = NULL, modellingParams = NULL, characteristics = NULL,
                      minTierSize = 10, seed=1, verbose=F){

  
#priceList has price for each sample of each model
#costBoundsList contains the maximum price 
#objective may be "cost" or "accuracy" -> permutations will be selected to achieve objective
#do we want a timelist eventually?
#fixedtier is whether to fix the first tier of permutations
  
MTDTresults = MTDT.algmClassifyR(MultiAssayExperiment, tierList, fixedTier = fixedTier, ssercutoffList = ssercutoffList, tierUnitCosts= tierUnitCosts,
                                 performanceType = performanceType, runtestorruntests = runtestorruntests , classes=classes , crossValParams=crossValParams, 
                                 modellingParams=modellingParams, characteristics=characteristics,minTierSize = minTierSize,  seed=seed, verbose=verbose)

return(MTDTresults)

}

