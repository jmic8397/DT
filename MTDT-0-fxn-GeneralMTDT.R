
GeneralMTDT <- function(MultiAssayExperiment,
                      modelList,
                      rsmpList,
                      tierList,
                      fixedTier = NULL,
                      ssercutoffList,
                      tierUnitCosts = c(100, 500, 1000)
                      ,
                      
                      costBound=0, objective=NULL,
                      resubstituteParams =ResubstituteParams(nFeatures = seq(5, 25, 5),
                                                             performanceType = "balanced error",
                                                             better = "lower"),
                      runtestorruntests = "runtest",
                      classes = NULL,
                      params = list(SelectParams(), TrainParams(), PredictParams()),
                      leave = 2,
                      percent=25,
                      minimumOverlapPercent = 80,
                      validation = c("permute", "leaveOut", "fold"),
                      parallelParams = bpparam(),
                      
                      easyDatasetID = "clinical",
                      hardDatasetID = names(dataList)[1],
                      featureSets = NULL,
                      metaFeatures = NULL,
                      datasetName = NULL,
                      classificationName = "Easy-Hard Classifier",
                      easyClassifierParams = list(minCardinality = 2, minPurity = 0.9),
                      hardClassifierParams = list(selParams, TrainParams(), PredictParams()),
                      
                      k=5, permutations=100, rsmp=100,
                      seed=1, verbose=F){



  
#priceList has price for each sample of each model
#costBoundsList contains the maximum price 
#objective may be "cost" or "accuracy" -> permutations will be selected to achieve objective
#do we want a timelist eventually?
#fixedtier is whether to fix the first tier of permutations

#determine inputs and call relevant function
if(class(MultiAssayExperiment) == "MultiAssayExperiment"){
  # print("Hi")
  
  #we pass the data into MTDT.algm which calls ClassifyR and returns the file containing the error
  
  
  
  
  MTDTresults = MTDT.algmClassifyR(MultiAssayExperiment,
                                   modelList,
                                   rsmpList,
                                   tierList,
                                   fixedTier,
                                   ssercutoffList,
                                   tierUnitCosts,
                                   resubstituteParams,
                                   runtestorruntests ,
                                   classes,
                                   params ,
                                   leave ,
                                   percent,
                                   minimumOverlapPercent,
                                   validation,
                                   parallelParams,
                                   easyDatasetID,
                                   hardDatasetID,
                                   featureSets,
                                   metaFeatures,
                                   datasetName,
                                   classificationName,
                                   easyClassifierParams,
                                   hardClassifierParams,
                                   k, permutations, rsmp,
                                   seed, verbose)

  
  
  }else{
    MTDTresults = MTDT.algm(MultiAssayExperiment,
                                     modelList,
                                     rsmpList,
                                     tierList,
                                      fixedTier = fixedTier,
                                     ssercutoffList = ssercutoffList,
                                     k=5, times=100, rsmp=100,
                                     seed=1, verbose=F)
  }
  


return(MTDTresults)

}

