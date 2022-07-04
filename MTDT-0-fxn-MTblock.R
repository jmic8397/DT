

#MTBlock for ClassifyR


####################################################################


#this is actually taking a model list not a single model
MTblockClassifyR <- function(data=NULL, id.retained=NULL,  
                     ssercutoff=0.2, tier = NULL, plotlabel = "",
                    runtestorruntests,classes = NULL, crossValParams = NULL, modellingParams = NULL, characteristics = NULL,
                    performanceType = "Sample Error", seed=1, verbose=F, finalTier,classIndex, z){
  
  errormessages = NULL
  
  
  #ok, we now run classifyR runtest and use the SelectResult object to append the error rates 
  
  #errormessages = capture.output(
    SSER <- sser_classifyR( MAEobject=data,
                     seed=seed, classes=classes, crossValParams = crossValParams, 
                     modellingParams =modellingParams, characteristics = characteristics, performanceType = performanceType,
                     verbose=verbose, classIndex, z)
  # )
  
  #tier specific error rate?
  TSER <- tser(SSER, sser)
    # print(TSER)
  tserplot = tser_plot(TSER, sser, ssercutoff)
  TSERcutoff <- tser_cutoff(SSER, sser, ssercutoff)
  tsercutoffplot = tsercutoff_plot(TSERcutoff$TSER_cutoff)
  stratplot = strat_plot(TSERcutoff, plotlabel)
  
  id.retained = TSERcutoff$Stratification %>% 
    dplyr::filter(strata=="Retained") %>% 
    dplyr::select(SampleID) %>% 
    as_vector()

  id.toprogress = TSERcutoff$Stratification %>% 
    dplyr::filter(strata=="To progress") %>% 
    dplyr::select(SampleID) %>% 
    as_vector()

  id.notprocessed = TSERcutoff$Stratification %>% 
    dplyr::filter(strata=="Not processed") %>% 
    as_vector()
  id.notprocessed = c(id.notprocessed)
  names(id.retained) <- NULL
  names(id.toprogress) <- NULL
  names(id.notprocessed) <- NULL

  
  if(verbose==T){
    TSERcutoff$TSER %>% print
    tibble(strata = c("Retained", "To progress", "Not processed"),
           n = c(length(id.retained), length(id.toprogress),
                 length(id.notprocessed))) %>% 
      print
    
  }

  
  return(list(SSER=SSER,
              TSER=TSER,
              TSERcutoff=TSERcutoff,
              id=list(id.retained=id.retained,
                      id.toprogress=id.toprogress,
                      id.notprocessed=id.notprocessed),
              plots=list(tserplot=tserplot,
                         tsercutoffplot = tsercutoffplot,
                         stratplot=stratplot),
              errormessages=errormessages)
  )
}






