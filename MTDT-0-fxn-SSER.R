

##########################################################################

#Sample specific error calculated using classifyR's predictions

sser_classifyR <- function(model=NULL, MAEobject=NULL, method=rsmpmethod, k=k, times=times, rsmp=rsmp,
                           seed=seed, classes, params, leave, percent, resubstituteParams, minimumOverlapPercent, validation, 
                           easyDatasetID, hardDatasetID,featureSets,metaFeatures, datasetName, classificationName, training, 
                           testing, verbose=verbose,classIndex, z){
  
  # if (is.null(model)) return(print("missing model)"))
  # if (is.null(data)) return(print("missing data"))
  # if (!("lda_diag") %in% class(model)) return(print("not a lda_diag object"))
  # 
  # 
  
  rsmp_method = data.frame(rsmpl_method="boot", rsmp=rsmp)
  
  library(sparsediscrim)
  library(tidyverse)
  library(dplyr)
  
  

  #removing NA's
  # dataAsFrame = as.data.frame(dataset)
  # data.ori = dataAsFrame
  # ID.ori = data.ori[, 1] %>% as_vector
  # dataAsFrame = dataAsFrame %>% dplyr::filter(complete.cases(.))
  # ID = dataAsFrame[, 1] %>% as_vector

  # ID.removed = setdiff(ID.ori, ID) %>% sort()
  # print(ID.removed)
  # print(MAEobject)
  # 
  # ID_ColName = colnames(colData(MAEobject)[1])
  # 
  # print(colData(MAEobject)[1])
  # IDnm = as.vector(colData(MAEobject)[['ID_ColName']])
  # ynm = as.vector(colData(MAEobject)[[classIndex]])
  # ID_nm = sym(IDnm)
  # y_nm = sym(ynm)
  # 
  # # print("ID")
  # # print(ID_nm)
  # 
  # 
  # 
  # removed = data.ori %>% 
  #   dplyr::filter(!!ID_nm %in% ID.removed) %>% 
  #   dplyr::select(!!ID_nm, !!y_nm)
  # 
  # y = dataAsFrame[, classIndex]
  # x = dataAsFrame[, -c(1, classIndex)]
  #L = predict(model, newdata=x)
  
  # phat = tibble(
  #   ID_nm := IDnm,
  #   y_nm :=  ynm,
  #   yhat := ynm , #This should be fixed? ****************************************
  #   # yhat = L$class %>% as.character() %>% as.numeric(),
  #   # phat = L$posterior[2],
  #   # yhat = L   %>% as.factor() %>% as.numeric() -1
  # )
  
  phat_rsmp = NULL
  
  set.seed(seed)
  

      #Insert cv params from user here
  
      # cv = CrossValParams()
      cv = CrossValParams(
        samplesSplits = "Permute k-Fold",
        permutations = 5,
        percentTest = 25,
        folds = 5,
        leave = 2,
        tuneMode = "none",
        parallelParams = bpparam()
      )
      
      #Select params
      # sp = SelectParams(featureSelection, characteristics = DataFrame(), minPresence = 1, intermediate = character(0),
      #                   subsetToSelections = TRUE, tuneParams = list(nFeatures = seq(10, 100, 10), performanceType = "Balanced Error"), ...)
      # 
      
      
      #insert model params from user here
      mp = ModellingParams(
        balancing = "none",
        transformParams = NULL,
        selectParams = SelectParams(),
        trainParams = TrainParams(),
        predictParams = PredictParams(),
        doImportance = FALSE
      )

      DMresults <- runTests(MAEobject, target=names(MAEobject)[z], outcomesColumns=classes, crossValParams = cv, modellingParams = mp)
      
      
      # show(DMresults)
      # print(predictions(DMresults))
      
      # This gives estimate error rate
      performance = calcCVperformance(DMresults, "Sample Error")
      error = performance(performance)
      # print(error)
      
      #Run prediction with classifyR
      # 
      # selParams <- SelectParams(featureSelection = differentMeansSelection, selectionName = "Difference in Means",
      #                           resubstituteParams = ResubstituteParams(1:10, "balanced error", "lower"))
      # 
      # easyHardCV <- runTestEasyHard(dataset, datasetName = "Test Data", classificationName = model, training, testing, easyClassifierParams = list(minCardinality = 2, minPurity = 0.9),
      #                               hardClassifierParams = list(selParams, TrainParams(), PredictParams()), easyDatasetID = "clinical", hardDatasetID = names(dataset)[1],
      #                               featureSets = NULL, metaFeatures = NULL, minimumOverlapPercent = 80,
      # 
      #                               verbose = 3)



      # 
      # selParams <- SelectParams(featureSelection = differentMeansSelection, selectionName = "Difference in Means",
      #                           resubstituteParams = ResubstituteParams(1:10, "balanced error", "lower"))
      # easyHardCV <- runTestEasyHard(dataset, datasetName = "Test Data", classificationName = "Easy-Hard", training = 1:10, testing = 1:10,
      #                               easyClassifierParams = list(minCardinality = 2, minPurity = 0.9),
      #                               hardClassifierParams = list(selParams, TrainParams(), PredictParams())
      # )
      # 
      # 
      # 
      
      
      
      # 
      # selParams <- SelectParams(featureSelection = differentMeansSelection, selectionName = "Difference in Means",
      #                           resubstituteParams = ResubstituteParams(1:10, "balanced error", "lower"))
      # 
      # easyHardCV <- runTestEasyHard(dataset, datasetName = "Test Data", classificationName = model, training = 1:10, testing = 1:10,
      #                               easyClassifierParams = list(minCardinality = 2, minPurity = 0.9),
      #                               hardClassifierParams = list(selParams, TrainParams(), PredictParams())
      # )      # 

      
      
      


    
    #
    # ClassifyRPredictions = df[[1]] %>% select(class) %>% mutate(class = ifelse(class=="Good", "1", "0"))
    # print(ClassifyRPredictions)
    # ClassifyRPredictions$class

# 
#       #append the repetition number and prediction to phat_rsmp
#       pred = dplyr::tibble(
#         !!IDnm := ID,
#         rpt = j,
#         !!ynm := y,
#         # yhatrsmp = predict %>% as.factor() %>%as.numeric()-1
#         yhatrsmp =ClassifyRPredictions$class
#         # yhatrsmp = sample(c(rep(1,100),rep(0,100)),length(y),replace = TRUE)
#         
#       )
#       phat_rsmp = dplyr::bind_rows(phat_rsmp, pred)
#     
  # }
  # print(phat_rsmp)

  # 
  # SSER = phat_rsmp %>% 
  #   dplyr::left_join(phat %>% dplyr::select(-class), by="Person_ID") %>% 
  #   dplyr::group_by(!!ID_nm) %>% 
  #   dplyr::summarise(
  #     sser=(1-mean(yhatrsmp==!!y_nm, na.rm=T)),
  #     yhatrsmperror=(1-mean(yhatrsmp==yhat, na.rm=T))) %>% 
  #   dplyr::left_join(phat, by=IDnm) %>% 
  #   dplyr::mutate(phat=NA, phat.SE=NA, phat.lower=NA, phat.upper=NA, 
  #                 phatrsmp_mean=NA, phatrsmp_median=NA, phatrsmp_GiniMd=NA, 
  #                 phatrsmp_SE=NA, phatrsmp_IQR=NA, absdiffphatphatrsmp_mean=NA) %>% 
  #   dplyr::select(!!ID_nm, !!y_nm, yhat, phat.SE, phat.lower, phat.upper, 
  #                 sser, yhatrsmperror,
  #                 starts_with("phatrsmp"), everything()) %>% 
  #   dplyr::arrange(!!ID_nm)
  # 
  
  # print(SSER)
  
  # 
  # SSER = phat_rsmp %>% 
  #   dplyr::left_join(phat %>% dplyr::select(-y.good), by="Person_ID") %>% 
  #   dplyr::group_by(!!ID_nm) %>% 
  #   dplyr::summarise(
  #     sser=(1-mean(yhatrsmp==!!y_nm, na.rm=T)),
  #     yhatrsmperror=(1-mean(yhatrsmp==yhat, na.rm=T))) %>% 
  #   dplyr::left_join(phat, by=IDnm) %>% 
  #   dplyr::mutate(phat=NA, phat.SE=NA, phat.lower=NA, phat.upper=NA, 
  #                 phatrsmp_mean=NA, phatrsmp_median=NA, phatrsmp_GiniMd=NA, 
  #                 phatrsmp_SE=NA, phatrsmp_IQR=NA, absdiffphatphatrsmp_mean=NA) %>% 
  #   dplyr::select(!!ID_nm, !!y_nm, yhat, phat.SE, phat.lower, phat.upper, 
  #                 sser, yhatrsmperror,
  #                 starts_with("phatrsmp"), everything()) %>% 
  #   dplyr::arrange(!!ID_nm)
  # 
  # 
  
  #Add this information eventually ****************************************************************
  
  # SSER = dplyr::left_join(
  #   phat,
  #   phat_rsmp %>% 
  #     dplyr::group_by(!!ID_nm) %>% 
  #     dplyr::summarise(
  #       sser=(1-mean(yhatrsmp==!!y_nm, na.rm=T)),
  #       yhatrsmperror=(1-mean(yhatrsmp==yhat, na.rm=T)),
  #       phatrsmp_mean=mean(phatrsmp, na.rm=T),
  #       phatrsmp_median=median(phatrsmp, na.rm=T),
  #       phatrsmp_GiniMd=Hmisc::GiniMd(phatrsmp, na.rm=T),
  #       phatrsmp_SE=sd(na.omit(phatrsmp))/sqrt(sum(!is.na(phatrsmp))),
  #       phatrsmp_IQR=IQR(phatrsmp, na.rm=T),
  #       absdiffphatphatrsmp_mean=mean(absdiffphatphatrsmp, na.rm=T)),
  #   by = IDnm) %>% 
  #   dplyr::arrange(!!ID_nm)
  
  # phat_rsmp = phat_rsmp %>% 
  #   dplyr::select(
  #     !!ID_nm, !!y_nm, rpt, 
  #     yhatrsmp, starts_with("phatrsmp"), 
  #     absdiffphatphatrsmp)  %>%
  #   dplyr::arrange(!!ID_nm, rpt)
  
  # phat_rsmp = phat_rsmp %>% 
  #   dplyr::select(
  #     !!ID_nm, !!y_nm, rpt, yhatrsmp) %>% 
  #   dplyr::mutate(phatrsmp=NA, phatrsmp.SE=NA, 
  #                 phatrsmp.lower=NA, phatrsmp.upper=NA, absdiffphatphatrsmp=NA) %>% 
  #   dplyr::arrange(!!ID_nm, rpt)
      
      SampleID = names(error$`Sample Error`)
      sser = unname(error$`Sample Error`)
      errorTable = tibble(SampleID,sser)
      view(errorTable)
  # 
  
  if (verbose==T){
    if (length(ID.removed)>0) {
      print(paste("Incomplete data.", length(ID.removed), "IDs removed:"))
      print(ID.removed)
    }
  }
  
  return(list(rsmp_method=rsmp_method,
              
              SSER=errorTable, 
              phat_rsmp=phat_rsmp))
}



#tier specific error
# TSE with different stability cutoffs
## col_name is the column name for the stability measure in the dataset in cvPSE
## col_name = sser = 1 - avg(yhatrsmp == y_nm) = error rate!
tser <- function(SSER=NULL, col_name=NULL){
  library(tidyverse)
  library(rms)
  
  # phat_rsmp = SSER$phat_rsmp
  SSER=SSER$SSER
  
  #print(SSER)
  
  # IDnm=colnames(SSER)[1]
  # ynm=colnames(phat_rsmp)[2]
  # ID_nm=sym(IDnm)
  # y_nm=sym(ynm)
  # 
  # y = SSER %>% dplyr::select(!!ID_nm, !!y_nm)
  
  #cutoffs are the average error rates for samples
  cutoffs = SSER %>% 
    dplyr::select({{col_name}}) %>% 
    distinct() %>% 
    as_vector() %>% 
    sort()

  names(cutoffs) <- NULL
  #print(cutoffs)
  cutoffs
  # phat_rsmp = dplyr::left_join(
  #   SSER, 
  #   phat_rsmp %>% dplyr::select(-!!y_nm),
  #   by=IDnm)
  
  err = dplyr::tibble(tser=NULL, cutoff=NULL)
  n = vector()
  cindices = vector()
  briers = vector()
  
  #for each cutoff 
  for (i in 1:length(cutoffs)){
    
    #threshold = sser (1-avg(yhatrsmp == y_nm))
    threshold = cutoffs[i]
    
    #filter for less than threshold?
    # df.phat_rsmp = phat_rsmp %>% dplyr::filter({{col_name}}<=threshold)
    df.sser = SSER %>% dplyr::filter({{col_name}}<=threshold)

    
    # #all below threshold
    # tser = df.phat_rsmp %>%
    #   #group by repetition number (several folds per repetition)
    #   dplyr::group_by(rpt) %>%
    #   #create new df for each repetition number containing mean accuracy
    #   dplyr::summarise(tser=(1-mean(!!y_nm==yhatrsmp, na.rm=T))) %>%
    #   #cutoff = threshold?
    #   dplyr::mutate(cutoff=threshold) %>%
    #   #select everything but repetitions, and tser
    #   dplyr::select(-rpt, tser)
    # 
    
    
    #all below threshold
    tser = df.sser %>%
      #group by repetition number (several folds per repetition)
      # dplyr::group_by(rpt) %>%
      #create new df for each repetition number containing mean accuracy
      dplyr::summarise(tser=(1-mean(sser, na.rm=T))) %>%
      #cutoff = threshold?
      dplyr::mutate(cutoff=threshold)
      #select everything but repetitions, and tser
      # dplyr::select(tser)
    
    
    print(tser)
    # if (is.error(with(df.sser, rms::val.prob(phat, !!y_nm, pl=F)))) {
    #   val = rep(NA, 18)
    # } else {
    #   #val = with sser, 
    #   #rms::val = validating predicted probabilities against binary events
    #   #returns a vector with the following named elements: Dxy, R2, D, D:Chi-sq,D:p, U, U:Chi-sq, U:p, Q, Brier, Intercept, Slope, S:z, S:p, Emax.
    #   val = with(df.sser, rms::val.prob(phat, !!y_nm, pl=F))
    #   names(val) = NULL
    # }
    # #validation metrics
    # cindices = c(cindices, val[2])
    # #validation metrics
    # briers = c(briers, val[11])
    #n retained = rows left in sser <= threshold error
    n = c(n, dim(df.sser)[1])
    err = dplyr::bind_rows(err, tser)
  }
  
  
  TSER = tibble(
    Index = rep(as_label(enquo(col_name)), length(cutoffs)),
    cutoff = cutoffs,
    tse = err %>%
      dplyr::group_by(cutoff) %>% 
      #tser = average tser
      dplyr::summarise(tser=mean(tser, na.rm=T)) %$% 
      tser,
    # cindices = cindices,
    # briers = briers,
    n_total = dim(SSER)[1],
    n_retained = n,
    prop_retained = n_retained/n_total,
    n_progress = n_total-n_retained,
    prop_progress = n_progress/n_total,
  )
  
  return(TSER)
}




# Plotting TSE - overall, included and excluded from the cutoff
tser_cutoff <- function(SSER=NULL, col_name=sser, mycutoff=0.5, finalTier){
  
  library(tidyverse)
  # Removed = SSER$removed
  # Removed = none()
  # phat_rsmp = SSER$phat_rsmp
  SSER=SSER$SSER
  # IDnm=colnames(phat_rsmp)[1]
  # ynm=colnames(phat_rsmp)[2]
  # ID_nm=sym(IDnm)
  # y_nm=sym(ynm)
  
  # phat_rsmp = dplyr::left_join(
  #   SSER,
  #   phat_rsmp %>% dplyr::select(-!!y_nm), 
  #   by=IDnm
  #   )
  
  print(SSER)
  
  id.retained = SSER %>% 
    dplyr::filter({{col_name}}<=mycutoff) %>% 
    dplyr::select(SampleID) %>%
    as_vector() %>% unique()

  
  id.toprogress = setdiff(
    SSER %>% dplyr::select(SampleID) %>% as_vector(), 
    id.retained
    )
  n=dim(SSER)[1]
  n_retained=length(id.retained)
  n_progress=n-n_retained
  
  
  
  TSER_overall = SSER %>% 
    
    dplyr::summarise(tser=(1-mean(sser, na.rm=T))) %>% 
    dplyr::mutate(strata="Overall",
                  n=n)
  
  TSER_retained = SSER %>% 
    dplyr::filter({SampleID} %in% id.retained) %>% 
    # dplyr::group_by(rpt) %>% 
    dplyr::summarise(tser=(1-mean(sser, na.rm=T))) %>% 
    dplyr::mutate(strata="Retained",
                  n=n_retained)

  
  TSER_progress = SSER %>% 
    dplyr::filter(!{SampleID} %in% id.retained) %>% 
    # dplyr::group_by(rpt) %>% 
    dplyr::summarise(tser=(1-mean(sser, na.rm=T))) %>% 
    dplyr::mutate(strata="To progress",
                  n=n_progress)
  
  
  TSER_cutoff = bind_rows(TSER_overall, TSER_retained, TSER_progress) %>% 
    dplyr::mutate(
      strata=factor(strata, levels=c("Overall", "Retained", "To progress"))
    )
  
  
  strata = dplyr::bind_rows(
    SSER %>% 
      dplyr::filter(SampleID %in% id.retained) %>% 
      dplyr::select(SampleID, sser) %>% 
      dplyr::mutate(strata = "Retained"),
    SSER %>% 
      dplyr::filter(SampleID %in% id.toprogress) %>% 
      dplyr::select(SampleID, sser) %>% 
      dplyr::mutate(strata = "To progress"),
    # Removed %>%
    #   dplyr::mutate(strata = "Not processed")
    ) %>%
    dplyr::mutate(strata=factor(strata, levels=c("Retained", "To progress", "Not processed")))
  
  #prev
  # strata = dplyr::bind_rows(
  #   SSER %>% 
  #     dplyr::filter(!!ID_nm %in% id.retained) %>% 
  #     dplyr::select(!!ID_nm, !!y_nm, yhat, sser) %>% 
  #     dplyr::mutate(strata = "Retained"),
  #   SSER %>% 
  #     dplyr::filter(!!ID_nm %in% id.toprogress) %>% 
  #     dplyr::select(!!ID_nm, !!y_nm, yhat, sser) %>% 
  #     dplyr::mutate(strata = "To progress"),
  #   Removed %>% 
  #     dplyr::mutate(strata = "Not processed")
  # ) %>% 
  #   dplyr::mutate(strata=factor(strata, levels=c("Retained", "To progress", "Not processed")))
  # 
  # 
  
  # tse_boxplot = ggplot(data=TSER) +
  #   geom_boxplot(aes(y=tser, x=Strata, fill=Strata)) +
  #   geom_line(aes(y=n, x=Strata, colour="red")) +
  #   ylim(0, 1) +
  #   ylab("TSE") + xlab("") +
  #   theme_bw() +
  #   theme(legend.position = "none")
  
  TSER_table = TSER_cutoff %>% 
    dplyr::group_by(strata) %>% 
    dplyr::summarise(mean(tser, na.rm=T), n=mean(n))
  
  return(list(TSER_cutoff=TSER_cutoff,
              Stratification=strata,
              TSER = TSER_table))
}



# TSER vs increasing cutoffs
tser_plot <- function(TSER, col_name, mycutoff) {
  library(tidyverse)
  library(grid)
  library(gridExtra)
  
  pTSER = TSER %>% dplyr::select(cutoff, tse) %>% 
    ggplot(aes(y=tse, x=cutoff)) +
    geom_point(colour="grey", size=2) +
    geom_line(colour="black", alpha=0.8) +
    geom_line(aes(y=..y..), stat='ecdf', colour="darkorange", alpha=0.5) +
    geom_vline(xintercept=mycutoff, colour="red", linetype="longdash") +
    geom_text(data=data.frame(x=mycutoff, y=-0.025), 
              aes(x=x, y=y, label=as.character(x)), colour="red") +
    scale_y_continuous(name="Tier-Specific Accuracy Rate", 
                       sec.axis=sec_axis(~., name="Proportion retained")) +
    xlab(paste(toupper(as_label(enquo(col_name))), "Cut-offs")) +
    ggtitle("Mis-classification Error Rate") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          # axis.title.y.right = element_text(colour="darkorange"),
          axis.text.y.right = element_text(colour="darkorange"),
          axis.ticks.y.right = element_line(colour="darkorange"))
  
  return(pTSER)
}


tsercutoff_plot <- function(TSERcutoff){
  
  # n_total=TSERcutoff[TSERcutoff$Strata=="Overall",][["n"]] %>% unique()
  # N = TSERcutoff %>% 
  #   dplyr::select(-rpt, -tser) %>% 
  #   dplyr::distinct() %>% 
  #   dplyr::filter(Strata!="Overall") %>% 
  #   dplyr::mutate(total=n_total,
  #                 prop=n/total,
  #                 cumprop=cumsum(n)/n_total)
  
  # print(TSERcutoff)

  TSERcutoff = TSERcutoff %>% filter(TSERcutoff$strata == "To progress" | TSERcutoff$strata == "Retained")
  p = TSERcutoff %>%
    ggplot() +
    geom_boxplot(aes(y=1-tser, x=strata, fill=strata)) +
    # geom_point(data=N, aes(y=cumprop, x=Strata)) +
    # geom_line(data=N, aes(y=..y.., x=Strata), stat='ecdf', colour="darkorange", alpha=0.5)
    # geom_line(data=N, aes(y=cumprop, x=Strata, group=1), colour="darkorange")
    ylim(0, 1) +
    scale_x_discrete(guide = guide_axis(n.dodge=3))+
    ylab("Tier-Specific Accuracy Rate") + xlab("") +
    theme_bw() +
    theme(legend.position = "none")
  
  return(p)
}




strat_plot <- function(TSER_cutoff=NULL, tier=""){
  
  stra = TSER_cutoff$Stratification
  IDnm=colnames(stra)[1]
  ID_nm=IDnm
  ID=stra[, 1]
  StraMeasure=colnames(stra)[4]
  StraMeasure_nm=StraMeasure
  
  stra = stra %>% 
    dplyr::mutate(tier=tier) %>%  
    dplyr::arrange(tier, strata, StraMeasure_nm) %>% 
    dplyr::mutate(ID_nm:=factor(ID_nm, levels=ID_nm))
  
  p = stra %>% 
    dplyr::arrange(ID_nm) %>% 
    ggplot() +
    geom_tile(aes(y=strata, x=ID_nm, fill=StraMeasure_nm), colour="black") +
    ylab("") + xlab("") + ggtitle(tier) +
    scale_fill_gradient(low="black", high="white") +
    theme(panel.background = element_blank(),
          axis.text.x = element_blank(),
          aspect.ratio=1/10)
  return(p)
}

  


  


