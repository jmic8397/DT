MTblock <- function(data=NULL, formula=NULL, model=NULL,
                    method="lrm", rsmpmethod="rcv",
                    k=5, times=100, rsmp=100,
                    id.include=NULL, ssercutoff=0.2,
                    tier = NULL, plotlabel = "",
                    seed=1, verbose=T){
  
  errormessages = NULL
  
  #convert data to dataframe
  d = as.data.frame(data)
  #convert ID and y to vector
  ID.ori = d[, 1] %>% as_vector
  y.ori = d[, 2] %>% as_vector
  
  IDnm = colnames(d)[1]
  ynm = colnames(d)[2]
  ID_nm = sym(IDnm)
  y_nm = sym(ynm)
  
  #filter data for those to be included from previous layer
  d = d %>% dplyr::filter(ID_nm %in% id.include)
  id.notincluded = setdiff(id.include, d %>% dplyr::select(ID_nm) %>% as_vector)
  names(id.notincluded) <- NULL
  
  
  #model logic
  
  #LRM
  if (method=="lrm"){
    m = model
    if (is.null(model)){
      m = rms::lrm(formula = formula, data = d, x=T, y=T, se.fit=T)
    }
    # print(m, coefs=F)
    errormessages = capture.output(
      SSER <- sser_lrm(model=m, data=d,
                       method=rsmpmethod, k=k, times=times, rsmp=rsmp, 
                       seed=seed, verbose=verbose)
    )
  }
  
  #GLM
  if (method=="glm"){
    m = model
    if (is.null(model)){
      m = glm(formula = formula, data = d, family=binomial(link="logit"))
    }
    # summary(m, coefs=F)
    errormessages = capture.output(
      SSER <- sser_lrm(model=m, data=d,
                       method=rsmpmethod, k=k, times=times, rsmp=rsmp, 
                       seed=seed, verbose=verbose)
    )
  }
  
  #DLDA
  if (method=="dlda"){
    m = model
    if (is.null(model)){
      y=d %>% dplyr::select(!!y_nm) %>% as_vector() %>% as.factor()
      names(y) = NULL
      m = sparsediscrim::lda_diag(
        y=y,
        x=d %>% dplyr::select(-ID_nm, -!!y_nm)
      )
      rm(y)
    }
    SSER <- sser_dlda(m, d, "boot", rsmp=rsmp, seed=seed)
  }
  
  
  #########################
  #KNN
  if (method=="knn"){
    #non parametric?
    m = model
    #if (is.null(model)){
     # m = rms::lrm(formula = formula, data = d, x=T, y=T, se.fit=T)
    #}
    # print(m, coefs=F)
    
    errormessages = capture.output(
      SSER <- sser_knn(model=m, data=d,
                       method=method, k, times, rsmp=rsmp, 
                       seed=seed, verbose=verbose)
    )
    
  }
  ##################################
  
  #########################
  #  #########################
  if (method=="svm"){

    m = model
    #if (is.null(model)){
    # m = rms::lrm(formula = formula, data = d, x=T, y=T, se.fit=T)
    #}
    # print(m, coefs=F)
    
    errormessages = capture.output(
      SSER <- sser_svm(model=m, data=d,
                       method=method, k, times, rsmp=rsmp, 
                       seed=seed, verbose=verbose)
    )
    
  }
  #########################
  #print("Test2")
  #print(SSER)
  #flush.console()

  #tier specific error rate?
  TSER <- tser(SSER, sser)
  tserplot = tser_plot(TSER, sser, ssercutoff)
  TSERcutoff <- tser_cutoff(SSER, sser, ssercutoff)
  tsercutoffplot = tsercutoff_plot(TSERcutoff$TSER_cutoff)
  stratplot = strat_plot(TSERcutoff, plotlabel)
  
  
  id.retained = TSERcutoff$Stratification %>% 
    dplyr::filter(strata=="Retained") %>% 
    dplyr::select(ID_nm) %>% 
    as_vector()
  
  # print("TSERRET")
  # print(id.retained)
  # 
  
  id.toprogress = TSERcutoff$Stratification %>% 
    dplyr::filter(strata=="To progress") %>% 
    dplyr::select(ID_nm) %>% 
    as_vector()
  id.notprocessed = TSERcutoff$Stratification %>% 
    dplyr::filter(strata=="Not processed") %>% 
    dplyr::select(ID_nm) %>% 
    as_vector()
  id.notprocessed = c(id.notprocessed, id.notincluded)
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








#MTBlock for ClassifyR


####################################################################


#this is actually taking a model list not a single model
MTblockClassifyR <- function(data=NULL, id.retained=NULL, formula=NULL, model=NULL, method="lrm", rsmpmethod="rcv",
                    k=5, permutations=100, rsmp=1, ssercutoff=0.2, tier = NULL, plotlabel = "",
                    seed=1, verbose=T,runtestorruntests,classes,  params,leave,percent,  resubstituteParams,     minimumOverlapPercent,       validation,
                    parallelParams,easyDatasetID ,hardDatasetID,featureSets , metaFeatures,
                    datasetName ,classificationName ,   training, testing, finalTier,classIndex, z){
  
  errormessages = NULL
  
  
  #ok, we now run classifyR runtest and use the SelectResult object to append the error rates 
  
  # #convert data to dataframe
  # MAEdataAsFrame = as.data.frame(colData(MAEobject))
  # 
  # print("Data")
  # print(MAEdataAsFrame)
  # 
  # 
  # #convert ID and y to vector
  # ID.ori = MAEdataAsFrame[1] %>% as_vector
  # y.ori = MAEdataAsFrame[classIndex] %>% as_vector
  # 
  # IDnm = colData(MAEdataAsFrame)[1]
  # ynm = colData(MAEdataAsFrame)[classIndex]
  # 
  
  
  # ID_nm = sym(IDnm)
  # y_nm = sym(ynm)
  # 
  
  # print("NAME INC ?")
  #print(ID_nm)
  #print(id.include)
  
  
  # 
  # MAEwithoutRetained = MAEobject
  # 
  # if(length(id.retained)>0){
  #   for(person in 1:length(id.retained)){
  #     retainList = MAEwithoutRetained$Person_ID != id.retained[person]
  #     # print(retainList)
  #     MAEwithoutRetained = MAEwithoutRetained[,retainList & !is.na(retainList),]
  #   }
  # }
  
  
  
  
  # print(d)
  
  #ok, we want to get our SSER object with classifyR, maybe we don;t need a seperate function here? 
  #depends on how much the input params will vary based on the model
  #for now, lets assume it will work with a single sser_classifyR
  
  
  # data=NULL, formula=NULL, model=NULL,
  # method="lrm", rsmpmethod="rcv",
  # k=5, times=100, rsmp=100,
  # id.include=NULL, ssercutoff=0.2,
  # tier = NULL, plotlabel = "",
  # seed=1, verbose=T,runtestorruntests,classes,  params,leave,percent,  resubstituteParams,     minimumOverlapPercent,       validation,
  # parallelParams,        easyDatasetID ,       hardDatasetID,       featureSets ,       metaFeatures ,
  # datasetName ,       classificationName ,   training, testing
  
  #errormessages = capture.output(
    SSER <- sser_classifyR(model=model, MAEobject=data,
                     method=rsmpmethod, k=k, times=times, rsmp=rsmp,
                     seed=seed, classes=classes, params, leave, percent, resubstituteParams, 
                     minimumOverlapPercent, validation, easyDatasetID, hardDatasetID,
                     featureSets,metaFeatures, datasetName, classificationName, training, 
                     testing, verbose=verbose, classIndex, z)
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
    dplyr::select(SampleID) %>% 
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
  
  # print("Cutoff")
  # print(TSERcutoff$TSER)
  
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






