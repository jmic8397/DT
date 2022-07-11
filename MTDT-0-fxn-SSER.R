

##########################################################################

#Sample specific error calculated using classifyR's predictions

sser_classifyR <- function( MAEobject=NULL,
                            classes, tier = NULL, crossValParams = NULL, 
                           modellingParams =NULL, characteristics = characteristics, seed=1, 
                           performanceType = "Sample Error", 
                           verbose=verbose,classIndex, minTierSize = 10, z){


  library(sparsediscrim)
  library(tidyverse)
  library(dplyr)
  library(e1071)
  
  
  set.seed(seed)
  

  #remove Na's for classifyR
  #this has already been done - does nothing unless samples with missing data exist
  if(tier != "sampleInfo"){
  MAECompleteCases = MAEobject[,complete.cases(MAEobject[,,tier]),tier]
  }
  else{
    MAECompleteCases= MAEobject
  }
  
  #errorTable
  errorTable = tibble(SampleID=character(),sser=double())
  ID.removed = list()

  

  # Check CV is valid
  if(length(unlist(colnames(MAECompleteCases[1]))) <= minTierSize){
    print(length(unlist(colnames(MAECompleteCases[1]))))

    print("Error: Low row count, unable to use CV")

  }else{

    DMresults <- runTests(MAECompleteCases, target=tier, outcomesColumns=classes, 
                          crossValParams = crossValParams, modellingParams = modellingParams[[z]], 
                          characteristics = characteristics, verbose =verbose)
    

    # This gives estimate error rate
    performance = calcCVperformance(DMresults, performanceType)
    error = performance(performance)
    
    #find NA's
    SampleIDComplete = names(error$`Sample Error`)
    SampleIDTotal = rownames(colData(MAEobject))
    ID.removed = setdiff(SampleIDTotal,SampleIDComplete )

    #construct Table
    SampleID = names(error$`Sample Error`)
    sser = unname(error$`Sample Error`)
    errorTable = tibble(SampleID,sser)
    view(errorTable)
    
  }
  
  if (verbose==T){
    if (length(ID.removed)>0) {
      print(paste("Incomplete data.", length(ID.removed), "IDs removed:"))
      # print(ID.removed)
    }
  }
  
  return(list(
              removed = ID.removed,
              SSER=errorTable))
}



#tier specific error
tser <- function(SSER=NULL, col_name=NULL){
  
  library(tidyverse)
  library(rms)
  
  SSER=SSER$SSER
  
  #Check non-empty tier
  if(length(dim(SSER)[1])!=0){
  
  #cutoffs are the average error rates for samples
  cutoffs = SSER %>% 
    dplyr::select({{col_name}}) %>% 
    distinct() %>% 
    as_vector() %>% 
    sort()
  
  #Initialize data structures
  names(cutoffs) <- NULL
  err = dplyr::tibble(tser=NULL, cutoff=NULL)
  n = vector()
  cindices = vector()
  briers = vector()

  }
  
  #Check non-empty rows
  if(dim(SSER)[1]!=0 ){
    
    #for each cutoff 
    for (i in 1:length(cutoffs)){
      
      #threshold = tser (1-avg(sser))
      threshold = cutoffs[i]
      
      #filter for less than threshold
      df.sser = SSER %>% dplyr::filter({{col_name}}<=threshold)
  
      #all below threshold
      tser = df.sser %>%
        #create new mean error for each threshold
        dplyr::summarise(tser=(1-mean(sser, na.rm=T))) %>%
        
        dplyr::mutate(cutoff=threshold)
  
      #n retained = rows left in sser <= threshold error *************** *&*
      if(dim(SSER)[1] - dim(df.sser)[1] > minTierSize){
      n = c(n, dim(df.sser)[1])
      err = dplyr::bind_rows(err, tser)
      }
      else{
        n = c(n, 0)
        err = dplyr::bind_rows(err, tser)
      }
    }
    }
  
    #Check nonempty rows
    if(dim(SSER)[1]!=0){
        
        TSER = tibble(
          Index = rep(as_label(enquo(col_name)), length(cutoffs)),
          cutoff = cutoffs,
          tse = err %>%
            dplyr::group_by(cutoff) %>% 
            dplyr::summarise(tser=mean(tser, na.rm=T)) %$% 
            tser,
          n_total = dim(SSER)[1],
          n_retained = n,
          prop_retained = n_retained/n_total,
          n_progress = n_total-n_retained,
          prop_progress = n_progress/n_total,
          )
        
        }
    
    
  #Check zero rows
  if(dim(SSER)[1]==0){
      TSER = tibble(
        Index = as_label(enquo(col_name)),
        cutoff = 1,
        tse = 1,
        n_total = 0,
        n_retained = 0,
        prop_retained = 0,
        n_progress = 0,
        prop_progress = 0)
    }
    

  return(TSER)
  }


# Plotting TSE - overall, included and excluded from the cutoff
tser_cutoff <- function(SSER=NULL, col_name=sser, mycutoff=0.5, finalTier, minTierSize=10){
  
  library(tidyverse)
  Removed = as.data.frame(SSER$removed)
  SSER=SSER$SSER
  
  #If this tier was processed
  if(dim(SSER)[1]>0){
  
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
  
  TSER_progress = SSER %>% 
    dplyr::filter(!{SampleID} %in% id.retained) %>% 
    dplyr::summarise(tser=(1-mean(sser, na.rm=T))) %>% 
    dplyr::mutate(strata="To progress",
                  n=n_progress)
  
  TSER_retained = SSER %>% 
    dplyr::filter({SampleID} %in% id.retained) %>% 
    dplyr::summarise(tser=(1-mean(sser, na.rm=T))) %>% 
    dplyr::mutate(strata="Retained",
                  n=n_retained)
  
  
  #if finalTier of tree, add all samples to retained
  #If less than min to progress, add all samples to retained
  if(finalTier == TRUE || ( n_progress < minTierSize) ){
    print("Too small")
    
    TSER_retained = SSER %>% 
      dplyr::summarise(tser=(1-mean(sser, na.rm=T))) %>% 
      dplyr::mutate(strata="Retained",
                    n=n)
    TSER_progress = tibble(strata="To progress",mean = 0,n=0)
    
    id.retained = SSER %>% select(SampleID) %>% as_vector() %>% unique()
    n_retained=length(id.retained)
    id.toprogress = NULL
    n_progress = 0
    
  }
  

  #Tier not processed, creating empty results object
  }else{
    id.retained = NULL
    id.toprogress = NULL
    n=0
    n_retained=0
    n_progress=0
    TSER_overall =  tibble(strata="Overall", mean = 0, n=0)
    TSER_retained = tibble(strata="Retained",mean = 0,n=0)
    TSER_progress = tibble(strata="To progress",mean = 0,n=0)
  }
  
  TSER_cutoff = bind_rows(TSER_overall, TSER_retained, TSER_progress) %>% 
    dplyr::mutate(
      strata=factor(strata, levels=c("Overall", "Retained", "To progress"))
    )
  
  print(TSER_cutoff)
  
  #If tier was processed
  if(dim(SSER)[1]>0){
  
  strata = dplyr::bind_rows(
    SSER %>% 
      dplyr::filter(SampleID %in% id.retained) %>% 
      dplyr::select(SampleID, sser) %>% 
      dplyr::mutate(strata = "Retained"),
    SSER %>% 
      dplyr::filter(SampleID %in% id.toprogress) %>% 
      dplyr::select(SampleID, sser) %>% 
      dplyr::mutate(strata = "To progress"),
    Removed %>%
      dplyr::mutate(strata = "Not processed")
    ) %>%
    dplyr::mutate(strata=factor(strata, levels=c("Retained", "To progress", "Not processed")))
  

  TSER_table = TSER_cutoff %>% 
    dplyr::group_by(strata) %>% 
    dplyr::summarise(mean(tser, na.rm=T), n=mean(n))
  
  #Else empty results object
  }else{
    
    strata = tibble(SampleID=character(),sser=double(),strata=factor(),Removed=character())
    TSER_table = TSER_cutoff %>% 
      dplyr::group_by(strata) %>% 
      dplyr::summarise(0, n=0)
    
  }
  
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

#Box plot
tsercutoff_plot <- function(TSERcutoff){
  print("PLOT")
  print("TSERcutoff")
  #Remove this if retaining entries with missing values
  TSERcutoff = TSERcutoff %>% filter(TSERcutoff$strata == "To progress" | TSERcutoff$strata == "Retained")
  p = TSERcutoff %>%
    ggplot() +
    #CHECING THIS &*&*
    geom_boxplot(aes(y=tser, x=strata, fill=strata)) +
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

  


  


