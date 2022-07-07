

#classifyR version
#we will use the function in producting the MTDT summary
#includes all variables needed for advanced runtests
#Assumes MAE being passed in
#maybe classificationname should be referenced from modelslist???? ******************
MTDT.algmClassifyR <-   function(MultiAssayExperiment,
                                 tierList,fixedTier = NULL,
                                 ssercutoffList,tierUnitCosts = c(100, 500, 1000), 
                                 performanceType = "Sample Error", runtestorruntests = "runtest",
                                 classes = NULL, crossValParams = NULL, modellingParams = NULL, 
                                 characteristics = NULL, 
                                 seed=1, verbose=F){
  
  
  
  

  #how many different datasets
  nblocks = length(tierList)
  
  
  #returns matrix of permutations of models
  myperms = perm(v=1:nblocks)
  
  
  #Fix tier if specified (e.g. Clinical first)
  if((!(is.null(fixedTier))) & nblocks>1){
    #adjust perms for clinical as first tier
    for(tier in 1:length(tierList)){
      if(tierList[tier] == fixedTier ){
        #swap fixed tier and tier 1, fixedTier is now tier 1
        temp = tierList[1]
        tierList[1] = tierList[tier]
        tierList[tier] = temp
        #concat a column of tier 1 with the remaining permutations to enforce this new fixing of 'fixedTier'
        myperms = data.frame(rep(1,nblocks-1),(perm(1:(nblocks-1))+1))
        nblocks = nblocks-1
        
      }
    }
  }
  
  
  MTlist = list()
  
  
  #errorchecking for inputs here
  

  #RuntestS

  #OK,  runtests was selected
  #also need alternative function for runtest() (no cv) and runtesteasyhard() ***********************************
  if(runtestorruntests=="runtests"){
    
    #for each permutation (row)
    for (nperm in 1:dim(myperms)[1]) {
      
      #create a new list for results
      MTunits = list()
    
      #create list for # retained for this permutation
      id.retained = list()
      
      #create sequence for current tier
      tierseq = NULL
      
      #for each tier of each permutation
      for (ntier in 1:dim(myperms)[2]){
        #set tierseq and tier as per myperms #
        tierseq = c(tierseq, tierList[[myperms[nperm, ntier]]])
      }
      
      #create tier-sequence string
      tierseq = paste0(tierseq, collapse="-")
      
      #for each tier of each perm
      for (ntier in 1:dim(myperms)[2]){
        
        #Final tier check
        # this will later enforce that all samples are allocated to a leaf node regardless of error
        finalTier = FALSE
        if(ntier == dim(myperms)[2]){
          finalTier = TRUE
        }
        
        #z = which layer of the MAE as per permutations (e.g. 1=clinical, 2=histo, perm 1-2 = clin-histo)
        z = myperms[nperm, ntier]
        

        
        #Can probably use names(experiments(measurements)) + clinical for tiers
        tier = tierList[[z]]
        ssercutoff = ssercutoffList[[z]]
        
        if(finalTier){
          ssercutoff = 1
        }
        
        
        #get the column index of the class from the MAE
        classIndex = grep(classes,colnames(colData(melanomaAssaysNorm)))
        
        #Subsetting NA's
        MAEwithoutRetained = MultiAssayExperiment
        retainList = unlist(id.retained)
        retainListLogical = !(rownames(colData(MAEwithoutRetained)) %in% retainList )
        
        if(length(id.retained)>0){
            MAEwithoutRetained = MAEwithoutRetained[,retainListLogical ,]
        }
        
        #Create tierstring
        #e.g. (1) Clin-Histo-Nano: <Clin>
        print(paste0("(", nperm, ") ", tierseq, ": <", tier, ">"))
        # print(paste0("nperms = ", nperm, " ; ntiers = ", ntier, " ; tier = ", tier)) #############<-------------------------
        
        
        #pass model and data to MTblock which returns retained status of current sequence layer
        MTunits[[ntier]] = MTblockClassifyR(
          data=MAEwithoutRetained, id.retained=id.retained,
          ssercutoff=ssercutoff, tier=tier, plotlabel=tier,
          seed=seed, verbose, runtestorruntests=runtestorruntests,
          classes=classes, crossValParams=crossValParams, modellingParams=modellingParams, 
          characteristics=characteristics, performanceType=performanceType, finalTier=finalTier, 
           classIndex, z)

        
        #retained is now current + new retained
        id.retained = c(id.retained, MTunits[[ntier]]$id$id.retained)

        retained = MTunits[[ntier]]$id$id.retained %>% unique() %>% length()
        toprogress = MTunits[[ntier]]$id$id.toprogress %>% unique() %>% length()
        notprocessed = MTunits[[ntier]]$id$id.notprocessed %>% unique() %>% length()
        total = retained+toprogress+notprocessed
        processed=retained+toprogress
        
        
        #not enough samples - terminate early
        print(paste0("Size Check"))
        print(paste0(total - retained))
        
        # if(total - retained < 10){
        #   
        # }
        # 

        
        print(paste0("    Total = ", total,  ""))
        print(paste0("    Processed =", processed, 
                     " (", retained, " retained, ",
                     toprogress, " to progress to next tier)"))
        print(paste0("    Not processed = ", notprocessed))
        
        #not enough samples for another tier
        #to set user defined threshold
        # if(total - retained < 10 ){
        #   print(paste0("Less than 10 to progress, Tree terminating early"))
        #   for(remainingTier in ntier:dim(myperms)[2]){
        #     print(remainingTier)
        #   MTunits[[remainingTier]]$id$id.retained = c(MTunits[[remainingTier]]$id$id.retained,MTunits[[remainingTier]]$id$id.toprogress,MTunits[[remainingTier]]$id$id.notprocessed)
        #   MTunits[[remainingTier]]$id$id.toprogress = tibble(SampleID=NULL)
        #   }
        #   break
        # }

        
      }
      
      MTlist[[nperm]] = MTunits
      
    }
    
    return(MTDTobject = list(MTlist=MTlist, myperms=myperms, dataList=MultiAssayExperiment,
                             tierList=tierList, ssercutoffList=ssercutoffList))
  }
}


MTDT.algmCost <- function(dataList, rsmpList, tierList,
  ssercutoffList, k=5, times=100, rsmp=100, seed=1, verbose=F){
  
  #how many different datasets
  nblocks = length(dataList)
  
  if(class(dataList) == "MultiAssayExperiment"){

        nblocks = length(experiments(measurements))
  }
  
  #returns matrix of permutations of models
  myperms = perm(1:nblocks)
  
  MTlist = list()

  #create a new list
  MTunits = list()
  
  #create # retained
  id.retained = NULL
  
  #create sequence
  tierseq = NULL
  
  #for each col of each row 
  for (ntier in 1:dim(myperms)[2]){
    #set tierseq and tier as per myperms #
    tierseq = c(tierseq, tierList[[myperms[nperm, ntier]]])
  }
  
  tierseq = paste0(tierseq, collapse="-")
  
  #for each col of each row
  for (ntier in 1:dim(myperms)[2]){
    
    #z = which # permutation
    z = myperms[nperm, ntier]
    data = dataList[[z]]
    # model = modelList[[z]]
    rsmpmethod = rsmpList[[z]]
    tier = tierList[[z]]
    ssercutoff = ssercutoffList[[z]]
    
    
    #e.g. (1) Clin-Histo-Nano: <Clin>
    print(paste0("(", nperm, ") ", tierseq, ": <", tier, ">"))
    # print(paste0("nperms = ", nperm, " ; ntiers = ", ntier, " ; tier = ", tier)) #############<-------------------------
    
    
    #pass model and data to MTblock which returns retained status of current layer
    MTunits[[ntier]] = MTblockClassifyR( data=data, id.retained=id.include,
      ssercutoff=ssercutoff, tier=tier, plotlabel=tier, runtestorruntests, classes = classes,
      crossValParams = crossValParams, modellingParams = modellingParams, characteristics = characteristics,
      performanceType = performanceType,seed=seed, verbose=verbose,finalTier,classIndex, z)

    
        
    id.retained = c(id.retained, MTunits[[ntier]]$id$id.retained)
    

    
    retained = MTunits[[ntier]]$id$id.retained %>% unique() %>% length()
    toprogress = MTunits[[ntier]]$id$id.toprogress %>% unique() %>% length()
    notprocessed = MTunits[[ntier]]$id$id.notprocessed %>% unique() %>% length()
    total = retained+toprogress+notprocessed
    processed=retained+toprogress
    print(paste0("    Total = ", total,  ""))
    print(paste0("    Processed = ", processed, 
                 " (", retained, " retained, ",
                 toprogress, " to progress to next tier)"))
    print(paste0("    Not processed = ", notprocessed))
   
    
  }
  
  MTlist[[nperm]] = MTunits
  

  return(MTDTobject = list(MTlist=MTlist,
                           myperms=myperms,
                           dataList=dataList,
                           
                           tierList=tierList,
                           ssercutoffList=ssercutoffList))
}





MTDT.ClassifyR.summary <- function(MTDTobject){

  myperms=MTDTobject$myperms
  MTlist=MTDTobject$MTlist
  tierList=MTDTobject$tierList
  
  dataList=MTDTobject$dataList
  
  if(class(dataList) == "MultiAssayExperiment"){
    dataList = colData(dataList)
  }
  
  ntiers=dim(myperms)[2]
  nperms=dim(myperms)[1]
  
  tserplots = list()
  stratplots = list()
  treeplots = list()
  stratplots.prop = list()
  TSERcutoffplots = list()
  TSERcutoffplots.tiers = list()
  
  TSER.overall <- list()
  strat.overall <- list()
  id <- list()
  tiersequence <- vector()
  TSERoverall.retained <- vector()
  TSERoverall.notretained <- vector()
  N.retained <- vector()
  N.notretained <- vector()
  
  for (nperm in 1:nperms){
    
    id.retained = NULL
    id.all = NULL
    id.notretained = NULL
    tierlabel = NULL
    seqname = paste0("Sequence.", nperm)
    TSERcutoffplots[[nperm]] = seqname
    
    TSER.retained = NULL
    TSER.notretained = NULL
    
    for (ntier in 1:ntiers){
      

      id.retained = c(id.retained,
                      MTlist[[nperm]][[ntier]]$id$id.retained)
   
      id.retained = id.retained %>% unique()

      id.all = c(id.all,
                 MTlist[[nperm]][[ntier]]$id$id.retained,
                 MTlist[[nperm]][[ntier]]$id$id.toprogress,
                 MTlist[[nperm]][[ntier]]$id$id.notprocessed)
      id.all = id.all %>% unique()
      id.notretained = setdiff(id.all, id.retained)
      tierlabel = c(tierlabel, tierList[[myperms[nperm, ntier]]])
      
      
      IDnm = colnames(dataList)[1]
      ID_nm = sym(IDnm)


      TSER.retained = dplyr::bind_rows(
        TSER.retained,
        MTlist[[nperm]][[ntier]]$SSER$SSER %>% 
          dplyr::filter(SampleID %in% id.retained)
      )
      TSER.notretained = dplyr::bind_rows(
        TSER.notretained,
        MTlist[[nperm]][[ntier]]$SSER$SSER  %>%
          dplyr::filter(SampleID %in% id.notretained)
      )
    }
    
    # print(TSER.retained)
    # print(TSER.notretained)
    

    tierlabel = paste(tierlabel, collapse = '-')
    
    id[[nperm]] = list(id.retained=id.retained, 
                       id.notretained=id.notretained)
    tiersequence[[nperm]] = tierlabel
    
    ynm = colnames(TSER.retained)[2]
    y_nm = ynm
  
    
    TSER.overall[[nperm]] = dplyr::bind_rows(
      TSER.retained %>% 
        dplyr::summarise(tser=(1-mean(sser, na.rm=T))) %>% 
        dplyr::mutate(strata=factor("Retained", levels=c("Retained", "Not retained"))),
      TSER.notretained %>% 
        dplyr::summarise(tser=(1-mean(sser, na.rm=T))) %>% 
        dplyr::mutate(strata=factor("Not retained", levels=c("Retained", "Not retained")))
    )
    print(TSER.overall[[nperm]])
    
    tserplots[[nperm]] = tsercutoff_plot(TSER.overall[[nperm]]) + 
      ggtitle(paste0("Overall: ", tierlabel))
    
    for (ntier in 1:ntiers){
      TSERcutoffplots.tiers[[ntier]] = MTlist[[nperm]][[ntier]]$plots$tsercutoffplot +
        ggtitle(tierList[[myperms[nperm, ntier]]])
      
    }
    
    TSERcutoffplots[[nperm]] <- TSERcutoffplots.tiers
    
    n.retained=0
    lvls = NULL
    stra = NULL

    for (ntier in 1:ntiers){
      print(ntier)
      
      n.from=n.retained + 1
      n.end = n.retained +
        dim(MTlist[[nperm]][[ntier]]$TSERcutoff$Stratification)[1]

      if(dim(MTlist[[nperm]][[ntier]]$TSERcutoff$Stratification)[1] == 0){
        print("Empty, skipping tier")
        next
      }
      
      stra.MTunit <- MTlist[[nperm]][[ntier]]$TSERcutoff$Stratification %>% 
        dplyr::arrange(strata, sser) %>% 
        dplyr::mutate(tier=tierList[[myperms[nperm, ntier]]]) %>% 
        tibble::tibble(., ID = n.from:n.end)
      
      print(stra.MTunit)
      
      n.retained = n.retained +
        MTlist[[nperm]][[ntier]]$TSERcutoff$Stratification %>% 
        dplyr::filter(strata=="Retained") %>% nrow()
      
      
      stra = dplyr::bind_rows(stra,  stra.MTunit)
      
      lvls = c(lvls,
               paste0("Retained (", tierList[[myperms[nperm, ntier]]], ")"), 
               paste0("To progress (", tierList[[myperms[nperm, ntier]]], ")"),
               paste0("Not processed (", tierList[[myperms[nperm, ntier]]], ")"))
    }
    
    stra = stra %>% 
      dplyr::mutate(
        tierstrata=paste0(strata, " (", tier, ")"),
        tierstrata=factor(tierstrata, levels=lvls)) %>% 
      dplyr::arrange(ID, tierstrata, sser)
    
    stratplots.prop[[nperm]] = stra %>% 
      ggplot() +
      geom_tile(aes(y=tierstrata, x=ID, fill=sser), colour="black") +
      ylab("") + xlab("") + ggtitle(paste("Sequence:", tierlabel, "\n Sample Error Cutoff: ", MTDTobject$ssercutoffList[1])) +
      scale_fill_gradient(low="black", high="white") +
      theme(panel.background = element_blank(),
            axis.text.x = element_blank(),
            aspect.ratio=1/4)
    
    library(dplyr)
    
    
    #############################################################################################
    #create tree 
    
    treeRoot <- Node$new(stra$tier[1])
    currNode = treeRoot
    
    for(ntier in 1:ntiers+1){
      
      #Retained
      num1 = stra %>% filter(tierstrata == levels(stra$tierstrata)[(((ntier-1)-1)*3)+1]) %>% nrow()
      tier = strsplit(levels(stra$tierstrata)[(((ntier-1)-1)*3)+1], '[()]')[[1]][2]
      if(num1>0){
        #finalTier Check -> retain all if final tier
        if(ntier == ntiers+1 || ((dim(MTlist[[nperm]][[ntier-1]]$TSERcutoff$Stratification)[1] == 0))){
          num1 = num1 + (stra %>% filter(tierstrata == levels(stra$tierstrata)[(((ntier-1)-1)*3)+3]) %>% nrow() )
        }
      retained = currNode$AddChild(levels(stra$tierstrata)[  (((ntier-1)-1)*3)+1], counter = num1, tier = tier, tierstrata = levels(stra$tierstrata)[(((ntier-1)-1)*3)+1], tierOrder = ntier-1)
      }
      
      #To Progress
      num2 = stra %>% filter(tierstrata == levels(stra$tierstrata)[(((ntier-1)-1)*3)+2]) %>% nrow()
      if(num2>0){
      progress = currNode$AddChild(levels(stra$tierstrata)[  (((ntier-1)-1)*3)+2], counter = num2,tier = tier, tierstrata = levels(stra$tierstrata)[(((ntier-1)-1)*3)+2], tierOrder = ntier-1)
      }else{
        # currNode$AddChild(levels(stra$tierstrata)[  (((ntier-1)-1)*3)+2], counter = 0,tier = tier, tierstrata = levels(stra$tierstrata)[(((ntier-1)-1)*3)+2], tierOrder = ntier-1)
        break;
      }
      
      #Not Processed
      num3 = stra %>% filter(tierstrata == levels(stra$tierstrata)[(((ntier-1)-1)*3)+3]) %>% nrow() 
      if(num3>0){
      notprog = currNode$AddChild(levels(stra$tierstrata)[  (((ntier-1)-1)*3)+3], counter = num3,tier = tier, tierstrata = levels(stra$tierstrata)[(((ntier-1)-1)*3)+3], tierOrder = ntier-1)
      }
      
      currNode = progress 
    }

    treeplots[[nperm]] = treeRoot
    
    ####################################################################################################
    
    IDnm = colnames(stra)[1]
    ID_nm = sym(IDnm)
    stra = stra %>% 
      dplyr::arrange(strata)
    idsort = stra[[1]] %>% unique()
    stratplots[[nperm]] = stra %>% 
      dplyr::arrange(strata) %>% 
      dplyr::mutate(ID=factor(!!ID_nm, levels=idsort)) %>% 
      ggplot() +
      geom_tile(aes(y=tierstrata, x=ID, fill=sser), colour="black") +
      ylab("") + xlab("") + ggtitle(paste("Sequence:", tierlabel)) +
      scale_fill_gradient(low="black", high="white") +
      theme(panel.background = element_blank(),
            axis.text.x = element_blank(),
            aspect.ratio=1/4)
    strat.overall[[nperm]] = stra
    
    TSERoverall.retained[[nperm]] = TSER.overall[[nperm]] %>% 
      dplyr::filter(strata=="Retained") %>% 
      dplyr::group_by(strata) %>% 
      dplyr::summarise(mean=mean(tser, na.rm=T)) %$% mean
    
    TSERoverall.notretained[[nperm]] = TSER.overall[[nperm]] %>% 
      dplyr::filter(strata=="Not retained") %>% 
      dplyr::group_by(strata) %>% 
      dplyr::summarise(mean=mean(tser, na.rm=T)) %$% mean
    N.retained[[nperm]] = id[[nperm]]$id.retained %>% length()
    N.notretained[[nperm]] = id[[nperm]]$id.notretained %>% length()
    
  }
  
  TSER.summary <- tibble::tibble(
    Tier.Sequence = tiersequence,
    TSA.retained = (1-TSERoverall.retained),
    N.retained = N.retained,
    TSA.notretained = (1-TSERoverall.notretained),
    N.notretained = N.notretained,
    Threshold = MTDTobject$ssercutoffList[1]
  )
  

  return(list(TSER.summary = TSER.summary,
              IDs = id,
              Results.TSER.overall = TSER.overall,
              Results.stratification = strat.overall,
              Plots.TSER = tserplots,
              Plots.TSERcutoff = TSERcutoffplots,
              Plots.Tree = treeplots,
              Plots.stratification = stratplots,
              Plots.strat.prop = stratplots.prop))
}


MTperm.summary <- function(MTDTobject){
  
  myperms = MTDTobject$myperms
  MTlist = MTDTobject$MTlist
  tierList = MTDTobject$tierList %>% unlist
  ssercutoffList = MTDTobject$ssercutoffList
  tierUnitCosts = tierUnitCosts

  
  nperms = dim(myperms)[1]
  ntiers = dim(myperms)[2]
  N.in = NULL
  N.processed = NULL
  N.notprocessed = NULL
  tierseq = NULL
  MTperm.summary = NULL
  
  for (nperm in 1:nperms) {
    
    seq = NULL
    
    for (ntier in 1:ntiers) {

      coordinate = myperms[nperm, ntier]
      
      retained = MTDTobject$MTlist[[nperm]][[ntier]]$id$id.retained %>% unique() %>% length()
      toprogress = MTDTobject$MTlist[[nperm]][[ntier]]$id$id.toprogress %>% unique() %>% length()
      notprocessed = MTDTobject$MTlist[[nperm]][[ntier]]$id$id.notprocessed %>% unique() %>% length()
      total = retained+toprogress+notprocessed
      processed=retained+toprogress
      
      MTperm.summary = dplyr::bind_rows(
        MTperm.summary,
        tibble::tibble(
          nperm = nperm,
          ntier = ntier,
          tier = tierList[[coordinate]],
          # model = modelclass[[coordinate]],
          sser = ssercutoffList[[coordinate]],
          Cost.perUnit = tierUnitCosts[[coordinate]],
          N.in = total,
          N.retained = retained,
          N.toprogress = toprogress,
          N.processed = processed,
          N.notprocessed = notprocessed,
          Cost.subtotal = Cost.perUnit*N.processed
        )
      )
      seq = c(seq, tierList[[coordinate]])
    }
    Seq = str_c(seq, collapse="-")
    tierseq = c(tierseq, rep(Seq, ntiers))
  }
  
  MTperm.summary$Tier.Sequence = tierseq
  MTperm.summary = MTperm.summary %>% 
    dplyr::select(Tier.Sequence, everything(.))
  
  
  return(MTperm.summary)
}



MTDT.ClassifyR.cost.summary <- function(MTDTobject, tierUnitCosts){
  
  MTDTsummary = MTDT.ClassifyR.summary(MTDTobject)
  MTperm.summary = MTperm.summary(MTDTobject)
  strat.overall = MTDTsummary$Results.stratification
  
  
  myperms = MTDTobject$myperms
  tierList = MTDTobject$tierList
  ntiers = dim(myperms)[2]
  nperms = MTDTsummary$Results.stratification %>% length()
  Costs = NULL
  
  if(length(tierUnitCosts)!=ntiers) return(print("length of tierUnitCost not equal to number of tiers"))
  
  for (nperm in 1:nperms){
    tierseq = MTDTsummary$TSER.summary$Tier.Sequence[nperm] %>% str_split("-") %>% as_vector()
    
    # print(strat.overall[[nperm]])
    # print("perm")
    # print(dim(myperms)[1])
    # print(MTDTobject$myperms)
    
    numTiers = length(unique(strat.overall[[nperm]]$tier))
    # print("numtiers")
    # print(c(1:numTiers))
    # print(strat.overall[[nperm]]$tier)
    # print(unique(strat.overall[[nperm]]$tier))
    c = strat.overall[[nperm]] %>% 
      dplyr::filter(!strata=="Not processed") %>% 
      dplyr::count(tier) %>% 
      dplyr::mutate(
        cost=n*tierUnitCosts[c(myperms[nperm,][c(1:numTiers)])],
        tier=factor(tier, levels=tierseq)) %>% 
      dplyr::arrange(tier)
    
    print(c)

    Costs = dplyr::bind_rows(
      Costs,
      tibble(Tier.Sequence = str_c(MTDTsummary$TSER.summary$Tier.Sequence[nperm] %>% str_split("-") %>% as_vector(), collapse="-"),
             Cost.byTier = str_c(paste0("$", c$cost), collapse="-"),
             Cost.Total=sum(c$cost, na.rm=T))
    )
    print(Costs)
    
  }

  TSER.summary = MTDTsummary$TSER.summary %>% 
    dplyr::mutate(
      N.Total = N.retained + N.notretained,
      Prop.notretained = N.notretained/N.Total) %>% 
    dplyr::select(Tier.Sequence, 
                  TSA.retained, TSA.notretained,
                  N.Total, N.retained, N.notretained, 
                  Prop.notretained, Threshold) %>% 
    dplyr::left_join(Costs, by="Tier.Sequence")
  
  print(TSER.summary)

  return(list(TSER.summary = TSER.summary,
              MT.summary = MTperm.summary,
              IDs = MTDTsummary$IDs,
              Results.TSA.retained.overall = MTDTsummary$Results.TSER.overall,
              Results.stratification = MTDTsummary$Results.stratification,
              Plots.TSER = MTDTsummary$Plots.TSER,
              Plots.TSERcutoff = MTDTsummary$Plots.TSERcutoff,
              Plots.Tree = MTDTsummary$Plots.Tree, 
              Plots.stratification = MTDTsummary$Plots.stratification,
              Plots.strat.prop = MTDTsummary$Plots.strat.prop))
  
}

