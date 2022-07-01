MTDT.algmCost <- function(dataList,
                      modelList,
                      rsmpList,
                      tierList,
                      ssercutoffList,
                      k=5, times=100, rsmp=100,
                      seed=1, verbose=F){
  
  #how many different datasets
  nblocks = length(dataList)
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
      
      #z = which # as per permutation
      z = myperms[nperm, ntier]
      data = dataList[[z]]
      model = modelList[[z]]
      rsmpmethod = rsmpList[[z]]
      tier = tierList[[z]]
      ssercutoff = ssercutoffList[[z]]
      modelclass = class(model)
      
      #assume ID first column, y second column
      IDnm = colnames(data)[1]
      ynm = colnames(data)[2]
      
      #sym takes strings and turns them into symbols
      ID_nm = sym(IDnm)
      y_nm = sym(ynm)
      
      #because id is a symbol now, !! tells dplyr to ignore its quotes and read its literal value instead
      #if id not in retained then we include it
      id.include = data %>% 
        dplyr::filter(!(!!ID_nm) %in% id.retained) %>% 
        dplyr::select(!!ID_nm) %>% 
        as_vector 
      names(id.include) <- NULL #set name of object to null? no col headings
      
      #model logic
      #print(modelclass)
      #modelclass
      #flush.console()
      #print("test")
      #flush.console()
      sys
      
      
      #linear model / linear regression
      if ("glm" %in% modelclass) {
        if ("lrm" %in% modelclass) {
          formula = model$sformula
          method = "lrm"
        } else {
          formula = model$formula
          method = "glm"
        }
      }
      
      #linear disc analysis
      if ("lda_diag" %in% modelclass) {
        formula = NULL
        method = "dlda"
        rsmpmethod = "boot"
      }
      
      #KNN
      ##############################
      if ("knn" %in% modelclass) {
        formula = NULL
        method = "knn"
        rsmpmethod = "rcv"
      }
      ############################
      
      #e.g. (1) Clin-Histo-Nano: <Clin>
      print(paste0("(", nperm, ") ", tierseq, ": <", tier, ">"))
      # print(paste0("nperms = ", nperm, " ; ntiers = ", ntier, " ; tier = ", tier)) #############<-------------------------
      
      
      #pass model and data to MTblock which returns retained status of current sequence layer
      MTunits[[ntier]] = MTblock(
        data=data, id.include=id.include,
        formula=formula, model=model, method=method, 
        rsmpmethod=rsmpmethod, k=k, times=times, rsmp=rsmp,
        ssercutoff=ssercutoff, tier=tier, plotlabel=tier,
        seed=seed, verbose=verbose
      )
      
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
                           modelList=modelList,
                           tierList=tierList,
                           ssercutoffList=ssercutoffList))
}



MTDT.summary <- function(MTDTobject){
  
  myperms=MTDTobject$myperms
  MTlist=MTDTobject$MTlist
  tierList=MTDTobject$tierList
  dataList=MTDTobject$dataList
  
  ntiers=dim(myperms)[2]
  nperms=dim(myperms)[1]
  # MTlist[[nperms]][[ntiers]]
  
  tserplots = list()
  stratplots = list()
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
    # id.retained.bytier = NULL
    # id.notretained.bytier = NULL
    tierlabel = NULL
    seqname = paste0("Sequence.", nperm)
    TSERcutoffplots[[nperm]] = seqname

    TSER.retained = NULL
    TSER.notretained = NULL

    for (ntier in 1:ntiers){
      
      id.retained = c(id.retained,
                      MTlist[[nperm]][[ntier]]$id$id.retained)
      id.retained = id.retained %>% unique()
      # id.retained.bytier[[nperm]][[ntier]] <- MTlist[[nperm]][[ntier]]$id$id.retained
      # id.notretained.bytier[[nperm]][[ntier]] <-c(MTlist[[nperm]][[ntier]]$id$id.toprogress,
      #                                             MTlist[[nperm]][[ntier]]$id$id.notprocessed)
      id.all = c(id.all,
                 MTlist[[nperm]][[ntier]]$id$id.retained,
                 MTlist[[nperm]][[ntier]]$id$id.toprogress,
                 MTlist[[nperm]][[ntier]]$id$id.notprocessed)
      id.all = id.all %>% unique()
      id.notretained = setdiff(id.all, id.retained)
      tierlabel = c(tierlabel, tierList[[myperms[nperm, ntier]]])
      
      
      #print(colnames(dataList[[ntier]]))
      IDnm = colnames(dataList[[ntier]])[1]
      ID_nm = sym(IDnm)

      #print(IDnm)
      #print("ID?")
      
      
      
      TSER.retained$y.good = TSER.retained$y.good %>% as.factor() %>% as.numeric()-1
      MTlist[[nperm]][[ntier]]$SSER$phat_rsmp$y.good = MTlist[[nperm]][[ntier]]$SSER$phat_rsmp$y.good %>% as.factor() %>% as.numeric()-1
      
      #print(TSER.retained$y.good)
      #print(MTlist[[nperm]][[ntier]]$SSER$phat_rsmp$y.good)
      
      TSER.retained = dplyr::bind_rows(
        TSER.retained,
        MTlist[[nperm]][[ntier]]$SSER$phat_rsmp %>% 
          dplyr::filter(!!ID_nm %in% id.retained)
      )
      TSER.notretained = dplyr::bind_rows(
        TSER.notretained,
        MTlist[[nperm]][[ntier]]$SSER$phat_rsmp %>%
          dplyr::filter(!!ID_nm %in% id.notretained)
      )
    }
    
    #print(TSER.retained)
    #print(id.retained)
    
    tierlabel = paste(tierlabel, collapse = '-')
    
    id[[nperm]] = list(id.retained=id.retained, 
                       id.notretained=id.notretained)
    tiersequence[[nperm]] = tierlabel
    
    
    #print(colnames(TSER.retained))
    
    ynm = colnames(TSER.retained)[2]
    y_nm = sym(ynm)
    
    #print("y?")
    #print(y_nm)
    
    ##CHECK THIS - GETTING 1 for retained and not retained?
    #print(colnames(TSER.retained))
    #print("Check")

    #print(TSER.retained$yhatrsmp)
    #print("Hi")

        
    TSER.overall[[nperm]] = dplyr::bind_rows(
      TSER.retained %>% 
        dplyr::group_by(rpt) %>% 
        dplyr::summarise(tser=(1-mean(yhatrsmp==y.good, na.rm=T))) %>% 
        dplyr::mutate(strata=factor("Retained", levels=c("Retained", "Not retained"))),
      TSER.notretained %>% 
        dplyr::group_by(rpt) %>% 
        dplyr::summarise(tser=(1-mean(yhatrsmp==y.good, na.rm=T))) %>% 
        dplyr::mutate(strata=factor("Not retained", levels=c("Retained", "Not retained")))
    )
    
    print(TSER.overall[[nperm]])

    tserplots[[nperm]] = tsercutoff_plot(TSER.overall[[nperm]]) + 
      ggtitle(paste0("Overall: ", tierlabel))
    
    for (ntier in 1:ntiers){
      TSERcutoffplots.tiers[[ntier]] = MTlist[[nperm]][[ntier]]$plots$tsercutoffplot +
          # ggtitle(paste0(tierlabel, " (", tierList[[myperms[nperm, ntier]]], ")"))
          ggtitle(tierList[[myperms[nperm, ntier]]])
        
    }
    TSERcutoffplots[[nperm]] <- TSERcutoffplots.tiers
    
    n.retained=0
    lvls = NULL
    stra = NULL
    ##
    for (ntier in 1:ntiers){
      n.from=n.retained + 1
      n.end = n.retained +
        dim(MTlist[[nperm]][[ntier]]$TSERcutoff$Stratification)[1]
      
      stra.MTunit <- MTlist[[nperm]][[ntier]]$TSERcutoff$Stratification %>% 
        dplyr::arrange(strata, sser) %>% 
        dplyr::mutate(tier=tierList[[myperms[nperm, ntier]]]) %>% 
        tibble::tibble(., ID = n.from:n.end)
      n.retained = n.retained +
        MTlist[[nperm]][[ntier]]$TSERcutoff$Stratification %>% 
        dplyr::filter(strata=="Retained") %>% nrow()
      

      stra$y.good = stra$y.good  %>%  as.factor() %>% as.numeric()-1
      stra$yhat = stra$yhat %>%  as.factor() %>% as.numeric()-1
      
      stra.MTunit$y.good = stra.MTunit$y.good %>%  as.factor() %>%as.numeric()-1
      stra.MTunit$yhat = stra.MTunit$yhat %>% as.factor() %>%as.numeric()-1
      
      stra = dplyr::bind_rows(stra,  stra.MTunit)
      
      lvls = c(lvls,
               paste0("Retained (", tierList[[myperms[nperm, ntier]]], ")"), 
               paste0("To progress (", tierList[[myperms[nperm, ntier]]], ")"),
               paste0("Not processed (", tierList[[myperms[nperm, ntier]]], ")"))
    }
    
    #print(stra.MTunit$tier)
    
    stra = stra %>% 
      dplyr::mutate(
        tierstrata=paste0(strata, " (", tier, ")"),
        tierstrata=factor(tierstrata, levels=lvls)) %>% 
      dplyr::arrange(ID, tierstrata, sser)
    
    stratplots.prop[[nperm]] = stra %>% 
      ggplot() +
      geom_tile(aes(y=tierstrata, x=ID, fill=sser), colour="black") +
      ylab("") + xlab("") + ggtitle(paste("Sequence:", tierlabel)) +
      scale_fill_gradient(low="black", high="white") +
      theme(panel.background = element_blank(),
            axis.text.x = element_blank(),
            aspect.ratio=1/4)
    
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
  

  # for (nperm in 1:nperms){
  #   for (ntier in  1:ntiers){
  #     tierList[[myperms[nperm, ntier]]] %>% print
  #     id.retained.bytier[[nperm]][[ntier]] %>% length %>% print
  #     id.notretained.bytier[[nperm]][[ntier]] %>% length %>% print
  #     # TSER.summary <- NULL
  #     # TSER.summary <- c(
  #     #   assign(tierList[[myperms[nperm, ntier]]], )
  #     # )    
  #     MTlist[[nperm]][[ntier]]$TSERcutoff$TSER %>% print
  #   }
  # }
  
  
  TSER.summary <- tibble::tibble(
    Tier.Sequence = tiersequence,
    TSER.retained = TSERoverall.retained,
    N.retained = N.retained,
    TSER.notretained = TSERoverall.notretained,
    N.notretained = N.notretained
  )
  
  # TSER.summary %>% print()
  
  return(list(TSER.summary = TSER.summary,
              IDs = id,
              Results.TSER.overall = TSER.overall,
              Results.stratification = strat.overall,
              Plots.TSER = tserplots,
              Plots.TSERcutoff = TSERcutoffplots,
              Plots.stratification = stratplots,
              Plots.strat.prop = stratplots.prop))
}


MTperm.summary <- function(MTDTobject){
  
  myperms = MTDTobject$myperms
  MTlist = MTDTobject$MTlist
  tierList = MTDTobject$tierList %>% unlist
  ssercutoffList = MTDTobject$ssercutoffList
  tierUnitCosts = tierUnitCosts
  modelclass = NULL
  for (i in 1:length(MTDTobject$modelList)){
    modelclass = c(modelclass, class(MTDTobject$modelList[[i]])[1])
  }  
  
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
          model = modelclass[[coordinate]],
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


MTDT.cost.summary <- function(MTDTobject, tierUnitCosts){
  
  MTDTsummary = MTDT.summary(MTDTobject)
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
    c = strat.overall[[nperm]] %>% 
      dplyr::filter(!strata=="Not processed") %>% 
      dplyr::count(tier) %>% 
      dplyr::mutate(
        cost=n*tierUnitCosts,
        tier=factor(tier, levels=tierseq)) %>% 
      dplyr::arrange(tier)
    
    Costs = dplyr::bind_rows(
      Costs,
      tibble(Tier.Sequence = str_c(c$tier, collapse="-"),
             Cost.byTier = str_c(paste0("$", c$cost), collapse="-"),
             Cost.Total=sum(c$cost, na.rm=T))
    )
  }
  
  TSER.summary = MTDTsummary$TSER.summary %>% 
    dplyr::mutate(
      N.Total = N.retained + N.notretained,
      Prop.notretained = N.notretained/N.Total) %>% 
    dplyr::select(Tier.Sequence, 
                  TSER.retained, TSER.notretained,
                  N.Total, N.retained, N.notretained, 
                  Prop.notretained) %>% 
    dplyr::left_join(Costs, by="Tier.Sequence")
  
  # TSER.summary %>% print()
  
  return(list(TSER.summary = TSER.summary,
              MT.summary = MTperm.summary,
              IDs = MTDTsummary$IDs,
              Results.TSER.overall = MTDTsummary$Results.TSER.overall,
              Results.stratification = MTDTsummary$Results.stratification,
              Plots.TSER = MTDTsummary$Plots.TSER,
              Plots.TSERcutoff = MTDTsummary$Plots.TSERcutoff,
              Plots.stratification = MTDTsummary$Plots.stratification,
              Plots.strat.prop = MTDTsummary$Plots.strat.prop))
         
}

