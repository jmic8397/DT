
MTDT.Plot <- Function(Summary, Results){

  MTDTsummary = Summary
  
  
#Modify for multiple cutoffs
#if multiple cutoffs...

# Scatterplot with total model error
theme_set(theme_bw())  # pre-set the bw theme.
g <- ggplot(MTDTsummary$TSER.summary, aes(((TSA.retained * N.retained ) + (TSA.notretained * N.notretained))/N.Total , Cost.Total)) + 
  labs(subtitle="Accuracy and Cost", title="Bubble Plot", x = "Total Sequence Accuracy", y = "Total Cost")
g + geom_jitter(aes(col=Tier.Sequence))
#add extra layer of jitter for each tier??



#INCLUDED CUTOFF IN TITLE
cap1 = paste("Table 1, Cutoff:", ssercutoffList[1], sep = " ", collapse = NULL)
MTDTsummary$TSER.summary %>% 
  DT::datatable(option=list(
    columnDefs=list(list(className="dt-center", targets=c(1, 8))),
    pageLength=20),caption =cap1 ) %>% 
  DT::formatRound(columns=c(2, 3, 7, 9), digits=2)



nperms = dim(MTDTresults$myperms)[1]



for (nperm in 1:nperms){
  
  
  #Print tree plot
  
  treeRoot = MTDTsummary$Plots.Tree[[nperm]]
  # plot(treeRoot)
  
  print(treeRoot, "counter")
  
  
  ############################
  
  n = nrow(colData(MultiAssayExperiment))
  
  GetEdgeLabel <- function(node) {
    if (!node$isRoot ) {
      label = paste0( ((node$counter/n)*100),  "% (", node$counter, ")")
    } else {
      label = node$name
    }
    return (label)
    
  }
  
  
  GetFillColour <- function(node){
    if(!is.null(node$tier)){
      name = node$tier
      
      name = name %>% charToRaw %>% as.numeric %>%paste(collapse = '')
      
      return(paste0('#',((strtoi(name)*255)%%100),'FF89'))}
    else{
      return("LightBlue")
    }
    
    
  }
  
  SetGraphStyle(treeRoot, rankdir = "LR")
  SetEdgeStyle(treeRoot, fontname = 'helvetica', label = GetEdgeLabel)
  SetNodeStyle(treeRoot, style="filled", shape = "box", fillcolor = GetFillColour, fontname = 'helvetica')
  Graph = ToDiagrammeRGraph(treeRoot)
  graphname = paste0("graph",nperm, ".png")
  export_graph(Graph, file_name = graphname, file_type = "png")
  graphpath = paste0(getwd(), "/",graphname)
  cat("![](",graphname,")")
  
  
  
  
  #strat plot of all tiers
  p1 = MTDTsummary$Plots.strat.prop[[nperm]] +
    theme(plot.title=element_text(face="bold"))
  #Tier specific error box plots
  p2 = MTDTsummary$Plots.TSERcutoff[[nperm]][[1]]
  p3 = MTDTsummary$Plots.TSERcutoff[[nperm]][[2]] + ylab("")
  p4 = MTDTsummary$Plots.TSERcutoff[[nperm]][[3]] + ylab("")
  #overall error
  p5 = MTDTsummary$Plots.TSER[[nperm]] + 
    ggtitle("Overall") + ylab("") +
    theme(plot.title=element_text(hjust=0.5, face="bold"))
  
  
  
  #MTDTsummary$Plots.Tree[[nperm]]
  
  gridExtra::grid.arrange(
    grobs = list( p1, p2, p3,p4,  p5),
    layout_matrix = rbind(c(1,1,1,1),
                          c(1,1,1,1),
                          c(2,3,4,5),
                          c(2,3,4,5),
                          c(2,3,4,5))
  )
  
  #rpart?
  library(rpart)
  
}
}