A.d = tierClin.d %>% 
  dplyr::mutate(Prim_SiteBinary=factor(Prim_SiteBinary))
B.d = tierHisto.d
C.d = tierNano.d

A.f = tierClin.formula
# A.f = y.good ~ 0 + Age + Female + Prim_SiteBinary + Prim_SunExpSite
B.f = tierHisto.formula

tierA = "Clin"
tierB = "Histo"
tierC = "Nano"

A.m = rms::lrm(A.f, data=A.d, x=T, y=T, se.fit=T)
# A.m = glm(A.f, data=A.d, family=binomial(logit))
B.m = rms::lrm(B.f, data=B.d, x=T, y=T, se.fit=T)
# B.m = glm(B.f, data=B.d, family=binomial(link="logit"))
y=C.d[,2] %>% as.factor; x=C.d[, -c(1,2)]
# C.m = sparsediscrim::dlda(y=y, x=x)
C.m = sparsediscrim::lda_diag(y=y, x=x)
rm(list=c("x", "y"))


dataList = list(A.d, B.d, C.d)
modelList = list(A.m, B.m, C.m)
rsmpList = list("rcv", "rcv", "boot")
tierList = list(tierA, tierB, tierC)
# ssercutoffList = c(0.2, 0.2, 0.2)
# tierUnitCosts = c(0, 500, 1000)
ssercutoffList = c(0.2, 0.3, 0.4)
tierUnitCosts = c(100, 500, 1000)


MTDTresults = MTDT.algm(dataList,
                        modelList,
                        rsmpList,
                        tierList,
                        ssercutoffList = ssercutoffList,
                        k=5, times=100, rsmp=100,
                        seed=1, verbose=F)
  
MTDTsummary = MTDT.cost.summary(MTDTresults, tierUnitCosts)


nperms = dim(MTDTresults$myperms)[1]

MTDTsummary$TSER.summary %>% 
  DT::datatable(option=list(
    columnDefs=list(list(className="dt-center", targets=c(1, 8))),
    pageLength=20)) %>% 
  DT::formatRound(columns=c(2, 3, 7, 9), digits=2)

for (nperm in 1:nperms){
  p1 = MTDTsummary$Plots.stratification[[nperm]] +
    theme(plot.title=element_text(face="bold"))
  p2 = MTDTsummary$Plots.TSERcutoff[[nperm]][[1]]
  p3 = MTDTsummary$Plots.TSERcutoff[[nperm]][[2]] + ylab("")
  p4 = MTDTsummary$Plots.TSERcutoff[[nperm]][[3]] + ylab("")
  p5 = MTDTsummary$Plots.TSER[[nperm]] + 
    ggtitle("Overall") + ylab("") +
    theme(plot.title=element_text(hjust=0.5, face="bold"))
  
  gridExtra::grid.arrange(
    grobs = list(p1, p2, p3, p4, p5),
    layout_matrix = rbind(c(1,1,1,1),
                          c(1,1,1,1),
                          c(2,3,4,5),
                          c(2,3,4,5),
                          c(2,3,4,5))
    )
  # ggsave(
  #   gridExtra::grid.arrange(
  #     grobs = list(p1, p2, p3, p4, p5),
  #     layout_matrix = rbind(c(1,1,1,1),
  #                           c(2,3,4,5),
  #                           c(2,3,4,5),
  #                           c(2,3,4,5))),
  #   filename=paste0("MTDTplots-seq", nperm, ".png")
  # )
}
  






for (nperm in 1:nperms){
  p1 = MTDTsummary$Plots.strat.prop[[nperm]] +
    theme(plot.title=element_text(face="bold"))
  p2 = MTDTsummary$Plots.TSERcutoff[[nperm]][[1]]
  p3 = MTDTsummary$Plots.TSERcutoff[[nperm]][[2]] + ylab("")
  p4 = MTDTsummary$Plots.TSERcutoff[[nperm]][[3]] + ylab("")
  p5 = MTDTsummary$Plots.TSER[[nperm]] + 
    ggtitle("Overall") + ylab("") +
    theme(plot.title=element_text(hjust=0.5, face="bold"))
  
  gridExtra::grid.arrange(
    grobs = list(p1, p2, p3, p4, p5),
    layout_matrix = rbind(c(1,1,1,1),
                          c(2,3,4,5),
                          c(2,3,4,5),
                          c(2,3,4,5))
  )
  # ggsave(
  #   gridExtra::grid.arrange(
  #     grobs = list(p1, p2, p3, p4, p5),
  #     layout_matrix = rbind(c(1,1,1,1),
  #                           c(2,3,4,5),
  #                           c(2,3,4,5),
  #                           c(2,3,4,5))),
  #   filename=paste0("MTDTplots(proportion)-seq", nperm, ".png")
  # )
}

save.image(file="MTDT.20210623.RData")
