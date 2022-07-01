set.seed(1)

# Load libraries
# library(SummarizedExperiment)
library(MultiAssayExperiment)
library(tidyverse)
library(magrittr)
library(rms)
library(Hmisc)
library(skimr)
library(knitr)
library(caret)
library(grid)
library(gridExtra)


# Load data
# setwd("/Users/andywang/Dropbox (Sydney Uni)/aw2021PhD/melanoma/")
# load("./melanomaAssays.RData")
load("./melanomaAssaysNorm.RData")
source("./MTDT-0-fxn-SSER.R")
source("./MTDT-0-fxn-MTblock.R")
source("./MTDT-0-fxn-others.R")
source("./MTDT-0-fxn-MTDT.R")
source("./MTDT-1-datachecking.R")

# setting data
dat = dat %>% 
  dplyr::filter(!is.na(class2vs2years)) %>% 
  dplyr::mutate(
    y.good = ifelse(class2vs2years=="Good", 1, 0),
    Female = ifelse(Person_Sex=="Female", 1, 0))


tierClin.d = dat %>% 
  dplyr::select(Person_ID, y.good, Age, Female, 
                Prim_SiteBinary, Prim_SunExpSite) %>% 
  as_tibble()
tierClin.formula = (y.good ~ Age+Female+Prim_SiteBinary+Prim_SunExpSite)


tierHisto.d = dat %>% 
  dplyr::select(Person_ID, y.good, 
                Prim_Stage, Prim_Nodular, 
                Tum_BRAFmut, Tum_NRASmut, 
                Tum_LargeCellSize, Tum_Pigment) %>% 
  as_tibble()
tierHisto.formula = (y.good ~ 
                    Prim_Stage + Prim_Nodular + Tum_BRAFmut + Tum_NRASmut +
                    Tum_LargeCellSize + Tum_Pigment)


tierNano.d <- t(melanomaAssaysNorm@ExperimentList@listData$NanoString) %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var="Tumour_ID") %>% 
  dplyr::left_join(dat %>% dplyr::select(Person_ID, y.good, Tumour_ID),
                   by = "Tumour_ID") %>% 
  dplyr::select(Person_ID, y.good, everything(.), -Tumour_ID)
tierNano.d %>%  
  tidyr::pivot_longer(c(-Person_ID, -y.good), names_to = "gene", values_to = "expr") %>% 
  dplyr::mutate(y.good=factor(y.good, labels=c("No", "Yes"))) %>% 
  dplyr::filter(!is.na(y.good)) %>% 
  ggplot(aes(x = expr,
             group = Person_ID %>% as.factor,
             colour = y.good)) + 
  stat_density(position = "identity", geom = "line") +
  labs(title = "Expression density plot overlaid by samples") +
  facet_grid(y.good~.) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

#######






  
  

