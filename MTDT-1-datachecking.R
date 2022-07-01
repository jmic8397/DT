# dat=SummarizedExperiment::colData(melanomaAssays) %>% as.data.frame
# saveRDS(dat, file="MelanomaClinicoPath.RDS")
# dat <- readRDS("./MelanomaClinicoPath.RDS")
dat=SummarizedExperiment::colData(melanomaAssaysNorm) %>% as.data.frame
describe(dat)

# dat$c = ifelse(dat$Person_FUStatus==levels(dat$Person_FUStatus)[3:5], 1, 0)

dat <- dat %>% 
  dplyr::mutate(
    time = Prognosis_TimeSinceLNMet,
    status = case_when(
      Person_FUStatus %in% c("Alive NSR", 
                             "Alive with Melanoma", 
                             "Alive, status unknown", 
                             "Dead, cause unknown", 
                             "Dead, not melanoma") ~ 0,
      Person_FUStatus %in% c("Dead, melanoma") ~ 1,
      TRUE ~ NA_real_
    ),
    Age = Age_Analysis/365.25,
    Prim_Nodular = case_when(
      is.na(Prim_Histol) ~ NA_character_,
      Prim_Histol %in% c("Acral Lentiginous", 
                         "Desmoplastic",
                         "Desmoplastic with Neurotropia", 
                         "Lentigo Maligna Tumour",
                         "Malignant Blue Naevus",
                         "SSM") ~ "No",
      Prim_Histol %in% c("NM", 
                         "SSM with NM", 
                         "NM (Spitz Naevus)",
                         "NM (Minimal deviated melanoma)") ~ "Yes"
    ),
    Prim_Nodular = factor(Prim_Nodular),
    Prim_SunExpSite = case_when(
      Prim_Site_SunExp == "Nil" ~ 0,
      Prim_Site_SunExp == "Intermittent" ~ 0,
      Prim_Site_SunExp == "Chronic" ~ 1,
      TRUE ~ NA_real_
    ),
    Prim_SunExpSite = factor(Prim_SunExpSite, 
                             levels = c(0, 1), 
                             labels = c("Nil or intermittent", 
                                        "Chronic")),
    Prim_SiteBinary = case_when(
      dat$Prim_Site %in% c("Acral - Toe",
                           "Foot",
                           "Hand",
                           "Lower Arm",
                           "Lower Leg",
                           "Upper Arm",
                           "Upper Leg") ~ "extremity",
      dat$Prim_Site %in% c("Abdomen",
                           "Back",
                           "Buttock",
                           "Chest",
                           "Ear",
                           "Head",
                           "Mucosal",
                           "Occult",
                           "Acral - Sole") ~ "axis",
      is.na(dat$Prim_Site) ~ NA_character_
    ),
    Prim_Site=factor(Prim_Site),
    Prim_ThicknessCat=case_when(
      (dat$Prim_Breslow<0.76) ~ "<0.76mm",
      (dat$Prim_Breslow>=0.76 & dat$Prim_Breslow<1.70) ~ "0.76-1.69mm",
      (dat$Prim_Breslow>=1.70) ~ ">=1.70mm",
      is.na(dat$Prim_Breslow) ~ NA_character_
    ),
    Prim_ThicknessCat=factor(Prim_ThicknessCat),
    Tum_LargeCellSize = case_when(
      Tum_CellSize %in% c("0", "1", "1M") ~ "No",
      Tum_CellSize %in% c("2", "2M") ~ "Yes", 
      TRUE ~ NA_character_
    ),
    Tum_LargeCellSize = factor(Tum_LargeCellSize)
  )


dat$Prim_StageOri <- dat$Prim_Stage
levels(dat$Prim_Stage) <- c(levels(dat$Prim_Stage)[1:2], 
                            "Stage III/IV",
                            "Stage III/IV")

dat$Prim_RegressOri <- dat$Prim_Regress
levels(dat$Prim_Regress) <- c("Absent", 
                              "Present",
                              "Present",
                              "Present")



# density(dat$Age, na.rm=T) %>% plot
# dat$Prim_Site_SunExp %>% table(useNA="ifany") %>% kable
# dat$Prim_SunExpSite %>% table(useNA="ifany") %>% kable
# dat$Prim_StageOri %>% table(useNA="ifany") %>% kable
# dat$Prim_Stage %>% table(useNA="ifany") %>% kable
# density(dat$Prim_Breslow, na.rm=T) %>% plot
# dat$Prim_Histol %>% table(useNA="ifany") %>% kable
# dat$Prim_Nodular %>% table(useNA="ifany") %>% kable
# dat$Prim_Regress %>% table(useNA="ifany") %>% kable
# dat$Prim_Ulc %>% table(useNA="ifany") %>% kable
# density(dat$Prim_Mitos, na.rm=T) %>% plot
# density(dat$Tum_MetSize, na.rm=T) %>% plot
# dat$Tum_Extranodal %>% table(useNA="ifany") %>% kable
# dat$Tum_CellType %>% table(useNA="ifany") %>% kable
# dat$Tum_CellSize %>% table(useNA="ifany") %>% kable
# dat$Tum_LargeCellSize %>% table(useNA="ifany") %>% kable
# dat$Tum_Pigment %>% table(useNA="ifany") %>% kable
# 
# dat$Tum_BRAFmut %>% table(useNA="ifany") %>% kable
# dat$Tum_NRASmut %>% table(useNA="ifany") %>% kable
# dat$Tum_FLT3mut %>% table(useNA="ifany") %>% kable
# dat$Tum_METmut %>% table(useNA="ifany") %>% kable
# dat$Tum_PIK3CAmut %>% table(useNA="ifany") %>% kable 
# dat$Tum_PDGFRAmut %>% table(useNA="ifany") %>% kable 
# dat$Tum_EGFRmut %>% table(useNA="ifany") %>% kable 
# dat$Tum_CKITmut %>% table(useNA="ifany") %>% kable 


