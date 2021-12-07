library("sas7bdat")
library("dplyr")
library("mice")

load_imputed_data <- function(fname, seed = 944) {
  dat_cc <- NULL
  
  if (file.exists(fname)) {
    message("Loading previously saved object ",fname)
    load(fname)
    return(dat_cc)
  } 
  
  dsout <- NULL
  if (file.exists("F:/Projects/Datasets/Weights disease risk score/IDB_control.rda")) {
    load("F:/Projects/Datasets/Weights disease risk score/IDB_control.rda")
  } else {
    rawdata <- sas7bdat::read.sas7bdat("F:/Projects/Datasets/Weights disease risk score/ades.sas7bdat") 
    patinfo <- sas7bdat::read.sas7bdat("F:/Projects/Datasets/Weights disease risk score/adsl.sas7bdat")
    #progdata <- sas7bdat::read.sas7bdat("F:/Projects/Datasets/Weights disease risk score/adprog.sas7bdat")
    
    dsout <- subset(rawdata, STUDYID %in% c("C-1801","105MS301", "109MS301", "109MS302") & AVISITN >= 0)
    dsout$TRIAL <- factor(dsout$STUDYID, levels = c("C-1801","105MS301", "109MS301", "109MS302"), 
                          labels = c("AFFIRM", "ADVANCE",  "DEFINE", "CONFIRM"))
    dsout <- left_join(dsout, patinfo, by = "USUBJID", keep = TRUE)
    save(dsout, file = "F:/Projects/Datasets/Weights disease risk score/IDB_control.rda")
  }
  
  # Convert variables
  dsout <- dsout %>% transform(WhiteRace = ifelse(dsout$RACE == "WHITE", 1, 0),
                               MaleGender = ifelse(dsout$SEX == "MALE", 1, 0))

  

  
  set.seed(seed)
  
  # Start multiple imputation
  impvars <- c("STUDYID.x", "USUBJID.x", "SUBJID.x", "TRIAL", "TRTA", "TRTAN", 
               "VISIT", "AVISIT", "AVISITN", "AVAL", 
               "AGE", "MaleGender", "WhiteRace", "BASE", "HEIGHTBL", "WEIGHTBL", "ONSYRS", 
               "DIAGYRS", "PRDMTGR", "PRMSGR", "RLPS1YR", "RLPS3YR", "GDLESBL", "T1LESBL", "T2LESBL",
               "NHPTMBL", "PASATABL", "T25FWABL", "EDSSBL", "TRELMOS")
  mice.prep <- mice(dsout[,impvars], maxit = 0)
  
  pM <- mice.prep$predictorMatrix
  pM[, "SUBJID.x"] <- -2
  pM[, c("STUDYID.x", "USUBJID.x", "AVISIT")] <- 0
  
  fit <- mice(dsout[,impvars], method = mice.prep$method, predictorMatrix = pM, m = 1)
  
  dat_cc <- complete(fit,1)
  save(dat_cc, file = fname)
  
  return(dat_cc)
}

prepare_nrs <- function(data, # List with data sets
                        STUDYID_control = c("DEFINE", "CONFIRM"),  # Merge DEINFE + CONFIRM
                        STUDYID_treat = c("AFFIRM"), # AFFIRM
                        SUBGROUP = FALSE,
                        END_VISIT = c("WEEK 36", "VISIT 9 WK 36")
) { 
  
  # Set up control group (DEINFE + CONFIRM)
  ds <- subset(data, TRIAL %in% c(STUDYID_control, STUDYID_treat) & TRTA == "Placebo" & AVISITN == 0)
  EDSS_END <-  subset(data, TRIAL %in% c(STUDYID_control, STUDYID_treat) & TRTA == "Placebo" & VISIT %in% END_VISIT)[,c("SUBJID.x", "AVAL")]
  
  ds_full <- (left_join(ds, EDSS_END, "SUBJID.x"))
  ds_full <- na.omit(ds_full)
  
  ds_full$Trt <- NA
  ds_full$Trt[which(ds_full$TRIAL %in% STUDYID_control)] <- 0
  ds_full$Trt[which(ds_full$TRIAL %in% STUDYID_treat)] <- 1
  
  ds_full$y <- ds_full$AVAL.y
  
  # DO we need to take a subgroup?
  if (SUBGROUP) {
    ds_full <- subset(ds_full, (BASE > 2.50 & Trt == 1) | Trt == 0)
  }
  
  
  ds_full
}