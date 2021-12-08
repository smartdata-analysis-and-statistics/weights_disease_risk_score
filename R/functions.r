library("sas7bdat")
library("dplyr")
library("mice")
library(MatchIt)

# target distribution importance weighting
tdw <- function(a, x, data, estimand="ate",
                method = "linear", rule = 1, f = 0, ties = mean)
{
  d <- density(data[,x])
  f <- approxfun(d, method=method, rule=rule, f=f, ties=ties)
  d0 <- density(data[data[,a]==0,x])
  f0 <- approxfun(d0, method=method, rule=rule, f=f, ties=ties)
  d1 <- density(data[data[,a]==1,x])
  f1 <- approxfun(d1, method=method, rule=rule, f=f, ties=ties)
  
  if(estimand=="ate") {
    data[data[,a]==1,"w"] <- f(data[data[,a]==1,x])/f1(data[data[,a]==1,x])
    data[data[,a]==0,"w"] <- f(data[data[,a]==0,x])/f0(data[data[,a]==0,x])
  } else if(estimand=="att") {
    data[data[,a]==1,"w"] <- 1
    data[data[,a]==0,"w"] <- f1(data[data[,a]==0,x])/f0(data[data[,a]==0,x])
  } else if(estimand=="atu") {
    data[data[,a]==1,"w"] <- f0(data[data[,a]==1,x])/f1(data[data[,a]==1,x])
    data[data[,a]==0,"w"] <- 1
  }
  TDW <- data[,"w"]
  return(TDW)
}

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
    
    dsout <- subset(rawdata, STUDYID %in% c("C-1801","105MS301", "109MS301", "109MS302") & AVISITN >= 0)
    dsout$TRIAL <- factor(dsout$STUDYID, levels = c("C-1801","105MS301", "109MS301", "109MS302"), 
                          labels = c("AFFIRM", "ADVANCE",  "DEFINE", "CONFIRM"))
    dsout <- left_join(dsout, patinfo, by = "USUBJID", keep = TRUE)
    save(dsout, file = "F:/Projects/Datasets/Weights disease risk score/IDB_control.rda")
  }
  
  # Convert variables
  dsout <- dsout %>% transform(WhiteRace = ifelse(dsout$RACE == "WHITE", 1, 0),
                               MaleGender = ifelse(dsout$SEX == "M", 1, 0))

  

  
  set.seed(seed)
  
  # Start multiple imputation
  impvars <- c("STUDYID.x", "USUBJID.x", "SUBJID.x", "TRIAL", "TRTA", "TRTAN", 
               "VISIT", "AVISIT", "AVISITN", "AVAL", 
               "AGE", "MaleGender", "WhiteRace", "BASE", "HEIGHTBL", "WEIGHTBL", "ONSYRS", 
               "DIAGYRS", "PRDMTGR", "PRMSGR", "RLPS1YR", "RLPS3YR", "GDLESBL", "T1LESBL", "T2LESBL",
               "NHPTMBL", "PASATABL", "T25FWABL", "EDSSBL", "TRELMOS", "SFPCSBL", "SFMCSBL", "BVZBL", "VFT25BL")
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
  
  ds <- data %>% 
    group_by(USUBJID.x, TRIAL) %>% 
    filter(AVISITN == 0 & TRIAL %in% c(STUDYID_control, STUDYID_treat) & TRTA == "Placebo") %>%
    summarize(
      AGE = first(AGE),
      SEX = first(MaleGender),
      RACE = first(WhiteRace),
      WEIGHTBL = first(WEIGHTBL),
      HEIGHTBL = first(HEIGHTBL),
      ONSYRS = first(ONSYRS),
      DIAGYRS = first(DIAGYRS),
      RLPS1YR = first(RLPS1YR),
      RLPS3YR = first(RLPS3YR),
      EDSSBL = first(EDSSBL),
      GDLESBL = first(GDLESBL),
      T1LESBL = first(T1LESBL),
      T2LESBL = first(T2LESBL),
      PRMSGR = first(PRMSGR),
      NHPTMBL = first(NHPTMBL),
      PASATABL = first(PASATABL),
      T25FWABL = first(T25FWABL),
      TRELMOS = first(TRELMOS),
      SFMCSBL = first(SFMCSBL),
      SFPCSBL = first(SFPCSBL),
      BVZBL = first(BVZBL),
      VFT25BL = first(VFT25BL)
    )
  
  EDSS_END <-  subset(data, TRIAL %in% c(STUDYID_control, STUDYID_treat) & TRTA == "Placebo" & VISIT %in% END_VISIT)[,c("USUBJID.x", "AVAL")]
  
  ds_full <- (left_join(ds, EDSS_END, "USUBJID.x"))
  ds_full <- na.omit(ds_full)
  
  ds_full$Trt <- NA
  ds_full$Trt[which(ds_full$TRIAL %in% STUDYID_control)] <- 0
  ds_full$Trt[which(ds_full$TRIAL %in% STUDYID_treat)] <- 1
  ds_full$y <- ds_full$AVAL
  ds_full$TrtGroup <- ifelse(ds_full$Trt == 0, "Control", "Active")
  
  
  # DO we need to take a subgroup?
  if (SUBGROUP) {
    ds_full <- subset(ds_full, (EDSSBL > 2.50 & Trt == 1) | Trt == 0)
  }
  
  
  ds_full
}

simulate_nrs <- function(data, seed = 944) {
  set.seed(seed)
  
  # Remove grouping variable
  data <- ungroup(data)
  data <- data %>% dplyr::select(-c(USUBJID.x, TRIAL)) # Remove key columns
  data <- data %>% mutate(imputed = 0)
  
  # Set baseline characteristics of the trial
  mean_EDSS_trt0 <- mean(data$EDSSBL)
    
  # Add empty rows to data
  data <- data %>% add_row(
    EDSSBL = rpois(500, mean_EDSS_trt0),
    Trt = rep(0,500),
    imputed = 1)
  data <- data %>% add_row(
    EDSSBL = rpois(2000, mean_EDSS_trt0),
    Trt = rep(1,2000),
    imputed = 1
    )
  
  # Impute!
  impvars <- c("AGE", "SEX", "RACE", "WEIGHTBL", "ONSYRS", 
               "DIAGYRS", "PRMSGR", "RLPS1YR", "RLPS3YR", "GDLESBL", "T1LESBL", "T2LESBL",
               "NHPTMBL", "PASATABL", "T25FWABL", "EDSSBL", "TRELMOS", "SFPCSBL", "SFMCSBL", "BVZBL", "VFT25BL",
               "Trt", "y", "imputed")
  mice.prep <- mice(data[,impvars], maxit = 0)
  
  pM <- mice.prep$predictorMatrix
  pM[, "imputed"] <- 0
  
  imeth <- mice.prep$method
  imeth["AGE"] <- "pmm"
  imeth["WEIGHTBL"] <- "rf"
  imeth["NHPTMBL"] <- "pmm"
  imeth["T25FWABL"] <- "pmm"
  imeth["BVZBL"] <- "norm"
  imeth["ONSYRS"] <- "pmm"
  imeth["DIAGYRS"] <- "pmm"
  imeth["RLPS1YR"] <- "pmm" 
  imeth["RLPS3YR"] <- "pmm" 


  

  fit <- mice(data[,impvars], method = imeth, predictorMatrix = pM, m = 1, maxit = 10, printFlag = FALSE)
  
  dat_cc <- complete(fit,1)
  
  # Remove patients with baseline EDSS <= 2.50
  dat_cc <- subset(dat_cc, ((EDSSBL > 2.50 & Trt == 1) | Trt == 0) & imputed == 1)
  dat_cc <- dat_cc %>% dplyr::select(-imputed)
  
  dat_cc$TrtGroup <- ifelse(dat_cc$Trt == 0, "Control", "Active")

  return(dat_cc)
}

plot_density  <- function(data, facet = "TRIAL") {
  
  npat <-  nrow(data)
  
  ggdat <- data.frame(x = c(data$AGE, 
                            data$WEIGHTBL, 
                            data$ONSYRS,
                            data$DIAGYRS,
                            data$EDSSBL,
                            data$RLPS1YR,
                            data$RLPS3YR,
                            log(data$GDLESBL),
                            data$T1LESBL,
                            data$T2LESBL,
                            log(data$NHPTMBL),
                            data$PASATABL,
                            log(data$T25FWABL),
                            data$TRELMOS,
                            data$SFPCSBL,
                            data$SFMCSBL,
                            data$BVZBL,
                            data$VFT25BL),
                      var = c(rep("Age (years)", npat), 
                              rep("Weight (kg)", npat), 
                              rep("Time since first symptoms (years)", npat),
          rep("Time since diagnosis (years)", npat),
          rep("Baseline EDSS", npat),
          rep("No. of Relapses (12mo)", npat),
          rep("No. of Relapses (36mo)", npat),
          rep("Log No. of Gd+ lesions", npat),
          rep("No. of T1 lesions", npat),
          rep("No. of T2 lesions",npat),
          rep("Log 9 Hole Peg Test Average Score", npat),
          rep("PASAT 3", npat), 
          rep("Log Timed 25 Foot Walk", npat),
          rep("No. of months Since Recent Pre-Study Relapse",npat),
          rep("SF-36 PCS", npat),
          rep("SF-36 MCS", npat),
          rep("Brain Volume Z-Score", npat),
          rep("VFT 2.5%", npat)),
          Treatment = rep(data[[facet]], 18))

  
  ggplot(ggdat, aes(x = x, color = Treatment)) + 
    geom_density() + 
    facet_wrap( ~ var, scales = "free", ncol = 2) +
    xlab("") +
    theme(legend.position = "bottom") +
    ggtitle("Baseline characteristics")
}

analyze_nrs <- function(data) {
  
  results <- data.frame(method = c("naive", "DSR NNM", "DSR OFM", "DSR IPW", "DSR TDW" ), 
                        est_beta = NA, 
                        est_se = NA,
                        est_beta_CIl = NA, 
                        est_beta_CIu = NA,
                        est_time = NA)
  
  
  
  ####################################################################################################
  # Naive treatment effect
  ####################################################################################################
  start.naive <- proc.time()
  fit <- glm(y ~ Trt, data = data)
  t.naive <- (proc.time() - start.naive)[1]
  
  results$est_beta[1] <- coef(fit)["Trt"]
  results$est_se[1] <- sqrt(vcov(fit)["Trt", "Trt"])
  results$est_beta_CIl[1] <- confint(fit)["Trt",]["2.5 %"]
  results$est_beta_CIu[1] <- confint(fit)["Trt",]["97.5 %"]
  results$est_time[1] <- t.naive
  
  ####################################################################################################
  # Baseline adjusted treatment effect
  ####################################################################################################
  
  ### a) fit a DRS using the unexposed group only
  start.pgs <- proc.time()
  pgs <- glm(y ~ AGE + 
               SEX + 
               RACE + 
               WEIGHTBL +
               ONSYRS +
               DIAGYRS + 
               PRMSGR + 
               RLPS1YR + # No. of Relapses 1 Yr Prior to Study
               RLPS3YR +
               GDLESBL  +
               T1LESBL  +
               T2LESBL +
               NHPTMBL +
               PASATABL +
               T25FWABL +
               EDSSBL + # Baseline EDSS
               TRELMOS  +
               SFPCSBL  +
               SFMCSBL +
               BVZBL  +
               VFT25BL 
             , data = data[data$Trt == 0,])
  
  data$pgs <- predict(pgs, newdata = data, "response")  
  t.pgs <- (proc.time() - start.pgs)[1]
  
  #data$Treatment <- factor(data$Trt, levels = c(0,1), labels = c("Control", "Comparator"))
  #g1 <- ggplot(data, aes(x = pgs, group = Treatment, color = Treatment)) + geom_density()
  
  # NNM
  start.nnm <- proc.time()
  nnmatch <- matchit(Trt ~ pgs, data = data, caliper = 0.025, estimand = "ATT")
  data$nnm <- nnmatch$weights
  fit <- glm(y ~ Trt, weights = nnm, data = data)
  t.nnm <- (proc.time() - start.nnm)[1]
  
  results$est_beta[2] <- fit$coefficient["Trt"]
  results$est_se[2] <- sqrt(vcov(fit)["Trt", "Trt"])
  results$est_beta_CIl[2] <- confint(fit)["Trt",]["2.5 %"]
  results$est_beta_CIu[2] <- confint(fit)["Trt",]["97.5 %"]
  results$est_time[2] <- t.nnm
  
  
  # Full matching
  start.ofm <- proc.time()
  options("optmatch_max_problem_size" = Inf)
  fullmatch <- matchit(Trt ~ pgs, data = data, method = "full", estimand = "ATT")     #if it crashes due to large data set, re-run the above option code and re-run the full matching code
  data$ofm <- fullmatch$weights
  fit <- glm(y ~ Trt, weights = ofm, data = data)
  t.ofm <- (proc.time() - start.ofm)[1]

  results$est_beta[3] <- fit$coefficient["Trt"]
  results$est_se[3] <- sqrt(vcov(fit)["Trt", "Trt"])
  results$est_beta_CIl[3] <- confint(fit)["Trt",]["2.5 %"]
  results$est_beta_CIu[3] <- confint(fit)["Trt",]["97.5 %"]
  results$est_time[3] <- t.ofm
  
  
  # IPW
  start.ipw <- proc.time()
  ps <- glm(Trt ~ pgs, data = data, family = binomial())         #fit a prognostic propensity score
  data$ps <- ps$fitted
  data$ipw <- data$Trt*1 + (1 - data$Trt)*(data$ps/(1 - data$ps))
  fit <- glm(y ~ Trt, weights = ipw, data = data)
  t.ipw <- (proc.time() - start.ipw)[1]
  #g2 <- ggplot(data, aes(x = ps, group = Treatment, color = Treatment)) + geom_density() + xlab("Propensity score")
  
  results$est_beta[4] <- fit$coefficient["Trt"]
  results$est_se[4] <- sqrt(vcov(fit)["Trt", "Trt"])
  results$est_beta_CIl[4] <- confint(fit)["Trt",]["2.5 %"]
  results$est_beta_CIu[4] <- confint(fit)["Trt",]["97.5 %"]
  results$est_time[4] <- t.ipw
  
  # Target distribution weighting
  start.tdw <- proc.time()
  data$tdw <- tdw("Trt", "pgs", as.data.frame(data), estimand = "att")
  fit <- glm(y ~ Trt, weights = tdw, data = data)
  t.tdw <- (proc.time() - start.tdw)[1]

  results$est_beta[5] <- fit$coefficient["Trt"]
  results$est_se[5] <- sqrt(vcov(fit)["Trt", "Trt"])
  results$est_beta_CIl[5] <- confint(fit)["Trt",]["2.5 %"]
  results$est_beta_CIu[5] <- confint(fit)["Trt",]["97.5 %"]
  results$est_time[5] <- t.tdw
  
  return(results)
}