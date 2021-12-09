library("sas7bdat")
library("dplyr")
library("mice")
library(MatchIt)
library(tidyr)
library(mvtnorm)
library(splines)

# target distribution importance weighting
tdw <- function(a, x, data, estimand = "ate",
                method = "linear", rule = 1, f = 0, ties = mean)
{
  d <- density(data[,x])
  f <- approxfun(d, method = method, rule = rule, f = f, ties = ties)
  d0 <- density(data[data[,a] == 0, x])
  f0 <- approxfun(d0, method = method, rule = rule, f = f, ties = ties)
  d1 <- density(data[data[,a] == 1, x])
  f1 <- approxfun(d1, method = method, rule = rule, f = f, ties = ties)
  
  if (estimand == "ate") {
    data[data[,a] == 1, "w"] <- f(data[data[,a] == 1, x])/f1(data[data[,a] == 1, x])
    data[data[,a] == 0, "w"] <- f(data[data[,a] == 0, x])/f0(data[data[,a] == 0, x])
  } else if (estimand == "att") {
    data[data[,a] == 1, "w"] <- 1
    data[data[,a] == 0, "w"] <- f1(data[data[,a] == 0, x])/f0(data[data[,a] == 0, x])
  } else if (estimand == "atu") {
    data[data[,a] == 1, "w"] <- f0(data[data[,a] == 1, x])/f1(data[data[,a] == 1, x])
    data[data[,a] == 0, "w"] <- 1
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
  
  ds_full
}

simulate_nrs <- function(data, load = TRUE, seed = 944, dir = "../Data/") {
  set.seed(seed)
  
  dat_cc <- NULL
  if (load & file.exists(paste(dir, "simdat.rda", sep = "" ))) {
    load(paste(dir, "simdat.rda", sep = "" ))
    return(dat_cc)
  }
  
  # Outcome model
  data$yinteger <- data$y * 2
  fit_om <- glm(yinteger ~ 
                  AGE + 
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
                  VFT25BL , data = data, family = quasipoisson())
  
  # Remove grouping variable
  data <- ungroup(data)
  data <- data %>% dplyr::select(-c(USUBJID.x, TRIAL)) # Remove key columns
  data <- data %>% mutate(imputed = 0)
  
  ncontrol <- 500
  ntreated <- 2000
  
  data$logT25FWABL <- log(data$T25FWABL + 1)
  data$logNHPTMBL <- log(data$NHPTMBL + 1)
  data$sqrtWEIGHTBL <- sqrt(data$WEIGHTBL)
  
  # SIgnificant variables:  GDLESBL, T1LESBL, RLPS3YR 
  
  simvars <- c("SFPCSBL", "SFMCSBL", "AGE", "logT25FWABL", "logNHPTMBL", "sqrtWEIGHTBL", "BVZBL")
  
  meanControl <- colMeans(subset(data, Trt == 0)[,simvars]) - c(3.5, 2.5, -2.5, 0.25, 0.20, -1, 0.3)
  vcovControl <- cov(subset(data, Trt == 0)[,simvars])
  rControl <- data.frame(rmvnorm(ncontrol, mean = meanControl, sigma = vcovControl))
  
  meanTreated <- colMeans(subset(data, Trt == 1)[,simvars]) + c(2.5, 3, -2, 0.25, 0.1, 1, 0.3)
  vcovTreated <- cov(subset(data, Trt == 1)[,simvars])
  rTreated <- data.frame(rmvnorm(ntreated, mean = meanTreated, sigma = vcovTreated))

  # Add empty rows to data
  data <- data %>% add_row(
    #SFPCSBL = rControl$SFPCSBL,
    #SFMCSBL = rControl$SFMCSBL,
    #AGE = rControl$AGE,
    #T25FWABL = exp(rControl$logT25FWABL) - 1,
    #NHPTMBL = exp(rControl$logNHPTMBL) - 1,
    #WEIGHTBL = rControl$sqrtWEIGHTBL**2,
    #BVZBL = rControl$BVZBL,
    #PASATABL = rbeta(ncontrol, 5.1, 1.15)*60,
    #GDLESBL = rbeta(ncontrol, 1.5, 130)*75,
    Trt = rep(0, ncontrol),
    imputed = 1)
  data <- data %>% add_row(
    #SFPCSBL = rTreated$SFPCSBL,
    #SFMCSBL = rTreated$SFMCSBL,
    #AGE  = rTreated$AGE,
    #T25FWABL = exp(rTreated$logT25FWABL) - 1,
    #NHPTMBL = exp(rTreated$logNHPTMBL) - 1,
    #WEIGHTBL = rTreated$sqrtWEIGHTBL**2,
    #BVZBL = rTreated$BVZBL,
    #PASATABL = rbeta(ntreated, 4.8, 1.00)*60,
    #GDLESBL = rbeta(ntreated, 1.8, 250)*200,
    Trt = rep(1, ntreated),
    imputed = 1
    )
  
  # Impute!
  impvars <- c("AGE", "SEX", "RACE", "WEIGHTBL", "ONSYRS", 
               "DIAGYRS", "PRMSGR", "RLPS1YR", "RLPS3YR", "GDLESBL", "T1LESBL", "T2LESBL",
               "NHPTMBL", "PASATABL", "T25FWABL", "EDSSBL", "TRELMOS", "SFPCSBL", "SFMCSBL", "BVZBL", "VFT25BL",
               "Trt", "imputed", "y")
  mice.prep <- mice(data[,impvars], maxit = 0)
  
  pM <- mice.prep$predictorMatrix
  pM[, "imputed"] <- 0
  pM["y", "Trt"] <- 0 #treatment does not affect the outcome
  
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
  imeth["TRELMOS"] <- "rf" 
  imeth["SFMCSBL"] <- "rf"
  imeth["SFMCSBL"] <- "rf"

  fit <- mice(data[,impvars], method = imeth, predictorMatrix = pM, m = 1, maxit = 10, printFlag = FALSE)
  
  dat_cc <- subset(complete(fit,1), imputed == 1)
  

  #lp_y <- stats::predict(fit_om, newdata = dat_cc, type = "response")
  #dat_cc$y <- (rpois(nrow(dat_cc), lambda = lp_y))/2
  
  # Selectively remove patients from the control group
  lp_include <- -log(dat_cc$T25FWABL + 1) - 0.09 * dat_cc$SFPCSBL - 0.05 * dat_cc$SFMCSBL - 0.7 * dat_cc$GDLESBL - 0.5 * log(dat_cc$T1LESBL + 1) - 0.05 * dat_cc$AGE - 0.8 * dat_cc$EDSSBL - 1 * dat_cc$RLPS3YR  - 1 * dat_cc$SEX
  
  p_include <- 1/(1 + exp(-(-mean(lp_include) + lp_include)))

  dat_cc$include <- rbinom(n = nrow(dat_cc), size = 1, prob = p_include)
  dat_cc$include[dat_cc$Trt == 0] <- 1 # Include all control patients
  
  
  dat_cc <- subset(dat_cc, include == 1)

  
  # Regenerate outcome y
  dat_cc <- dat_cc %>% dplyr::select(-imputed)
  dat_cc$TrtGroup <- ifelse(dat_cc$Trt == 0, "Control", "Active")

  save(dat_cc, file = paste(dir, "simdat.rda", sep = "" ))

  return(dat_cc)
}

plot_density  <- function(data, facet = "TRIAL", palette = "Set1") {
  
  npat <-  nrow(data)
  
  ggdat <- data.frame(x = c(data$AGE, 
                            data$WEIGHTBL, 
                            data$ONSYRS,
                            data$DIAGYRS,
                            data$EDSSBL,
                            data$RLPS1YR,
                            data$RLPS3YR,
                            log(data$GDLESBL + 1),
                            data$T1LESBL,
                            data$T2LESBL,
                            log(data$NHPTMBL + 1),
                            data$PASATABL,
                            log(data$T25FWABL + 1),
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
    ggtitle("Baseline characteristics") +
    scale_color_brewer(palette = palette) + 
    scale_fill_brewer(palette = palette)
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
  
  data$pgs <- predict(pgs, newdata = data, type = "response")  
  t.pgs <- (proc.time() - start.pgs)[1]
  

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

analyze_nrs_bs <- function(data, iter = 1000, seed = 944, load = TRUE, dir = "../Data/") {
  B.RESULTS <- NULL
  
  if (load & file.exists(paste(dir, "bootstrap.rda", sep = ""))) {
    load(paste(dir, "bootstrap.rda", sep = ""))
    return(B.RESULTS)
  }
  
  
  pb <- txtProgressBar(min = 0, max = iter, initial = 0) 
  for (b in seq(iter))
  {
    set.seed(seed + b)
    
    data.b <- data[sample(rownames(data), nrow(data), replace = TRUE),]
    b.results <- analyze_nrs(data.b)
    
    B.RESULTS <- rbind(B.RESULTS, data.frame(iter = b, b.results))
    setTxtProgressBar(pb,b)
  }
  close(pb)
  
  B.RESULTS$method <- factor(B.RESULTS$method, labels = unique(B.RESULTS$method),
                             levels = unique(B.RESULTS$method))
  
  save(B.RESULTS, file = paste(dir, "bootstrap.rda", sep = ""))
  
  return(B.RESULTS)
}


plot_distr_ps <- function(data, palette = "Set1") 
{
  # Calculate the PS
  fit <- glm(Trt ~ AGE + 
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
        VFT25BL, data = data, family = binomial())
  
  data$pps <- stats::predict(fit, type = "response")
  
  data <- data %>%
    spread(TrtGroup, "pps", sep = "_p")
  
  ggplot(data) + 
    geom_histogram(breaks = seq(0, 1, length = 50), aes(x = TrtGroup_pActive, fill = "Patients treated with active placebo"), color = "red", alpha = 0.5) + 
    geom_histogram(breaks = seq(0, 1, length = 50), aes(x = TrtGroup_pControl, y = -..count.., fill = "Patients treated with control placebo"), color = "blue", alpha = 0.5) + 
    ylab("Number of patients") + 
    xlab("Probability of receiving active placebo") +
    geom_hline(yintercept = 0, lwd = 0.5) +
    scale_y_continuous(labels = abs) +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    xlim(c(0,1)) +
    scale_color_brewer(palette = palette) + 
    scale_fill_brewer(palette = palette)
}