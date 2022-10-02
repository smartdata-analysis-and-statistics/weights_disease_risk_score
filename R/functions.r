library(sas7bdat)
library(broom)
library(dplyr)
library(mice)
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
  } else {
    stop("Estimand not supported!")
  }
  TDW <- data[,"w"]
  return(TDW)
}


plot_density <- function(data, facet = "TRIAL", palette = "Set1") {
  
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

estimate_pgs <- function(data) {
  fit <- glm(y ~ AGE + 
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
  
  data$pgs <- predict(fit, newdata = data, type = "response")  
  return(data)
}

analyze_nrs <- function(data, seed = 1234) {
  set.seed(seed)
  
  results <- data.frame("method" = character(),
                        "est_beta" = numeric(),
                        "est_se" = numeric(),
                        "est_beta_CIl" = numeric(),
                        "est_beta_CIu" = numeric(),
                        "est_time" = numeric(),
                        "seed" = numeric())
  
  ##############################################################################
  # Naive treatment effect
  ##############################################################################
  start.naive <- proc.time()
  fit <- lm(y ~ Trt,data = data)
  fitci <- tidy(fit, conf.int = TRUE)
  t.naive <- (proc.time() - start.naive)[1]
  
  results <- results %>% add_row(
    data.frame("method"  = "naive",
               "est_beta" = (fitci %>% filter(term == "Trt"))$estimate,
               "est_se" = (fitci %>% filter(term == "Trt"))$std.error,
               "est_beta_CIl" = (fitci %>% filter(term == "Trt"))$conf.low,
               "est_beta_CIu" = (fitci %>% filter(term == "Trt"))$conf.high,
               "est_time" = t.naive))
  rm(fit, fitci)
  
  ##############################################################################
  # NNM
  ##############################################################################
  start.nnm <- proc.time()
  nnmatch <- matchit(Trt ~ pgs, data = data, caliper = 0.025)
  data$nnm <- nnmatch$weights
  fit <- glm(y ~ Trt,weights = nnm,data = data)
  fitci <- tidy(fit, conf.int = TRUE)
  t.nnm <- (proc.time() - start.nnm)[1]
  
  results <- results %>% add_row(data.frame("method"  = "DSR NNM",
                                            "est_beta" = (fitci %>% filter(term == "Trt"))$estimate,
                                            "est_se" = (fitci %>% filter(term == "Trt"))$std.error,
                                            "est_beta_CIl" = (fitci %>% filter(term == "Trt"))$conf.low,
                                            "est_beta_CIu" = (fitci %>% filter(term == "Trt"))$conf.high,
                                            "est_time" = t.nnm))
  rm(fit, fitci)
  
  ##############################################################################
  # OFM
  ############################################################################## 
  start.ofm <- proc.time()
  options("optmatch_max_problem_size" = Inf)
  fullmatch <- matchit(Trt ~ pgs, data = data, method = "full")
  data$ofm <- fullmatch$weights
  fit <- glm(y ~ Trt,weights = ofm, data = data)
  fitci <- tidy(fit, conf.int = TRUE)
  t.ofm <- (proc.time() - start.ofm)[1]
  
  results <- results %>% add_row(data.frame("method"  = "DSR OFM",
                                            "est_beta" = (fitci %>% filter(term == "Trt"))$estimate,
                                            "est_se" = (fitci %>% filter(term == "Trt"))$std.error,
                                            "est_beta_CIl" = (fitci %>% filter(term == "Trt"))$conf.low,
                                            "est_beta_CIu" = (fitci %>% filter(term == "Trt"))$conf.high,
                                            "est_time" = t.ofm))
  rm(fit, fitci)
  
  ##############################################################################
  # IPW
  ############################################################################## 
  start.ipw <- proc.time()
  ps <- glm(Trt~pgs,data = data, family = "binomial")
  data$ps <- ps$fitted
  data$ipw <- data$Trt*1 + (1 - data$Trt)*(data$ps/(1 - data$ps))
  fit <- glm(y ~ Trt, weights = ipw, data = data)
  fitci <- tidy(fit, conf.int = TRUE)
  t.ipw <- (proc.time() - start.ipw)[1]
  
  results <- results %>% add_row(data.frame("method"  = "DSR IPW",
                                            "est_beta" = (fitci %>% filter(term == "Trt"))$estimate,
                                            "est_se" = (fitci %>% filter(term == "Trt"))$std.error,
                                            "est_beta_CIl" = (fitci %>% filter(term == "Trt"))$conf.low,
                                            "est_beta_CIu" = (fitci %>% filter(term == "Trt"))$conf.high,
                                            "est_time" = t.ipw))
  rm(fit, fitci)
  
  ##############################################################################
  # TDW
  ##############################################################################  
  start.tdw <- proc.time()
  data$tdw <- tdw("Trt","pgs", data, estimand = "att")
  fit <- glm(y ~ Trt, weights = tdw, data = data)
  fitci <- tidy(fit, conf.int = TRUE)
  t.tdw <- (proc.time() - start.tdw)[1]
  
  results <- results %>% add_row(data.frame("method"  = "DSR TDW",
                                            "est_beta" = (fitci %>% filter(term == "Trt"))$estimate,
                                            "est_se" = (fitci %>% filter(term == "Trt"))$std.error,
                                            "est_beta_CIl" = (fitci %>% filter(term == "Trt"))$conf.low,
                                            "est_beta_CIu" = (fitci %>% filter(term == "Trt"))$conf.high,
                                            "est_time" = t.tdw))
  
  results$seed <- seed
  
  return(results)
}

analyze_nrs_bs <- function(data, iter = 1000, seed = 944, load = TRUE, dir = "../Data", fn = "bootstrap.rda") {
  
  B.RESULTS <- data.frame(iter = numeric(),
                          method = character(),
                          est_beta = numeric(),
                          est_se = numeric(),
                          est_beta_CIl = numeric(),
                          est_beta_CIu = numeric(),
                          est_time = numeric(),
                          seed = numeric()
                          )
  
  if (load & file.exists(file.path(dir, fn))) {
    message("Loading saved results")
    load(file.path(dir, fn))
    return(B.RESULTS)
  }
  
  if (!file.exists(dir)) {
    message("Creating new directory to save results")
    dir.create(dir)
  }
  
  pb <- txtProgressBar(min = 0, max = iter, initial = 0) 
  for (b in seq(iter))
  {
    data.b <- data[sample(rownames(data), nrow(data), replace = TRUE),]
    b.results <- analyze_nrs(data.b, seed = seed + b)
    
    B.RESULTS <- B.RESULTS %>% add_row(data.frame(iter = b, b.results))
    setTxtProgressBar(pb,b)
  }
  close(pb)
  
  B.RESULTS$method <- factor(B.RESULTS$method, labels = unique(B.RESULTS$method),
                             levels = unique(B.RESULTS$method))
  
  save(B.RESULTS, file = file.path(dir, fn))
  
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
    ggtitle("Distribution of the propensity score") +
    scale_color_brewer(palette = palette) + 
    scale_fill_brewer(palette = palette)
}

plot_distr_pgs <- function(data, palette = "Set1") {
  data <- data %>%
    spread(TrtGroup, "pgs", sep = "_p")
  
  ggplot(data) + 
    geom_histogram(breaks = seq(0, 7, 0.5), aes(x = TrtGroup_pActive, fill = "Patients treated with active placebo"), color = "red", alpha = 0.5) + 
    geom_histogram(breaks = seq(0, 7, 0.5), aes(x = TrtGroup_pControl, y = -..count.., fill = "Patients treated with control placebo"), color = "blue", alpha = 0.5) + 
    ylab("Number of patients") + 
    xlab("Predicted EDSS after 36 weeks if treated with control placebo") +
    geom_hline(yintercept = 0, lwd = 0.5) +
    scale_y_continuous(labels = abs) +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    ggtitle("Distribution of the disease risk score") +
    scale_color_brewer(palette = palette) + 
    scale_fill_brewer(palette = palette)
}

plot_balance <- function(data, palette = "Set1") # show results for matching with replacement?
{
  
  nnmatch <- matchit(Trt ~ AGE + 
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
                       VFT25BL, data = data, caliper = 0.025, estimand = "ATT")
  
    # Calculate SMD
    sfit <- summary(nnmatch)
    
    # Unmatched data
    ggdat <- data.frame("SMD" = abs(sfit$sum.all[, "Std. Mean Diff."]),
                        "covariate" = rownames(sfit$sum.all))
    ggdat <- ggdat %>% arrange(desc(SMD))
    ggdat$covariate <- factor(ggdat$covariate, levels = ggdat$covariate, labels = ggdat$covariate)
    
    
    ggplot(data = ggdat,
           mapping = aes(x = covariate, y = SMD)) +
      #geom_line() +
      xlab("") + 
      ylab("Absolute standardized mean difference") +
      geom_point() +
      geom_hline(yintercept = 0.1, color = "black", size = 0.1, lty = 2) +
      coord_flip() +
      theme_bw() + 
      theme(legend.key = element_blank())  + 
      theme(legend.position = "bottom") + 
      theme(legend.title = element_blank())
  
}