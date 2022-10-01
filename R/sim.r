################################################################################
# *** PGS weigthing methods -  T.L. Nguyen 2021 *** #
################################################################################
require(MatchIt)

# From Wyss et al. Pharmacoepi and Drug Safety 2015: 10.1002/pds.3810
Alpha <- function(d,bounds) {runif(d,bounds[1],bounds[2])}
Beta <- function(d,bounds) {runif(d,bounds[1],bounds[2])}

simData <- function(n, alpha, beta, hte = "no")
{
  x.1 <- replicate(96,rbinom(n,1,prob = 0.5))
  x.97 <- replicate(4,rnorm(n,mean = 0,sd = 1))
  x <- cbind(x.1,x.97)
  
  a <- rbinom(n,1,prob = plogis(-0.847 + rowSums(t(alpha*t(x)))))
  int <- ifelse(hte == "no", 0, 0.7)

  
  y.1 <- rbinom(n,1,prob = plogis(-0.847 + x[,1]*int + rowSums(t(beta*t(x)))))
  y.0 <- rbinom(n,1,prob = plogis(-0.847 + rowSums(t(beta*t(x)))))
  y <- a*y.1 + (1 - a)*y.0
  
  data.frame(x,a,y.1,y.0,y)
}

simDataHist <- function(n,beta)
{
  x.1 <- replicate(96,rbinom(n,1,prob = 0.5))
  x.97 <- replicate(4,rnorm(n,mean = 0,sd = 1))
  x <- cbind(x.1,x.97)
  y <- rbinom(n,1,prob = plogis(-0.847 + rowSums(t(beta*t(x)))))
  data.frame(x,y)
}

runSim <- function(n = 1000, nh = 10000, iter = 1000, seed = 93828, load = TRUE, 
                   dir = "../Data", fn = "simulation_results.rda") {
  
  results <- NULL
  
  if (load & file.exists(file.path(dir, fn))) {
    message("Loading saved results")
    load(file.path(dir, fn))
    return(results)
  }
  
  scenarios <- data.frame("A" = c(-0.182,0.182), 
                          "B" = c(-0.405,0.405), 
                          "C" = c(-0.7,0.7) )
  TIME <- NULL
  ESTIMATE <- NULL
  
  
  for (s in colnames(scenarios))
  {
    set.seed(seed)
    alpha <- Alpha(100,scenarios[,s])
    set.seed(709530)
    beta <- Beta(100,scenarios[,s])
    
    pb <- txtProgressBar(min = 0, max = iter, initial = 0) 
    for (i in seq(iter))
    {
      setTxtProgressBar(pb,i)
      for (hte in c("no","yes"))
      {
        set.seed(4309 + i)
        data <- simData(n,alpha,beta,hte)
        data.h <- simDataHist(nh,beta)
        
        att <- mean(data[data$a == 1,]$y.1) - mean(data[data$a == 1,]$y.0)
        naive <- coef(lm(y~a,data = data))["a"]
        
        fmla.pgs <- as.formula(paste("y~",paste0(colnames(data[,1:100]), collapse = "+")))
        pgs <- glm(fmla.pgs,data = data.h,family = "binomial")
        data$pgs <- predict(pgs,newdata = data, "response")
        
        # NNM
        start.nnm <- proc.time()
        nnmatch <- matchit(a~pgs, data = data, caliper = 0.025)
        data$nnm <- nnmatch$weights
        NNM <- coef(glm(y~a,weights = nnm,data = data))["a"]
        t.nnm <- (proc.time() - start.nnm)[1]
        
        # Full matching
        start.ofm <- proc.time()
        fullmatch <- matchit(a~pgs, data = data, method = "full")
        data$ofm <- fullmatch$weights
        OFM <- coef(glm(y~a,weights = ofm,data = data))["a"]
        t.ofm <- (proc.time() - start.ofm)[1]
        
        # IPW
        start.ipw <- proc.time()
        ps <- glm(a~pgs,data = data,family="binomial")
        data$ps <- ps$fitted
        data$ipw <- data$a*1 + (1-data$a)*(data$ps/(1-data$ps))
        IPW <- coef(glm(y~a,weights=ipw,data=data))["a"]
        t.ipw <- (proc.time() - start.ipw)[1]
        
        # Target distribution weighting
        start.tdw <- proc.time()
        data$tdw <- tdw("a","pgs",data,estimand="att")
        TDW <- coef(glm(y~a,weights=tdw,data=data))["a"]
        t.tdw <- (proc.time() - start.tdw)[1]
        
        time <- data.frame(NNM = t.nnm, OFM = t.ofm, 
                           IPW = t.ipw, TDW = t.tdw)
        TIME <- rbind(TIME,time)
        
        estimate <- data.frame( scenario=s, hte, i, att, naive, NNM, OFM, IPW, TDW,
                                prev.a = mean(data$a), prev.y = mean(data$y)
        )
        ESTIMATE <- rbind(ESTIMATE,estimate)
      }
    }
    close(pb)
  }
  
  colnames(TIME) <- c("t.nnm", "t.ofm", "t.ipw", "t.tdw")
  results <- cbind(ESTIMATE,TIME)
  
  save(results, file = file.path(dir, fn))
  return(results)
}





