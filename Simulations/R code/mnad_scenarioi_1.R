install.packages("DTRreg")
library(DTRreg)
source("choose_alpha.R")

expit <- function(x) exp(x)/(1+exp(x))
set.seed(543) # each scenario used a different seed. See README file.

# gamma parameters following Chakraborty et al (2013) to control for irregularity in the generated data
g <- matrix(NA, nrow = 9, ncol = 7)
g[1,] <- c(0,0,0,0,0,0,0)
g[2,] <- c(0,0,0,0,0.01,0,0)
g[3,] <- c(0,0,-0.5,0,0.5,0,-0.5)
g[4,] <- c(0,0,-0.5,0,0.99,0,-0.98)
g[5,] <- c(0,0,-0.5,0,1,0.5,-0.5)
g[6,] <- c(0,0,-0.5,0,0.25,0.5,0.5)
g[7,] <- c(0,0,-0.25,0,0.75,0.5,0.5)
g[8,] <- c(0,0,0,0,1,0,-1)
g[9,] <- c(0,0,0,0,0.25,0,-0.24)

# delta parameters following Chakraborty et al (2013) to control for irregularity in the generated data
d <- matrix(NA, nrow = 9, ncol = 2)
d[1,] <- c(0.5,0.5)
d[2,] <- c(0.5,0.5)
d[3,] <- c(0.5,0.5)
d[4,] <- c(0.5,0.5)
d[5,] <- c(1,0)
d[6,] <- c(0.1,0.1)
d[7,] <- c(0.1,0.1)
d[8,] <- c(0,0)
d[9,] <- c(0,0)

######################### m-out-of-n bootstrap : adpative alpha #############################


sc <- seq(1,9)
# number of simulated dataset
Nsimul <- 500 
# number of boostrap samples
Nboot <- 1000
# sample size
n <- 300

# model specification
blip.model <- list(~ O1, ~ O2 + A1)
treat.model <- list(A1~1, A2~1) 
tf.model <- list(~ O1, ~ O1 + A1 + O1*A1)

# allocate space and predefined quantities
estm <- vector(mode = "list", length = 1)

# specify scenario with i between 1 and 9. Here scenario 1
i <- 1
  # reset estimates to NA for new scenario
  estm[[1]] <- matrix(NA, nrow = Nsimul , ncol = Nboot + 4)
  for(s in 1:Nsimul) # loop over number of simulations
  {
    # treatment A1, A2: P(Aj = 1) = P(Aj = 0) = 0.5
    A1 <- rbinom(n, size = 1, prob = 0.5)
    A2 <- rbinom(n, size = 1, prob = 0.5)
  
    # covariates O1, O2: coded as -1, 1, where O2 depends on A1, O1 and (delta_1,delta_2)
    O1 <- 2*rbinom(n, size = 1, prob = 0.5) - 1
    O2 <- 2*rbinom(n, size = 1, prob = expit(d[sc[i],1]*O1 + d[sc[i],2]*(2*A1-1))) - 1
    
    # generated outcome Y2 (Y1 set to 0), using parameters (gamma_1,...,gamma_7)
    Y2 <- g[sc[i],1] + g[sc[i],2]*O1 + g[sc[i],3]*A1 + g[sc[i],4]*O1*A1 + g[sc[i],5]*A2 + g[sc[i],6]*O2*A2 + g[sc[i],7]*A1*A2 + rnorm(n)
    
    # generated dataset
    complete <- cbind(A1, A2, O1, O2, Y2)
    
    # fit dWOLS to the generated dataset, using all n=300 observations
    proba <- list(as.vector(rep(0.5,n)))
    res.n <- try(DTRreg(outcome = Y2, blip.mod = blip.model, treat.mod = treat.model, tf.mod = tf.model, treat.mod.man = rep(proba,2), method = "dwols", var.estim = "bootstrap", data = as.data.frame(complete)))
    es <- try(extract(res.n))
    
    # save estimates using all observations in the first column
    estm[[1]][s,1] <- es
  
    # estimate of nonregularity
    phat <- res.n$nonreg[2]
    estm[[1]][s,2] <- phat
    
    # choice of alpha
    alpha <- dbalpha(data = complete, psin = es, blip.model = blip.model, treat.model = treat.model, tf.model = tf.model)
    estm[[1]][s,3] <- alpha
    print(c(s,alpha))
    
    # resampling size
    m <- n^((1 + alpha*(1-phat))/(1 + alpha))
    estm[[1]][s,4] <- m

    # probability treatment with m
    proba <- list(as.vector(rep(0.5,floor(m))))
    
    # bootstrap resampling + estimate
    for(b in 1:Nboot) # loop over number of bootstrap samples
    {
      # resample with replacement 
      index <- sample(1:n, floor(m), replace = TRUE)
      boot <- complete[index,]
      
      # fit the model to bootstrap sample
      res <- try(DTRreg(outcome = Y2, blip.mod = blip.model, treat.mod = treat.model, tf.mod = tf.model, treat.mod.man = rep(proba,2), method = "dwols", data = as.data.frame(boot)))
      esb <- try(extract(res))
      
      # save bootstrap estimates i in the (i+1) column
      estm[[1]][s, b + 4] <- esb
    }
  }
  # linux command to save results of the simulations in a CSV file
  colnames(estm[[1]])[1:4] <- c("psi10","phat","alpha","m")
  name1 <- paste("mnad_psi10_scenario", paste(sc[i]),"_1.csv",sep ="")
  write.csv(estm[[1]], file = name1, row.names = FALSE)










