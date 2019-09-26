extract <-function(out)
{
  psi10 <-out["psi"][[1]][[1]][1] 
  # psi11 <-out["psi"][[1]][[1]][2]
  return(psi10)
}

dbalpha <- function(data, psin, B1 = 500, B2 = 500, blip.model, treat.model, tf.model)
{
  data <- as.data.frame(data)
  
  # initialize
  coverage <- 0
  probaN <- list(as.vector(rep(0.5,n)))
  n <- nrow(data)
  alpha <- 0.025

  while(alpha <= 0.15)
  {
    # reset matrix to save estimates
    est <- matrix(NA, ncol = B2 + 3, nrow = B1)
    for(j in 1:B1) # loop over B1 first stage bootstrap samples
    {
      # draw a n-out-of-n bootstrap sample
      index <- sample(1:n, n, replace = TRUE)
      boot1 <- as.data.frame(data[index,])
      
      # fit the model to b1-th bootstrap sample
      res1 <- try(DTRreg(outcome = Y2, blip.mod = blip.model, treat.mod = treat.model, tf.mod = tf.model, treat.mod.man = rep(probaN,2), method = "dwols", var.estim = "bootstrap", data = boot1))
      esb1 <- try(extract(res1))
      est[j,1] <- esb1
      
      # estimate m for each b1 bootstrap sample
      phat <- res1$nonreg[2]
      est[j,2] <- phat
      
      # resampling size
      m <- n^((1 + alpha*(1-phat))/(1 + alpha)) 
      est[j,3] <- m
      
      # probability treatment with m
      proba <- list(as.vector(rep(0.5,floor(m))))
      
      for(k in 1:B2) # loop over B2 second stage bootstrap samples
      {
        # resample with replacement 
        index <- sample(1:n, floor(m), replace = TRUE)
        boot2 <- boot1[index,]
        
        # fit the model to bootstrap sample
        res2 <- try(DTRreg(outcome = Y2, blip.mod = blip.model, treat.mod = treat.model, tf.mod = tf.model, treat.mod.man = rep(proba,2), method = "dwols", data = as.data.frame(boot2)))
        esb2 <- try(extract(res2))
        
        # save the (b1,b2) bootstrap estimates j in the (k+3) column
        est[j,k + 3] <- esb2
      }
      quan <- quantile(sqrt(m)*(est[j,4:(503)]-esb1), probs = c(0.025,0.975))
      coverage <- coverage + ((esb1 - quan[2]/sqrt(m) <= psin) & (esb1 - quan[1]/sqrt(m) >= psin))
    }
    # for each B1 first stage bootstrap sample, calculate 0.025, 0.975 percentiles
    coverage <- coverage/B1
    if(coverage >= 0.95)
    {
      print(alpha)
      break
    }
    else
    {
      print(alpha)
      alpha <- alpha + 0.025
    }
  }
  return(as.numeric(alpha))
}