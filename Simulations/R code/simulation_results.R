# reproduce result tables in Chakraborty et al. (2013) with dWOLS as method of analysis

# Table 2: average bootstrap resample size
# Table 3: coverage rate for 1000 simulated dataset
# Table 4: average width for 1000 simulated dataset

expit <- function(x) exp(x)/(1+exp(x))

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

# true paramater psi_10
true <- c(0, 0, -1, -1.48, -1, -0.5 +0.25 + (1/2*expit(0.2) + 1/2*expit(0)) - (1/2*expit(0) + 1/2*expit(-0.2))*0.75, -0.25 + 0.75 - 0.25 + (1/2*expit(0.2) + 1/2*expit(0)) - (1/2*expit(0) + 1/2*expit(-0.2)), -1, -0.24)

# percentile_nn for regular bootstrap
percentile_nn <- function(x)
{
  phi1n <- x[1] 
  dis <- x[2:1001]-phi1n
  quan <- sqrt(300) * quantile(dis, probs = c(0.025, 0.975))
  return(quan)
}

# construct CI for n-out-of-n regular bootstrap
consCI_nn <- function(x)
{
  ph1n <- x[1]
  l <- x[2]
  u <- x[3]
  CI <- c(ph1n - u/sqrt(300), ph1n - l/sqrt(300))
  return(CI)
}

# percentile_mn for regular bootstrap
percentile_mn <- function(x)
{
  phi1n <- x[1] 
  m <- x[3]
  dis <- x[4:1003]-phi1n
  quan <- sqrt(m) * quantile(dis, probs = c(0.025, 0.975))
  return(quan)
}

# construct CI for n-out-of-n regular bootstrap
consCI_mn <- function(x)
{
  ph1n <- x[1]
  m <- x[2]
  l <- x[3]
  u <- x[4]
  CI <- c(ph1n - u/sqrt(m), ph1n - l/sqrt(m))
  return(CI)
}

# scenarios
sc <- seq(1,9)

######################### for psi10 #############################

# Table 3 and Table 4 summarize the information in Figure 1
table2 <- matrix(NA, nrow = 4, ncol = 9) # average m
table3 <- matrix(NA, nrow = 4, ncol = 9) # coverage rate
table4 <- matrix(NA, nrow = 4, ncol = 9) # mean width
rownames(table2) <- rownames(table3) <- rownames(table4) <- c("nn","mn0.05","mn0.1","mnAD")
colnames(table2) <- colnames(table3) <- colnames(table4) <- paste("sc",paste(sc),sep = "")

# m-out-of-n bootstrap alpha=0.1
di <- c("mn0.1_psi10_scenario")

for(i in 1:9)
{
  name <- paste(di, paste(sc[i]), paste(".csv"), sep = "")
  dat <- read.csv(name, header = TRUE) 
  table2[3,i] <- mean(dat$m) # average m
  ul <- t(apply(dat, 1, percentile_mn))
  dat2 <- cbind(dat$psi10, dat$m, ul[,1], ul[,2])
  CI <- t(apply(dat2, 1, consCI_mn))
  cov <- ifelse(CI[,1] < true[i] & CI[,2] > true[i], 1, 0)
  print(prop.test(sum(cov), 1000, p=0.95)$p.value) # test H0: prop=0.5 vs HA: prop!=0.95
  table3[3,i] <- mean(cov)
  table4[3,i] <- mean((CI[,2]-CI[,1]))
}

# n-out-of-n bootstrap
di <- c("nn_psi10_scenario")

for(i in 1:9)
{
  name <- paste(di, paste(sc[i]), paste(".csv"), sep = "")
  dat <- read.csv(name, header = TRUE)
  table1[1,i] <- 300 # resampling size always 300
  ul <- t(apply(dat, 1, percentile_nn))
  dat2 <- cbind(dat$psi10, ul[,1], ul[,2])
  CI <- t(apply(dat2, 1, consCI_nn))
  cov <- ifelse(CI[,1] < true[i] & CI[,2] > true[i], 1, 0)
  print(prop.test(sum(cov), 1000, p=0.95)$p.value)
  table2[1,i] <- mean(cov)
  table3[1,i] <- mean((CI[,2]-CI[,1]))
}

# m-out-of-n bootstrap alpha=0.05
di <- c("mn0.05_psi10_scenario")

for(i in 1:9)
{
  name <- paste(di, paste(sc[i]), paste(".csv"), sep = "")
  dat <- read.csv(name, header = TRUE)
  table2[2,i] <- mean(dat$m) # average m
  ul <- t(apply(dat, 1, percentile_mn))
  dat2 <- cbind(dat$psi10, dat$m, ul[,1], ul[,2])
  CI <- t(apply(dat2, 1, consCI_mn))
  cov <- ifelse(CI[,1] < true[i] & CI[,2] > true[i], 1, 0)
  print(prop.test(sum(cov), 1000, p=0.95)$p.value)
  table3[2,i] <- mean(cov)
  table4[2,i] <- mean((CI[,2]-CI[,1]))
}


# m-out-of-n bootstrap adaptive alpha
di <- c("mnad_psi10_scenario")
mean_alpha <- rep(NA, 9)

for(i in 1:9)
{
  name <- paste(di, paste(sc[i]), paste(".csv"), sep = "")
  dat <- read.csv(name, header = TRUE)
  table2[4,i] <- mean(dat$m) # average m
  mean_alpha[i] <- mean(dat$alpha) # average alpha
  dat <- dat[,-which(colnames(dat)=="alpha")] # delete column with alpha to apply percentile_mn
  ul <- t(apply(dat, 1, percentile_mn))
  dat2 <- cbind(dat$psi10, dat$m, ul[,1], ul[,2])
  CI <- t(apply(dat2, 1, consCI_mn))
  cov <- ifelse(CI[,1] < true[i] & CI[,2] > true[i], 1, 0)
  print(prop.test(sum(cov), 1000, p=0.95)$p.value)
  table3[4,i] <- mean(cov)
  table4[4,i] <- mean((CI[,2]-CI[,1]))
}

