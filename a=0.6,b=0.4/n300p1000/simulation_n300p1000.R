library(stats)
library(iterators)
library(parallel)
library(foreach)
library(doParallel)
library(MASS)
library(ncvreg)
library(bda)
library(rJava)
library("xlsx")
library("matrixStats")


simdata <- function(n, p, alpha, beta, binaryOutcome = FALSE, seed) {
  set.seed(seed)
  Z1<-matrix(rbinom(5*n, 1, 0.3),n,5)
  Sigma <- matrix(c(1,0.3,0.3,0.3,0.3,0.3,1,0.3,0.3,0.3,0.3,0.3,1,0.3,0.3,0.3,0.3,0.3,1,0.3,0.3,0.3,0.3,0.3,1),5,5)
  Z2<- mvrnorm(n, rep(0, 5), Sigma)
  Z<-cbind(Z1,Z2)
  #对于不同的人，不同甲基化，混杂变量对中介变量的影响相同�?
  phi<-c(.1,.3,.4,.4,.6,.1,.3,.4,.4,.6)#phi: confounders --> exposure
  e0<-rnorm(100, 0, 1)
  lopr<-Z%*%phi+e0
  pr<-1/(1+exp(-lopr))  # 混杂变量使得X (0,1) 的取值概�?
  X<-matrix(rbinom(n,1,pr),n,1)
  ck <- t(runif(p, 0, 2))
  M <- matrix(0, n, p)
  phik<-c(0.2,0.2,0.3,0.5,0.6,0.2,0.2,0.3,0.5,0.6) #混杂变量对中介变量的影响
  # M(连续�?) #
  for (i in 1:n) {
    e <- rnorm(p, 0, 1.2)
    M[i, ] <- ck + X[i] * alpha + Z[i,] %*% phik *rep(1,p) + e #confounders --> mediators : exp(0.4)
  } 
  colnames(M) <- paste0("M", 1:ncol(M))
  XMZ <- cbind(X, M, Z)  #  [X M]
  B <- c(0.5, beta, phi)  # (p+1) times 1 beta固定是因为对于不同的样本，Mi都不是中间变量，因此与中介变量相关的alpha和beta固定
  E <- rnorm(n, 0, 1)
  Y <- 0.5 + XMZ %*% B + E  #  the response  n times 1
  
  if(binaryOutcome)
    Y <- matrix(rbinom(n, 1, 1/(1+exp(-Y))), nrow = n)
  #计算PS�?
  fit <- glm(X~Z[,1]+Z[,2]+Z[,3]+Z[,4]+Z[,5]+Z[,6]+Z[,7]+Z[,8]+Z[,9]+Z[,10],
             family = binomial(link='logit') )
  PS = predict(fit, type='response') # 选择预测后的输出结果用在binomial数据，response表示输出结果预测响应变量�?1的概�?
  PS <- matrix(PS, n, 1)
  return(list(Y = Y, M = M, X = X, Z = Z, pr=PS, n = n, p = p))
}

results_Pscw_n300p1000<-list() 
results_Iptw_n300p1000<-list() 
results_Pscov_n300p1000<-list() 
results_Non_n300p1000<-list() 
results_Cov_n300p1000<-list() 
for (v in 1:500) {
  n <- 300   # sample size
  p <- 1000 # the dimension of mediators                                                                                                                           
  alpha <- rep(0, p)  # the regression coefficients (exposure --> mediators)
  beta1 <- rep(0, p) # regression coefficients beta (mediators --> outcome);continuous outcomethe 
  #beta2 <- rep(0, p) # for binary outcome
  # the first four markers are true mediators
  alpha[1:8] <- c(0.30,0.36,0.45,0.48,0.60,0.72,0.90,1.20)
  beta1[1:8] <- c(0.20,0.24,0.30,0.32,0.40,0.48,0.60,0.80)
  #beta2[1:4] <- c(1.45,1.5,1.55,1.6)
  # these are not true mediators 
  alpha[11:12] <- 0.90
  beta1[9:10] <- 0.60
  sim_data<-simdata(n,p,alpha,beta1,seed=v)
  bresults_Pscw_n300p1000[[v]]=myhima(sim_data,method='Psuni',v)
  results_Iptw_n300p1000[[v]]=myhima(sim_data,method='Iptw',v)
  results_Pscov_n300p1000[[v]]=myhima(sim_data,method='Pscov',v)
  results_Non_n300p1000[[v]]=myhima(sim_data,method='Non',v)
  results_Cov_n300p1000[[v]]=myhima(sim_data,method='Cov',v)
}


