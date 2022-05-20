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
  #对于不同的人，不同甲基化，混杂变量对中介变量的影响相同。
  phi<-c(.1,.3,.4,.4,.6,.1,.3,.4,.4,.6)#phi: confounders --> exposure
  e0<-rnorm(100, 0, 1)
  lopr<-Z%*%phi+e0
  pr<-1/(1+exp(-lopr))  # 混杂变量使得X (0,1) 的取值概率
  X<-matrix(rbinom(n,1,pr),n,1)
  ck <- t(runif(p, 0, 2))
  M <- matrix(0, n, p)
  phik<-c(0.2,0.2,0.3,0.5,0.6,0.2,0.2,0.3,0.5,0.6) #混杂变量对中介变量的影响
  # M(连续型) #
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
  #计算PS值
  fit <- glm(X~Z[,1]+Z[,2]+Z[,3]+Z[,4]+Z[,5]+Z[,6]+Z[,7]+Z[,8]+Z[,9]+Z[,10],
             family = binomial(link='logit') )
  PS = predict(fit, type='response') # 选择预测后的输出结果用在binomial数据，response表示输出结果预测响应变量为1的概率
  PS <- matrix(PS, n, 1)
  return(list(Y = Y, M = M, X = X, Z = Z, pr=PS, n = n, p = p))
}


############### Pre-settings ###############
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
results_Pscw_n300p1000<-list() 


loopcount<-1
v=1 ### seeds
  for (v in 1:500) 
{ 

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
  
############### 第一步：SIS 估计α 和 β ############### 
    
e_a<- rep(0,p) ### estimation of a 
p_a<-rep(0,p)  ### p-value of a

## weights calculate
Nmw<-rep(0,n) ### weights 
Nmw=ifelse(sim_data$X==1,1/sim_data$pr,1/(1-sim_data$pr)) 

## α by taking propensity score as covariate
    message("estimating α  ..." , "(", Sys.time(), ")","(", loopcount, ")")
for(i in 1:p) 
{
  glmfit<-summary(glm(sim_data$M[,i] ~ sim_data$X+sim_data$pr))$coefficient
  e_a[i]<-glmfit[2,1] # α估计值
  p_a[i]<-glmfit[2,4] #α估计的p-value
} 
# p_a_BH<-p.adjust(p_a,method='BH')

    
## β by applying weights
    
p_b<-rep(0,p) ### p-value of b
coef_x<-summary(glm(sim_data$Y~sim_data$X,weights=Nmw))$coefficient[2,1]

      message("estimating β ...", "(", Sys.time(), ")","(", loopcount, ")")
for(i in 1:p) 
{ coef_x<-summary(glm(sim_data$Y~sim_data$X+sim_data$M[,i],weights=Nmw))$coefficient[2,1]
  glmfit<-summary(glm((sim_data$Y-coef_x*sim_data$X)~sim_data$M[,i]))$coefficient
  p_b[i]<-glmfit[2,4] #β估计的p-value
}  

## 根据p-value 筛选β
b_sort<-sort(p_b)
d<-ceiling(2*n/log(n))
ID <- which(p_b<= b_sort[d])  # the index of top mediators
M_SIS <- sim_data$M[, ID] 
XM <- cbind(M_SIS, sim_data$X)


############### 第二步：MCP 进一步筛选β ###############
MCPfit <- ncvreg(XM, sim_data$Y, 
                 penalty.factor = c(rep(1, ncol(M_SIS)),0))


#plot(MCPfit)
lam <- MCPfit$lambda[which.min(BIC(MCPfit))]
MCPcoefficients <- coef(MCPfit, lambda = lam)
est <- MCPcoefficients[2:(d + 1)]
ID_1_non <- which(est != 0)
beta_est <- est[ID_1_non]  # The non-zero MCP estimators of beta
ID_test <- ID[ID_1_non]  # The index of the ID of non-zero beta in Y ~ M

############### 联合显著性检验 ###############
 
## α
a_N <- e_a[ID_test] # 正常权重的α

 
alpha_est_ID_test <- as.numeric(a_N)  #  the estimator for alpha
alpha <- p_a[ID_test, drop = FALSE]
P_adjust_alpha <- length(ID_test) * alpha  # the adjusted p-value for alpha (bonferroni)
P_adjust_alpha[P_adjust_alpha > 1] <- 1
P_fdr_alpha <- p.adjust(alpha, "fdr")  # the adjusted p-value for alpha (FDR)
alpha_est <- alpha_est_ID_test
P_BH_alpha<- p.adjust(alpha,'BH')
P_BY_alpha<- p.adjust(alpha,'BY')

# β
YMX <- data.frame(Y = sim_data$Y, sim_data$M[, ID_test, drop = FALSE], X = sim_data$X)
res <- summary(glm(Y ~ ., data = YMX))$coefficients
est <- res[2:(length(ID_test) + 1), 1]  # the estimator for beta
pvalue<-res[2:(length(ID_test) + 1), 4]
P_adjust_beta <- length(ID_test) * pvalue  # the adjused p-value for beta (bonferroni)
P_adjust_beta[P_adjust_beta > 1] <- 1
P_fdr_beta <- p.adjust(pvalue, "fdr")  # the adjusted p-value for beta (FDR)
ab_est <- alpha_est * beta_est
P_BH_beta<- p.adjust(pvalue,'BH')
P_BY_beta<- p.adjust(pvalue,'BY')

## Use the maximum value as p value 
PA <- rbind(P_adjust_beta, P_adjust_alpha)
P_value <- apply(PA, 2, max)
FDRA <- rbind(P_fdr_beta, P_fdr_alpha)
FDR <- apply(FDRA, 2, max)
BYA <- rbind(P_BY_beta, P_BY_alpha)
BY <- apply(BYA, 2, max)


# Total effect
gamma_est <- coef_x
results_Pscw_n300p1000[[v]] <- data.frame(alpha = alpha_est, beta = beta_est, gamma = gamma_est, 
                      `alpha*beta` = ab_est, `% total effect` = ab_est/gamma_est * 100, 
                      `BF.P` = P_value, `BH.FDR` = FDR,`BY.p` = BY, check.names = FALSE)
loopcount<-loopcount + 1
}


# 估计值
  estimate_a<-matrix(c(rep(0,4000)),500,8,dimnames=list(NULL,c("M1", "M2","M3", "M4","M5","M6","M7","M8")))
  estimate_b<-matrix(c(rep(0,4000)),500,8,dimnames=list(NULL,c("M1", "M2","M3", "M4","M5","M6","M7","M8")))
  estimate_ab<-matrix(c(rep(0,4000)),500,8,dimnames=list(NULL,c("M1", "M2","M3", "M4","M5","M6","M7","M8")))
  for (i in 1:500){
    estimate_a[i,1]<-results_Pscw_n300p1000[[i]]["M1",1]
    estimate_a[i,2]<-results_Pscw_n300p1000[[i]]["M2",1]
    estimate_a[i,3]<-results_Pscw_n300p1000[[i]]["M3",1]
    estimate_a[i,4]<-results_Pscw_n300p1000[[i]]["M4",1]
    estimate_a[i,5]<-results_Pscw_n300p1000[[i]]["M5",1]
    estimate_a[i,6]<-results_Pscw_n300p1000[[i]]["M6",1]
    estimate_a[i,7]<-results_Pscw_n300p1000[[i]]["M7",1]
    estimate_a[i,8]<-results_Pscw_n300p1000[[i]]["M8",1]
  }
  for (i in 1:500){
    estimate_b[i,1]<-results_Pscw_n300p1000[[i]]["M1",2]
    estimate_b[i,2]<-results_Pscw_n300p1000[[i]]["M2",2]
    estimate_b[i,3]<-results_Pscw_n300p1000[[i]]["M3",2]
    estimate_b[i,4]<-results_Pscw_n300p1000[[i]]["M4",2]
    estimate_b[i,5]<-results_Pscw_n300p1000[[i]]["M5",2]
    estimate_b[i,6]<-results_Pscw_n300p1000[[i]]["M6",2]
    estimate_b[i,7]<-results_Pscw_n300p1000[[i]]["M7",2]
    estimate_b[i,8]<-results_Pscw_n300p1000[[i]]["M8",2]
  }  
  
  for (i in 1:500){
    estimate_ab[i,1]<-results_Pscw_n300p1000[[i]]["M1",4]
    estimate_ab[i,2]<-results_Pscw_n300p1000[[i]]["M2",4]
    estimate_ab[i,3]<-results_Pscw_n300p1000[[i]]["M3",4]
    estimate_ab[i,4]<-results_Pscw_n300p1000[[i]]["M4",4]
    estimate_ab[i,5]<-results_Pscw_n300p1000[[i]]["M5",4]
    estimate_ab[i,6]<-results_Pscw_n300p1000[[i]]["M6",4]
    estimate_ab[i,7]<-results_Pscw_n300p1000[[i]]["M7",4]
    estimate_ab[i,8]<-results_Pscw_n300p1000[[i]]["M8",4]
  } 
  a_mean_std<-rbind(colMeans(estimate_a,na.rm=TRUE),colVars(estimate_a,na.rm=TRUE))
  a_mean_std
  b_mean_std<-rbind(colMeans(estimate_b,na.rm=TRUE),colVars(estimate_b,na.rm=TRUE))
  b_mean_std
  ab_mean_std<-rbind(colMeans(estimate_ab,na.rm=TRUE),colVars(estimate_ab,na.rm=TRUE))
  ab_mean_std
  est_Psuni<-ab_mean_std
  
  # MCP筛选结果
  M<-matrix(c(rep(0,4000)),500,8,dimnames = list(NULL,c("M1", "M2","M3", "M4","M5","M6","M7","M8")))
  for(i in 1:500){
    M[i,1]<-"M1"%in%rownames(results_Pscw_n300p1000[[i]])
    M[i,2]<-"M2"%in%rownames(results_Pscw_n300p1000[[i]])
    M[i,3]<-"M3"%in%rownames(results_Pscw_n300p1000[[i]])
    M[i,4]<-"M4"%in%rownames(results_Pscw_n300p1000[[i]])
    M[i,5]<-"M5"%in%rownames(results_Pscw_n300p1000[[i]])
    M[i,6]<-"M6"%in%rownames(results_Pscw_n300p1000[[i]])
    M[i,7]<-"M7"%in%rownames(results_Pscw_n300p1000[[i]])
    M[i,8]<-"M8"%in%rownames(results_Pscw_n300p1000[[i]])
  }
  sum_M<-matrix(c(rep(0,16)),2,8,dimnames = list(c("0","1"),c("M1", "M2","M3", "M4","M5","M6","M7","M8")))
  for(i in 1:8){
    sum_M[1,i]<-table(M[,i])["0"]
    sum_M[2,i]<-table(M[,i])["1"]
  }
  sum_M
  mcp_Psuni<-sum_M
  # 真阳性
  M<-matrix(c(rep(0,4000)),500,8,dimnames = list(NULL,c("M1", "M2","M3", "M4","M5","M6","M7","M8")))
  for(i in 1:500){
    M[i,1]<-"M1"%in%rownames(results_Pscw_n300p1000[[i]]) & results_Pscw_n300p1000[[i]]['M1','BH.FDR']<0.05
    M[i,2]<-"M2"%in%rownames(results_Pscw_n300p1000[[i]]) & results_Pscw_n300p1000[[i]]['M2','BH.FDR']<0.05
    M[i,3]<-"M3"%in%rownames(results_Pscw_n300p1000[[i]]) & results_Pscw_n300p1000[[i]]['M3','BH.FDR']<0.05
    M[i,4]<-"M4"%in%rownames(results_Pscw_n300p1000[[i]]) & results_Pscw_n300p1000[[i]]['M4','BH.FDR']<0.05
    M[i,5]<-"M5"%in%rownames(results_Pscw_n300p1000[[i]]) & results_Pscw_n300p1000[[i]]['M5','BH.FDR']<0.05
    M[i,6]<-"M6"%in%rownames(results_Pscw_n300p1000[[i]]) & results_Pscw_n300p1000[[i]]['M6','BH.FDR']<0.05
    M[i,7]<-"M7"%in%rownames(results_Pscw_n300p1000[[i]]) & results_Pscw_n300p1000[[i]]['M7','BH.FDR']<0.05
    M[i,8]<-"M8"%in%rownames(results_Pscw_n300p1000[[i]]) & results_Pscw_n300p1000[[i]]['M8','BH.FDR']<0.05
  }
  sum_Mt<-matrix(c(rep(0,16)),2,8,dimnames = list(c("0","1"),c("M1", "M2","M3", "M4","M5","M6","M7","M8")))
  for(i in 1:8){
    sum_Mt[1,i]<-table(M[,i])["0"]
    sum_Mt[2,i]<-table(M[,i])["1"]
  }
  sum_Mt
  test_Psuni<-sum_Mt
  # 假阳性
  miss<-c()
  for(i in 1:500){ 
    miss[i]<-length(setdiff(rownames(results_Pscw_n300p1000[[i]]),c("M1", "M2","M3", "M4","M5","M6")))
  }
  round(mean(miss,na.rm = FALSE),4)
  
  fpr<-c()
  for (j in 1:500)
  { 
    miss<-c()
    fp<-setdiff(rownames(results_Pscw_n300p1000[[j]]),c("M1", "M2","M3", "M4")) 
    for (i in 1:length(fp))
    {
      miss[i]= fp[i] %in% rownames(results_Pscw_n300p1000[[j]]) & results_Pscw_n300p1000[[j]][fp[i],'BH.FDR']<0.05
      fpr[j]<-sum(miss)
    }
    
  }
  mean(fpr)
  fp_Psuni<-mean(fpr)
  
  # 功效图
  # 功效图
  power_Psuni<-(sum_Mt[2,]/500)[1:8]
  a<- c(0.30,0.36,0.45,0.48,0.60,0.72,0.90,1.20)

  b<- c(0.20,0.24,0.30,0.32,0.40,0.48,0.60,0.80)
  t<-a/0.6
  plot(t,power_Psuni)
  
  