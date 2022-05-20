
myhima=function(sim_data,method,loop)
{
  ############### 缁楊兛绔村銉窗SIS 娴兼媽顓镐吉 閸?? 灏? ############### 
  p=dim(sim_data[['M']])[2]
  n=dim(sim_data[['M']])[1]
  e_a<- rep(0,p) ### estimation of a 
  p_a<-rep(0,p)  ### p-value of a
  
  ## weights calculate
  Nmw<-rep(0,n) ### weights 
  Nmw=ifelse(sim_data$X==1,1/sim_data$pr,1/(1-sim_data$pr)) 
  
  ## 浼? by taking propensity score as covariate
  message("estimating alpha  ..." , "(", Sys.time(), ")",'loop count:',v,',','method=',method)
  
  if(method %in% c('Pscov','Psuni')){
      for(i in 1:p) {
        glmfit<-summary(glm(sim_data[['M']][,i] ~ sim_data$X+sim_data$pr))$coefficient
        e_a[i]<-glmfit[2,1] # 浼导鎷岊吀閸??
        p_a[i]<-glmfit[2,4] #浼导鎷岊吀閻ㄥ埦-value
        } 
  } else if (method %in% c('Iptw')) {
      for(i in 1:p) {
       glmfit<-summary(glm(sim_data[['M']][,i] ~ sim_data$X,weights=Nmw))$coefficient
       e_a[i]<-glmfit[2,1] # 浼导鎷岊吀閸??
       p_a[i]<-glmfit[2,4] #浼导鎷岊吀閻ㄥ埦-value
      } 
  } else if (method %in% c('Non','None')){
      for(i in 1:p) {
      glmfit<-summary(glm(sim_data[['M']][,i] ~ sim_data$X))$coefficient
      e_a[i]<-glmfit[2,1] # 浼导鎷岊吀閸??
      p_a[i]<-glmfit[2,4] #浼导鎷岊吀閻ㄥ埦-value
    } 
  } else if (method %in% c('Cov')){
      for(i in 1:p) {
      glmfit<-summary(glm(sim_data$M[,i] ~ sim_data$X+sim_data$Z[,1]+sim_data$Z[,2]+sim_data$Z[,3]+sim_data$Z[,4]
                          +sim_data$Z[,5]+sim_data$Z[,6]+sim_data$Z[,7]+
                            sim_data$Z[,8]+sim_data$Z[,9]+sim_data$Z[,10]))$coefficient
      e_a[i]<-glmfit[2,1] # α估计值
      p_a[i]<-glmfit[2,4] #α估计的p-value
    } 
  } else{print('method input error')}
  

  ## 灏? by applying weights
  p_b<-rep(0,p) ### p-value of b
  
  if (method %in% c('Iptw','Psuni')) {
    totaleffect<-summary(glm(sim_data$Y~sim_data$X,weights=Nmw))$coefficient[2,1]
  } else if (method %in% c('Pscov')) {
    totaleffect<-summary(glm(sim_data$Y~sim_data$X+sim_data$pr))$coefficient[2,1]
  } else if (method %in% c('Non','None')) {
    totaleffect<-summary(glm(sim_data$Y~sim_data$X))$coefficient[2,1]
  } else if (method %in% c('Cov')) {
    totaleffect <- summary(glm(sim_data$Y~sim_data$X+sim_data$Z[,1]+sim_data$Z[,2]+sim_data$Z[,3]+sim_data$Z[,4]
                             +sim_data$Z[,5]+sim_data$Z[,6]+sim_data$Z[,7]+
                               sim_data$Z[,8]+sim_data$Z[,9]+sim_data$Z[,10]))$coefficient[2,1]
  }
    else{ print('Method input error')
  }
  
  
  message("estimating beta  ..." , "(", Sys.time(), ")",'loop count:',v,',','method=',method)
  
  if(method %in% c('Pscov')){
      for(i in 1:p) {   
        glmfit<-summary(glm(sim_data$Y~sim_data[['M']][,i]+sim_data$X+sim_data$pr))$coefficient
        p_b[i]<-glmfit[2,4] #灏炬导鎷岊吀閻ㄥ埦-value
      }  
  } else if (method %in% c('Iptw','Psuni')){
      for(i in 1:p) { 
        coef_x<-summary(glm(sim_data$Y~sim_data$X+sim_data[['M']][,i],weights=Nmw))$coefficient[2,1]
        glmfit<-summary(glm((sim_data$Y-coef_x*sim_data$X)~sim_data[['M']][,i]))$coefficient
        p_b[i]<-glmfit[2,4] #灏炬导鎷岊吀閻ㄥ埦-value
      }  
  } else if (method %in% c('Non','None')){
      for(i in 1:p) 
      { 
        glmfit<-summary(glm(sim_data$Y~sim_data$X+sim_data[['M']][,i]))$coefficient
        p_b[i]<-glmfit[3,4] #灏炬导鎷岊吀閻ㄥ埦-value
      }
  } else if (method %in% c('Cov')){
      for(i in 1:p) 
      {
        glmfit<-summary(glm(sim_data$Y~sim_data$M[,i]+sim_data$X+sim_data$Z[,1]+sim_data$Z[,2]+sim_data$Z[,3]+sim_data$Z[,4]
                          +sim_data$Z[,5]+sim_data$Z[,6]+sim_data$Z[,7]+
                            sim_data$Z[,8]+sim_data$Z[,9]+sim_data$Z[,10]))$coefficient
        p_b[i]<-glmfit[2,4] #β估计的p-value
      }  
  } else{print('method input error')} 
  
  ## 閺嶈宓乸-value 缁涙盯鈧）?
  b_sort<-sort(p_b)
  d<-ceiling(2*n/log(n))
  ID <- which(p_b<= b_sort[d])  # the index of top mediators
  M_SIS <- sim_data[['M']][, ID] 
  XM <- cbind(M_SIS, sim_data$X)
  
  library(ncvreg)
  ############### 缁楊兛绨╁銉窗MCP 鏉╂稐绔村銉х摣闁）? ###############
  
  if (method %in% c('None','Non','Psuni','Iptw')){
    MCPfit <- ncvreg(XM, sim_data$Y, 
                  penalty.factor = c(rep(1, ncol(M_SIS)),0))
  } else if (method %in% c('Pscov')) {
    XM <- cbind(M_SIS, sim_data$X,sim_data$pr)
    MCPfit <- ncvreg(XM, sim_data$Y, 
                   penalty.factor = c(rep(1, ncol(M_SIS)),0,0))
  } else if (method %in% c('Cov')) {
    XM <- cbind(M_SIS, sim_data$X,sim_data$Z[,1],sim_data$Z[,2],sim_data$Z[,3],sim_data$Z[,4]
                ,sim_data$Z[,5],sim_data$Z[,6],sim_data$Z[,7],
                sim_data$Z[,8],sim_data$Z[,9],sim_data$Z[,10])
    MCPfit <- ncvreg(XM, sim_data$Y, 
                     penalty.factor = c(rep(1, ncol(M_SIS)), 0,0,0,0,0,0,0,0,0,0,0))
  } else {print('method input error')}
  
  
  #plot(MCPfit)
  lam <- MCPfit$lambda[which.min(BIC(MCPfit))]
  MCPcoefficients <- coef(MCPfit, lambda = lam)
  est <- MCPcoefficients[2:(d + 1)]
  ID_1_non <- which(est != 0)
  beta_est <- est[ID_1_non]  # The non-zero MCP estimators of beta
  ID_test <- ID[ID_1_non]  # The index of the ID of non-zero beta in Y ~ M
  
  ############### 閼辨柨鎮庨弰鎹愭啿閹勵梾妤?? ###############
  
  ## 浼?
  a_N <- e_a[ID_test] # 濮濓絽鐖堕弶鍐櫢閻ㄥ嫀?
  
  
  alpha_est_ID_test <- as.numeric(a_N)  #  the estimator for alpha
  alpha <- p_a[ID_test, drop = FALSE]
  P_adjust_alpha <- length(ID_test) * alpha  # the adjusted p-value for alpha (bonferroni)
  P_adjust_alpha[P_adjust_alpha > 1] <- 1
  P_fdr_alpha <- p.adjust(alpha, "fdr")  # the adjusted p-value for alpha (FDR)
  alpha_est <- alpha_est_ID_test
  P_BH_alpha<- p.adjust(alpha,'BH')
  P_BY_alpha<- p.adjust(alpha,'BY')
  
  # 灏?
  YMX <- data.frame(Y = sim_data$Y, sim_data[['M']][, ID_test, drop = FALSE], X = sim_data$X)
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
  gamma_est <- totaleffect
  results<- data.frame(alpha = alpha_est, beta = beta_est, gamma = gamma_est, 
                              `alpha*beta` = ab_est, `% total effect` = ab_est/gamma_est * 100, 
                              `BF.P` = P_value, `BH.FDR` = FDR,`BY.p` = BY, check.names = FALSE)
  
  return(results)
}

