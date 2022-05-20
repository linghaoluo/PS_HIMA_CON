

## a*b estimation  
{
estimate_ab_Pscw<-matrix(c(rep(0,4000)),500,8,
                    dimnames=list(NULL,c("M1", "M2","M3", "M4","M5","M6","M7","M8")))

estimate_ab_Iptw<-matrix(c(rep(0,4000)),500,8,
                         dimnames=list(NULL,c("M1", "M2","M3", "M4","M5","M6","M7","M8")))

estimate_ab_Non<-matrix(c(rep(0,4000)),500,8,
                         dimnames=list(NULL,c("M1", "M2","M3", "M4","M5","M6","M7","M8")))

estimate_ab_Pscov<-matrix(c(rep(0,4000)),500,8,
                         dimnames=list(NULL,c("M1", "M2","M3", "M4","M5","M6","M7","M8")))

estimate_ab_Cov<-matrix(c(rep(0,4000)),500,8,
                          dimnames=list(NULL,c("M1", "M2","M3", "M4","M5","M6","M7","M8")))

for (i in 1:500){
  estimate_ab_Pscw[i,1]<-results_Pscw_n500p1000[[i]]["M1",4]  
  estimate_ab_Pscw[i,2]<-results_Pscw_n500p1000[[i]]["M2",4]
  estimate_ab_Pscw[i,3]<-results_Pscw_n500p1000[[i]]["M3",4]
  estimate_ab_Pscw[i,4]<-results_Pscw_n500p1000[[i]]["M4",4]
  estimate_ab_Pscw[i,5]<-results_Pscw_n500p1000[[i]]["M5",4]
  estimate_ab_Pscw[i,6]<-results_Pscw_n500p1000[[i]]["M6",4]
  estimate_ab_Pscw[i,7]<-results_Pscw_n500p1000[[i]]["M7",4]
  estimate_ab_Pscw[i,8]<-results_Pscw_n500p1000[[i]]["M8",4]
  
  estimate_ab_Iptw[i,1]<-results_Iptw_n500p1000[[i]]["M1",4]
  estimate_ab_Iptw[i,2]<-results_Iptw_n500p1000[[i]]["M2",4]
  estimate_ab_Iptw[i,3]<-results_Iptw_n500p1000[[i]]["M3",4]
  estimate_ab_Iptw[i,4]<-results_Iptw_n500p1000[[i]]["M4",4]
  estimate_ab_Iptw[i,5]<-results_Iptw_n500p1000[[i]]["M5",4]
  estimate_ab_Iptw[i,6]<-results_Iptw_n500p1000[[i]]["M6",4]
  estimate_ab_Iptw[i,7]<-results_Iptw_n500p1000[[i]]["M7",4]
  estimate_ab_Iptw[i,8]<-results_Iptw_n500p1000[[i]]["M8",4]
  
  estimate_ab_Non[i,1]<-results_Non_n500p1000[[i]]["M1",4]
  estimate_ab_Non[i,2]<-results_Non_n500p1000[[i]]["M2",4]
  estimate_ab_Non[i,3]<-results_Non_n500p1000[[i]]["M3",4]
  estimate_ab_Non[i,4]<-results_Non_n500p1000[[i]]["M4",4]
  estimate_ab_Non[i,5]<-results_Non_n500p1000[[i]]["M5",4]
  estimate_ab_Non[i,6]<-results_Non_n500p1000[[i]]["M6",4]
  estimate_ab_Non[i,7]<-results_Non_n500p1000[[i]]["M7",4]
  estimate_ab_Non[i,8]<-results_Non_n500p1000[[i]]["M8",4]
  
  estimate_ab_Pscov[i,1]<-results_Pscov_n500p1000[[i]]["M1",4]
  estimate_ab_Pscov[i,2]<-results_Pscov_n500p1000[[i]]["M2",4]
  estimate_ab_Pscov[i,3]<-results_Pscov_n500p1000[[i]]["M3",4]
  estimate_ab_Pscov[i,4]<-results_Pscov_n500p1000[[i]]["M4",4]
  estimate_ab_Pscov[i,5]<-results_Pscov_n500p1000[[i]]["M5",4]
  estimate_ab_Pscov[i,6]<-results_Pscov_n500p1000[[i]]["M6",4]
  estimate_ab_Pscov[i,7]<-results_Pscov_n500p1000[[i]]["M7",4]
  estimate_ab_Pscov[i,8]<-results_Pscov_n500p1000[[i]]["M8",4]
  
  estimate_ab_Cov[i,1]<-results_Cov_n500p1000[[i]]["M1",4]
  estimate_ab_Cov[i,2]<-results_Cov_n500p1000[[i]]["M2",4]
  estimate_ab_Cov[i,3]<-results_Cov_n500p1000[[i]]["M3",4]
  estimate_ab_Cov[i,4]<-results_Cov_n500p1000[[i]]["M4",4]
  estimate_ab_Cov[i,5]<-results_Cov_n500p1000[[i]]["M5",4]
  estimate_ab_Cov[i,6]<-results_Cov_n500p1000[[i]]["M6",4]
  estimate_ab_Cov[i,7]<-results_Cov_n500p1000[[i]]["M7",4]
  estimate_ab_Cov[i,8]<-results_Cov_n500p1000[[i]]["M8",4]
} 


est_Pscw=data.frame(rbind(colMeans(estimate_ab_Pscw,na.rm=T),
                 (colSds(estimate_ab_Pscw,na.rm=T))**2),row.names = c('mean','sd'))
est_Pscw['method']='Pscw'

est_Pscov=data.frame(rbind(colMeans(estimate_ab_Pscov,na.rm=T),
                          (colSds(estimate_ab_Pscov,na.rm=T))**2),row.names = c('mean','sd'))
est_Pscov['method']='Pscov'

est_Iptw=data.frame(rbind(colMeans(estimate_ab_Iptw,na.rm=T),
                          (colSds(estimate_ab_Iptw,na.rm=T))**2),row.names = c('mean','sd'))
est_Iptw['method']='Iptw'

est_Cov=data.frame(rbind(colMeans(estimate_ab_Cov,na.rm=T),
                         (colSds(estimate_ab_Cov,na.rm=T))**2),row.names = c('mean','sd'))
est_Cov['method']='Cov'

est_Non=data.frame(rbind(colMeans(estimate_ab_Non,na.rm=T),
                          (colSds(estimate_ab_Non,na.rm=T))**2),row.names = c('mean','sd'))
est_Non['method']='Non'

est=rbind(est_Pscw,est_Pscov,est_Iptw,est_Cov,est_Non)
rm(est_Pscw,est_Pscov,est_Iptw,est_Cov,est_Non,estimate_ab_Cov,estimate_ab_Iptw,estimate_ab_Non,
   estimate_ab_Pscov,estimate_ab_Pscw)
}

## MCP&testing estimation
choose=function(results)
{
  M<-matrix(c(rep(0,4000)),500,8,dimnames = list(NULL,c("M1", "M2","M3", "M4","M5","M6","M7","M8")))
  Mt<-matrix(c(rep(0,4000)),500,8,dimnames = list(NULL,c("M1", "M2","M3", "M4","M5","M6","M7","M8")))
  for(i in 1:500){
    M[i,1]<-"M1"%in%rownames(results[[i]])
    M[i,2]<-"M2"%in%rownames(results[[i]])
    M[i,3]<-"M3"%in%rownames(results[[i]])
    M[i,4]<-"M4"%in%rownames(results[[i]])
    M[i,5]<-"M5"%in%rownames(results[[i]])
    M[i,6]<-"M6"%in%rownames(results[[i]])
    M[i,7]<-"M7"%in%rownames(results[[i]])
    M[i,8]<-"M8"%in%rownames(results[[i]])
    
    Mt[i,1]<-"M1"%in%rownames(results[[i]]) & results[[i]]['M1','BH.FDR']<0.05
    Mt[i,2]<-"M2"%in%rownames(results[[i]]) & results[[i]]['M2','BH.FDR']<0.05
    Mt[i,3]<-"M3"%in%rownames(results[[i]]) & results[[i]]['M3','BH.FDR']<0.05
    Mt[i,4]<-"M4"%in%rownames(results[[i]]) & results[[i]]['M4','BH.FDR']<0.05
    Mt[i,5]<-"M5"%in%rownames(results[[i]]) & results[[i]]['M5','BH.FDR']<0.05
    Mt[i,6]<-"M6"%in%rownames(results[[i]]) & results[[i]]['M6','BH.FDR']<0.05
    Mt[i,7]<-"M7"%in%rownames(results[[i]]) & results[[i]]['M7','BH.FDR']<0.05
    Mt[i,8]<-"M8"%in%rownames(results[[i]]) & results[[i]]['M8','BH.FDR']<0.05
  }
  
  return(data.frame(rbind(colSums(M),colSums(Mt)),row.names=c('MCPselection','Testing'))) 
}

{
  Choose_Cov=data.frame(choose(results_Cov_n500p1000))
  Choose_Cov['Method']='Cov'
 
  Choose_Iptw=data.frame(choose(results_Iptw_n500p1000))
  Choose_Iptw['Method']='Iptw'
  
  Choose_Pscov=data.frame(choose(results_Pscov_n500p1000))
  Choose_Pscov['Method']='Pscov'
  
  Choose_Non=data.frame(choose(results_Non_n500p1000))
  Choose_Non['Method']='Non'
  
  Choose_Psuni=data.frame(choose(results_Pscw_n500p1000))
  Choose_Psuni['Method']='Psuni'
  
  MCPandtesting=data.frame(rbind(Choose_Cov,Choose_Iptw,Choose_Pscov,Choose_Non,Choose_Psuni))
  rm(Choose_Cov,Choose_Iptw,Choose_Pscov,Choose_Non,Choose_Psuni)
}


# 假阳性-平均每一次假阳性个数
fpr_num=function(results){
  fpr<-c()
  for (j in 1:500)
  { 
    # results = results_Cov_n500p10000
    miss<-c()
    fp<-setdiff(rownames(results[[j]]),c("M1", "M2","M3", "M4",'M5','M6','M7','M8')) 
    for (i in 1:length(fp))
    {
      miss[i]=results[[j]][fp[i],'BH.FDR']<0.05
      fpr[j]<-sum(miss)
    }
  }
  return(mean(fpr))
}    

# FDR (假阳性/总阳性)
FDR=function(results){
  fpr<-c()
  tp<-c()
  fdr<-c()
  for (j in 1:500)
  { 
    #results = results_Iptw_n500p10000
    miss<-c()
    fp<-setdiff(rownames(results[[j]]),c("M1", "M2","M3", "M4",'M5','M6','M7','M8')) 
    for (i in 1:length(fp))
    {
      miss[i]=results[[j]][fp[i],'BH.FDR']<0.05
      fpr[j]<-sum(miss)
    }
    tp[j] = sum(results[[j]][,'BH.FDR']<0.05)
  }
  return(mean(fpr/tp,na.rm = TRUE))
}    


Fpr=rbind(fpr_num(results_Cov_n500p1000),fpr_num(results_Pscov_n500p1000),
          fpr_num(results_Iptw_n500p1000),fpr_num(results_Pscw_n500p1000),fpr_num(results_Non_n500p1000))
Fpr=data.frame(Fpr,
               row.names=c('Cov','Pscov','Iptw','Pscw','Non'))
Fpr                                                                                                               

# FDR - 每一次假阳性比例
Fdr=rbind(FDR(results_Cov_n500p1000),FDR(results_Pscov_n500p1000),
          FDR(results_Iptw_n500p1000),FDR(results_Pscw_n500p1000),FDR(results_Non_n500p1000))
Fdr=data.frame(Fdr,
               row.names=c('Cov','Pscov','Iptw','Pscw','Non'))
Fdr
