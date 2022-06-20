
### Adjust the data
intsec=intersect(LUNG$submitter_id.samples,Meth$ID)
LUNG=LUNG[which(LUNG$submitter_id.samples %in% intsec),]
Meth=Meth[which(Meth$ID %in% intsec),]
dim(Meth)
rm(Clinical.choose,fit,site,PS,Bind,Clinical,ID)

loc=match(Meth$ID,LUNG$submitter_id.samples)
LUNG=LUNG[loc,]

####
results_PSW<-list()
results_PSR<-list()
results_PSW<-list()
results_Non<-list()

LUNG=data.frame(LUNG)
Y = as.numeric(matrix(LUNG$dlco_predictive_percent,length(LUNG$dlco_predictive_percent),1))
X = as.numeric(matrix(LUNG$X,length(LUNG$X),1))
X
Z = cbind((apply(as.matrix(LUNG[,2]),2,as.numeric)),LUNG[,c(8,10)])
pr=as.numeric(matrix(LUNG$PS,length(LUNG$PS),1))

sim_data<-list(Y = scale(Y),
               X = as.numeric(matrix(LUNG$X,length(LUNG$X),1)), 
               M = Meth[,-1],
               Z = Z,
               pr=as.numeric(matrix(LUNG$PS,length(LUNG$PS),1)))

