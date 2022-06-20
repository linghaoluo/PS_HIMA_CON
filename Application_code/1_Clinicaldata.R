library(data.table)

# Working directory
setwd('D:/Personal working files/Lab/Application_final_v')

############# Clinical data reading ################
LUAD.GDC<-(fread('TCGA-LUAD.GDC_phenotype.tsv',header=T))
LUSC.GDC<-(fread('TCGA-LUSC.GDC_phenotype.tsv',header=T))

vari_choose<-c("dlco_predictive_percent","age_at_initial_pathologic_diagnosis",
               "gender.demographic","submitter_id",'submitter_id.samples','tobacco_smoking_history','race.demographic')

LUAD.choose<-LUAD.GDC[,c("dlco_predictive_percent","age_at_initial_pathologic_diagnosis",'race.demographic',
                         "gender.demographic","submitter_id",'submitter_id.samples','tobacco_smoking_history')]
LUSC.choose<-LUSC.GDC[,c("dlco_predictive_percent","age_at_initial_pathologic_diagnosis",'race.demographic',
                         "gender.demographic","submitter_id",'submitter_id.samples','tobacco_smoking_history')]

LUNG<-rbind(LUAD.choose,LUSC.choose)
sex<-ifelse(LUNG$gender.demographic=='male',1,0 )
LUNG$gender=sex

## Clean the NAs
LUNG<-na.omit(LUNG)
rm(LUAD.choose,LUAD.GDC,LUSC.choose,LUSC.GDC)

## adjust the exposure 
X<-ifelse(LUNG[,'tobacco_smoking_history']==2|LUNG[,'tobacco_smoking_history']==4,1,0) #take 1 as current smoke or quit<15
X<-as.numeric(X)
LUNG<-cbind(LUNG,X)
LUNG=data.frame(LUNG)

table(LUNG$race.demographic)
LUNG$race=ifelse(LUNG$race.demographic=='white',1,0)

# Calculate PS
fit <- glm(as.numeric(X) ~ as.numeric(age_at_initial_pathologic_diagnosis)+gender+race,data=LUNG,
           family = binomial(link='logit') )

summary(fit)
PS = predict(fit, type='response') # predict the PS
PS <- matrix(PS, length(PS), 1)
hist(PS)
LUNG<-cbind(LUNG,PS)


