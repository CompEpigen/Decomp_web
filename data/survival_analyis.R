library(survival)
library(survminer)
library(MeDeCom)
load("FactorViz_outputs/medecom_set.RData")
load("FactorViz_outputs/ann_S.RData")
rem.samples <- grepl("11A",ann.S$tcga_barcodes)
ann.S <- ann.S[!rem.samples,]
props.LMC5 <- getProportions(medecom.set,K=9,lambda=0.001)[6,!rem.samples]
dat <- data.frame(time=as.numeric(as.character(ann.S$days_to_death)),status=as.numeric(ann.S$vital_status),LMC5=as.numeric(as.character(props.LMC5)),age=as.numeric(as.character(ann.S$age_at_initial_pathologic_diagnosis)),sex=ann.S$sex)
dat <- dat[!is.na(dat$status),]
dat$time[dat$status==1] <- as.numeric(as.character(ann.S$days_to_last_followup))[dat$status==1]
dat$status <- (dat$status %% 2)+1
form <- as.formula(Surv(time,status)~LMC5+age+sex)
mod <- coxph(form,dat)
summary(mod)
#                      coef  exp(coef)   se(coef)      z Pr(>|z|)
#LMC5            -6.224e-01  5.366e-01  4.636e-01 -1.343    0.179
#age              2.020e-05  1.000e+00  1.806e-05  1.118    0.263
#sexmale          1.014e-02  1.010e+00  1.317e-01  0.077    0.939
#stagestage i    -4.939e-01  6.102e-01  7.919e-01 -0.624    0.533
#stagestage ia   -3.025e-01  7.390e-01  5.922e-01 -0.511    0.610
#stagestage ib   -5.245e-01  5.919e-01  5.988e-01 -0.876    0.381
#stagestage ii   -1.440e+01  5.562e-07  1.470e+03 -0.010    0.992
#stagestage iia  -2.456e-01  7.822e-01  6.158e-01 -0.399    0.690
#stagestage iib  -4.921e-01  6.114e-01  6.135e-01 -0.802    0.422
#stagestage iiia -4.377e-01  6.455e-01  6.147e-01 -0.712    0.476
#stagestage iiib -8.435e-01  4.302e-01  8.316e-01 -1.014    0.310
#stagestage iv   -1.076e+00  3.410e-01  7.187e-01 -1.497    0.134

library(survival)
library(survminer)
library(MeDeCom)
load("FactorViz_outputs/medecom_set.RData")
load("FactorViz_outputs/ann_S.RData")
rem.samples <- grepl("11A",ann.S$Comment..TCGA.Barcode.)
ann.S <- ann.S[!rem.samples,]
props.LMC6 <- getProportions(medecom.set,K=7,lambda=0.001)[6,!rem.samples]
dat <- data.frame(time=as.numeric(as.character(ann.S$days_to_death))/365,status=as.numeric(ann.S$vital_status),LMC6=as.numeric(as.character(props.LMC6)),age=as.numeric(as.character(ann.S$age_at_diagnosis)),sex=ann.S$gender,stage=ann.S$tumor_stage)
dat$time[dat$status==1] <- as.numeric(as.character(ann.S$days_to_last_follow_up))[dat$status==1]/365
dat$status <- (dat$status %% 2)+1
form <- as.formula(Surv(time,status)~LMC6+age+sex+stage)
mod <- coxph(form,dat)
summary(mod)
#                      coef  exp(coef)   se(coef)      z Pr(>|z|)  
#LMC6             1.057e+00  2.878e+00  4.946e-01  2.137   0.0326 *
#age              1.523e-05  1.000e+00  1.814e-05  0.840   0.4009  
#sexmale          1.532e-02  1.015e+00  1.303e-01  0.118   0.9064  
#stagestage i    -5.592e-01  5.717e-01  7.914e-01 -0.707   0.4799  
#stagestage ia   -2.925e-01  7.464e-01  5.908e-01 -0.495   0.6205  
#stagestage ib   -5.232e-01  5.926e-01  5.944e-01 -0.880   0.3787  
#stagestage ii   -1.431e+01  6.116e-07  1.474e+03 -0.010   0.9923  
#stagestage iia  -2.128e-01  8.083e-01  6.120e-01 -0.348   0.7280  
#stagestage iib  -5.126e-01  5.989e-01  6.110e-01 -0.839   0.4015  
#stagestage iiia -4.505e-01  6.373e-01  6.137e-01 -0.734   0.4629  
#stagestage iiib -7.248e-01  4.844e-01  8.241e-01 -0.880   0.3791  
#stagestage iv   -1.012e+00  3.636e-01  7.128e-01 -1.419   0.1558  

props.LMC1 <- getProportions(medecom.set,K=7,lambda=0.001)[1,!rem.samples]
dat <- data.frame(time=as.numeric(as.character(ann.S$days_to_death))/365,status=as.numeric(ann.S$vital_status),LMC1=as.numeric(as.character(props.LMC1)),age=as.numeric(as.character(ann.S$age_at_diagnosis)),sex=ann.S$gender,stage=ann.S$tumor_stage)
dat$time[dat$status==1] <- as.numeric(as.character(ann.S$days_to_last_follow_up))[dat$status==1]/365
dat$status <- (dat$status %% 2)+1
form <- as.formula(Surv(time,status)~LMC1+age+sex+stage)
mod <- coxph(form,dat)
summary(mod)
#                      coef  exp(coef)   se(coef)      z Pr(>|z|)
#LMC1             4.091e-01  1.505e+00  7.003e-01  0.584    0.559
#age              2.192e-05  1.000e+00  1.865e-05  1.175    0.240
#sexmale         -2.030e-02  9.799e-01  1.294e-01 -0.157    0.875
#stagestage i    -3.693e-01  6.912e-01  7.834e-01 -0.471    0.637
#stagestage ia   -2.745e-01  7.600e-01  5.923e-01 -0.463    0.643
#stagestage ib   -4.370e-01  6.460e-01  5.947e-01 -0.735    0.462
#stagestage ii   -1.438e+01  5.686e-07  1.468e+03 -0.010    0.992
#stagestage iia  -2.014e-01  8.176e-01  6.177e-01 -0.326    0.744
#stagestage iib  -4.143e-01  6.608e-01  6.101e-01 -0.679    0.497
#stagestage iiia -4.068e-01  6.658e-01  6.151e-01 -0.661    0.508
#stagestage iiib -7.148e-01  4.893e-01  8.247e-01 -0.867    0.386
#stagestage iv   -9.911e-01  3.712e-01  7.151e-01 -1.386    0.166

props.LMC2 <- getProportions(medecom.set,K=7,lambda=0.001)[2,!rem.samples]
dat <- data.frame(time=as.numeric(as.character(ann.S$days_to_death))/365,status=as.numeric(ann.S$vital_status),LMC2=as.numeric(as.character(props.LMC2)),age=as.numeric(as.character(ann.S$age_at_diagnosis)),sex=ann.S$gender,stage=ann.S$tumor_stage)
dat$time[dat$status==1] <- as.numeric(as.character(ann.S$days_to_last_follow_up))[dat$status==1]/365
dat$status <- (dat$status %% 2)+1
form <- as.formula(Surv(time,status)~LMC2+age+sex+stage)
mod <- coxph(form,dat)
summary(mod)
#                      coef  exp(coef)   se(coef)      z Pr(>|z|)
#LMC2            -1.367e-01  8.722e-01  5.143e-01 -0.266    0.790
#age              2.011e-05  1.000e+00  1.833e-05  1.097    0.273
#sexmale         -2.847e-02  9.719e-01  1.299e-01 -0.219    0.827
#stagestage i    -3.136e-01  7.308e-01  7.855e-01 -0.399    0.690
#stagestage ia   -2.460e-01  7.820e-01  5.914e-01 -0.416    0.677
#stagestage ib   -3.969e-01  6.724e-01  5.932e-01 -0.669    0.503
#stagestage ii   -1.436e+01  5.792e-07  1.467e+03 -0.010    0.992
#stagestage iia  -1.508e-01  8.601e-01  6.117e-01 -0.246    0.805
#stagestage iib  -3.900e-01  6.770e-01  6.087e-01 -0.641    0.522
#stagestage iiia -3.775e-01  6.856e-01  6.128e-01 -0.616    0.538
#stagestage iiib -6.808e-01  5.062e-01  8.310e-01 -0.819    0.413
#stagestage iv   -9.403e-01  3.905e-01  7.158e-01 -1.314    0.189

props.LMC3 <- getProportions(medecom.set,K=7,lambda=0.001)[3,!rem.samples]
dat <- data.frame(time=as.numeric(as.character(ann.S$days_to_death))/365,status=as.numeric(ann.S$vital_status),LMC3=as.numeric(as.character(props.LMC3)),age=as.numeric(as.character(ann.S$age_at_diagnosis)),sex=ann.S$gender,stage=ann.S$tumor_stage)
dat$time[dat$status==1] <- as.numeric(as.character(ann.S$days_to_last_follow_up))[dat$status==1]/365
dat$status <- (dat$status %% 2)+1
form <- as.formula(Surv(time,status)~LMC3+age+sex+stage)
mod <- coxph(form,dat)
summary(mod)
#                      coef  exp(coef)   se(coef)      z Pr(>|z|)
#LMC3             1.343e-01  1.144e+00  5.777e-01  0.233    0.816
#age              1.950e-05  1.000e+00  1.811e-05  1.077    0.282
#sexmale         -2.184e-02  9.784e-01  1.297e-01 -0.168    0.866
#stagestage i    -3.197e-01  7.264e-01  7.841e-01 -0.408    0.683
#stagestage ia   -2.418e-01  7.852e-01  5.925e-01 -0.408    0.683
#stagestage ib   -3.916e-01  6.760e-01  5.954e-01 -0.658    0.511
#stagestage ii   -1.438e+01  5.695e-07  1.467e+03 -0.010    0.992
#stagestage iia  -1.420e-01  8.677e-01  6.132e-01 -0.232    0.817
#stagestage iib  -3.824e-01  6.822e-01  6.096e-01 -0.627    0.530
#stagestage iiia -3.620e-01  6.963e-01  6.157e-01 -0.588    0.557
#stagestage iiib -6.808e-01  5.062e-01  8.331e-01 -0.817    0.414
#stagestage iv   -9.369e-01  3.918e-01  7.186e-01 -1.304    0.192

props.LMC4 <- getProportions(medecom.set,K=7,lambda=0.001)[4,!rem.samples]
dat <- data.frame(time=as.numeric(as.character(ann.S$days_to_death))/365,status=as.numeric(ann.S$vital_status),LMC4=as.numeric(as.character(props.LMC4)),age=as.numeric(as.character(ann.S$age_at_diagnosis)),sex=ann.S$gender,stage=ann.S$tumor_stage)
dat$time[dat$status==1] <- as.numeric(as.character(ann.S$days_to_last_follow_up))[dat$status==1]/365
dat$status <- (dat$status %% 2)+1
form <- as.formula(Surv(time,status)~LMC4+age+sex+stage)
mod <- coxph(form,dat)
summary(mod)
#                      coef  exp(coef)   se(coef)      z Pr(>|z|)
#LMC4            -1.060e+00  3.466e-01  6.846e-01 -1.548    0.122
#age              1.721e-05  1.000e+00  1.812e-05  0.950    0.342
#sexmale         -8.806e-03  9.912e-01  1.296e-01 -0.068    0.946
#stagestage i    -3.859e-01  6.798e-01  7.824e-01 -0.493    0.622
#stagestage ia   -2.537e-01  7.759e-01  5.907e-01 -0.430    0.668
#stagestage ib   -3.976e-01  6.719e-01  5.919e-01 -0.672    0.502
#stagestage ii   -1.437e+01  5.762e-07  1.470e+03 -0.010    0.992
#stagestage iia  -1.308e-01  8.774e-01  6.113e-01 -0.214    0.831
#stagestage iib  -3.724e-01  6.891e-01  6.083e-01 -0.612    0.540
#stagestage iiia -3.317e-01  7.177e-01  6.131e-01 -0.541    0.589
#stagestage iiib -6.275e-01  5.339e-01  8.258e-01 -0.760    0.447
#stagestage iv   -9.191e-01  3.989e-01  7.132e-01 -1.289    0.197

props.LMC5 <- getProportions(medecom.set,K=7,lambda=0.001)[5,!rem.samples]
dat <- data.frame(time=as.numeric(as.character(ann.S$days_to_death))/365,status=as.numeric(ann.S$vital_status),LMC5=as.numeric(as.character(props.LMC5)),age=as.numeric(as.character(ann.S$age_at_diagnosis)),sex=ann.S$gender,stage=ann.S$tumor_stage)
dat$time[dat$status==1] <- as.numeric(as.character(ann.S$days_to_last_follow_up))[dat$status==1]/365
dat$status <- (dat$status %% 2)+1
form <- as.formula(Surv(time,status)~LMC5+age+sex+stage)
mod <- coxph(form,dat)
summary(mod)
#                      coef  exp(coef)   se(coef)      z Pr(>|z|)
#LMC5            -6.224e-01  5.366e-01  4.636e-01 -1.343    0.179
#age              2.020e-05  1.000e+00  1.806e-05  1.118    0.263
#sexmale          1.014e-02  1.010e+00  1.317e-01  0.077    0.939
#stagestage i    -4.939e-01  6.102e-01  7.919e-01 -0.624    0.533
#stagestage ia   -3.025e-01  7.390e-01  5.922e-01 -0.511    0.610
#stagestage ib   -5.245e-01  5.919e-01  5.988e-01 -0.876    0.381
#stagestage ii   -1.440e+01  5.562e-07  1.470e+03 -0.010    0.992
#stagestage iia  -2.456e-01  7.822e-01  6.158e-01 -0.399    0.690
#stagestage iib  -4.921e-01  6.114e-01  6.135e-01 -0.802    0.422
#stagestage iiia -4.377e-01  6.455e-01  6.147e-01 -0.712    0.476
#stagestage iiib -8.435e-01  4.302e-01  8.316e-01 -1.014    0.310
#stagestage iv   -1.076e+00  3.410e-01  7.187e-01 -1.497    0.134

props.LMC7 <- getProportions(medecom.set,K=7,lambda=0.001)[7,!rem.samples]
dat <- data.frame(time=as.numeric(as.character(ann.S$days_to_death))/365,status=as.numeric(ann.S$vital_status),LMC7=as.numeric(as.character(props.LMC7)),age=as.numeric(as.character(ann.S$age_at_diagnosis)),sex=ann.S$gender,stage=ann.S$tumor_stage)
dat$time[dat$status==1] <- as.numeric(as.character(ann.S$days_to_last_follow_up))[dat$status==1]/365
dat$status <- (dat$status %% 2)+1
form <- as.formula(Surv(time,status)~LMC7+age+sex+stage)
mod <- coxph(form,dat)
summary(mod)
#                      coef  exp(coef)   se(coef)      z Pr(>|z|)
#LMC7             3.158e-01  1.371e+00  7.763e-01  0.407    0.684
#age              1.968e-05  1.000e+00  1.812e-05  1.086    0.277
#sexmale         -3.319e-02  9.674e-01  1.309e-01 -0.254    0.800
#stagestage i    -3.246e-01  7.228e-01  7.818e-01 -0.415    0.678
#stagestage ia   -2.565e-01  7.738e-01  5.910e-01 -0.434    0.664
#stagestage ib   -4.175e-01  6.587e-01  5.929e-01 -0.704    0.481
#stagestage ii   -1.435e+01  5.878e-07  1.468e+03 -0.010    0.992
#stagestage iia  -1.571e-01  8.546e-01  6.118e-01 -0.257    0.797
#stagestage iib  -3.979e-01  6.717e-01  6.089e-01 -0.654    0.513
#stagestage iiia -3.830e-01  6.818e-01  6.129e-01 -0.625    0.532
#stagestage iiib -7.423e-01  4.760e-01  8.294e-01 -0.895    0.371
#stagestage iv   -9.710e-01  3.787e-01  7.135e-01 -1.361    0.174


dat <- data.frame(time=as.numeric(as.character(ann.S$days_to_death))/365,status=as.numeric(ann.S$vital_status),LMC5=factor(ifelse(as.numeric(as.character(props.LMC5))<median(as.numeric(as.character(props.LMC5))),"low","high")))
dat$time[dat$status==1] <- as.numeric(as.character(ann.S$days_to_last_follow_up))[dat$status==1]/365
dat$status <- (dat$status %% 2)+1
mod <- survfit(Surv(time,status)~LMC5,data=dat)
ggsurvplot(mod,data=dat,pval=T)

props.LMC1 <- getProportions(medecom.set,K=7,lambda=0.001)[1,!rem.samples]
dat <- data.frame(time=as.numeric(as.character(ann.S$days_to_death))/365,status=as.numeric(ann.S$vital_status),LMC1=factor(ifelse(as.numeric(as.character(props.LMC1))<median(as.numeric(as.character(props.LMC1))),"low","high")))
dat$time[dat$status==1] <- as.numeric(as.character(ann.S$days_to_last_follow_up))[dat$status==1]/365
dat$status <- (dat$status %% 2)+1
mod <- survfit(Surv(time,status)~LMC1,data=dat)
ggsurvplot(mod,data=dat,pval=T)
ggsave("KM_curve_high_low_LMC1.pdf")
props.LMC2 <- getProportions(medecom.set,K=7,lambda=0.001)[2,!rem.samples]
dat <- data.frame(time=as.numeric(as.character(ann.S$days_to_death))/365,status=as.numeric(ann.S$vital_status),LMC2=factor(ifelse(as.numeric(as.character(props.LMC2))<median(as.numeric(as.character(props.LMC2))),"low","high")))
dat$time[dat$status==1] <- as.numeric(as.character(ann.S$days_to_last_follow_up))[dat$status==1]/365
dat$status <- (dat$status %% 2)+1
mod <- survfit(Surv(time,status)~LMC2,data=dat)
ggsurvplot(mod,data=dat,pval=T)
ggsave("KM_curve_high_low_LMC2.pdf")
props.LMC3 <- getProportions(medecom.set,K=7,lambda=0.001)[3,!rem.samples]
dat <- data.frame(time=as.numeric(as.character(ann.S$days_to_death))/365,status=as.numeric(ann.S$vital_status),LMC3=factor(ifelse(as.numeric(as.character(props.LMC3))<median(as.numeric(as.character(props.LMC3))),"low","high")))
dat$time[dat$status==1] <- as.numeric(as.character(ann.S$days_to_last_follow_up))[dat$status==1]/365
dat$status <- (dat$status %% 2)+1
mod <- survfit(Surv(time,status)~LMC3,data=dat)
ggsurvplot(mod,data=dat,pval=T)
ggsave("KM_curve_high_low_LMC3.pdf")
props.LMC4 <- getProportions(medecom.set,K=7,lambda=0.001)[4,!rem.samples]
dat <- data.frame(time=as.numeric(as.character(ann.S$days_to_death))/365,status=as.numeric(ann.S$vital_status),LMC4=factor(ifelse(as.numeric(as.character(props.LMC4))<median(as.numeric(as.character(props.LMC4))),"low","high")))
dat$time[dat$status==1] <- as.numeric(as.character(ann.S$days_to_last_follow_up))[dat$status==1]/365
dat$status <- (dat$status %% 2)+1
mod <- survfit(Surv(time,status)~LMC4,data=dat)
ggsurvplot(mod,data=dat,pval=T)
ggsave("KM_curve_high_low_LMC4.pdf")
props.LMC5 <- getProportions(medecom.set,K=7,lambda=0.001)[5,!rem.samples]
dat <- data.frame(time=as.numeric(as.character(ann.S$days_to_death))/365,status=as.numeric(ann.S$vital_status),LMC4=factor(ifelse(as.numeric(as.character(props.LMC5))<median(as.numeric(as.character(props.LMC5))),"low","high")))
dat$time[dat$status==1] <- as.numeric(as.character(ann.S$days_to_last_follow_up))[dat$status==1]/365
dat$status <- (dat$status %% 2)+1
mod <- survfit(Surv(time,status)~LMC4,data=dat)
ggsurvplot(mod,data=dat,pval=T)
ggsave("KM_curve_high_low_LMC5.pdf")
props.LMC6 <- getProportions(medecom.set,K=7,lambda=0.001)[6,!rem.samples]
dat <- data.frame(time=as.numeric(as.character(ann.S$days_to_death))/365,status=as.numeric(ann.S$vital_status),LMC6=factor(ifelse(as.numeric(as.character(props.LMC6))<median(as.numeric(as.character(props.LMC6))),"low","high")))
dat$time[dat$status==1] <- as.numeric(as.character(ann.S$days_to_last_follow_up))[dat$status==1]/365
dat$status <- (dat$status %% 2)+1
mod <- survfit(Surv(time,status)~LMC6,data=dat)
ggsurvplot(mod,data=dat,pval=T)
ggsave("KM_curve_high_low_LMC6.pdf")
props.LMC7 <- getProportions(medecom.set,K=7,lambda=0.001)[7,!rem.samples]
dat <- data.frame(time=as.numeric(as.character(ann.S$days_to_death))/365,status=as.numeric(ann.S$vital_status),LMC7=factor(ifelse(as.numeric(as.character(props.LMC7))<median(as.numeric(as.character(props.LMC7))),"low","high")))
dat$time[dat$status==1] <- as.numeric(as.character(ann.S$days_to_last_follow_up))[dat$status==1]/365
dat$status <- (dat$status %% 2)+1
mod <- survfit(Surv(time,status)~LMC7,data=dat)
ggsurvplot(mod,data=dat,pval=T)
ggsave("KM_curve_high_low_LMC7.pdf")

#0.022

props <- as.numeric(as.character(props.LMC6))
q.33 <- quantile(props,.33)
q.66 <- quantile(props,.66)
dat <- data.frame(time=as.numeric(as.character(ann.S$days_to_death))/365,status=as.numeric(ann.S$vital_status),LMC6=cut(props.LMC6,breaks=c(0,q.33,q.66,1),labels=c("low","medium","high")))
dat$time[dat$status==1] <- as.numeric(as.character(ann.S$days_to_last_follow_up))[dat$status==1]/365
dat$status <- (dat$status %% 2)+1
mod <- survfit(Surv(time,status)~LMC6,data=dat)
ggsurvplot(mod,data=dat,pval=T)

#0.19
props <- as.numeric(as.character(props.LMC5))
q.25 <- quantile(props,.25)
q.5 <- quantile(props,.5)
q.75 <- quantile(props,.75)
dat <- data.frame(time=as.numeric(as.character(ann.S$days_to_death))/365,status=as.numeric(ann.S$vital_status),LMC5=cut(props.LMC5,breaks=c(0,q.25,q.5,q.75,1),labels=c("low","q25","q50","high")))
dat$time[dat$status==1] <- as.numeric(as.character(ann.S$days_to_last_follow_up))[dat$status==1]/365
dat$status <- (dat$status %% 2)+1
mod <- survfit(Surv(time,status)~LMC5,data=dat)
ggsurvplot(mod,data=dat,pval=T)

