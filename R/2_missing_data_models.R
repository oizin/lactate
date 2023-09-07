###############################################################################
# Analysis of missing data
#
# Contents:
# - comparison of missing vs. non missing lactate pops
# - missingness logistic regression model
# - analysis of weighting incl. table 1
#
#
## IMPORT PACKAGES  ============================================================
# packages
library(ggplot2)
library(data.table)
library(ggExtra)
library(tableone)
library(survey)
library(mgcv)
library(xgboost)
library(fastDummies)
source("R/functions/cross_validation.R")

SEED <- 111

## IMPORT EICU =================================================================
# data
icu <- fread("data/cohort.csv")

# choice of lactate and glucose vars
icu[,glucose := glucose_lact]
icu[,lactate := lactate_gluc]

# descriptives
icu[,.N]
icu[,.N,by=missing_lactate]

## MISSING WEIGHTS =============================================================

#### Logistic regression missingness model #####
# missing lactate var
icu[,missing_lactate1 := fifelse(is.na(lactate),1,0,0)]

# missingness model
m <- glm(missing_lactate1 ~ 
           apacheivascore*sepsis_apache +
           operative + 
           oobventday1 +
           oobintubday1 +
           organ_system + 
           age +
           unittype_f,
         data=icu,family=binomial())
summary(m)
m_res <- broom::tidy(m,conf.int = FALSE, conf.level = 0.95, exponentiate = FALSE)
setDT(m_res)
m_res[,lower := estimate - 1.96*std.error]
m_res[,higher := estimate + 1.96*std.error]
fwrite(m_res,file = "results/missing_model_glm.csv")
# save the model
saveRDS(m,file = "models/missingness_model_glm.rds")

# add weights to data
icu[!is.na(lactate),missing_weight_glm := 1/(1-predict(m,newdata=icu[!is.na(lactate)],
                                                   type = "response"))]
icu[is.na(lactate),missing_weight_glm := 1/(predict(m,newdata=icu[is.na(lactate)],
                                                type = "response"))]
icu[,missing_prob_glm := predict(m,newdata=icu,type = "response")]
MLmetrics::AUC(icu$missing_prob_glm,icu$missing_lactate1)
res_glm <- evaluate_predictions(icu$missing_prob_glm,icu$missing_lactate1,rep(1,nrow(icu)))
res_glm
fwrite(data.frame(res_glm),"results/missing_model_glm_eval.csv",row.names = TRUE)

#fwrite(icu,"data/eicu_1patient_weight.csv")
# distribution of prob missing
icu[,missing_prob_glm_grp := cut(missing_prob_glm,seq(0,1,0.1),ordered=TRUE)]
tab <- icu[,.(missing_lactate1=mean(missing_lactate1),missing_prob_glm=mean(missing_prob_glm),.N),by=missing_prob_glm_grp]
tab <- tab[order(missing_prob_glm_grp)]
tab[,perfect := seq(0.05,0.95,0.1)]
tab[,se := sqrt(missing_lactate1*(1-missing_lactate1)/N)]
tab[,labels:=c("0.0-0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5",
               "0.5-0.6","0.6-0.7","0.7-0.8","0.8-0.9","0.9-1.0")]
p0 <- ggplot(icu,aes(x=missing_prob_glm,fill=missing_lactate)) + 
  geom_density(alpha=0.5) +
  theme_bw(base_size=14)+
  coord_cartesian(xlim=c(0,1)) +
  scale_fill_discrete(name="Blood lactate",labels=c("Measured","Missing")) +
  labs(x="Probability of missing lactate",title="(A)") +
  theme(legend.position = "bottom")
p0
p1 <- ggplot(tab,aes(x=missing_prob_glm,y=missing_lactate1)) +
  geom_abline(slope=1,intercept=0,linetype=2,col="grey50") +
  geom_point() +
  geom_pointrange(aes(ymax=missing_lactate1+2*se,ymin=missing_lactate1-2*se),size=0.1) +
  scale_y_continuous(breaks=seq(0,1,0.1),limit=c(0,1)) +
  scale_x_continuous(breaks=seq(0,1,0.1),limit=c(0,1)) +
  labs(x = "Model estimated probability\nof missingness",y="Observed proportion\nof missingness",
       title="(B)") +
  theme_bw(base_size = 14)+
  theme(panel.grid.minor = element_blank())
p1
ggsave("graphs/missing_model_glm_calibration.png")
p2 <- ggplot(icu[!is.na(lactate)],aes(x=missing_weight_glm)) + 
  stat_ecdf(geom = "step") +
  coord_cartesian(xlim=c(1,30)) +
  theme_bw(base_size=14) +
  scale_x_log10(n.breaks=8) +
  labs(x="Missingness weight",y="Cumulative probability",title="(D)")
p2
p3 <- ggplot(icu[!is.na(lactate)],aes(y=missing_prob_glm,x=apacheivascore,col=factor(sepsis_apache))) + 
  geom_point(alpha=0.05) +
  geom_smooth(se=FALSE,size=1.5) +
  geom_smooth(aes(group=factor(sepsis_apache),col=NULL),size=0.1,se=FALSE,col="black") +
  scale_color_discrete(name="Sepsis on\nadmission",labels=c("No","Yes")) +
  coord_cartesian(ylim=c(0,1)) +
  theme_bw(base_size=14)+
  labs(y="Probability of missing lactate",x="APACHE score",title="(C)")+
  theme(legend.position = "bottom")
p3
gridExtra::grid.arrange(p0,p1,p3,p2,ncol=2)
pp <- gridExtra::arrangeGrob(p0,p1,p3,p2,ncol=2)
ggsave(plot = pp,filename = "graphs/missing_glm_model_summary.png",height = 7,width = 9)

#### XGBoost missingness model #####
icu[,apache_dx1 := apache_dx]
to_other <- icu[,.N,by=apache_dx1][N < 100]$apache_dx1
icu[apache_dx1 %in% to_other,apache_dx1 := "other"]
icu[,apache_dx1 := make.names(apache_dx1)]
vars <- c("potassium","sodium","chloride",
          "creatinine","bun","bicarbonate","calcium","gcs_verbal",
          "gcs_eyes","diabetes","bmi","gcs",
          "oobventday1","oobintubday1","apacheivascore","operative","sepsis_apache",
          "age","gender","unittype_f",
          "apache_dx1")
icu1 <- icu[,..vars]
icu1 <- dummy_cols(icu1,select_columns = c("apache_dx1","unittype_f"),remove_selected_columns = TRUE)
icu1 <- as.matrix(icu1)
# missingness model
set.seed(SEED)
xgb_mods <- xgb.cv(params = list(objective = "binary:logistic"),
       data=as.matrix(icu1),
       label=icu$missing_lactate1,
       prediction=TRUE,
       callbacks = list(cb.cv.predict(save_models = TRUE)),
       nfold=10,
       nrounds=2000,
       early_stopping_rounds=20,
       metrics=list("logloss","auc"))
saveRDS(xgb_mods,file = "models/missingness_model_xgb.rds")
icu[,missing_prob_xgb := xgb_mods$pred]
icu[,missing_weight_xgb := 1/(1-missing_prob_xgb)]
MLmetrics::AUC(icu$missing_prob_xgb,icu$missing_lactate1)
res_xgb <- evaluate_predictions(icu$missing_prob_xgb,icu$missing_lactate1,rep(1,nrow(icu1)))
res_xgb
fwrite(data.frame(res_xgb),"results/missing_model_xgb_eval.csv",row.names = TRUE)

#fwrite(icu,"data/eicu_1patient_weight.csv")
# distribution of prob missing
icu[,missing_prob_xgb_grp := cut(missing_prob_xgb,seq(0,1,0.1),ordered=TRUE)]
tab <- icu[,.(missing_lactate1=mean(missing_lactate1),missing_prob_xgb=mean(missing_prob_xgb),.N),by=missing_prob_xgb_grp]
tab <- tab[order(missing_prob_xgb_grp)]
tab[,perfect := seq(0.05,0.95,0.1)]
tab[,se := sqrt(missing_lactate1*(1-missing_lactate1)/N)]
tab[,labels:=c("0.0-0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5",
               "0.5-0.6","0.6-0.7","0.7-0.8","0.8-0.9","0.9-1.0")]
p0 <- ggplot(icu,aes(x=missing_prob_xgb,fill=missing_lactate)) + 
  geom_density(alpha=0.5) +
  theme_bw(base_size=14)+
  coord_cartesian(xlim=c(0,1)) +
  scale_fill_discrete(name="Blood lactate",labels=c("Measured","Missing")) +
  labs(x="Probability of missing lactate",title="(A)") +
  theme(legend.position = "bottom")
p0
p1 <- ggplot(tab,aes(x=missing_prob_xgb,y=missing_lactate1)) +
  geom_abline(slope=1,intercept=0,linetype=2,col="grey50") +
  geom_point() +
  geom_pointrange(aes(ymax=missing_lactate1+2*se,ymin=missing_lactate1-2*se),size=0.1) +
  scale_y_continuous(breaks=seq(0,1,0.1),limit=c(0,1)) +
  scale_x_continuous(breaks=seq(0,1,0.1),limit=c(0,1)) +
  labs(x = "Model estimated probability\nof missingness",y="Observed proportion\nof missingness",
       title="(B)") +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank())
p1
ggsave("graphs/missing_model_xgb_calibration.png")
p2 <- ggplot(icu[!is.na(lactate)],aes(x=missing_weight_xgb)) + 
  stat_ecdf(geom = "step") +
  coord_cartesian(xlim=c(1,30)) +
  theme_bw(base_size=14) +
  scale_x_log10(n.breaks=8) +
  labs(x="Missingness weight",y="Cumulative probability",title="(D)")
p2
p3 <- ggplot(icu[!is.na(lactate)],aes(y=missing_prob_xgb,x=apacheivascore,col=factor(sepsis_apache))) + 
  geom_point(alpha=0.05) +
  geom_smooth(se=FALSE,size=1.5) +
  geom_smooth(aes(group=factor(sepsis_apache),col=NULL),size=0.1,se=FALSE,col="black") +
  scale_color_discrete(name="Sepsis on\nadmission",labels=c("No","Yes")) +
  coord_cartesian(ylim=c(0,1)) +
  theme_bw(base_size=14)+
  labs(y="Probability of missing lactate",x="APACHE score",title="(C)")+
  theme(legend.position = "bottom")
p3
gridExtra::grid.arrange(p0,p1,p3,p2,ncol=2)
pp <- gridExtra::arrangeGrob(p0,p1,p3,p2,ncol=2)
ggsave(plot = pp,filename = "graphs/missing_xgb_model_summary.png",height = 7,width = 9)

# variable importance
imprt <- lapply(xgb_mods$models,function(x) xgb.importance(model=x))
imprt <- do.call("rbind",imprt)
imprt <- imprt[,.(Gain = mean(Gain),Gain_sd = sd(Gain)),by=Feature]
imprt[,Feature1 := as.numeric((gsub(pattern = "f",replacement = "",x = Feature)))]
imprt[,name := colnames(icu1)[Feature1+1]]
imprt[order(-Gain)][,.(name,Gain,Gain_sd)][1:20]
fwrite(imprt,"results/missing_model_xgb_import.csv")

## FURTHER ASSESS WEIGHTS ======================================================
# This section includes:
# the distribution of the weights
# an analysis of who has large weights
# construction of table1

# distribution of the weights --------------------------------------------------
quantile(icu$missing_weight_glm,na.rm=TRUE,probs = c(seq(0,0.8,0.1),seq(0.9,1.0,0.01)))
icu[,.N,by=.(missing_weight_glm > 10)]
ggplot(icu[!is.na(lactate)],aes(x=missing_weight_glm)) + 
  stat_ecdf(geom = "step") +
  coord_cartesian(xlim=c(1,20)) +
  scale_y_continuous(n.breaks=6) +
  theme_bw(base_size=14) +
  labs(x="Missingness weight",y="Cumulative proportion")
ggsave(filename = "graphs/ecdf_glm_weights.png")
ggsave(filename = "graphs/ecdf_glm_weights_log10.png")
ggplot(icu,aes(x=missing_weight_glm)) + 
  geom_density(fill="lightblue") +
  scale_x_continuous(limits = c(1,15)) +
  theme_bw(base_size=14)+
  labs(x="Missingness weight")
ggsave(filename = "graphs/dist_glm_weights.png")
ggplot(icu,aes(x=missing_weight_glm)) + 
  geom_density(fill="lightblue") +
  theme_bw(base_size=14)+
  scale_x_log10(n.breaks=8) +
  labs(x="Missingness weight")
ggsave(filename = "graphs/dist_glm_weights_log10.png")


## why large weights? ----------------------------------------------------------
icu[!is.na(lactate) & missing_weight_glm > 10,apacheivascore]
icu[!is.na(lactate) & missing_weight_glm > 10,age]
icu[!is.na(lactate) & missing_weight_glm > 10,apache_dx]
# apache score
ggplot(icu[!is.na(missing_weight_glm) & !is.na(lactate)],
       aes(x=apacheivascore,fill=missing_weight_glm > 10)) +
  geom_density(alpha=0.5) +
  scale_fill_discrete(name="Weight > 10") +
  theme_bw(base_size=14) +
  labs(x = "Apache IV score")
ggsave(filename = "graphs/dist_glm_large_weights.png")
# apache score and diagnosis
large_weight_diag <- icu[!is.na(missing_weight_glm) & !is.na(lactate),
                         .(apacheivascore=mean(apacheivascore),.N),by=.(apache_dx,
                         large_weight=missing_weight_glm > 10)]
large_weight_diag_f <- large_weight_diag[large_weight == FALSE]
large_weight_diag_t <- large_weight_diag[large_weight == TRUE]
large_weight_diag <- merge(large_weight_diag_f[,.(apache_dx,apacheivascore=apacheivascore,N=N)],
                           large_weight_diag_t[,.(apache_dx,apacheivascore_large=apacheivascore,N_large=N)],
                           by="apache_dx")
large_weight_diag[order(-N_large)]
large_weight_diag[,diff := apacheivascore_large - apacheivascore]
large_weight_diag[order(diff)]
# operative
icu[!is.na(lactate) & missing_weight_glm > 10,.N,by=operative][order(operative)]
icu[!is.na(lactate),.N,by=operative][order(operative)]
# unit
tab1 <- icu[!is.na(lactate) & missing_weight_glm > 10,.N,by=unittype]
tab1[,p:= N/sum(N)]
tab1[order(unittype)]
tab2 <- icu[!is.na(lactate),.N,by=unittype][order(unittype)]
tab2[,p:= N/sum(N)]
tab2[order(unittype)]
