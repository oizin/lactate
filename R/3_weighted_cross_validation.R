## IMPORT PACKAGES  ============================================================
# packages
library(ggplot2)
library(data.table)
library(mgcv)
library(xgboost)
# library(mgcViz)
# library(gratia)
# library(MLmetrics)
source("R/functions/cohort.R")
source("R/functions/cross_validation.R")

SEED <- 111

missing_model_glm <- readRDS("models/missingness_model_glm.rds")
missing_model_xgb <- readRDS("models/missingness_model_xgb.rds")

## IMPORT EICU =================================================================
## ONE MEASURE PER PATIENT (refine)
icu <- fread("data/cohort.csv")
icu[,.N]

# choice of lactate and glucose vars
icu[,glucose := glucose_lact]
icu[,lactate := lactate_gluc]

# vars
icu[,insulin_f := factor(insulin_f)]
icu[,missing_lactate1 := fifelse(is.na(lactate),1,0,0)]

## MODELLING DATASET(S) ========================================================

# descriptives
icu[,.N]
icu[,.N,by=diabetes]

# model dataset
model_data <- icu
model_data[lactate <= 0.1,lactate := median(lactate,na.rm=TRUE)]
model_data[glucose <= 0.1,glucose := median(glucose,na.rm=TRUE)]
model_data[is.na(hosp_mort),hosp_mort := 0]

# new vars
model_data[,insulin_diabetes_f := interaction(insulin_f,diabetes_f)]
lactate_bins <- c(min(model_data$lactate,na.rm=TRUE)-0.1,1,1.3,1.7,2.3,max(model_data$lactate,na.rm=TRUE)+0.1)
glucose_bins <- c(min(model_data$glucose,na.rm=TRUE)-0.1,7*18,7.6*18,8.2*18,9.0*18,max(model_data$glucose,na.rm=TRUE)+0.1)
model_data[,lactate_binned := cut(lactate,lactate_bins,right=TRUE)]
model_data[,glucose_binned := cut(glucose,glucose_bins,right=TRUE)]

# descriptives
model_data[,.N]

# weights
# glm
model_data[,prob_missing_glm := predict(missing_model_glm,newdata=model_data,type="response")]
model_data[,weight_glm := 1/(1-prob_missing_glm) * mean(model_data$missing_lactate1)]
#model_data[weight_glm > 8.5,weight_glm := 8.5]
model_data[is.na(weight_glm),weight_glm := median(model_data$weight_glm,na.rm=TRUE)]
# xgb
model_data[,prob_missing_xgb := missing_model_xgb$pred]
model_data[,weight_xgb := 1/(1-prob_missing_xgb) * mean(model_data$missing_lactate1)]
model_data[weight_xgb > 10,weight_xgb := 10]
model_data[is.na(weight_xgb),weight_xgb := median(model_data$weight_xgb,na.rm=TRUE)]

# no missings
model_data <- model_data[!is.na(lactate) & !is.na(glucose)]
model_data[,.N]

## CROSS VALIDATION FUNCTIONS ==================================================

# see R/_cross_validation.R

## MODELS ======================================================================

lbls <- c("hosp_mort ~ glucose_binned" = "LR: glucose",
          "hosp_mort ~ s(log(glucose))" = "GAM: glucose",
          "hosp_mort ~ lactate_binned" = "LR: lactate",
          "hosp_mort ~ s(log(lactate))" = "GAM: lactate",
          "hosp_mort ~ glucose_binned + lactate_binned" = 
            "LR: glucose + lactate",
          "hosp_mort ~ s(log(glucose)) + s(log(lactate))" = 
            "GAM: glucose + lactate",
          # "hosp_mort ~ s(log(glucose), by = insulin_f) + s(log(lactate))"=
          #   "GAM: (glucose | insulin) + lactate",
          # "hosp_mort ~ glucose_binned * insulin_f + lactate_binned"=
          #   "LR: (glucose : insulin) + lactate",
          "hosp_mort ~ s(log(glucose), log(lactate))"="GAM: glucose * lactate",
          "hosp_mort ~ glucose_binned * lactate_binned"="LR: glucose * lactate",
          "xgb ~ 1" = "XGBoost: glucose, lactate")
models <- unname(lapply(names(lbls),formula))
print(models)
model_names <- as_labeller(lbls)
lbsl_order <- unlist(model_names(names(lbls)))

## PATIENT GROUPED CROSS VALIDATION (GLM WEIGHT) ===============================

### ALL PATIENT #######
## split data
NFOLDS <- 10
# reproducibility
set.seed(SEED)
ids <- unique(model_data$patienthealthsystemstayid)
cv_folds <- kFold(NFOLDS,"patienthealthsystemstayid",ids,1/model_data$weight_glm)
cv_folds$splits
lapply(cv_folds$weights,sum)
lapply(cv_folds$splits,length)

## train/test
cv_results <- train_test(cv_splits = cv_folds,data = model_data,
                         models = models,weight_var = "weight_glm")
cv_results[,model := as.character(rep(models,NFOLDS))]
cv_results_long <- melt(cv_results,
                        id.vars = c("fold","model"),
                        variable.name = "metric")
head(model_data)

## tables/graphs - raw
## tables
tab_all <- dcast(cv_results_long[,.(value=mean(value)),by=.(model,metric)],model ~ metric, value.var = "value")
tab_all[,group := "overall"]

## plot CV results
cv_results_long[,model_disp := factor(unlist(model_names(model)),levels = lbsl_order,ordered=TRUE)]
p1 <- ggplot(cv_results_long[metric %in% c("wAUC","wBrier","wPRAUC","wLogLoss")],
       aes(x = metric,y = value,fill=model_disp)) +
  #geom_point(position = position_dodge(0.5)) +
  geom_boxplot() +
  scale_fill_discrete() +
  labs(y = "",x = "") +
  scale_x_discrete(breaks=NULL)+
  coord_flip() +
  facet_wrap(~metric,scale="free",
             labeller = as_labeller(c("wAUC"="AUC-ROC","wPRAUC"="AUC-PR",
                                      "wBrier"="Brier","wLogLoss"="Log loss"))) +
  theme_bw(base_size=14) +
  theme(legend.title = element_blank(),legend.position="right")
p1

## APACHE
# overall
sel <- model_data$predictedhospitalmortality != -1
evaluate_predictions(model_data[sel]$predictedhospitalmortality,
                     model_data[sel]$hosp_mort,
                     model_data[sel]$weight_glm)


### DIABETICS #######
model_data_dm <- model_data[diabetes == 1]
ids <- unique(model_data_dm$patienthealthsystemstayid)
set.seed(SEED)
cv_folds <- kFold(NFOLDS,"patienthealthsystemstayid",ids,1/model_data_dm$weight_glm)
cv_folds$splits
lapply(cv_folds$weights,sum)
lapply(cv_folds$splits,length)

## train/test
cv_results_dm <- train_test(cv_splits = cv_folds,data = model_data_dm,
                         models = models,weight_var = "weight_glm")
cv_results_dm[,model := as.character(rep(models,NFOLDS))]
cv_results_dm_long <- melt(cv_results_dm,
                        id.vars = c("fold","model"),
                        variable.name = "metric")
head(model_data_dm)

## tables/graphs - raw
## tables
tab_dm <- dcast(cv_results_dm_long[,.(value=mean(value)),by=.(model,metric)],model ~ metric, value.var = "value")
tab_dm[,group := "dm"]

## plot CV results
cv_results_dm_long[,model_disp := factor(unlist(model_names(model)),levels = lbsl_order,ordered=TRUE)]
p2 <- ggplot(cv_results_dm_long[metric %in% c("wAUC","wBrier","wPRAUC","wLogLoss")],
       aes(x = metric,y = value,fill=model_disp)) +
  #geom_point(position = position_dodge(0.5)) +
  geom_boxplot() +
  scale_fill_discrete() +
  labs(y = "",x = "") +
  scale_x_discrete(breaks=NULL)+
  coord_flip() +
  facet_wrap(~metric,scale="free",
             labeller = as_labeller(c("wAUC"="AUC-ROC","wPRAUC"="AUC-PR",
                                      "wBrier"="Brier","wLogLoss"="Log loss"))) +
  theme_bw(base_size=14) +
  theme(legend.title = element_blank(),legend.position="right")
p2

## APACHE
# diabetic
sel <- model_data_dm$predictedhospitalmortality != -1 & model_data_dm$diabetes == 1
evaluate_predictions(model_data_dm[sel]$predictedhospitalmortality,
                     model_data_dm[sel]$hosp_mort,
                     model_data_dm[sel]$weight_glm)

### NON-DIABETICS #######
model_data_nd <- model_data[diabetes == 0]
ids <- unique(model_data_nd$patienthealthsystemstayid)
set.seed(SEED)
cv_folds <- kFold(NFOLDS,"patienthealthsystemstayid",ids,1/model_data_nd$weight_glm)
cv_folds$splits
lapply(cv_folds$weights,sum)
lapply(cv_folds$splits,length)

## train/test
cv_results_nd <- train_test(cv_splits = cv_folds,data = model_data_nd,
                         models = models,weight_var = "weight_glm")
cv_results_nd[,model := as.character(rep(models,NFOLDS))]
cv_results_nd_long <- melt(cv_results_nd,
                        id.vars = c("fold","model"),
                        variable.name = "metric")
head(model_data_nd)

## tables/graphs - raw
## tables
tab_nd <- dcast(cv_results_nd_long[,.(value=mean(value)),by=.(model,metric)],model ~ metric, value.var = "value")
tab_nd[,group := "nd"]

## plot CV results
cv_results_nd_long[,model_disp := factor(unlist(model_names(model)),levels = lbsl_order,ordered=TRUE)]
p3 <- ggplot(cv_results_nd_long[metric %in% c("wAUC","wBrier","wPRAUC","wLogLoss")],
       aes(x = metric,y = value,fill=model_disp)) +
  #geom_point(position = position_dodge(0.5)) +
  geom_boxplot() +
  scale_fill_discrete() +
  labs(y = "",x = "") +
  scale_x_discrete(breaks=NULL)+
  coord_flip() +
  facet_wrap(~metric,scale="free",
             labeller = as_labeller(c("wAUC"="AUC-ROC","wPRAUC"="AUC-PR",
                                      "wBrier"="Brier","wLogLoss"="Log loss"))) +
  theme_bw(base_size=14) +
  theme(legend.title = element_blank(),legend.position="right")
p3
gridExtra::grid.arrange(p1,p2,p3,ncol=1)
pp <- gridExtra::arrangeGrob(p1 + labs(title="(A)"),
                             p2 + labs(title = "(B)"),
                             p3 + labs(title = "(C)"),ncol=1)
ggsave(plot = pp,filename = "graphs/cv_results_glm_weights.png",height = 16,width = 9)
ggsave(plot = pp,filename = "graphs/cv_results_glm_weights.eps",height = 16,width = 9,dpi=600)

## APACHE
# non-diabetic
sel <- model_data_nd$predictedhospitalmortality != -1 & model_data_nd$diabetes == 0
evaluate_predictions(model_data_nd[sel]$predictedhospitalmortality,
                     model_data_nd[sel]$hosp_mort,
                     model_data_nd[sel]$weight_glm)

### SAVE ####
fwrite(x = rbind(tab_all,tab_dm,tab_nd),file = "results/cross_validation_results_glm_weights.csv")

## PATIENT GROUPED CROSS VALIDATION (XGB WEIGHT) ===============================

### ALL PATIENT #######
## split data
NFOLDS <- 10
# reproducibility
set.seed(SEED)
ids <- unique(model_data$patienthealthsystemstayid)
set.seed(SEED)
cv_folds <- kFold(NFOLDS,"patienthealthsystemstayid",ids,1/model_data$weight_xgb)
cv_folds$splits
lapply(cv_folds$weights,sum)
lapply(cv_folds$splits,length)

## train/test
cv_results <- train_test(cv_splits = cv_folds,data = model_data,
                         models = models,weight_var = "weight_xgb")
cv_results[,model := as.character(rep(models,NFOLDS))]
cv_results_long <- melt(cv_results,
                        id.vars = c("fold","model"),
                        variable.name = "metric")
head(model_data)

## tables/graphs - raw
## tables
tab_all <- dcast(cv_results_long[,.(value=mean(value)),by=.(model,metric)],model ~ metric, value.var = "value")
tab_all[,group := "overall"]

## plot CV results
cv_results_long[,model_disp := factor(unlist(model_names(model)),levels = lbsl_order,ordered=TRUE)]
p1 <- ggplot(cv_results_long[metric %in% c("wAUC","wBrier","wPRAUC","wLogLoss")],
       aes(x = metric,y = value,fill=model_disp)) +
  #geom_point(position = position_dodge(0.5)) +
  geom_boxplot() +
  scale_fill_discrete() +
  labs(y = "",x = "") +
  scale_x_discrete(breaks=NULL)+
  coord_flip() +
  facet_wrap(~metric,scale="free",
             labeller = as_labeller(c("wAUC"="AUC-ROC","wPRAUC"="AUC-PR",
                                      "wBrier"="Brier","wLogLoss"="Log loss"))) +
  theme_bw(base_size=14) +
  theme(legend.title = element_blank(),legend.position="right")

### DIABETICS #######
model_data_dm <- model_data[diabetes == 1]
ids <- unique(model_data_dm$patienthealthsystemstayid)
set.seed(SEED)
cv_folds <- kFold(NFOLDS,"patienthealthsystemstayid",ids,1/model_data_dm$weight_glm)
cv_folds$splits
lapply(cv_folds$weights,sum)
lapply(cv_folds$splits,length)

## train/test
cv_results_dm <- train_test(cv_splits = cv_folds,data = model_data_dm,
                            models = models,weight_var = "weight_glm")
cv_results_dm[,model := as.character(rep(models,NFOLDS))]
cv_results_dm_long <- melt(cv_results_dm,
                           id.vars = c("fold","model"),
                           variable.name = "metric")
head(model_data_dm)

## tables/graphs - raw
## tables
tab_dm <- dcast(cv_results_dm_long[,.(value=mean(value)),by=.(model,metric)],model ~ metric, value.var = "value")
tab_dm[,group := "dm"]

## plot CV results
cv_results_dm_long[,model_disp := factor(unlist(model_names(model)),levels = lbsl_order,ordered=TRUE)]
p2 <- ggplot(cv_results_dm_long[metric %in% c("wAUC","wBrier","wPRAUC","wLogLoss")],
       aes(x = metric,y = value,fill=model_disp)) +
  #geom_point(position = position_dodge(0.5)) +
  geom_boxplot() +
  scale_fill_discrete() +
  labs(y = "",x = "") +
  scale_x_discrete(breaks=NULL)+
  coord_flip() +
  facet_wrap(~metric,scale="free",
             labeller = as_labeller(c("wAUC"="AUC-ROC","wPRAUC"="AUC-PR",
                                      "wBrier"="Brier","wLogLoss"="Log loss"))) +
  theme_bw(base_size=14) +
  theme(legend.title = element_blank(),legend.position="right")

## APACHE
# diabetic
sel <- model_data_dm$predictedhospitalmortality != -1 & model_data_dm$diabetes == 1
evaluate_predictions(model_data_dm[sel]$predictedhospitalmortality,
                     model_data_dm[sel]$hosp_mort,
                     model_data_dm[sel]$weight_glm)

### NON-DIABETICS #######
model_data_nd <- model_data[diabetes == 0]
ids <- unique(model_data_nd$patienthealthsystemstayid)
set.seed(SEED)
cv_folds <- kFold(NFOLDS,"patienthealthsystemstayid",ids,1/model_data_nd$weight_glm)
cv_folds$splits
lapply(cv_folds$weights,sum)
lapply(cv_folds$splits,length)

## train/test
cv_results_nd <- train_test(cv_splits = cv_folds,data = model_data_nd,
                            models = models,weight_var = "weight_glm")
cv_results_nd[,model := as.character(rep(models,NFOLDS))]
cv_results_nd_long <- melt(cv_results_nd,
                           id.vars = c("fold","model"),
                           variable.name = "metric")
head(model_data_nd)

## tables/graphs - raw
## tables
tab_nd <- dcast(cv_results_nd_long[,.(value=mean(value)),by=.(model,metric)],model ~ metric, value.var = "value")
tab_nd[,group := "nd"]

## plot CV results
cv_results_nd_long[,model_disp := factor(unlist(model_names(model)),levels = lbsl_order,ordered=TRUE)]
p3 <- ggplot(cv_results_nd_long[metric %in% c("wAUC","wBrier","wPRAUC","wLogLoss")],
       aes(x = metric,y = value,fill=model_disp)) +
  #geom_point(position = position_dodge(0.5)) +
  geom_boxplot() +
  scale_fill_discrete() +
  labs(y = "",x = "") +
  scale_x_discrete(breaks=NULL)+
  coord_flip() +
  facet_wrap(~metric,scale="free",
             labeller = as_labeller(c("wAUC"="AUC-ROC","wPRAUC"="AUC-PR",
                                      "wBrier"="Brier","wLogLoss"="Log loss"))) +
  theme_bw(base_size=14) +
  theme(legend.title = element_blank(),legend.position="right")

## APACHE
# non-diabetic
sel <- model_data_nd$predictedhospitalmortality != -1 & model_data_nd$diabetes == 0
evaluate_predictions(model_data_nd[sel]$predictedhospitalmortality,
                     model_data_nd[sel]$hosp_mort,
                     model_data_nd[sel]$weight_glm)

gridExtra::grid.arrange(p1,p2,p3,ncol=1)
pp <- gridExtra::arrangeGrob(p1 + labs(title="(A)"),
                             p2 + labs(title = "(B)"),
                             p3 + labs(title = "(C)"),ncol=1)
ggsave(plot = pp,filename = "graphs/cv_results_xgb_weights.png",height = 16,width = 9)


### SAVE ####
fwrite(x = rbind(tab_all,tab_dm,tab_nd),file = "results/cross_validation_results_xgb_weights.csv")
