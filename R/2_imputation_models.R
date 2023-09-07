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

icu[lactate <= 0.1,lactate := median(lactate,na.rm=TRUE)]

## MISSING WEIGHTS =============================================================
#### GLM imputation model ####
m <- gam(lactate~ 
           s(apacheivascore) +
           sepsis_apache +
           operative + 
           oobventday1 +
           oobintubday1 +
           organ_system + 
           age +
           unittype_f,
         data=icu)
summary(m)
m_res <- broom::tidy(m,conf.int = FALSE, conf.level = 0.95,parametric=TRUE)
setDT(m_res)
m_res[,lower := estimate - 1.96*std.error]
m_res[,higher := estimate + 1.96*std.error]
fwrite(m_res,file = "results/imputation_model_glm.csv")
m_res1 <- broom::tidy(m,conf.int = FALSE, conf.level = 0.95,parametric=FALSE)
fwrite(m_res1,file = "results/imputation_model_glm_spline.csv")
# save the model
saveRDS(m,file = "models/imputation_model_glm.rds")
icu[,pred_glm := predict(m,newdata = icu)]
res_glm <- icu[,.(rmse=sqrt(mean((pred_glm - lactate)^2,na.rm=TRUE)))]
fwrite(data.frame(res_glm),"results/imputation_model_glm_eval.csv",row.names = TRUE)

#### XGBoost imputation model #####
icu[,apache_dx1 := apache_dx]
to_other <- icu[,.N,by=apache_dx1][N < 100]$apache_dx1
icu[apache_dx1 %in% to_other,apache_dx1 := "other"]
icu[,apache_dx1 := make.names(apache_dx1)]
vars <- c("potassium","sodium","chloride","glucose_mean",
          "creatinine","bun","bicarbonate","calcium","gcs_verbal",
          "gcs_eyes","diabetes","bmi","gcs",
          "oobventday1","oobintubday1","apacheivascore","operative","sepsis_apache",
          "age","gender","unittype_f",
          "apache_dx1")
icu1 <- icu[,..vars]
icu1 <- dummy_cols(icu1,select_columns = c("apache_dx1","unittype_f"),remove_selected_columns = TRUE)
icu1 <- as.matrix(icu1)

# imputation model
missing_outcome <- icu[,!is.na(lactate)]
set.seed(SEED)
xgb_mods1 <- xgb.cv(params = list(objective = "reg:squarederror"),
                    data=as.matrix(icu1[missing_outcome,]),
                    label=icu[missing_outcome]$lactate,
                    prediction=TRUE,
                    nfold=10,
                    nround=2000,
                    callbacks = list(cb.cv.predict(save_models = TRUE)),
                    early_stopping_rounds=20,
                    metrics=list("rmse","mae"))
imprt <- lapply(xgb_mods1$models,function(x) xgb.importance(model=x))
imprt <- do.call("rbind",imprt)
imprt <- imprt[,.(Gain = mean(Gain),Gain_sd = sd(Gain)),by=Feature]
imprt[,Feature1 := as.numeric((gsub(pattern = "f",replacement = "",x = Feature)))]
imprt[,name := colnames(icu1)[Feature1+1]]
imprt[order(-Gain)][,.(name,Gain,Gain_sd)][1:20]
fwrite(imprt,"results/impute_model_xgb_import.csv")
# save the model
saveRDS(xgb_mods1,file = "models/imputation_model_xgb.rds")
icu[missing_outcome,pred_xgb := xgb_mods1$pred]
res_xgb <- icu[,.(rmse=sqrt(mean((pred_xgb - lactate)^2,na.rm=TRUE)))]
fwrite(data.frame(res_xgb),"results/imputation_model_xgb_eval.csv",row.names = TRUE)
