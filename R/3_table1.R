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

## IMPORT EICU =================================================================
# data
icu <- fread("data/cohort.csv")

## VARS
# choice of lactate and glucose vars
icu[,glucose := glucose_lact]
icu[,lactate := lactate_gluc]

# factors
icu[,unittype_f := factor(unittype)]
icu[,hosp_mort_f := factor(hosp_mort,levels=c(0,1),labels=c("survived","died"))]
icu[,diabetes_f := factor(diabetes,levels=c(0,1),labels=c("non-diabetic","diabetic"))]
icu[,insulin_f := factor(insulin,levels=c(0,1),labels=c("non-treated","treated"))]
icu[,missing_lactate := fifelse(is.na(lactate),"missing","measured")]
icu[,missing_lactate1 := fifelse(is.na(lactate),1,0,0)]
icu[,missing_outcome := fifelse(is.na(hosp_mort),"missing","measured")]
icu[is.infinite(lactate_max),lactate_max := NA]
icu[,gcs := gcs_verbal + gcs_eyes + gcs_motor]
icu[,bmi := admissionweight/((admissionheight/100)^2)]
icu[bmi < 10 | bmi > 60,bmi := NA]

## models
# glm
m_glm <- readRDS(file = "models/missingness_model_glm.rds")
icu[,missing_prob_glm := predict(m_glm,newdata = icu,type="response")]
evaluate_predictions(icu$missing_prob_glm,icu$missing_lactate1,rep(1,nrow(icu)))
# xgb
m_xgb <- readRDS(file = "models/missingness_model_xgb.rds")
icu[,missing_prob_xgb := m_xgb$pred]
evaluate_predictions(icu$missing_prob_xgb,icu$missing_lactate1,rep(1,nrow(icu)))

# probs
icu[is.na(lactate),sample_prob_glm := missing_prob_glm]
icu[!is.na(lactate),sample_prob_glm := 1-missing_prob_glm]
icu[is.na(lactate),sample_prob_xgb := missing_prob_xgb]
icu[!is.na(lactate),sample_prob_xgb := 1-missing_prob_xgb]

## table 1 (generic) -----------------------------------------------------------

# lactate numbers:
l_vars <- c("lactate_gluc","lactate_mean","lactate_max","lactate0")
CreateTableOne(vars=l_vars,
               data=icu,test = FALSE,addOverall = TRUE)
icu_measured_w <- svydesign(~1,
                            probs=~sample_prob_glm,
                            data=icu[!is.na(sample_prob_glm) & !is.na(lactate)])
svyCreateTableOne(vars=l_vars,data=icu_measured_w)
# number of patients
length(unique(icu$patienthealthsystemstayid))
length(unique(icu[!is.na(lactate)]$patienthealthsystemstayid))
length(unique(icu[is.na(lactate)]$patienthealthsystemstayid))
# number hospitals
length(unique(icu$hospitalid))

# missing glucose
icu[is.na(lactate),.N,by=.(is.na(lactate))]
length(icu$patienthealthsystemstayid)
length(unique(icu[!is.na(lactate)]$patienthealthsystemstayid))
length(unique(icu[is.na(lactate)]$patienthealthsystemstayid))

## compare weighted populations - table 1 (GLM ) -------------------------------
# table 1 missing vs. non missing
icu_all_w <- svydesign(~1,
                       probs=~sample_prob_glm,
                       data=icu[!is.na(sample_prob_glm)])
t1_vars <- c("apacheivascore",
             "gender","bmi","gcs",
             "glucose_lact","glucose0","glucose_mean","glucose_max","glucose80","glucose50",
             "insulin_f","oobventday1","oobintubday1",
             "hosp_mort_f","icu_los_hours","icu_mort","hosp_mort",
             "unittype_f","ethnicity",
             "sepsis_apache",
             "creatinine","bilirubin",
             "operative",
             "age","diabetes_f")
t1_fvars <- c("operative","sepsis_apache",
              "icu_mort","hosp_mort")
svyCreateTableOne(vars=t1_vars,factorVars = t1_fvars,
                  strata=c("missing_lactate"),
                  data=icu_all_w)
# table 1 non missing vs. rest of data
icu_measured_w <- svydesign(~1,
                            probs=~sample_prob_glm,
                            data=icu[!is.na(sample_prob_glm) & !is.na(lactate)])
tab1 <- print(svyCreateTableOne(vars=t1_vars,factorVars = t1_fvars,
                                data=icu_measured_w))
tab2 <- print(CreateTableOne(vars=t1_vars,factorVars = t1_fvars,
                             data=icu,
                             strata="missing_lactate",test = FALSE,addOverall = TRUE))
tab <- as.data.table(cbind(tab1,tab2),keep.rownames = TRUE)
names(tab) <- c("var","Weighted","Overall","Measured","Missing")
tab
fwrite(tab,"results/table1_missing_glm.csv")

## compare weighted populations - table 1 (xgboost) ----------------------------
# table 1 missing vs. non missing
icu_all_w <- svydesign(~1,
                       probs=~sample_prob_xgb,
                       data=icu[!is.na(sample_prob_xgb)])
t1_vars <- c("apacheivascore",
             "gender","bmi","gcs",
             "glucose_lact","glucose0","glucose_mean","glucose_max","glucose80","glucose50",
             "insulin_f","oobventday1","oobintubday1",
             "hosp_mort_f","icu_los_hours","icu_mort","hosp_mort",
             "unittype_f","ethnicity",
             "sepsis_apache",
             "creatinine","bilirubin",
             "operative",
             "age","diabetes_f")
t1_fvars <- c("operative","sepsis_apache",
              "icu_mort","hosp_mort")
svyCreateTableOne(vars=t1_vars,factorVars = t1_fvars,
                  strata=c("missing_lactate"),
                  data=icu_all_w)
# table 1 non missing vs. rest of data
icu_measured_w <- svydesign(~1,
                            probs=~sample_prob_xgb,
                            data=icu[!is.na(sample_prob_xgb) & !is.na(lactate)])
tab1 <- print(svyCreateTableOne(vars=t1_vars,factorVars = t1_fvars,
                                data=icu_measured_w))
tab2 <- print(CreateTableOne(vars=t1_vars,factorVars = t1_fvars,
                             data=icu,
                             strata="missing_lactate",test = FALSE,addOverall = TRUE))
tab <- as.data.table(cbind(tab1,tab2),keep.rownames = TRUE)
names(tab) <- c("var","Weighted","Overall","Measured","Missing")
tab
fwrite(tab,"results/table1_missing_xgb.csv")
