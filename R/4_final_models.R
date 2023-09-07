## IMPORT PACKAGES  ============================================================
# packages
library(ggplot2)
library(data.table)
library(mgcv)
# library(mgcViz)
library(gratia)
# library(MLmetrics)
source("R/functions/cohort.R")
source("R/functions/cross_validation.R")

missing_model_glm <- readRDS("models/missingness_model_glm.rds")
missing_model_xgb <- readRDS("models/missingness_model_xgb.rds")

## IMPORT EICU =================================================================
icu <- fread("data/cohort.csv")
icu[,.N]

# choice of lactate and glucose vars
icu[,glucose := glucose_lact]
icu[,lactate := lactate_gluc]

# vars
icu[,insulin_f := factor(insulin_f)]
icu[,missing_lactate1 := fifelse(is.na(lactate),1,0,0)]

## ADD WEIGHTS =================================================================

icu[,prob_missing_glm := predict(missing_model_glm,newdata=icu,type="response")]
icu[,missing_weight_glm := 1/(1-prob_missing_glm) * mean(icu$missing_lactate1)]
icu[missing_weight_glm > 8.5,missing_weight_glm := 8.5]
icu[,prob_missing_xgb := missing_model_xgb$pred]
icu[,missing_weight_xgb := 1/(1-prob_missing_xgb) * mean(icu$missing_lactate1)]
icu[missing_weight_xgb > 10,missing_weight_xgb := 10]

## MODELLING DATASET(S) ========================================================

# analysis dataset
model_data <- icu[!is.na(lactate) & lactate > 0 & !is.na(missing_weight_glm)]

# descriptives
model_data[,.N]
model_data[,.N,by=diabetes]

# model dataset
model_data[is.na(hosp_mort),hosp_mort := 0]

# new vars
model_data[,insulin_diabetes_f := interaction(insulin_f,diabetes_f)]
lactate_bins <- c(min(model_data$lactate,na.rm=TRUE)-0.1,1,1.3,1.7,2.3,max(model_data$lactate,na.rm=TRUE)+0.1)
glucose_bins <- c(min(model_data$glucose,na.rm=TRUE)-0.1,7*18,7.6*18,8.2*18,9.0*18,max(model_data$glucose,na.rm=TRUE)+0.1)
model_data[,lactate_binned := cut(lactate,lactate_bins,right=TRUE)]
model_data[,glucose_binned := cut(glucose,glucose_bins,right=TRUE)]

# potenital exclusions
#model_data <- model_data[insulin_pres_init_dy >0 | is.na(insulin_pres_init_dy)]
#model_data <- model_data[insulin_infuse_init_dy >0 | is.na(insulin_infuse_init_dy)]

# log vars
model_data[,ln_glucose:= log(glucose)]
model_data[,ln_lactate := log(lactate)]

## MODELS (GLM WEIGHTS) =========================================================

GLUC_BREAKS <- c(50,80,100,150,200,300,400)
GLUC_LABEL <- function(x) exp(x)
LACT_BREAKS <- c(0.5,1.0,2.0,5.0,10.0,30.0)
LACT_LABEL <- function(x) exp(x)

## Binned models
summary(gam(hosp_mort ~ glucose_binned + lactate_binned + s(hospitalid,bs="re"),
            family=binomial(),
            weights = missing_weight_glm,
            data=model_data[diabetes == 0]))
summary(mm <- gam(hosp_mort ~ glucose_binned*lactate_binned + s(hospitalid,bs="re"),
            family=binomial(),
            weights = missing_weight_glm,
            data=model_data[diabetes == 0]))
fwrite(x = broom::tidy(mm,parametric=TRUE),file = "results/interaction_binned.csv")

## (A) The impact of adjusting for lactate on blood glucose and outcome
m1 <- gam(hosp_mort ~ s(ln_glucose) + s(hospitalid,bs="re"),
          family=binomial(),
          weights = missing_weight_glm,
          data=model_data)
sm1 <- add_confint(smooth_estimates(m1))
setDT(sm1)
m2 <- gam(hosp_mort ~ s(ln_lactate) + s(ln_glucose) + s(hospitalid,bs="re"),
          family=binomial(),
          weights = missing_weight_glm,
          data=model_data)
sm2 <- add_confint(smooth_estimates(m2))
setDT(sm2)
sm2[,adjusted := "GAM: glucose +\nlactate"]
sm1[,adjusted := "GAM: glucose"]
sm1_2 <- rbind(sm1,sm2,fill=TRUE)
p1_2 <- ggplot(sm1_2[smooth == "s(ln_glucose)"],
               aes(x=ln_glucose,y=est,col=adjusted)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci,
                  ymax = upper_ci,fill=adjusted,col=NULL),alpha=0.1) +
  geom_point(data=model_data[1:5000],
    aes(x=ln_glucose,y=-1,col=NULL,fill=NULL),alpha=0.1,shape="|") +
  coord_cartesian(xlim=c(4,6),ylim=c(-1,2.5)) +
  scale_x_continuous(breaks = log(c(GLUC_BREAKS)),labels=GLUC_LABEL) +
  labs(x = "Blood glucose (mg/dL)",
       y = "Log odds\n(hospital mortality)",
       title = "(A)") +
  geom_hline(yintercept = 0,linetype=2,alpha=0.8) +
  theme_bw(base_size=14) +
  theme(legend.title=element_blank())
p1_2

## Diabetic vs. non-diabetic
m3a <- gam(hosp_mort ~ s(ln_lactate) + s(ln_glucose) + s(hospitalid,bs="re"),
           family=binomial(),
           weights = missing_weight_glm,
           data=model_data[diabetes == 1])
sm3a <- add_confint(smooth_estimates(m3a))
m4a <- gam(hosp_mort ~ s(ln_lactate) + s(ln_glucose) + s(hospitalid,bs="re"),
           family=binomial(),
           weights = missing_weight_glm,
           data=model_data[diabetes == 0])
sm4a <- add_confint(smooth_estimates(m4a))
setDT(sm3a)
setDT(sm4a)
sm3a[,group := "Diabetic"]
sm4a[,group := "Non-diabetic"]
sm3a_4a <- rbind(sm3a,sm4a,fill=TRUE)
p3a_4a <- ggplot(sm3a_4a[!is.na(ln_glucose)],
       aes(x=ln_glucose,y=est,col=group)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci,
                  ymax = upper_ci,fill=group,col=NULL),alpha=0.1) +
  geom_point(data=model_data[diabetes == 1],
             aes(x=ln_glucose,y=-1,col=NULL,fill=NULL),alpha=0.1,shape="|",
             col="#00BFC4") +
  geom_point(data=model_data[diabetes == 0],
             aes(x=ln_glucose,y=-0.9,col=NULL,fill=NULL),alpha=0.1,shape="|",
             col="#F8766D") +
  coord_cartesian(xlim=c(4,6),ylim=c(-1,2.5)) +
  scale_x_continuous(breaks = log(c(GLUC_BREAKS)),labels=GLUC_LABEL) +
  labs(x = "Blood glucose (mg/dL)",
       y = "Log odds\n(hospital mortality)",
       title = "(B)") +
  geom_hline(yintercept = 0,linetype=2,alpha=0.8) +
  theme_bw(base_size=14) +
  theme(legend.title=element_blank())

## (B) The impact of insulin on diabetics
# m3 <- gam(hosp_mort ~ s(ln_lactate) + s(ln_glucose,by=insulin_f) + s(hospitalid,bs="re"),
#           family=binomial(),
#           weights = missing_weight_glm,
#           data=model_data[diabetes == 1])
# sm3 <- add_confint(smooth_estimates(m3))
# setDT(sm3)
# sm3[,insulin := fifelse(insulin_f == "treated","Treated\nwith insulin","No treatment")]
# p3 <- ggplot(sm3[!is.na(ln_glucose)],
#        aes(x=ln_glucose,y=est,col=insulin)) +
#   geom_line() +
#   geom_ribbon(aes(ymin = lower_ci,
#                   ymax = upper_ci,fill=insulin,col=NULL),alpha=0.1) +
#   geom_point(data=model_data[diabetes == 1 & insulin == 0],
#              aes(x=ln_glucose,y=-1,col=NULL,fill=NULL),alpha=0.1,shape="|",
#              col="#00BFC4") +
#   geom_point(data=model_data[diabetes == 1 & insulin == 1],
#              aes(x=ln_glucose,y=-0.9,col=NULL,fill=NULL),alpha=0.1,shape="|",
#              col="#F8766D") +
#   coord_cartesian(xlim=c(4,6),ylim=c(-1,2.5)) +
#   scale_x_continuous(breaks = log(c(GLUC_BREAKS)),labels=GLUC_LABEL) +
#   labs(x = "Blood glucose (mg/dL)",
#        y = "Log odds\n(hospital mortality)",
#        title = "(C)") +
#   geom_hline(yintercept = 0,linetype=2,alpha=0.8) +
#   theme_bw(base_size=14) +
#   theme(legend.title=element_blank())
# 
# ## (C) The impact of insulin on non-diabetics
# m4 <- gam(hosp_mort ~ s(ln_lactate) + s(ln_glucose,by=insulin_f) + s(hospitalid,bs="re"),
#           family=binomial(),
#           weights = missing_weight_glm,
#           data=model_data[diabetes == 0])
# sm4 <- add_confint(smooth_estimates(m4))
# setDT(sm4)
# sm4[,insulin := fifelse(insulin_f == "treated","Treated\nwith insulin","No treatment")]
# p4 <- ggplot(sm4[!is.na(ln_glucose)],
#        aes(x=ln_glucose,y=est,col=insulin)) +
#   geom_line() +
#   geom_ribbon(aes(ymin = lower_ci,
#                   ymax = upper_ci,fill=insulin,col=NULL),alpha=0.1) +
#   geom_point(data=model_data[diabetes == 0 & insulin == 0],
#              aes(x=ln_glucose,y=-1,col=NULL,fill=NULL),alpha=0.1,shape="|",
#              col="#00BFC4") +
#   geom_point(data=model_data[diabetes == 0 & insulin == 1],
#              aes(x=ln_glucose,y=-0.9,col=NULL,fill=NULL),alpha=0.1,shape="|",
#              col="#F8766D") +
#   coord_cartesian(xlim=c(4,6),ylim=c(-1,2.5)) +
#   scale_x_continuous(breaks = log(c(GLUC_BREAKS)),labels=GLUC_LABEL) +
#   labs(x = "Blood glucose (mg/dL)",
#        y = "Log odds\n(hospital mortality)",
#        title = "(D)") +
#   geom_hline(yintercept = 0,linetype=2,alpha=0.8) +
#   theme_bw(base_size=14) +
#   theme(legend.title=element_blank())
# interaction
m5 <- gam(hosp_mort ~ s(ln_lactate,ln_glucose),
         family=binomial(),
         weights = missing_weight_glm,
         data=model_data[diabetes == 0])
sm5 <- add_confint(smooth_estimates(m5))
setDT(sm5)
p5 <- ggplot(sm5[ln_glucose > 4 & ln_glucose < 6 & 
                   ln_lactate >= log(0.5) & ln_lactate < log(30)]
             ,aes(x=exp(ln_glucose),y=exp(ln_lactate))) +
  scale_x_log10(expand=c(0, 0)) +
  scale_y_log10(expand=c(0, 0)) +
  geom_raster(aes(fill=est),interpolate=TRUE) +
  scale_fill_viridis_b(n.breaks=10,name="Log odds\n(hosp. mortality)") +
  labs(x = "Blood glucose (mg/dL)",
       y = "Blood lactate (mmol/L)",
       title = "(C)") +
  #coord_cartesian(ylim=c(log(0.5),log(30)),xlim=c(4,6)) +
  theme_bw(base_size=14) 

m6 <- gam(hosp_mort ~ s(ln_lactate,ln_glucose),
          family=binomial(),
          weights = missing_weight_glm,
          data=model_data[diabetes == 1])
sm6 <- add_confint(smooth_estimates(m6))
setDT(sm6)
p6 <- ggplot(sm6[ln_glucose > 4 & ln_glucose < 6 & 
                   ln_lactate >= log(0.5) & ln_lactate < log(30)]
             ,aes(x=exp(ln_glucose),y=exp(ln_lactate))) +
  scale_x_log10(expand=c(0, 0)) +
  scale_y_log10(expand=c(0, 0)) +
  geom_raster(aes(fill=est),interpolate=TRUE) +
  scale_fill_viridis_b(n.breaks=10,name="Log odds\n(hosp. mortality)") +
  labs(x = "Blood glucose (mg/dL)",
       y = "Blood lactate (mmol/L)",
       title = "(D)") +
  #coord_cartesian(ylim=c(log(0.5),log(30)),xlim=c(4,6)) +
  theme_bw(base_size=14) 

gridExtra::grid.arrange(p1_2,p3a_4a,p5,p6,ncol=1)
p_all <- gridExtra::arrangeGrob(p1_2,p3a_4a,p5,p6,ncol=1)
ggsave(plot = p_all,filename = "graphs/gam_effects_glm_weights.png",height=12,width=8) 
ggsave(plot = p_all,filename = "graphs/gam_effects_glm_weights.eps",height=12,width=8,dpi=600) 

## MODELS (XGB WEIGHTS) =========================================================

## Binned models
summary(gam(hosp_mort ~ glucose_binned + lactate_binned + s(hospitalid,bs="re"),
            family=binomial(),
            weights = missing_weight_xgb,
            data=model_data))

## (A) The impact of adjusting for lactate on blood glucose and outcome
m1 <- gam(hosp_mort ~ s(ln_glucose) + s(hospitalid,bs="re"),
          family=binomial(),
          weights = missing_weight_xgb,
          data=model_data)
sm1 <- add_confint(smooth_estimates(m1))
setDT(sm1)
m2 <- gam(hosp_mort ~ s(ln_lactate) + s(ln_glucose) + s(hospitalid,bs="re"),
          family=binomial(),
          weights = missing_weight_xgb,
          data=model_data)
sm2 <- add_confint(smooth_estimates(m2))
setDT(sm2)
sm2[,adjusted := "GAM: glucose +\nlactate"]
sm1[,adjusted := "GAM: glucose"]
sm1_2 <- rbind(sm1,sm2,fill=TRUE)
p1_2 <- ggplot(sm1_2[smooth == "s(ln_glucose)"],
               aes(x=ln_glucose,y=est,col=adjusted)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci,
                  ymax = upper_ci,fill=adjusted,col=NULL),alpha=0.1) +
  geom_point(data=model_data[1:5000],
             aes(x=ln_glucose,y=-1,col=NULL,fill=NULL),alpha=0.1,shape="|") +
  coord_cartesian(xlim=c(4,6),ylim=c(-1,2.5)) +
  scale_x_continuous(breaks = log(c(GLUC_BREAKS)),labels=GLUC_LABEL) +
  labs(x = "Blood glucose (mg/dL)",
       y = "Log odds\n(hospital mortality)",
       title = "(A)") +
  geom_hline(yintercept = 0,linetype=2,alpha=0.8) +
  theme_bw(base_size=14) +
  theme(legend.title=element_blank())
p1_2

## Diabetic vs. non-diabetic
m3a <- gam(hosp_mort ~ s(ln_lactate) + s(ln_glucose) + s(hospitalid,bs="re"),
           family=binomial(),
           weights = missing_weight_xgb,
           data=model_data[diabetes == 1])
sm3a <- add_confint(smooth_estimates(m3a))
m4a <- gam(hosp_mort ~ s(ln_lactate) + s(ln_glucose) + s(hospitalid,bs="re"),
           family=binomial(),
           weights = missing_weight_xgb,
           data=model_data[diabetes == 0])
sm4a <- add_confint(smooth_estimates(m4a))
setDT(sm3a)
setDT(sm4a)
sm3a[,group := "Diabetic"]
sm4a[,group := "Non-diabetic"]
sm3a_4a <- rbind(sm3a,sm4a,fill=TRUE)
p3a_4a <- ggplot(sm3a_4a[!is.na(ln_glucose)],
                 aes(x=ln_glucose,y=est,col=group)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci,
                  ymax = upper_ci,fill=group,col=NULL),alpha=0.1) +
  geom_point(data=model_data[diabetes == 1],
             aes(x=ln_glucose,y=-1,col=NULL,fill=NULL),alpha=0.1,shape="|",
             col="#00BFC4") +
  geom_point(data=model_data[diabetes == 0],
             aes(x=ln_glucose,y=-0.9,col=NULL,fill=NULL),alpha=0.1,shape="|",
             col="#F8766D") +
  coord_cartesian(xlim=c(4,6),ylim=c(-1,2.5)) +
  scale_x_continuous(breaks = log(c(GLUC_BREAKS)),labels=GLUC_LABEL) +
  labs(x = "Blood glucose (mg/dL)",
       y = "Log odds\n(hospital mortality)",
       title = "(B)") +
  geom_hline(yintercept = 0,linetype=2,alpha=0.8) +
  theme_bw(base_size=14) +
  theme(legend.title=element_blank())

# ## (B) The impact of insulin on diabetics
# m3 <- gam(hosp_mort ~ s(ln_lactate) + s(ln_glucose,by=insulin_f) + s(hospitalid,bs="re"),
#           family=binomial(),
#           weights = missing_weight_xgb,
#           data=model_data[diabetes == 1])
# sm3 <- add_confint(smooth_estimates(m3))
# setDT(sm3)
# sm3[,insulin := fifelse(insulin_f == "treated","Treated\nwith insulin","No treatment")]
# p3 <- ggplot(sm3[!is.na(ln_glucose)],
#              aes(x=ln_glucose,y=est,col=insulin)) +
#   geom_line() +
#   geom_ribbon(aes(ymin = lower_ci,
#                   ymax = upper_ci,fill=insulin,col=NULL),alpha=0.1) +
#   geom_point(data=model_data[diabetes == 1 & insulin == 0],
#              aes(x=ln_glucose,y=-1,col=NULL,fill=NULL),alpha=0.1,shape="|",
#              col="#00BFC4") +
#   geom_point(data=model_data[diabetes == 1 & insulin == 1],
#              aes(x=ln_glucose,y=-0.9,col=NULL,fill=NULL),alpha=0.1,shape="|",
#              col="#F8766D") +
#   coord_cartesian(xlim=c(4,6),ylim=c(-1,2.5)) +
#   scale_x_continuous(breaks = log(c(GLUC_BREAKS)),labels=GLUC_LABEL) +
#   labs(x = "Blood glucose (mg/dL)",
#        y = "Log odds\n(hospital mortality)",
#        title = "(C)") +
#   geom_hline(yintercept = 0,linetype=2,alpha=0.8) +
#   theme_bw(base_size=14) +
#   theme(legend.title=element_blank())
# 
# ## (C) The impact of insulin on non-diabetics
# m4 <- gam(hosp_mort ~ s(ln_lactate) + s(ln_glucose,by=insulin_f) + s(hospitalid,bs="re"),
#           family=binomial(),
#           weights = missing_weight_xgb,
#           data=model_data[diabetes == 0])
# sm4 <- add_confint(smooth_estimates(m4))
# setDT(sm4)
# sm4[,insulin := fifelse(insulin_f == "treated","Treated\nwith insulin","No treatment")]
# p4 <- ggplot(sm4[!is.na(ln_glucose)],
#              aes(x=ln_glucose,y=est,col=insulin)) +
#   geom_line() +
#   geom_ribbon(aes(ymin = lower_ci,
#                   ymax = upper_ci,fill=insulin,col=NULL),alpha=0.1) +
#   geom_point(data=model_data[diabetes == 0 & insulin == 0],
#              aes(x=ln_glucose,y=-1,col=NULL,fill=NULL),alpha=0.1,shape="|",
#              col="#00BFC4") +
#   geom_point(data=model_data[diabetes == 0 & insulin == 1],
#              aes(x=ln_glucose,y=-0.9,col=NULL,fill=NULL),alpha=0.1,shape="|",
#              col="#F8766D") +
#   coord_cartesian(xlim=c(4,6),ylim=c(-1,2.5)) +
#   scale_x_continuous(breaks = log(c(GLUC_BREAKS)),labels=GLUC_LABEL) +
#   labs(x = "Blood glucose (mg/dL)",
#        y = "Log odds\n(hospital mortality)",
#        title = "(D)") +
#   geom_hline(yintercept = 0,linetype=2,alpha=0.8) +
#   theme_bw(base_size=14) +
#   theme(legend.title=element_blank())

# interaction
m5 <- gam(hosp_mort ~ s(ln_lactate,ln_glucose) + s(hospitalid,bs="re"),
          family=binomial(),
          weights = missing_weight_glm,
          data=model_data[diabetes == 0])
sm5 <- add_confint(smooth_estimates(m5))
setDT(sm5)
p5 <- ggplot(sm5[ln_glucose > 4 & ln_glucose < 6 & 
                   ln_lactate >= log(0.5) & ln_lactate < log(30)]
             ,aes(x=exp(ln_glucose),y=exp(ln_lactate))) +
  scale_x_log10(expand=c(0, 0)) +
  scale_y_log10(expand=c(0, 0)) +
  geom_raster(aes(fill=est),interpolate=TRUE) +
  scale_fill_viridis_b(n.breaks=10,name="Log odds\n(hosp. mortality)") +
  labs(x = "Blood glucose (mg/dL)",
       y = "Blood lactate (mmol/L)",
       title = "(C)") +
  #coord_cartesian(ylim=c(log(0.5),log(30)),xlim=c(4,6)) +
  theme_bw(base_size=14) 

m6 <- gam(hosp_mort ~ s(ln_lactate,ln_glucose) + s(hospitalid,bs="re"),
          family=binomial(),
          weights = missing_weight_glm,
          data=model_data[diabetes == 1])
sm6 <- add_confint(smooth_estimates(m6))
setDT(sm6)
p6 <- ggplot(sm6[ln_glucose > 4 & ln_glucose < 6 & 
                   ln_lactate >= log(0.5) & ln_lactate < log(30)]
             ,aes(x=exp(ln_glucose),y=exp(ln_lactate))) +
  scale_x_log10(expand=c(0, 0)) +
  scale_y_log10(expand=c(0, 0)) +
  geom_raster(aes(fill=est),interpolate=TRUE) +
  scale_fill_viridis_b(n.breaks=10,name="Log odds\n(hosp. mortality)") +
  labs(x = "Blood glucose (mg/dL)",
       y = "Blood lactate (mmol/L)",
       title = "(D)") +
  #coord_cartesian(ylim=c(log(0.5),log(30)),xlim=c(4,6)) +
  theme_bw(base_size=14) 

gridExtra::grid.arrange(p1_2,p3a_4a,p5,p6,ncol=1)
p_all <- gridExtra::arrangeGrob(p1_2,p3a_4a,p5,p6,ncol=1)
ggsave(plot = p_all,filename = "graphs/gam_effects_xgb_weights.png",height=12,width=8) 
ggsave(plot = p_all,filename = "graphs/gam_effects_xgb_weights.eps",height=12,width=8,dpi=600) 

## MODELS XGB MODEL (MEAN WEIGHTS) =============================================

library(xgboost)
library(iml)

SEED <- 222
set.seed(SEED)
input_data <- model_data[,.(lactate,glucose,diabetes)]
valid_ids <- sample(1:nrow(input_data),floor(nrow(input_data)/10))
train_ids <- (1:nrow(input_data))[!1:nrow(input_data) %in% valid_ids]
train_data <- xgb.DMatrix(data=as.matrix(input_data[train_ids]),
                          label=model_data[train_ids]$hosp_mort)
valid_data <- xgb.DMatrix(data=as.matrix(input_data[valid_ids]),
                          label=model_data[valid_ids]$hosp_mort)

model_data[,missing_weight_all := 0.5*missing_weight_xgb + 0.5*missing_weight_glm]
xg1 <- xgb.train(data = train_data,
        watchlist = list(valid_data=valid_data),
        params=list(objective="binary:logistic",
                    eta=0.001),
        weights=model_data[train_ids]$missing_weight_all,
        early_stopping_rounds = 20,
        nrounds=5000)

g_qs <- quantile(model_data$glucose,probs=c(0.01,0.99),na.rm=TRUE)
l_qs <- quantile(model_data$lactate,probs=c(0.01,0.99),na.rm=TRUE)

df <- expand.grid(lactate = c(2,4,6),
           glucose = quantile(model_data$glucose,probs=seq(0.01,0.99,0.01),na.rm=TRUE),
           diabetes=c(0,1))
setDT(df)
df[,preds := predict(xg1,as.matrix(df))]
df1 <- df[,.(hosp_mort=mean(preds)),by=.(lactate,glucose,diabetes)]
df1[,diabetes_f := factor(diabetes,levels=c(0,1),labels=c("Non-diabetic","Diabetic"))]
ggplot(df1,aes(x=glucose,y=hosp_mort,col=factor(lactate))) +
  geom_step(alpha=0.5) +
  facet_wrap(~diabetes_f,ncol=1) +
  scale_x_log10(n.breaks=7) +
  scale_colour_manual(values=c("#9ECAE1","#6BAED6","#4292C6","#2171B5"),
                      name="Blood lactate\n(mmol/L)") +
  geom_smooth(se=FALSE,span=1) +
  theme_bw(base_size=14) +
  labs(x="Blood glucose (mg/dL)",y="Hospital mortality (%)") 
ggsave("graphs/xgb_effects_lact_all_weight.png",width=8,height=6)
ggsave("graphs/xgb_effects_lact_all_weight.eps",width=8,height=6,dpi=600)

################################################################################

