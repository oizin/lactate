library(rjson)
library(bigrquery)
library(data.table)
library(glue)
source("R/functions/cohort.R")

projectid <- fromJSON(file = ".projectid.json")$projectid
email <- fromJSON(file = ".projectid.json")$email
path <- fromJSON(file = ".projectid.json")$path

bigrquery::bq_auth(path=path,email=email)

read_sql <- function(projectid,table,page_size=NULL) {
  sql <- glue("select * from {projectid}.eicu_lactate.{table}")
  tab <- bq_project_query(projectid, sql)
  tab <- bq_table_download(tab, n_max = Inf,page_size=page_size)
  setDT(tab)
  tab
}

eicu <- read_sql(projectid,"eicu_stay_outcome")
fwrite(eicu,"data/eicu_stay_outcome.csv")
#eicu <- fread("data/eicu_stay_outcome.csv")

# cohort
eicu[!is.na(glucose),n_gluc := seq(1,.N),by=patientunitstayid]
eicu[!is.na(lactate),n_lact := seq(1,.N),by=patientunitstayid]
eicu[n_gluc == 1,glucose0 := glucose,by=patientunitstayid]
eicu[n_lact == 1,lactate0 := lactate,by=patientunitstayid]
eicu[n_gluc == 1,first_glucose_time := time]
eicu[n_lact == 1,first_lactate_time := time]

# matched glucose and lactate measurements
TIME_WINDOW_HR <- 12
eicu[!is.na(lactate),lact_time1 := min(time,na.rm=TRUE),by=patientunitstayid]
eicu[,lact_meas := fifelse(time == lact_time1,1,0)]
eicu[lact_meas == 1,lactate_gluc := lactate]
eicu[,glucose_lact := as.numeric(NA)]
eicu[,glucose_lact_win := as.numeric(NA)]
for (i in 1:TIME_WINDOW_HR) {
  print(i)
  eicu[lact_meas == 1,start_win := time - TIME_WINDOW_HR*60]
  eicu[lact_meas == 1,end_win := time + TIME_WINDOW_HR*60]
  eicu[,start_win := min(start_win,na.rm=TRUE),by=patientunitstayid]
  eicu[,end_win := max(end_win,na.rm=TRUE),by=patientunitstayid]
  eicu[is.infinite(start_win),start_win := NA]
  eicu[is.infinite(end_win),end_win := NA]
  eicu[time >= start_win & time <= end_win & is.na(glucose_lact),
       glucose_lact := mean(glucose,na.rm=TRUE),by=patientunitstayid]
  eicu[time >= start_win & time <= end_win & is.na(glucose_lact_win),
       glucose_lact_win := i,by=patientunitstayid]
}

# cohort
icu <- eicu_cohort(eicu,remove_missing=FALSE,sepsis=FALSE)
icu[is.na(glucose_lact),glucose_lact := glucose_mean]

# insulin vars
icu[insulin == 1,insulin_after_gluc1 := fifelse(first_glucose_time <= insulin_init_time,1,0)]
icu[insulin == 1,insulin_before_gluc1 := fifelse(first_glucose_time > insulin_init_time,1,0)]

## VARS
# factors
icu[,unittype_f := factor(unittype)]
icu[,hosp_mort_f := factor(hosp_mort,levels=c(0,1),labels=c("survived","died"))]
icu[,diabetes_f := factor(diabetes,levels=c(0,1),labels=c("non-diabetic","diabetic"))]
icu[,insulin_f := factor(insulin,levels=c(0,1),labels=c("non-treated","treated"))]

## vars
icu[,missing_lactate := fifelse(is.na(lactate_gluc),"missing","measured")]
icu[,missing_outcome := fifelse(is.na(hosp_mort),"missing","measured")]
icu[is.infinite(lactate_max),lactate_max := NA]
icu[,gcs := gcs_verbal + gcs_eyes + gcs_motor]
icu[,bmi := admissionweight/((admissionheight/100)^2)]
icu[bmi < 10 | bmi > 60,bmi := NA]
icu <- icu[order(patientunitstayid)]

fwrite(icu,"data/cohort.csv")
