
with tmp as (
  select patientunitstayid,unitdischargestatus 
  from `eicu_crd.patient`
)
select patientunitstayid,
    unitdischargestatus,
    (case when unitdischargestatus = "Expired" then 1 else 0 end) as icu_mort
from tmp