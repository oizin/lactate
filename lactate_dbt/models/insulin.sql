
with tmp as (
  select patientunitstayid,
        time,
        case when route = "sc" then 1 else 0 end as insulin_sc,
        case when route = "iv" then 1 else 0 end as insulin_iv
  from `eicu_crd_glucose.insulin_orders`
  where time > -12*60 and time < 60*24
  union distinct 
  select patientunitstayid,
        time,
        case when route = "sc" then 1 else 0 end as insulin_sc,
        case when route = "iv" then 1 else 0 end as insulin_iv
  from `eicu_crd_glucose.insulin_administrations`
  where time > -12*60 and time < 60*24
)

select patientunitstayid,
      min(time) as insulin_init_time,
      max(insulin_sc) as insulin_sc,
      max(insulin_iv) as insulin_iv
from tmp
group by patientunitstayid
