
with tmp as (
  select patientunitstayid,
        time,
        glucose_lab,
        glucose_bedside,
        lactate
  from {{ ref('_gluc_lact_labs') }}
  union distinct
  select patientunitstayid,
        time,
        CAST(NULL AS FLOAT64) as glucose_lab,
        CAST(glucose_bedside AS FLOAT64) as glucose_bedside,
        CAST(NULL AS FLOAT64) as lactate
  from {{ ref('_gluc_bedside') }}
)

select patientunitstayid,
        time,
        avg(glucose_lab) as glucose_lab,
        avg(glucose_bedside) as glucose_bedside,
        avg(lactate) as lactate
from tmp
group by patientunitstayid,
        time
order by patientunitstayid,
        time