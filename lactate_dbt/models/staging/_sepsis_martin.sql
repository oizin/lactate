
select patientunitstayid,STRING_AGG(concat(icd9code)) as icd9code,
      min(diagnosisoffset) as sepsis_time,
      1 as sepsis_martin
from `eicu_crd.diagnosis` 
where icd9code like '%038%' or icd9code like '%020.0%' or icd9code like '%790.7%'
or icd9code like '%117.9%' or icd9code like '%112.5%' or icd9code like '%112.81%'
group by patientunitstayid
order by patientunitstayid