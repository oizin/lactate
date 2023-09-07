
with tmp as (
select patientunitstayid,labresultoffset as time,
case when lower(labname) = 'total bilirubin' then labresult else null end as bilirubin,
case when lower(labname) = 'potassium' then labresult else null end as potassium,
case when lower(labname) = 'sodium' then labresult else null end as sodium,
case when lower(labname) = 'hgb' then labresult else null end as hgb,
case when lower(labname) = 'chloride' then labresult else null end as chloride,
case when lower(labname) = 'hct' then labresult else null end as hct,
case when lower(labname) = 'creatinine' then labresult else null end as creatinine,
case when lower(labname) = 'bun' then labresult else null end as bun,
case when lower(labname) = 'calcium' then labresult else null end as calcium,
case when lower(labname) in ('bicarbonate','hco3') then labresult else null end as bicarbonate,
case when lower(labname) = 'anion gap' then labresult else null end as anion_gap,
case when lower(labname) = 'base excess' then labresult else null end as base_excess,
case when lower(labname) = 'ph' then labresult else null end as ph,
case when lower(labname) = 'paco2' then labresult else null end as paco2,
case when lower(labname) = 'pao2' then labresult else null end as pao2,
from `eicu_crd.lab`
where lower(labname) in ('bilirubin','potassium','sodium','hgb','chloride',
                        'hct','creatinine','bun','calcium','bicarbonate',
                        'anion gap','total bilirubin','base excess','ph',
                        'paco2','pao2') and 
      labresultoffset > - 60*12 and
      labresultoffset < 60*24
order by patientunitstayid,labresultoffset
)
select patientunitstayid,
      time,
       avg(bilirubin) as bilirubin,
       avg(potassium) as potassium,
       avg(sodium) as sodium,
       avg(hgb) as hgb,
       avg(chloride) as chloride,
       avg(hct) as hct,
       avg(bicarbonate) as bicarbonate,
       avg(creatinine) as creatinine,
       avg(bun) as bun,
       avg(calcium) as calcium,
       avg(anion_gap) as anion_gap,
       avg(base_excess) as base_excess,
       avg(ph) as ph,
       avg(paco2) as paco2,
       avg(pao2) as pao2
from tmp
group by patientunitstayid,time