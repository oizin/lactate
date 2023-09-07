
select t1.patientunitstayid,
        t1.patienthealthsystemstayid,
        coalesce(t2.time,t3.time) as time,
        t2.glucose_lab,
        t2.glucose_bedside,
        coalesce(t2.glucose_lab,t2.glucose_bedside) as glucose,
        t2.lactate,
        t3.bilirubin,
        t3.potassium,
        t3.sodium,
        t3.hgb,
        t3.chloride,
        t3.hct,
        t3.bicarbonate,
        t3.creatinine,
        t3.bun,
        t3.calcium,
        t3.anion_gap,
        t3.base_excess,
        t3.ph,
        t3.paco2,
        t3.pao2,
        (case when insulin_sc = 1 then 1 else 0 end) as insulin_sc,
        (case when insulin_iv = 1 then 1 else 0 end) as insulin_iv,
        (case when insulin_sc = 1 or insulin_iv = 1 then 1 else 0 end) as insulin,
        insulin_init_time,
        icu_mort,
        hosp_mort,
        apacheivascore,
        predictedhospitalmortality,
        predictedicumortality,
        elective_admission,
        organ_system,
        apache_dx,
        operative,
        dka,
        hhs,
        acid_base_disturbance,
        hypoglycaemia,
        t6.sepsis as sepsis_apache,
        t9.diabetes,
        t1.unitvisitnumber,
        t1.hospitalid,
        t1.region,
        t1.unittype,
        t1.hospitaladmitoffset,
        t1.hospitaldischargeoffset,
        t1.unitadmitoffset,
        t1.unitdischargeoffset,
        t1.apache_iv,
        t1.hospitaldischargeyear,
        t1.age,
        t1.gender,
        t1.ethnicity,
        t1.admissionheight,
        t1.admissionweight,
        t1.dischargeweight,
        t1.icu_los_hours,
        saps3day1,
        saps3today,
        saps3yesterday,
        ventday1,
        oobventday1,
        oobintubday1,
        day1meds,
        day1verbal,
        day1motor,
        day1eyes
from 
{{ ref('_icustay_detail') }} t1
left join
{{ ref('glucose') }} t2
on t1.patientunitstayid = t2.patientunitstayid
left join
{{ ref('_labs') }} t3
on t1.patientunitstayid = t3.patientunitstayid and
  t2.time = t3.time
left join
{{ ref('insulin') }} t4
on t1.patientunitstayid = t4.patientunitstayid
left join
{{ ref('_apache') }} t5
on t1.patientunitstayid = t5.patientunitstayid
left join
{{ ref('_diagnosis') }} t6
on t1.patientunitstayid = t6.patientunitstayid
left join
{{ ref('_elective') }} t7
on t1.patientunitstayid = t7.patientunitstayid
left join
{{ ref('_organ_system') }} t8
on t1.patientunitstayid = t8.patientunitstayid
left join
`eicu_crd_glucose.diabetes` t9
on t1.patientunitstayid = t9.patientunitstayid
left join
{{ ref('_patient') }} t10
on t1.patientunitstayid = t10.patientunitstayid
left join
`eicu_crd.apachepredvar` t11
on t1.patientunitstayid = t11.patientunitstayid
order by patientunitstayid, time





