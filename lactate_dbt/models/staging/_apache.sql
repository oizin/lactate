
SELECT patientunitstayid,apachescore as apacheivascore,
predictedhospitalmortality,predictedicumortality
FROM `eicu_crd.apachepatientresult` 
WHERE apacheversion = 'IVa'