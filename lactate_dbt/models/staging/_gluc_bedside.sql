
SELECT patientunitstayid,
      nursingchartoffset as time,
      nursingchartvalue as glucose_bedside
FROM `eicu_crd.nursecharting`
WHERE lower(nursingchartcelltypevalname) like '%bedside glucose%' AND
    nursingchartoffset > -60*12 AND nursingchartoffset < 60*24