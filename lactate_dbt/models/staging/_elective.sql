
with tmp as (
  select patientunitstayid,admitdxenteredoffset,admitdxpath,admitdxname,admitdxtext 
  from `eicu_crd.admissiondx`
  where admitdxpath like "%Elective%"
)

select patientunitstayid,
      (case 
        when admitdxname = "Yes" then 1 else 0 
      end) as elective_admission
from tmp