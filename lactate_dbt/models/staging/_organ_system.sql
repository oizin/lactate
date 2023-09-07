
with tmp as (
  select patientunitstayid,admitdxenteredoffset,admitdxpath,admitdxname,admitdxtext 
  from `eicu_crd.admissiondx`
  where admitdxpath like "%Organ System%"
)

select patientunitstayid,
      admitdxtext as organ_system
from tmp