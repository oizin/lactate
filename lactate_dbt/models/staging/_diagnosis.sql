

with tmp as (
  select patientunitstayid,admitdxenteredoffset,admitdxpath,admitdxname,admitdxtext 
  from `eicu_crd.admissiondx`
  where admitdxpath like "%All Diagnosis%"
)

select patientunitstayid,
      (case 
        when admitdxpath like "%Operative%" then 1 else 0 
      end) as operative,
      admitdxtext as apache_dx,
      (case
        when lower(admitdxtext) like "%diabetic ketoacidosis%" then 1 else 0
      end) as dka,
      (case
        when lower(admitdxtext) like "%hyperglycemic hyperosmolar%" then 1 else 0
      end) as hhs,
      (case
        when lower(admitdxtext) like "%acid-base/electrolyte disturbance%" then 1 else 0
      end) as acid_base_disturbance,
      (case
        when lower(admitdxtext) like "%hypoglycemia%" then 1 else 0
      end) as hypoglycaemia,
      (case
        when lower(admitdxtext) like "%sepsis%" then 1 else 0
      end) as sepsis
from tmp