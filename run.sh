#!/usr/bin/env bash

cd lactate_dbt
dbt run
cd ..
Rscript ./R/0_data.R
Rscript ./R/2_imputation_models.R
Rscript ./R/2_missing_data_models.R
Rscript ./R/3_weighted_cross_validation.R
Rscript ./R/3_imputed_cross_validation.R
Rscript ./R/3_table1.R
Rscript ./R/3_additional_cv_sensitivity.R
Rscript ./R/4_final_models.R
