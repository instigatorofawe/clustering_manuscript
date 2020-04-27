library(pracma)
library(tictoc)
library(glmnet)
library(ROCR)
library(matrixStats)
library(parallel)
library(ggplot2)
library(xgboost)

tic("Total time to run")

clinical.data = readRDS("data/eicu/clinical_data_icd9_sofa_vent.rds")
patient.result = readRDS("data/eicu/patient_data.rds")
sofa.scores = readRDS("data/eicu/sofa_scores.rds")
comorbidities = readRDS("data/eicu/comorbidities.rds")
broad.spectrum = readRDS("data/eicu/has.broad.spectrum.rds")
significant.rx = readRDS("data/eicu/has.significant.rx.combined.rds")
comorbidities = rbind(comorbidities, t(as.matrix(broad.spectrum)), t(significant.rx))

source("src/R/eicu/functions/generate_sampling_rate_table.R")
source("src/R/eicu/functions/eval_carry_forward.R")
source("src/R/eicu/functions/eval_interval.R")
source("src/R/eicu/functions/eval_max_in_past_2.R")
source("src/R/eicu/functions/eval_sum_in_past.R")
source("src/R/eicu/functions/eval_early_prediction_timestamps_combined_rf_comorbidities.R")
source("src/R/eicu/functions/eval_table_with_sofa_2.R")
source("src/R/eicu/functions/eval_table_with_sofa_comorbidities.R")
source("src/R/eicu/functions/generate_table_with_sofa_timestamps.R")

lengths = sapply(sofa.scores, function(x) length(x$timestamps))

sepsis.labels = sapply(sofa.scores[lengths>0], function(x) rowSums(x[2:7])>=2)
has.sepsis = sapply(sepsis.labels, any)

sepsis.timestamps = mapply(function(x,y) x$timestamps[y], sofa.scores[lengths>0][has.sepsis],sepsis.labels[has.sepsis])

# Determine shock onsets
shock.labels = mapply(function(x,y) x&y$lactate&y$vasopressors, sepsis.labels, sofa.scores[lengths>0])
has.shock = sapply(shock.labels, function(x) any(x,na.rm=T)) #& has.dx[lengths>0]

sepsis.label.lengths = sapply(sepsis.labels,length)
shock.lengths=sapply(shock.labels,length)
sofa.timestamps.lengths = sapply(sofa.scores[lengths>0],function(x)length(x$timestamps))

shock.onsets = mapply(function(x,y) min(x$timestamps[y],na.rm=T),sofa.scores[lengths>0][has.shock],shock.labels[has.shock])
load("data/eicu/reference_dataset_rx_combined_2.rdata")

mortality.patients = patient.result$patientunitstayid[patient.result$hospitaldischargestatus=="Expired"]
has.mortality = sapply(clinical.data[lengths>0], function(x) is.element(x$subject.id,mortality.patients))
discharge.times = sapply(clinical.data[lengths>0][!has.mortality], function(x) patient.result$hospitaldischargeoffset[patient.result$patientunitstayid==x$subject.id])
icu.discharge.times = sapply(clinical.data[lengths>0][!has.mortality], function(x) patient.result$unitdischargeoffset[patient.result$patientunitstayid==x$subject.id])
