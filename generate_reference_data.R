rm(list=ls())

library(pracma)
library(tictoc)
library(glmnet)
library(ROCR)
library(matrixStats)
library(parallel)
library(ggplot2)
library(xgboost)
library(lubridate)

eicu.clinical.data = readRDS("data/eicu/clinical_data_icd9_sofa_vent.rds")
eicu.patient.result = readRDS("data/eicu/patient_data.rds")
eicu.sofa.scores = readRDS("data/eicu/sofa_scores.rds")

source("src/R/eicu/functions/generate_sampling_rate_table.R")
source("src/R/eicu/functions/eval_carry_forward.R")
source("src/R/eicu/functions/eval_interval.R")
source("src/R/eicu/functions/eval_max_in_past_2.R")
source("src/R/eicu/functions/eval_sum_in_past.R")
source("src/R/eicu/functions/eval_early_prediction_timestamps_combined_rf.R")
source("src/R/eicu/functions/eval_table_with_sofa_2.R")
source("src/R/eicu/functions/generate_table_with_sofa_timestamps.R")
source("src/R/eicu/functions/generate_table_with_sofa_timestamps_history.R")
source("src/R/eicu/functions/eval_table_with_sofa_timestamps_history.R")

lengths = sapply(eicu.sofa.scores, function(x) length(x$timestamps))
eicu.sepsis.labels = sapply(eicu.sofa.scores[lengths>0], function(x) rowSums(x[2:7])>=2)
eicu.has.sepsis = sapply(eicu.sepsis.labels, any)


# Determine shock onsets
eicu.shock.labels = mapply(function(x,y) x&y$lactate&y$vasopressors, eicu.sepsis.labels, eicu.sofa.scores[lengths>0])
eicu.has.shock = sapply(eicu.shock.labels, function(x) any(x,na.rm=T)) #& has.dx[lengths>0]

sepsis.label.lengths = sapply(eicu.sepsis.labels,length)
shock.lengths=sapply(eicu.shock.labels,length)

eicu.shock.onsets = mapply(function(x,y) min(x$timestamps[y],na.rm=T),eicu.sofa.scores[lengths>0][eicu.has.shock],eicu.shock.labels[eicu.has.shock])

### Training data
num.cores = 4
cluster = makeCluster(num.cores)
clusterExport(cluster,c("generate.table.with.sofa.timestamps.history","eval.table.with.sofa.timestamps.history","eval.carry.forward","eval.sum.in.past","eval.max.in.past","eval.interval","generate.table.with.sofa.timestamps.history","eicu.sofa.scores","eicu.clinical.data","eicu.shock.onsets","lengths","eicu.has.sepsis","eicu.has.shock"))
clusterEvalQ(cluster,library(lubridate))
tic("Generate data tables (parallel)")
nonsepsis.training.data = parLapply(cluster, 1:sum(!eicu.has.sepsis), function(x) generate.table.with.sofa.timestamps.history(min(eicu.sofa.scores[lengths>0][!eicu.has.sepsis][[x]]$timestamps),max(eicu.sofa.scores[lengths>0][!eicu.has.sepsis][[x]]$timestamps),10,eicu.clinical.data[lengths>0][!eicu.has.sepsis][[x]],60,12))
nonshock.training.data = parLapply(cluster, 1:sum(eicu.has.sepsis&!eicu.has.shock), function(x) generate.table.with.sofa.timestamps.history(min(eicu.sofa.scores[lengths>0][eicu.has.sepsis&!eicu.has.shock][[x]]$timestamps),max(eicu.sofa.scores[lengths>0][eicu.has.sepsis&!eicu.has.shock][[x]]$timestamps),10,eicu.clinical.data[lengths>0][eicu.has.sepsis&!eicu.has.shock][[x]],60,12))
preshock.training.data = parLapply(cluster, 1:sum(eicu.has.shock), function(x) generate.table.with.sofa.timestamps.history(min(eicu.sofa.scores[lengths>0][eicu.has.shock][[x]]$timestamps),eicu.shock.onsets[x]-60,10,eicu.clinical.data[lengths>0][eicu.has.shock][[x]],60,12))
toc()
stopCluster(cluster)

save(nonsepsis.training.data,nonshock.training.data,preshock.training.data,file="data/eicu/eicu.lstm.train.v2.rdata")

### Test data
num.cores = 4
cluster = makeCluster(num.cores)
clusterExport(cluster,c("eval.table.with.sofa.timestamps.history","eval.carry.forward","eval.sum.in.past","eval.max.in.past","eval.interval","eicu.sofa.scores","eicu.clinical.data","eicu.shock.onsets","lengths","eicu.has.sepsis","eicu.has.shock"))
clusterEvalQ(cluster,library(lubridate))
tic("Generate data tables (parallel)")
nonsepsis.data = parLapply(cluster, 1:sum(!eicu.has.sepsis), function(x) eval.table.with.sofa.timestamps.history(eicu.sofa.scores[lengths>0][!eicu.has.sepsis][[x]]$timestamps,eicu.clinical.data[lengths>0][!eicu.has.sepsis][[x]],60,12))
nonshock.data = parLapply(cluster, 1:sum(eicu.has.sepsis&!eicu.has.shock), function(x) eval.table.with.sofa.timestamps.history(eicu.sofa.scores[lengths>0][eicu.has.sepsis&!eicu.has.shock][[x]]$timestamps,eicu.clinical.data[lengths>0][eicu.has.sepsis&!eicu.has.shock][[x]],60,12))
preshock.data = parLapply(cluster, 1:sum(eicu.has.shock), function(x) eval.table.with.sofa.timestamps.history(eicu.sofa.scores[lengths>0][eicu.has.shock][[x]]$timestamps[eicu.sofa.scores[lengths>0][eicu.has.shock][[x]]$timestamps<=eicu.shock.onsets[x]],eicu.clinical.data[lengths>0][eicu.has.shock][[x]],60,12))
toc()
stopCluster(cluster)

save(nonsepsis.data,nonshock.data,preshock.data,file="data/eicu/eicu.lstm.test.v2.rdata")

# Note: takes about 9 hours to run on an i7-7740x
# Might be faster with more cores on a 3900x or better