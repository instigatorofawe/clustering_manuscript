rm(list=ls())

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
demographics = readRDS("data/eicu/demographics.rds")
comorbidities = rbind(comorbidities, t(as.matrix(broad.spectrum)), t(significant.rx), demographics)

source("src/R/eicu/functions/generate_sampling_rate_table.R")
source("src/R/eicu/functions/eval_carry_forward.R")
source("src/R/eicu/functions/eval_interval.R")
source("src/R/eicu/functions/eval_max_in_past_2.R")
source("src/R/eicu/functions/eval_sum_in_past.R")
source("src/R/eicu/functions/eval_early_prediction_timestamps_combined_rf_comorbidities.R")
source("src/R/eicu/functions/eval_table_with_sofa.R")
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

# Generate table of sampling rates
tic("Generating sampling rate data")
sampling.rate.data = lapply(clinical.data, generate.sampling.rate.table)
toc()

num.cores = detectCores()
cluster = makeCluster(num.cores)
clusterExport(cluster,c("eval.table.with.sofa","eval.carry.forward","eval.sum.in.past","eval.max.in.past","eval.interval","generate.table.with.sofa.timestamps","sofa.scores","sampling.rate.data","shock.onsets","lengths","has.sepsis","has.shock"))

tic("Generate data tables (parallel)")
nonsepsis.data.sampling.rate = parLapply(cluster, 1:sum(!has.sepsis), function(x) generate.table.with.sofa.timestamps(min(sofa.scores[lengths>0][!has.sepsis][[x]]$timestamps),max(sofa.scores[lengths>0][!has.sepsis][[x]]$timestamps),100,sampling.rate.data[lengths>0][!has.sepsis][[x]]))
nonshock.data.sampling.rate = parLapply(cluster, 1:sum(has.sepsis&!has.shock), function(x) generate.table.with.sofa.timestamps(min(sofa.scores[lengths>0][has.sepsis&!has.shock][[x]]$timestamps),max(sofa.scores[lengths>0][has.sepsis&!has.shock][[x]]$timestamps),100,sampling.rate.data[lengths>0][has.sepsis&!has.shock][[x]]))
preshock.data.sampling.rate = parLapply(cluster, 1:sum(has.shock), function(x) generate.table.with.sofa.timestamps(shock.onsets[x]-120,shock.onsets[x]-60,100,sampling.rate.data[lengths>0][has.shock][[x]]))
toc()
stopCluster(cluster)


num.cores = detectCores()
cluster = makeCluster(num.cores)
clusterExport(cluster,c("eval.table.with.sofa.comorbidities","eval.carry.forward","eval.sum.in.past","eval.max.in.past","eval.interval","generate.table.with.sofa.timestamps","sofa.scores","clinical.data","shock.onsets","lengths","has.sepsis","has.shock","nonsepsis.data.sampling.rate","nonshock.data.sampling.rate","preshock.data.sampling.rate","comorbidities"))
tic("Data tables (parallel)")
nonsepsis.data = parLapply(cluster, 1:sum(!has.sepsis), function(x) eval.table.with.sofa.comorbidities(nonsepsis.data.sampling.rate[[x]]$timestamps,clinical.data[lengths>0][!has.sepsis][[x]],comorbidities[,!has.sepsis][,x]))
nonshock.data = parLapply(cluster, 1:sum(has.sepsis&!has.shock), function(x) eval.table.with.sofa.comorbidities(nonshock.data.sampling.rate[[x]]$timestamps,clinical.data[lengths>0][has.sepsis&!has.shock][[x]],comorbidities[,has.sepsis&!has.shock][,x]))
preshock.data = parLapply(cluster, 1:sum(has.shock), function(x) eval.table.with.sofa.comorbidities(preshock.data.sampling.rate[[x]]$timestamps,clinical.data[lengths>0][has.shock][[x]],comorbidities[,has.shock][,x]))
toc()

stopCluster(cluster)

nonsepsis.sample = runif(n = length(nonsepsis.data)) < 0.7
nonshock.sample = runif(n = length(nonshock.data)) < 0.7
preshock.sample = runif(n = length(preshock.data)) < 0.7

a0 = do.call(rbind,nonsepsis.data[nonsepsis.sample])
a1 = do.call(rbind,nonshock.data[nonshock.sample])
a2 = do.call(rbind, preshock.data[preshock.sample])

b0 = do.call(rbind,nonsepsis.data.sampling.rate[nonsepsis.sample])[,-1]
b1 = do.call(rbind,nonshock.data.sampling.rate[nonshock.sample])[,-1]
b2 = do.call(rbind,preshock.data.sampling.rate[preshock.sample])[,-1]

x0 = cbind(a0,b0)
x1 = cbind(a1,b1)
x2 = cbind(a2,b2)

#rm(a0,a1,a2,b0,b1,b2)
#rm(nonsepsis.data.sampling.rate,nonshock.data.sampling.rate,preshock.data.sampling.rate)

x = rbind(x0,x1,x2)
x = as.matrix(x)
means = colMeans(x,na.rm=T)
for (i in 1:dim(x)[2]) {
    x[is.na(x[,i]),i] = means[i]
}

y = as.factor(c(rep(0,dim(x0)[1]),rep(0,dim(x1)[1]),rep(1,dim(x2)[1])))

rm(x0,x1,x2)

tic("XGBoost")
model = xgboost(x,y,nrounds=25)
toc()

#tic("Glmnet")
#model = glmnet(x,y,family="binomial")
#toc()

a0 = do.call(rbind,nonsepsis.data[!nonsepsis.sample])
a1 = do.call(rbind,nonshock.data[!nonshock.sample])
a2 = do.call(rbind, preshock.data[!preshock.sample])

b0 = do.call(rbind,nonsepsis.data.sampling.rate[!nonsepsis.sample])[,-1]
b1 = do.call(rbind,nonshock.data.sampling.rate[!nonshock.sample])[,-1]
b2 = do.call(rbind,preshock.data.sampling.rate[!preshock.sample])[,-1]

x0 = cbind(a0,b0)
x1 = cbind(a1,b1)
x2 = cbind(a2,b2)

x = as.matrix(rbind(x0,x1,x2))
for (i in 1:dim(x)[2]) {
    x[is.na(x[,i]),i] = means[i]
}
y.cts = c(rep(0,dim(x0)[1]),rep(0,dim(x1)[1]),rep(1,dim(x2)[1]))
y = as.factor(c(rep(0,dim(x0)[1]),rep(0,dim(x1)[1]),rep(1,dim(x2)[1])))

z = predict(model, x)

#index = max(which(model$df<=10))

#pred = prediction(z[,index],y.cts)
pred = prediction(z,y.cts)

perf = performance(pred,"auc")

fprintf("Point-by-point classification: %f AUC\n",perf@y.values[[1]])

num.cores = detectCores()
cluster = makeCluster(num.cores)
clusterExport(cluster,c("eval.early.prediction.timestamps.combined.rf.comorbidities","eval.table.with.sofa","eval.table.with.sofa.comorbidities","eval.carry.forward","eval.sum.in.past","eval.max.in.past","eval.interval","means","model","clinical.data","sampling.rate.data","lengths","has.sepsis","has.shock","preshock.sample","shock.onsets","comorbidities"))
clusterEvalQ(cluster,library(glmnet))

tic("Early prediction, parallel (all, with timestamps)")
nonsepsis.predictions = parLapply(cluster, 1:sum(!has.sepsis), function(x) eval.early.prediction.timestamps.combined.rf.comorbidities(model,clinical.data[lengths>0][!has.sepsis][[x]],sampling.rate.data[lengths>0][!has.sepsis][[x]],comorbidities[,!has.sepsis][,x],means))
nonshock.predictions = parLapply(cluster, 1:sum(has.sepsis&!has.shock), function(x) eval.early.prediction.timestamps.combined.rf.comorbidities(model,clinical.data[lengths>0][has.sepsis&!has.shock][[x]],sampling.rate.data[lengths>0][has.sepsis&!has.shock][[x]],comorbidities[,has.sepsis&!has.shock][,x],means))
preshock.predictions = parLapply(cluster, 1:sum(has.shock), function(x) eval.early.prediction.timestamps.combined.rf.comorbidities(model,clinical.data[lengths>0][has.shock][[x]],sampling.rate.data[lengths>0][has.shock][[x]],comorbidities[,has.shock][,x],means,-Inf,shock.onsets[x]))
# shock.predictions = parLapply(cluster, 1:sum(has.shock), function(x) eval.early.prediction.timestamps.combined.rf(model,clinical.data[lengths>0][has.shock][[x]],means))
toc()
stopCluster(cluster)

# Determine threshold on training set.
#nonsepsis.maxes = sapply(nonsepsis.predictions, function(x) max(x$predictions[,index]))[nonsepsis.sample]
#nonshock.maxes = sapply(nonshock.predictions, function(x) max(x$predictions[,index]))[nonshock.sample]
#shock.maxes = sapply(preshock.predictions, function(x) max(x$predictions[,index]))[preshock.sample]
nonsepsis.maxes = sapply(nonsepsis.predictions, function(x) max(x$predictions))[nonsepsis.sample]
nonshock.maxes = sapply(nonshock.predictions, function(x) max(x$predictions))[nonshock.sample]
shock.maxes = sapply(preshock.predictions, function(x) max(x$predictions))[preshock.sample]
shock.maxes = shock.maxes[!is.infinite(shock.maxes)]

pred = prediction(c(nonsepsis.maxes,nonshock.maxes,shock.maxes),c(rep(0,length(nonsepsis.maxes)+length(nonshock.maxes)),rep(1,length(shock.maxes))))
perf = performance(pred,"tpr","fpr")
losses = mapply(function(x,y) sqrt(x^2+(1-y)^2), perf@x.values[[1]], perf@y.values[[1]])
threshold = perf@alpha.values[[1]][which.min(losses)]

#nonsepsis.maxes = sapply(nonsepsis.predictions, function(x) max(x$predictions[,index]))[!nonsepsis.sample]
#nonshock.maxes = sapply(nonshock.predictions, function(x) max(x$predictions[,index]))[!nonshock.sample]
#shock.maxes = sapply(preshock.predictions, function(x) max(x$predictions[,index]))[!preshock.sample]
nonsepsis.maxes = sapply(nonsepsis.predictions, function(x) max(x$predictions))[!nonsepsis.sample]
nonshock.maxes = sapply(nonshock.predictions, function(x) max(x$predictions))[!nonshock.sample]
shock.maxes = sapply(preshock.predictions, function(x) max(x$predictions))[!preshock.sample]

onset.times = shock.onsets[!preshock.sample][!is.infinite(shock.maxes)]
shock.maxes = shock.maxes[!is.infinite(shock.maxes)]

has.detection = shock.maxes >= threshold
detection.times = sapply(preshock.predictions[!preshock.sample][has.detection], function(x) x$timestamps[which.min(x$predictions>=threshold)])
early.pred.times = onset.times[has.detection] - detection.times

pred = prediction(c(nonsepsis.maxes,nonshock.maxes,shock.maxes),c(rep(0,length(nonsepsis.maxes)+length(nonshock.maxes)),rep(1,length(shock.maxes))))
perf = performance(pred,"auc")

fprintf("Early prediction: %f AUC\n", perf@y.values[[1]])
perf = performance(pred,"ppv","sens")
perf = performance(pred,"tpr","fpr")

ppv = sum(shock.maxes>=threshold)/(sum(shock.maxes>=threshold)+sum(nonsepsis.maxes>=threshold)+sum(nonshock.maxes>=threshold))
sens = sum(shock.maxes>=threshold)/length(shock.maxes)
spec = (sum(nonsepsis.maxes<threshold)+sum(nonshock.maxes<threshold))/(length(nonsepsis.maxes)+length(nonshock.maxes))
# save(nonsepsis.sample, nonshock.sample, preshock.sample, nonsepsis.data.sampling.rate, nonshock.data.sampling.rate, preshock.data.sampling.rate, nonsepsis.data, nonshock.data, preshock.data, model, nonsepsis.predictions, nonshock.predictions, preshock.predictions, shock.predictions, threshold, file="reference_dataset_2.rdata")
toc()

save(nonsepsis.sample, nonshock.sample, preshock.sample, nonsepsis.data, nonshock.data, preshock.data, nonsepsis.data.sampling.rate, nonshock.data.sampling.rate, preshock.data.sampling.rate, model, comorbidities, nonsepsis.predictions, nonshock.predictions, preshock.predictions, threshold, file="data/eicu/reference_dataset_rx_combined_2.rdata")

#data = data.frame(x=perf@y.values[[1]],y=median.edts/60)
