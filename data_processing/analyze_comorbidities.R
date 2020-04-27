rm(list=ls())

library(ROCR)
library(glmnet)
library(pracma)

clinical.data = readRDS("processed/clinical_data_icd9_sofa_vent.rds")
patient.result = readRDS("processed/patient_data.rds")
sofa.scores = readRDS("processed/sofa_scores.rds")

past.history = readRDS("past.history.rds")

lengths = sapply(sofa.scores, function(x) length(x$timestamps))
sepsis.labels = sapply(sofa.scores[lengths>0], function(x) rowSums(x[2:7])>=2)
has.sepsis = sapply(sepsis.labels, any)

shock.labels = mapply(function(x,y) x&y$lactate&y$vasopressors, sepsis.labels, sofa.scores[lengths>0])
has.shock = sapply(shock.labels, function(x) any(x,na.rm=T)) #& has.dx[lengths>0]

shock.subjects = sapply(clinical.data[has.shock], function(x) x$subject.id)
clinical.subjects = sapply(clinical.data[lengths>0], function(x) x$subject.id)
mortality.subjects = patient.result$patientunitstayid[patient.result$unitdischargestatus=="Expired"]

# For each comorbidity, calculate a contingency table for mortality
unique.comorbidities = unique(past.history$pasthistorypath)
all.subjects = patient.result$patientunitstayid

has.mortality = is.element(clinical.subjects,mortality.subjects)

p.mortality=sapply(unique.comorbidities, function(x)  {
    comorbidity.subjects = unique(past.history$patientunitstayid[past.history$pasthistorypath==x])
    has.comorbidity = is.element(clinical.subjects,comorbidity.subjects)

    contingency = array(NA, dim=c(2,2))
    contingency[1,1] = sum(!has.mortality&!has.comorbidity)# No mortality, no comorbidity
    contingency[1,2] = sum(!has.mortality&has.comorbidity)# No mortality with comorbidity
    contingency[2,1] = sum(has.mortality&!has.comorbidity)
    contingency[2,2] = sum(has.mortality&has.comorbidity)
    return(fisher.test(contingency)$p.value)
})

sum(p.adjust(p.mortality,"bonferroni")<0.01)
sum(p.adjust(p.mortality,"BH")<0.01)

p.shock=sapply(unique.comorbidities, function(x)  {
    comorbidity.subjects = unique(past.history$patientunitstayid[past.history$pasthistorypath==x])
    has.comorbidity = is.element(clinical.subjects,comorbidity.subjects)
    
    contingency = array(NA, dim=c(2,2))
    contingency[1,1] = sum(!has.shock&!has.comorbidity)
    contingency[1,2] = sum(!has.shock&has.comorbidity)
    contingency[2,1] = sum(has.shock&!has.comorbidity)
    contingency[2,2] = sum(has.shock&has.comorbidity)
    return(fisher.test(contingency)$p.value)
})

sum(p.adjust(p.shock,"bonferroni")<0.01)
sum(p.adjust(p.shock,"BH")<0.01)

p.shock[p.adjust(p.shock,"bonferroni")<0.01]

indices = which(p.adjust(p.shock,"bonferroni")<0.01)
significant.comorbidities = unique.comorbidities[indices]

comorbidity.subjects = lapply(significant.comorbidities, function(x) unique(past.history$patientunitstayid[past.history$pasthistorypath==x]))
has.comorbidities = sapply(clinical.subjects, function(x) sapply(comorbidity.subjects, function(y) is.element(x,y)))
saveRDS(has.comorbidities,file="comorbidities_shock.rds")

load("reference_dataset.rdata")

# Determine threshold on training set.
index = max(which(model$df<=10))
lengths = sapply(sofa.scores, function(x) length(x$timestamps))
sepsis.labels = sapply(sofa.scores[lengths>0], function(x) rowSums(x[2:7])>=2)
has.sepsis = sapply(sepsis.labels, any)

# Determine shock onsets
shock.labels = mapply(function(x,y) x&y$lactate&y$vasopressors, sepsis.labels, sofa.scores[lengths>0])
has.shock = sapply(shock.labels, function(x) any(x,na.rm=T)) #& has.dx[lengths>0]

sepsis.label.lengths = sapply(sepsis.labels,length)
shock.lengths=sapply(shock.labels,length)
sofa.timestamps.lengths = sapply(sofa.scores[lengths>0],function(x)length(x$timestamps))

shock.onsets = mapply(function(x,y) min(x$timestamps[y],na.rm=T),sofa.scores[lengths>0][has.shock],shock.labels[has.shock])

nonsepsis.maxes = sapply(nonsepsis.predictions, function(x) max(x$predictions[,index]))[nonsepsis.sample]
nonshock.maxes = sapply(nonshock.predictions, function(x) max(x$predictions[,index]))[nonshock.sample]
shock.maxes = sapply(preshock.predictions, function(x) max(x$predictions[,index]))[preshock.sample]
shock.maxes = shock.maxes[!is.infinite(shock.maxes)]

pred = prediction(c(nonsepsis.maxes,nonshock.maxes,shock.maxes),c(rep(0,length(nonsepsis.maxes)+length(nonshock.maxes)),rep(1,length(shock.maxes))))
perf = performance(pred,"tpr","fpr")
losses = mapply(function(x,y) sqrt(x^2+(1-y)^2), perf@x.values[[1]], perf@y.values[[1]])
threshold = perf@alpha.values[[1]][which.min(losses)]

nonsepsis.maxes = sapply(nonsepsis.predictions, function(x) max(x$predictions[,index]))[!nonsepsis.sample]
nonshock.maxes = sapply(nonshock.predictions, function(x) max(x$predictions[,index]))[!nonshock.sample]
shock.maxes = sapply(preshock.predictions, function(x) max(x$predictions[,index]))[!preshock.sample]

onset.times = shock.onsets[!preshock.sample][!is.infinite(shock.maxes)]
shock.maxes = shock.maxes[!is.infinite(shock.maxes)]

has.detection = shock.maxes >= threshold
detection.times = sapply(preshock.predictions[!preshock.sample][has.detection], function(x) x$timestamps[which.min(x$predictions>=threshold)])
early.pred.times = onset.times[has.detection] - detection.times

pred = prediction(c(nonsepsis.maxes,nonshock.maxes,shock.maxes),c(rep(0,length(nonsepsis.maxes)+length(nonshock.maxes)),rep(1,length(shock.maxes))))
perf = performance(pred,"auc")
roc = performance(pred,"tpr","fpr")

fprintf("Early prediction: %f AUC\n", perf@y.values[[1]])

# True positive subjects
nonsepsis.maxes = sapply(nonsepsis.predictions, function(x) max(x$predictions[,index]))
nonshock.maxes = sapply(nonshock.predictions, function(x) max(x$predictions[,index]))
shock.maxes = sapply(preshock.predictions, function(x) max(x$predictions[,index]))

true.positive.subjects = clinical.subjects[has.shock][shock.maxes>=threshold]
# False positive subjects
false.positive.subjects = c(clinical.subjects[!has.sepsis][nonsepsis.maxes>=threshold],clinical.subjects[has.sepsis&!has.shock][nonshock.maxes>=threshold])

p.prediction = sapply(unique.comorbidities, function(x) {
    comorbidity.subjects = unique(past.history$patientunitstayid[past.history$pasthistorypath==x])
    contingency = array(NA, dim=c(2,2))
    contingency[1,1] = length(true.positive.subjects)-length(intersect(true.positive.subjects,comorbidity.subjects))# True positive & no comorbidity
    contingency[1,2] = length(false.positive.subjects)-length(intersect(false.positive.subjects,comorbidity.subjects))# False positive & no comorbidity
    contingency[2,1] = length(intersect(true.positive.subjects,comorbidity.subjects))# True positive & comorbidity
    contingency[2,2] = length(intersect(false.positive.subjects,comorbidity.subjects))# False positive & comorbidity
    return(fisher.test(contingency)$p.value)
})

sum(p.adjust(p.prediction,"bonferroni")<0.01,na.rm=T)
sum(p.adjust(p.prediction,"BH")<0.01,na.rm=T)

p.prediction[p.adjust(p.prediction,"BH")<0.01]

indices = which(p.adjust(p.prediction,"BH")<0.01)
significant.comorbidities = unique.comorbidities[indices]

comorbidity.subjects = lapply(significant.comorbidities, function(x) unique(past.history$patientunitstayid[past.history$pasthistorypath==x]))
has.comorbidities = sapply(clinical.subjects, function(x) sapply(comorbidity.subjects, function(y) is.element(x,y)))
saveRDS(has.comorbidities,file="comorbidities.rds")