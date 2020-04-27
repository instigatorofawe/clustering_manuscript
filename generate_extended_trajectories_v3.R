rm(list=ls())
source("load_reference_data_2.R")

library(kernlab)
library(ggplot2)
library(tidyr)
library(tictoc)
library(class)

# Get time of threshold crossing for each patient
shock.maxes = sapply(shock.predictions, function(x) max(x$predictions))
has.detection.shock = shock.maxes >= threshold
shock.detection.times = sapply(shock.predictions[has.detection.shock], function(x) x$timestamps[min(which(x$predictions>=threshold))])

#shock.detection.times = sapply(shock.predictions[has.detection.shock], function(x) x$timestamps[which.min(x$predictions>=threshold)])
shock.onset.times = shock.onsets[has.detection.shock]
shock.ewts = shock.detection.times - shock.onset.times

preshock.maxes = sapply(preshock.predictions, function(x) max(x$predictions))
has.detection.preshock = preshock.maxes >= threshold
preshock.detection.times = sapply(preshock.predictions[has.detection.preshock], function(x) x$timestamps[min(which(x$predictions>=threshold))])
preshock.dataset.length = sapply(preshock.predictions[has.detection.preshock], function(x) x$timestamps[length(x$timestamps)] - x$timestamps[1])

nonsepsis.maxes = sapply(nonsepsis.predictions, function(x) max(x$predictions))
has.detection.nonsepsis = nonsepsis.maxes >= threshold
nonsepsis.detection.times = sapply(nonsepsis.predictions[has.detection.nonsepsis], function(x) x$timestamps[min(which(x$predictions>=threshold))])

nonshock.maxes = sapply(nonshock.predictions, function(x) max(x$predictions))
has.detection.nonshock = nonshock.maxes >= threshold
nonshock.detection.times = sapply(nonshock.predictions[has.detection.nonshock], function(x) x$timestamps[min(which(x$predictions>=threshold))])

windows = seq(from=-12,to=30,by=1)

# Build trajectories: post-threshold crossing
nonsepsis.post.threshold.trajectories = mapply(function(x,y) {
    eval.timestamps = windows*60+y
    return(eval.carry.forward(eval.timestamps,x$timestamps,x$predictions))
},nonsepsis.predictions[has.detection.nonsepsis],nonsepsis.detection.times,SIMPLIFY = F)

nonshock.post.threshold.trajectories = mapply(function(x,y) {
    eval.timestamps = windows*60+y
    return(eval.carry.forward(eval.timestamps,x$timestamps,x$predictions))
},nonshock.predictions[has.detection.nonshock],nonshock.detection.times,SIMPLIFY = F)

shock.post.threshold.trajectories = mapply(function(x,y) {
    eval.timestamps = windows*60+y
    return(eval.carry.forward(eval.timestamps,x$timestamps,x$predictions))
},shock.predictions[has.detection.shock],shock.detection.times,SIMPLIFY = F)

nonsepsis.post.table = do.call(rbind,nonsepsis.post.threshold.trajectories)
nonshock.post.table = do.call(rbind,nonshock.post.threshold.trajectories)
shock.post.table = do.call(rbind,shock.post.threshold.trajectories)

ext.trajectories = rbind(nonsepsis.post.table,nonshock.post.table,shock.post.table)
labels = c(rep(1,dim(nonsepsis.post.table)[1]),rep(2,dim(nonshock.post.table)[1]),rep(3,dim(shock.post.table)[1]))

save(nonsepsis.post.table,nonshock.post.table,shock.post.table,ext.trajectories,labels,file="extended.detection.trajectories.v3.rdata")
