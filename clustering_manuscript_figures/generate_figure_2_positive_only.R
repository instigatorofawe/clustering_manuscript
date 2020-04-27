rm(list=ls())

library(pracma)
library(tidyr)
library(ggplot2)
library(matrixStats)
library(scales)

source("src/R/eicu/load_reference_data_2.R")
load("data/eicu/extended.detection.trajectories.v3.rdata")
load("data/eicu/spec.clustering.detection.rdata")

source("src/R/eicu/functions/generate_sampling_rate_table.R")
source("src/R/eicu/functions/eval_carry_forward.R")
source("src/R/eicu/functions/eval_interval.R")
source("src/R/eicu/functions/eval_max_in_past_2.R")
source("src/R/eicu/functions/eval_sum_in_past.R")
source("src/R/eicu/functions/eval_early_prediction_timestamps_combined_rf_comorbidities.R")
source("src/R/eicu/functions/eval_table_with_sofa_2.R")
source("src/R/eicu/functions/eval_table_with_sofa_comorbidities.R")
source("src/R/eicu/functions/generate_table_with_sofa_timestamps.R")

# Labels
labels = c(rep(1,dim(nonsepsis.post.table)[1]),rep(2,dim(nonshock.post.table)[1]),rep(3,dim(shock.post.table)[1]))

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


shock.post.detection.data = vector(mode="list",length=sum(has.detection.shock))
for (i in 1:sum(has.detection.shock)) {
    fprintf("%d of %d...\n",i,sum(has.detection.shock))
    shock.post.detection.data[[i]] = eval.table.with.sofa(shock.detection.times[i]+(-12:12)*60, clinical.data[lengths>0][has.shock][has.detection.shock][[i]])
}

nonshock.post.detection.data = vector(mode="list",length=sum(has.detection.nonshock))
for (i in 1:sum(has.detection.nonshock)) {
    fprintf("%d of %d...\n",i,sum(has.detection.nonshock))
    nonshock.post.detection.data[[i]] = eval.table.with.sofa(nonshock.detection.times[i]+(-12:12)*60, clinical.data[lengths>0][has.sepsis&!has.shock][has.detection.nonshock][[i]])
}

nonsepsis.post.detection.data = vector(mode="list",length=sum(has.detection.nonsepsis))
for (i in 1:sum(has.detection.nonsepsis)) {
    fprintf("%d of %d...\n",i,sum(has.detection.nonsepsis))
    nonsepsis.post.detection.data[[i]] = eval.table.with.sofa(nonsepsis.detection.times[i]+(-12:12)*60, clinical.data[lengths>0][!has.sepsis][has.detection.nonsepsis][[i]])
}

windows = -12:12

confusion = array(NA,dim=c(4,3))
for (i in 1:4) {
    for (j in 1:3) {
        confusion[i,j] = sum(sc@.Data==i&labels==j)
    }
}
tp.fraction = apply(confusion,MARGIN=1,function(x) x[3]/sum(x))
cluster.indices = order(tp.fraction,decreasing=T)



nonsepsis.hr = t(sapply(nonsepsis.post.detection.data, function(x) x$hr))
nonshock.hr = t(sapply(nonshock.post.detection.data, function(x) x$hr))
shock.hr = t(sapply(shock.post.detection.data, function(x) x$hr))
table.hr = rbind(nonsepsis.hr,nonshock.hr,shock.hr)


risk.score.centroids = do.call(rbind,lapply(1:4, function(x) {
    colMeans(table.hr[sc@.Data==x&labels==3,],na.rm=T)
}))
risk.score.sds = do.call(rbind,lapply(1:4, function(x) {
    colSds(table.hr[sc@.Data==x&labels==3,],na.rm=T)
}))

trajectory.confidence.bounds = data.frame(t=rep(NA,length(windows)*4),cluster=rep(NA,length(windows)*4),center=rep(NA,length(windows)*4),lower=rep(NA,length(windows)*4),upper=rep(NA,length(windows)*4))
for (i in 1:4) {
    indices = (i-1)*length(windows)+(1:length(windows))
    trajectory.confidence.bounds$t[indices] = windows
    trajectory.confidence.bounds$cluster[indices] = which(cluster.indices==i)
    trajectory.confidence.bounds$center[indices] = risk.score.centroids[i,]
    trajectory.confidence.bounds$lower[indices] = risk.score.centroids[i,] - risk.score.sds[i,]*1
    trajectory.confidence.bounds$upper[indices] = risk.score.centroids[i,] + risk.score.sds[i,]*1
}
trajectory.confidence.bounds$cluster=as.factor(trajectory.confidence.bounds$cluster)

# ggplot(trajectory.confidence.bounds,aes(group=cluster,x=t,y=center,color=cluster))+geom_line(size=2,alpha=1)+ylab("HR")+geom_ribbon(aes(ymin=lower,ymax=upper,fill=cluster),alpha=0.2)+
#     labs(fill="Cluster", color="Cluster")+guides(alpha=FALSE)+xlab("Time (hrs relative to detection)")

ggplot(trajectory.confidence.bounds,aes(group=cluster,x=t,y=center,color=cluster))+geom_line(size=2)+
    geom_ribbon(aes(ymin=lower,ymax=upper,fill=cluster),alpha=0.3)+
    scale_discrete_manual(aesthetics=c("color","fill"),name="Post-Prediction Clusters",breaks=c(1,2,3,4),labels=c("1 (High-risk)","2","3","4 (Low-risk)"),values=c("#f8766d","#7cae00","#00bfc4","#c77cff"))+
    labs(fill="Cluster",color="Cluster")+xlab("Time (hrs relative to early prediction)")+ylab("Heart Rate (bpm)")+
    theme(legend.justification=c(0,0),
          legend.position=c(0,0),
          legend.box.margin=margin(c(5,5,5,5)),
          legend.background = element_rect(fill=alpha('white',0.5)),
          axis.text.x = element_text(size=12),
          axis.title.x = element_text(size=14),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=14),
          legend.text = element_text(size=12),
          legend.title = element_text(size=14))

nonsepsis.lact = t(sapply(nonsepsis.post.detection.data, function(x) x$lact))
nonshock.lact = t(sapply(nonshock.post.detection.data, function(x) x$lact))
shock.lact = t(sapply(shock.post.detection.data, function(x) x$lact))
table.lact = rbind(nonsepsis.lact,nonshock.lact,shock.lact)


risk.score.centroids = do.call(rbind,lapply(1:4, function(x) {
    colMeans(table.lact[sc@.Data==x&labels==3,],na.rm=T)
}))
risk.score.sds = do.call(rbind,lapply(1:4, function(x) {
    colSds(table.lact[sc@.Data==x&labels==3,],na.rm=T)
}))

trajectory.confidence.bounds = data.frame(t=rep(NA,length(windows)*4),cluster=rep(NA,length(windows)*4),center=rep(NA,length(windows)*4),lower=rep(NA,length(windows)*4),upper=rep(NA,length(windows)*4))
for (i in 1:4) {
    indices = (i-1)*length(windows)+(1:length(windows))
    trajectory.confidence.bounds$t[indices] = windows
    trajectory.confidence.bounds$cluster[indices] = which(cluster.indices==i)
    trajectory.confidence.bounds$center[indices] = risk.score.centroids[i,]
    trajectory.confidence.bounds$lower[indices] = risk.score.centroids[i,] - risk.score.sds[i,]*1
    trajectory.confidence.bounds$upper[indices] = risk.score.centroids[i,] + risk.score.sds[i,]*1
}
trajectory.confidence.bounds$cluster=as.factor(trajectory.confidence.bounds$cluster)

# ggplot(trajectory.confidence.bounds,aes(group=cluster,x=t,y=center,color=cluster))+geom_line(size=2,alpha=1)+ylab("Lactate")+geom_ribbon(aes(ymin=lower,ymax=upper,fill=cluster),alpha=0.2)+
#     labs(fill="Cluster", color="Cluster")+guides(alpha=FALSE)+xlab("Time (hrs relative to detection)")

ggplot(trajectory.confidence.bounds,aes(group=cluster,x=t,y=center,color=cluster))+geom_line(size=2)+
    geom_ribbon(aes(ymin=lower,ymax=upper,fill=cluster),alpha=0.3)+
    scale_discrete_manual(aesthetics=c("color","fill"),name="Post-Prediction Clusters",breaks=c(1,2,3,4),labels=c("1 (High-risk)","2","3","4 (Low-risk)"),values=c("#f8766d","#7cae00","#00bfc4","#c77cff"))+
    labs(fill="Cluster",color="Cluster")+xlab("Time (hrs relative to early prediction)")+ylab("Lactate (mmol/L)")+
    theme(legend.justification=c(0,1),
          legend.position=c(0,1),
          legend.box.margin=margin(c(5,5,5,5)),
          legend.background = element_rect(fill=alpha('white',0.5)),
          axis.text.x = element_text(size=12),
          axis.title.x = element_text(size=14),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=14),
          legend.text = element_text(size=12),
          legend.title = element_text(size=14))


nonsepsis.sbp = t(sapply(nonsepsis.post.detection.data, function(x) x$sbp))
nonshock.sbp = t(sapply(nonshock.post.detection.data, function(x) x$sbp))
shock.sbp = t(sapply(shock.post.detection.data, function(x) x$sbp))
table.sbp = rbind(nonsepsis.sbp,nonshock.sbp,shock.sbp)


risk.score.centroids = do.call(rbind,lapply(1:4, function(x) {
    colMeans(table.sbp[sc@.Data==x&labels==3,],na.rm=T)
}))
risk.score.sds = do.call(rbind,lapply(1:4, function(x) {
    colSds(table.sbp[sc@.Data==x&labels==3,],na.rm=T)
}))

trajectory.confidence.bounds = data.frame(t=rep(NA,length(windows)*4),cluster=rep(NA,length(windows)*4),center=rep(NA,length(windows)*4),lower=rep(NA,length(windows)*4),upper=rep(NA,length(windows)*4))
for (i in 1:4) {
    indices = (i-1)*length(windows)+(1:length(windows))
    trajectory.confidence.bounds$t[indices] = windows
    trajectory.confidence.bounds$cluster[indices] = which(cluster.indices==i)
    trajectory.confidence.bounds$center[indices] = risk.score.centroids[i,]
    trajectory.confidence.bounds$lower[indices] = risk.score.centroids[i,] - risk.score.sds[i,]*1
    trajectory.confidence.bounds$upper[indices] = risk.score.centroids[i,] + risk.score.sds[i,]*1
}
trajectory.confidence.bounds$cluster=as.factor(trajectory.confidence.bounds$cluster)

#trajectory.confidence.bounds$line_a = rep(1,dim(trajectory.confidence.bounds)[1])
#trajectory.confidence.bounds$ribbon_a = rep(0.8,dim(trajectory.confidence.bounds)[1])

# ggplot(trajectory.confidence.bounds,aes(group=cluster,x=t,y=center,color=cluster))+geom_line(size=2,alpha=1)+ylab("SBP")+geom_ribbon(aes(ymin=lower,ymax=upper,fill=cluster),alpha=0.2)+
#     labs(fill="Cluster", color="Cluster")+guides(alpha=FALSE)+xlab("Time (hrs relative to detection)")


ggplot(trajectory.confidence.bounds,aes(group=cluster,x=t,y=center,color=cluster))+geom_line(size=2)+
    geom_ribbon(aes(ymin=lower,ymax=upper,fill=cluster),alpha=0.3)+
    scale_discrete_manual(aesthetics=c("color","fill"),name="Post-Prediction Clusters",breaks=c(1,2,3,4),labels=c("1 (High-risk)","2","3","4 (Low-risk)"),values=c("#f8766d","#7cae00","#00bfc4","#c77cff"))+
    labs(fill="Cluster",color="Cluster")+xlab("Time (hrs relative to early prediction)")+ylab("Systolic Blood Pressure (mmHg)")+
    theme(legend.justification=c(0,0),
          legend.position=c(0,0),
          legend.box.margin=margin(c(5,5,5,5)),
          legend.background = element_rect(fill=alpha('white',0.5)),
          axis.text.x = element_text(size=12),
          axis.title.x = element_text(size=14),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=14),
          legend.text = element_text(size=12),
          legend.title = element_text(size=14))