rm(list=ls())
source("load_reference_data_2.R")
load("spec.clustering.detection.rdata")

library(kernlab)
library(ggplot2)
library(tidyr)
library(tictoc)
library(class)
library(entropy)

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

windows = seq(from=-12,to=12,by=1)

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

trajectories = rbind(nonsepsis.post.table,nonshock.post.table,shock.post.table)
labels = c(rep(1,dim(nonsepsis.post.table)[1]),rep(2,dim(nonshock.post.table)[1]),rep(3,dim(shock.post.table)[1]))


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

confusion = array(NA,dim=c(4,3))
for (i in 1:4) {
    for (j in 1:3) {
        confusion[i,j] = sum(sc@.Data==i&labels==j)
    }
}
tp.fraction = apply(confusion,MARGIN=1,function(x) x[3]/sum(x))
cluster.indices = order(tp.fraction,decreasing=T)

###
risk.score.centroids = do.call(rbind,lapply(1:4, function(x) {
    colMeans(trajectories[sc@.Data==x,],na.rm=T)
}))
risk.score.sds = do.call(rbind,lapply(1:4, function(x) {
    colSds(trajectories[sc@.Data==x,],na.rm=T)
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

ggplot(trajectory.confidence.bounds,aes(group=cluster,x=t,y=center,color=cluster))+geom_line(size=2)+
    geom_ribbon(aes(ymin=lower,ymax=upper,fill=cluster),alpha=0.3)+
    geom_hline(yintercept=threshold)+
    scale_color_manual(breaks=c(1,2,3,4),labels=c("1 (High-risk)","2","3","4 (Low-risk)"),values=c("#f8766d","#7cae00","#00bfc4","#c77cff"))+
    scale_fill_manual(breaks=c(1,2,3,4),labels=c("1 (High-risk)","2","3","4 (Low-risk)"),values=c("#f8766d","#7cae00","#00bfc4","#c77cff"))+
    labs(fill="Cluster",color="Cluster")+xlab("Time (hrs relative to early prediction)")+ylab("Risk score")+guides(alpha=FALSE)+
    theme(legend.justification=c(0,1),
          legend.position=c(0,1),
          legend.box.margin=margin(c(5,5,5,5)))
###


nonsepsis.hr = t(sapply(nonsepsis.post.detection.data, function(x) x$hr))
nonshock.hr = t(sapply(nonshock.post.detection.data, function(x) x$hr))
shock.hr = t(sapply(shock.post.detection.data, function(x) x$hr))
table.hr = rbind(nonsepsis.hr,nonshock.hr,shock.hr)

nonsepsis.lact = t(sapply(nonsepsis.post.detection.data, function(x) x$lact))
nonshock.lact = t(sapply(nonshock.post.detection.data, function(x) x$lact))
shock.lact = t(sapply(shock.post.detection.data, function(x) x$lact))
table.lact = rbind(nonsepsis.lact,nonshock.lact,shock.lact)

nonsepsis.sbp = t(sapply(nonsepsis.post.detection.data, function(x) x$sbp))
nonshock.sbp = t(sapply(nonshock.post.detection.data, function(x) x$sbp))
shock.sbp = t(sapply(shock.post.detection.data, function(x) x$sbp))
table.sbp = rbind(nonsepsis.sbp,nonshock.sbp,shock.sbp)

high.risk = sc@.Data == cluster.indices[1]
cluster.2 = sc@.Data == cluster.indices[2]
cluster.3 = sc@.Data == cluster.indices[3]
low.risk = sc@.Data == cluster.indices[4]

# Compute KL divergence between high-risk and low-risk for each window for risk score, hr, lact, sbp
risk.kl.divergences = sapply(1:length(windows), function(x) {
    range = c(min(trajectories[,x],na.rm=T),max(trajectories[,x],na.rm=T))
    # Pad with epsilon 1
    low.risk.d = discretize(trajectories[low.risk&!is.na(trajectories[,x]),x],numBins=20,r=range)+1
    high.risk.d = discretize(trajectories[high.risk&!is.na(trajectories[,x]),x],numBins=20,r=range)+1
    
    # KL.empirical(low.risk.d,high.risk.d)
	KL.empirical(high.risk.d,low.risk.d)
})

# Repeat for each of the physiological variables
hr.kl.divergences = sapply(1:length(windows), function(x) {
    range = c(min(table.hr[,x],na.rm=T),max(table.hr[,x],na.rm=T))
    # Pad with epsilon 1
    low.risk.d = discretize(table.hr[low.risk&!is.na(table.hr[,x]),x],numBins=20,r=range)+1
    high.risk.d = discretize(table.hr[high.risk&!is.na(table.hr[,x]),x],numBins=20,r=range)+1
    
    # KL.empirical(low.risk.d,high.risk.d)
    KL.empirical(high.risk.d,low.risk.d)
})

lact.kl.divergences = sapply(1:length(windows), function(x) {
    range = c(min(table.lact[,x],na.rm=T),max(table.lact[,x],na.rm=T))
    # Pad with epsilon 1
    low.risk.d = discretize(table.lact[low.risk&!is.na(table.lact[,x]),x],numBins=20,r=range)+1
    high.risk.d = discretize(table.lact[high.risk&!is.na(table.lact[,x]),x],numBins=20,r=range)+1
    
    # KL.empirical(low.risk.d,high.risk.d)
    KL.empirical(high.risk.d,low.risk.d)
})

sbp.kl.divergences = sapply(1:length(windows), function(x) {
    range = c(min(table.sbp[,x],na.rm=T),max(table.sbp[,x],na.rm=T))
    # Pad with epsilon 1
    low.risk.d = discretize(table.sbp[low.risk&!is.na(table.sbp[,x]),x],numBins=20,r=range)+1
    high.risk.d = discretize(table.sbp[high.risk&!is.na(table.sbp[,x]),x],numBins=20,r=range)+1
    
    # KL.empirical(low.risk.d,high.risk.d)
    KL.empirical(high.risk.d,low.risk.d)
})

### Generate visualizations
qplot(windows,risk.kl.divergences,geom="line")+ylab("KL Divergence")+
    xlab("Time relative to early prediction (hrs)")+
    ggtitle("Risk score")

qplot(windows,hr.kl.divergences,geom="line")+ylab("KL Divergence")+
    xlab("Time relative to early prediction (hrs)")+
    ggtitle("Heart rate")

qplot(windows,lact.kl.divergences,geom="line")+ylab("KL Divergence")+
    xlab("Time relative to early prediction (hrs)")+
    ggtitle("Lactate")

qplot(windows,sbp.kl.divergences,geom="line")+ylab("KL Divergence")+
    xlab("Time relative to early prediction (hrs)")+
    ggtitle("Systolic blood pressure")

### Combined visualization
data = data.frame(t=windows,risk=risk.kl.divergences,hr=hr.kl.divergences,lact=lact.kl.divergences,sbp=sbp.kl.divergences)
data.formatted = gather(data,key="key",value="value",-t)

ggplot(data.formatted,aes(x=t,y=value,color=key))+geom_line(size=1)+ylab("KL Divergence")+xlab("Time (hrs relative to early prediction)")

#

lact.kl.divergences = sapply(1:length(windows), function(x) {
    range = c(min(table.lact[,x],na.rm=T),max(table.lact[,x],na.rm=T))
    # Pad with epsilon 1
    low.risk.d = discretize(table.lact[low.risk&!is.na(table.lact[,x]),x],numBins=20,r=range)+1
    high.risk.d = discretize(table.lact[high.risk&!is.na(table.lact[,x]),x],numBins=20,r=range)+1
    
    # KL.empirical(low.risk.d,high.risk.d)
    KL.empirical(high.risk.d,low.risk.d)
})

qplot(windows,lact.kl.divergences,geom="line")+ylab("KL Divergence")+
    xlab("Time relative to early prediction (hrs)")+
    ggtitle("Lactate (high-risk vs low-risk)")
