rm(list=ls())
source("src/R/eicu/load_reference_data_2.R")
load("data/eicu/extended.detection.trajectories.v3.rdata")
load("data/eicu/spec.clustering.detection.rdata")

# saveRDS(sc,file="post.detection.clusters.rds")

library(matrixStats)
library(ggplot2)

windows = seq(from=-12,to=30,by=1)

confusion = array(NA,dim=c(4,3))
for (i in 1:4) {
    for (j in 1:3) {
        confusion[i,j] = sum(sc@.Data==i&labels==j)
    }
}

ewts = c(10,13,15,30)

tp.fraction = apply(confusion,MARGIN=1,function(x) x[3]/sum(x))
cluster.indices = order(tp.fraction,decreasing=T)

# saveRDS(cluster.indices,file="post.detection.indices.rds")

trajectories = vector(mode="list",length=4)
for (i in 1:4) {
    # trajectories[[i]] = ext.trajectories[sc@.Data==cluster.indices[i],1:which(windows==ewts[i])]
    trajectories[[i]] = ext.trajectories[sc@.Data==cluster.indices[i],]
    
}

means = lapply(trajectories,function(x) colMeans(x,na.rm=T))
means = do.call(c,means)
sds = lapply(trajectories, function(x) colSds(x,na.rm=T))
sds = do.call(c,sds)

lowers = lapply(trajectories,function(x) apply(x, 2, function(y) quantile(y,0.16,na.rm=T)))
lowers = do.call(c,lowers)
uppers = lapply(trajectories,function(x) apply(x, 2, function(y) quantile(y,0.84,na.rm=T)))
uppers = do.call(c,uppers)

groups = lapply(1:4, function(x) rep(x,length(windows)))
groups = as.factor(do.call(c,groups))

pre.ewt = lapply(1:4, function(x) c(rep(T,which(windows==ewts[x])),rep(F,length(windows)-which(windows==ewts[x]))))
pre.ewt = do.call(c,pre.ewt)

cluster.windows = lapply(1:4, function(x) windows)
cluster.windows = do.call(c,cluster.windows)

data = data.frame(time=cluster.windows,mean=means,lower=lowers,upper=uppers,group=groups,pre.ewt=pre.ewt)

ggplot(data,aes(x=time,y=mean,color=group))+geom_line(size=2,aes(linetype=!pre.ewt))+
    geom_ribbon(aes(ymin=lower,ymax=upper,fill=group),alpha=0.3)+
    geom_hline(yintercept=threshold)+
    scale_discrete_manual(aesthetics=c("color","fill"),name="Post-Prediction Clusters",breaks=c(1,2,3,4),labels=c("1 (High-risk)","2","3","4 (Low-risk)"),values=c("#f8766d","#7cae00","#00bfc4","#c77cff"))+
    labs(fill="Cluster",color="Cluster")+xlab("Time (hrs relative to early prediction)")+ylab("Risk score")+guides(alpha=FALSE,linetype=FALSE)+
    theme(legend.justification=c(0,1),
          legend.position=c(0,1),
          legend.box.margin=margin(c(5,5,5,5)),
          axis.text.x = element_text(size=12),
          axis.title.x = element_text(size=14),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=14),
          legend.text = element_text(size=12),
          legend.title = element_text(size=14),
          legend.background = element_rect(alpha("white",0.5)))

