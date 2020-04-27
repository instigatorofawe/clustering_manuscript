rm(list=ls())
library(kernlab)
library(ggplot2)
library(tidyr)

source("src/R/eicu/load_reference_data_2.R")
source("src/R/eicu/functions/eval_carry_forward.R")

abx.data = readRDS("data/eicu/eicu.abx.data.rds")

abx.subjects = unique(abx.data$patientunitstayid)
abx.times = sapply(abx.subjects, function(x) {
    current.abx.data = abx.data[abx.data$patientunitstayid==x,]
    if (dim(current.abx.data)[1]==0) {
        return(NA)
    } else {
        return(min(current.abx.data$drugorderoffset))
    }
})

patient.weights = sapply(clinical.data, function(x) patient.result$admissionweight[which(patient.result$patientunitstayid==x$subject.id)])
fluid.data = readRDS("data/eicu/fluid.data.rds")
fluid.subjects = sapply(fluid.data, function(x) x$subject.id)

#windows = seq(from=0,to=12,by=1)
#ext.windows = seq(from=-12,to=48,by=1)
ext.windows = seq(from=-12,to=12,by=1)


nonsepsis.maxes = sapply(nonsepsis.predictions, function(x) max(x$predictions))
has.detection.nonsepsis = nonsepsis.maxes >= threshold
nonshock.maxes = sapply(nonshock.predictions, function(x) max(x$predictions))
has.detection.nonshock = nonshock.maxes >= threshold
shock.maxes = sapply(shock.predictions, function(x) max(x$predictions))
has.detection.shock = shock.maxes >= threshold

has.abx.nonsepsis = sapply(clinical.data[lengths>0][!has.sepsis][has.detection.nonsepsis], function(x) is.element(x$subject.id,abx.subjects))
abx.times.nonsepsis = sapply(clinical.data[lengths>0][!has.sepsis][has.detection.nonsepsis][has.abx.nonsepsis], function(x) {
    abx.times[which(abx.subjects==x$subject.id)]
})

has.abx.nonshock = sapply(clinical.data[lengths>0][has.sepsis&!has.shock][has.detection.nonshock], function(x) is.element(x$subject.id,abx.subjects))
abx.times.nonshock = sapply(clinical.data[lengths>0][has.sepsis&!has.shock][has.detection.nonshock][has.abx.nonshock], function(x) {
    abx.times[which(abx.subjects==x$subject.id)]
})

has.abx.shock = sapply(clinical.data[lengths>0][has.shock][has.detection.shock], function(x) is.element(x$subject.id,abx.subjects))
abx.times.shock = sapply(clinical.data[lengths>0][has.shock][has.detection.shock][has.abx.shock], function(x) {
    abx.times[which(abx.subjects==x$subject.id)]
})

## Reformat abx times

abx.times.nonsepsis.alt = rep(NA,length(has.abx.nonsepsis))
abx.times.nonsepsis.alt[has.abx.nonsepsis] = abx.times.nonsepsis

abx.times.nonshock.alt = rep(NA,length(has.abx.nonshock))
abx.times.nonshock.alt[has.abx.nonshock] = abx.times.nonshock

abx.times.shock.alt = rep(NA,length(has.abx.shock))
abx.times.shock.alt[has.abx.shock] = abx.times.shock

##

shock.vasopressor.times = sapply(1:sum(has.detection.shock), function(x) {
    current.data = sofa.scores[lengths>0][has.shock][has.detection.shock][[x]]
    return(current.data$timestamps[min(which(current.data$vasopressors))])
})

shock.fluids = sapply(1:sum(has.detection.shock), function(x) {
    current.subject = clinical.data[lengths>0][has.shock][has.detection.shock][[x]]$subject.id
    if (!any(fluid.subjects==current.subject)) {
        return(rep(NA,length(windows)))
    }
    
    fluid.index = which(fluid.subjects==current.subject)
    sapply(sofa.scores[lengths>0][has.shock][has.detection.shock][[x]]$timestamps, function(x) {
        sum(fluid.data[[fluid.index]]$values[fluid.data[[fluid.index]]$timestamps<=x&fluid.data[[fluid.index]]$timestamps>=x-180])
    })
},simplify=F)

shock.urine = sapply(1:sum(has.detection.shock), function(x) {
    current.data = clinical.data[lengths>0][has.shock][has.detection.shock][[x]]
    if (is.null(current.data$urine)) {
        return(NA)
    } else {
        return(eval.sum.in.past(sofa.scores[lengths>0][has.shock][has.detection.shock][[x]]$timestamps,current.data$urine$timestamps,current.data$urine$values,180))
    }
},simplify=F)
shock.cvp = sapply(1:sum(has.detection.shock), function(x) {
    current.data = clinical.data[lengths>0][has.shock][has.detection.shock][[x]]
    if (is.null(current.data$cvp)) {
        return(NA)
    } else {
        return(eval.carry.forward(sofa.scores[lengths>0][has.shock][has.detection.shock][[x]]$timestamps,current.data$cvp$timestamps,current.data$cvp$values))
    }
},simplify=F)

shock.fluid.resusc = mapply(function(fluids,urine,cvp,weight) (fluids/weight>=30) | (cvp>=8&cvp<=12) | (urine/weight>0.5), shock.fluids, shock.urine, shock.cvp, patient.weights[lengths>0][has.shock][has.detection.shock])
shock.fluid.resusc.times = mapply(function(a,b) a$timestamps[min(which(b))], sofa.scores[lengths>0][has.shock][has.detection.shock], shock.fluid.resusc)


nonshock.vasopressor.times = sapply(1:sum(has.detection.nonshock), function(x) {
    current.data = sofa.scores[lengths>0][has.sepsis&!has.shock][has.detection.nonshock][[x]]
    return(current.data$timestamps[min(which(current.data$vasopressors))])
})


nonshock.fluids = sapply(1:sum(has.detection.nonshock), function(x) {
    current.subject = clinical.data[lengths>0][has.sepsis&!has.shock][has.detection.nonshock][[x]]$subject.id
    if (!any(fluid.subjects==current.subject)) {
        return(rep(NA,length(windows)))
    }
    
    fluid.index = which(fluid.subjects==current.subject)
    sapply(sofa.scores[lengths>0][has.sepsis&!has.shock][has.detection.nonshock][[x]]$timestamps, function(x) {
        sum(fluid.data[[fluid.index]]$values[fluid.data[[fluid.index]]$timestamps<=x&fluid.data[[fluid.index]]$timestamps>=x-180])
    })
},simplify=F)

nonshock.urine = sapply(1:sum(has.detection.nonshock), function(x) {
    current.data = clinical.data[lengths>0][has.sepsis&!has.shock][has.detection.nonshock][[x]]
    if (is.null(current.data$urine)) {
        return(NA)
    } else {
        return(eval.sum.in.past(sofa.scores[lengths>0][has.sepsis&!has.shock][has.detection.nonshock][[x]]$timestamps,current.data$urine$timestamps,current.data$urine$values,180))
    }
},simplify=F)

nonshock.cvp = sapply(1:sum(has.detection.nonshock), function(x) {
    current.data = clinical.data[lengths>0][has.sepsis&!has.shock][has.detection.nonshock][[x]]
    if (is.null(current.data$cvp)) {
        return(NA)
    } else {
        return(eval.carry.forward(sofa.scores[lengths>0][has.sepsis&!has.shock][has.detection.nonshock][[x]]$timestamps,current.data$cvp$timestamps,current.data$cvp$values))
    }
},simplify=F)

nonshock.fluid.resusc = mapply(function(fluids,urine,cvp,weight) (fluids/weight>=30) | (cvp>=8&cvp<=12) | (urine/weight>0.5), nonshock.fluids, nonshock.urine, nonshock.cvp, patient.weights[lengths>0][has.sepsis&!has.shock][has.detection.nonshock])
nonshock.fluid.resusc.times = mapply(function(a,b) a$timestamps[min(which(b))], sofa.scores[lengths>0][has.sepsis&!has.shock][has.detection.nonshock], nonshock.fluid.resusc)


nonsepsis.vasopressor.times = sapply(1:sum(has.detection.nonsepsis), function(x) {
    current.data = sofa.scores[lengths>0][!has.sepsis][has.detection.nonsepsis][[x]]
    return(current.data$timestamps[min(which(current.data$vasopressors))])
})


nonsepsis.fluids = sapply(1:sum(has.detection.nonsepsis), function(x) {
    current.subject = clinical.data[lengths>0][!has.sepsis][has.detection.nonsepsis][[x]]$subject.id
    if (!any(fluid.subjects==current.subject)) {
        return(rep(NA,length(windows)))
    }
    
    fluid.index = which(fluid.subjects==current.subject)
    sapply(sofa.scores[lengths>0][!has.sepsis][has.detection.nonsepsis][[x]]$timestamps, function(x) {
        sum(fluid.data[[fluid.index]]$values[fluid.data[[fluid.index]]$timestamps<=x&fluid.data[[fluid.index]]$timestamps>=x-180])
    })
},simplify=F)

nonsepsis.urine = sapply(1:sum(has.detection.nonsepsis), function(x) {
    current.data = clinical.data[lengths>0][!has.sepsis][has.detection.nonsepsis][[x]]
    if (is.null(current.data$urine)) {
        return(NA)
    } else {
        return(eval.sum.in.past(sofa.scores[lengths>0][!has.sepsis][has.detection.nonsepsis][[x]]$timestamps,current.data$urine$timestamps,current.data$urine$values,180))
    }
},simplify=F)

nonsepsis.cvp = sapply(1:sum(has.detection.nonsepsis), function(x) {
    current.data = clinical.data[lengths>0][!has.sepsis][has.detection.nonsepsis][[x]]
    if (is.null(current.data$cvp)) {
        return(NA)
    } else {
        return(eval.carry.forward(sofa.scores[lengths>0][!has.sepsis][has.detection.nonsepsis][[x]]$timestamps,current.data$cvp$timestamps,current.data$cvp$values))
    }
},simplify=F)

nonsepsis.fluid.resusc = mapply(function(fluids,urine,cvp,weight) (fluids/weight>=30) | (cvp>=8&cvp<=12) | (urine/weight>0.5), nonsepsis.fluids, nonsepsis.urine, nonsepsis.cvp, patient.weights[lengths>0][!has.sepsis][has.detection.nonsepsis])
nonsepsis.fluid.resusc.times = mapply(function(a,b) a$timestamps[min(which(b))], sofa.scores[lengths>0][!has.sepsis][has.detection.nonsepsis], nonsepsis.fluid.resusc)

###

# Determine which comes first
nonsepsis.intervention.times = cbind(nonsepsis.fluid.resusc.times,nonsepsis.vasopressor.times,abx.times.nonsepsis.alt)
nonshock.intervention.times = cbind(nonshock.fluid.resusc.times,nonshock.vasopressor.times,abx.times.nonshock.alt)
shock.intervention.times = cbind(shock.fluid.resusc.times,shock.vasopressor.times,abx.times.shock.alt)

nonsepsis.first.intervention = apply(nonsepsis.intervention.times, 1, function(x) {
    if (any(!is.na(x))) {
        return(which.min(x))
    } else {
        return(NA)
    }
})
nonshock.first.intervention = apply(nonshock.intervention.times, 1, function(x) {
    if (any(!is.na(x))) {
        return(which.min(x))
    } else {
        return(NA)
    }
})
shock.first.intervention = apply(shock.intervention.times, 1, function(x) {
    if (any(!is.na(x))) {
        return(which.min(x))
    } else {
        return(NA)
    }
})
nonsepsis.has.intervention = !is.na(nonsepsis.first.intervention)
nonshock.has.intervention = !is.na(nonshock.first.intervention)
shock.has.intervention = !is.na(shock.first.intervention)

nonsepsis.first.intervention.time = apply(nonsepsis.intervention.times[nonsepsis.has.intervention,], 1, function(x) min(x,na.rm=T))
nonshock.first.intervention.time = apply(nonshock.intervention.times[nonshock.has.intervention,], 1, function(x) min(x,na.rm=T))
shock.first.intervention.time = apply(shock.intervention.times[shock.has.intervention,], 1, function(x) min(x,na.rm=T))

## Generate trajectories
intervention.trajectories.nonsepsis = mapply(function(x,y) {
    eval.timestamps = windows*60 + y
    return(eval.carry.forward(eval.timestamps,x$timestamps,x$predictions))
},nonsepsis.predictions[has.detection.nonsepsis][nonsepsis.has.intervention],nonsepsis.first.intervention.time,SIMPLIFY=F)
nonsepsis.trajectories = do.call(rbind,intervention.trajectories.nonsepsis)

intervention.trajectories.nonshock = mapply(function(x,y) {
    eval.timestamps = windows*60 + y
    return(eval.carry.forward(eval.timestamps,x$timestamps,x$predictions))
},nonshock.predictions[has.detection.nonshock][nonshock.has.intervention],nonshock.first.intervention.time,SIMPLIFY=F)
nonshock.trajectories = do.call(rbind,intervention.trajectories.nonshock)

intervention.trajectories.shock = mapply(function(x,y) {
    eval.timestamps = windows*60 + y
    return(eval.carry.forward(eval.timestamps,x$timestamps,x$predictions))
},shock.predictions[has.detection.shock][shock.has.intervention],shock.first.intervention.time,SIMPLIFY=F)
shock.trajectories = do.call(rbind,intervention.trajectories.shock)

trajectories = rbind(nonsepsis.trajectories,nonshock.trajectories,shock.trajectories)
means = colMeans(trajectories,na.rm=T)
for (i in 1:dim(trajectories)[2]) {
    trajectories[is.na(trajectories[,i]),i] = means[i]
}

## Extended trajectories
intervention.trajectories.ext.nonsepsis = mapply(function(x,y) {
    eval.timestamps = ext.windows*60 + y
    return(eval.carry.forward(eval.timestamps,x$timestamps,x$predictions))
},nonsepsis.predictions[has.detection.nonsepsis][nonsepsis.has.intervention],nonsepsis.first.intervention.time,SIMPLIFY=F)
nonsepsis.ext.trajectories = do.call(rbind,intervention.trajectories.ext.nonsepsis)

intervention.trajectories.ext.nonshock = mapply(function(x,y) {
    eval.timestamps = ext.windows*60 + y
    return(eval.carry.forward(eval.timestamps,x$timestamps,x$predictions))
},nonshock.predictions[has.detection.nonshock][nonshock.has.intervention],nonshock.first.intervention.time,SIMPLIFY=F)
nonshock.ext.trajectories = do.call(rbind,intervention.trajectories.ext.nonshock)

intervention.trajectories.ext.shock = mapply(function(x,y) {
    eval.timestamps = ext.windows*60 + y
    return(eval.carry.forward(eval.timestamps,x$timestamps,x$predictions))
},shock.predictions[has.detection.shock][shock.has.intervention],shock.first.intervention.time,SIMPLIFY=F)
shock.ext.trajectories = do.call(rbind,intervention.trajectories.ext.shock)

ext.trajectories = rbind(nonsepsis.ext.trajectories,nonshock.ext.trajectories,shock.ext.trajectories)
ext.means = colMeans(ext.trajectories,na.rm=T)
for (i in 1:dim(ext.trajectories)[2]) {
    ext.trajectories[is.na(ext.trajectories[,i]),i] = ext.means[i]
}

## Do clustering
# s <- function(x1, x2, alpha=1) {
#     exp(- alpha * norm(as.matrix(x1-x2), type="F"))
# }

# make.similarity <- function(my.data, similarity) {
#     N <- nrow(my.data)
#     S <- matrix(rep(NA,N^2), ncol=N)
#     diag(S) = 1
#     for(i in 1:(N-1)) {
#         for (j in (i+1):N) {
#             S[i,j] = similarity(my.data[i,],my.data[j,])
#             S[j,i] = S[i,j]
#         }
#     }
#     S
# }

# make.affinity <- function(S, n.neighboors=2) {
#     N <- length(S[,1])

#     if (n.neighboors >= N) {  # fully connected
#         A <- S
#     } else {
#         A <- matrix(rep(0,N^2), ncol=N)
#         for(i in 1:N) { # for each line
#             # only connect to those points with larger similarity 
#             best.similarities <- sort(S[i,], decreasing=TRUE)[1:n.neighboors]
#             for (s in best.similarities) {
#                 j <- which(S[i,] == s)
#                 A[i,j] <- S[i,j]
#                 A[j,i] <- S[i,j] # to make an undirected graph, ie, the matrix becomes symmetric
#             }
#         }
#     }
#     A  
# }
# tic()
# S <- make.similarity(trajectories, s)
# toc()
# A <- make.affinity(S, 3)  # use 3 neighboors (includes self)
# D <- diag(apply(A, 1, sum))
# U <- D - A
# evL <- eigen(U, symmetric=TRUE)
# plot(1:10, rev(evL$values)[1:10], log="y")
# abline(v=5.5, col="red", lty=2)

# tic()
# sc = specc(trajectories,centers=5)
# toc()

# saveRDS(sc,file="sc.first.intervention.rds")

sc = readRDS("data/eicu/sc.first.intervention.rds")

## Visualize
labels = c(rep(0,sum(nonsepsis.has.intervention)+sum(nonshock.has.intervention)),rep(1,sum(shock.has.intervention)))

has.mortality.nonsepsis = has.mortality[!has.sepsis][has.detection.nonsepsis][nonsepsis.has.intervention]
has.mortality.nonshock = has.mortality[has.sepsis&!has.shock][has.detection.nonshock][nonshock.has.intervention]
has.mortality.shock = has.mortality[has.shock][has.detection.shock][shock.has.intervention]
mortality.labels = c(has.mortality.nonsepsis,has.mortality.nonshock,has.mortality.shock)

sc.confusion = array(NA, dim=c(5,2))
for (i in 1:5) {
    sc.confusion[i,1] = sum(sc@.Data == i & labels)
    sc.confusion[i,2] = sum(sc@.Data == i & !labels)
}

# sc.confusion = array(NA, dim=c(5,2))
# for (i in 1:5) {
#     sc.confusion[i,1] = sum(sc@.Data == i & mortality.labels)
#     sc.confusion[i,2] = sum(sc@.Data == i & !mortality.labels)
# }


tp.fraction = apply(sc.confusion,MARGIN=1,function(x) x[1]/sum(x))
cluster.indices = order(tp.fraction,decreasing=T)

risk.score.centroids = do.call(rbind,lapply(1:5, function(y) {
    colMeans(ext.trajectories[sc@.Data==y,],na.rm=T)
}))
risk.score.sds = do.call(rbind,lapply(1:5, function(y) {
    colSds(ext.trajectories[sc@.Data==y,],na.rm=T)
}))

trajectory.confidence.bounds = data.frame(t=rep(NA,length(ext.windows)*5),cluster=rep(NA,length(ext.windows)*5),center=rep(NA,length(ext.windows)*5),lower=rep(NA,length(ext.windows)*5),upper=rep(NA,length(ext.windows)*5))
for (i in 1:5) {
    indices = (i-1)*length(ext.windows)+(1:length(ext.windows))
    trajectory.confidence.bounds$t[indices] = ext.windows
    trajectory.confidence.bounds$cluster[indices] = which(cluster.indices==i)
    trajectory.confidence.bounds$center[indices] = risk.score.centroids[i,]
    trajectory.confidence.bounds$lower[indices] = risk.score.centroids[i,] - risk.score.sds[i,]*1
    trajectory.confidence.bounds$upper[indices] = risk.score.centroids[i,] + risk.score.sds[i,]*1
}
trajectory.confidence.bounds$cluster=as.factor(trajectory.confidence.bounds$cluster)

ggplot(trajectory.confidence.bounds,aes(group=cluster,x=t,y=center,color=cluster))+geom_line(size=2)+
    geom_ribbon(aes(ymin=lower,ymax=upper,fill=cluster),alpha=0.3)+
    geom_hline(yintercept=threshold)+
    scale_discrete_manual(aesthetics=c("color","fill"),name="Post-Intervention Clusters",breaks=c(1,2,3,4,5),labels=c("1 (High-risk)","2","3","4","5 (Low-risk)"),values=c("#f8766d","#a3a500","#00bf7d","#00b0f6","#e76bf3"))+
    labs(fill="Cluster", color="Cluster")+guides(alpha=FALSE)+xlab("Time (hrs relative to first intervention)")+ylab("Risk Score")+
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

### Analyses
# first.interventions = c(nonsepsis.first.intervention[nonsepsis.has.intervention],nonshock.first.intervention[nonshock.has.intervention],shock.first.intervention[shock.has.intervention])

# # Examine concordance
# sc.post.detection = readRDS("post.detection.clusters.rdata")
# post.cluster.indices = readRDS("post.detection.indices.rdata")

# post.intervention.clusters = sc.post.detection@.Data[c(nonsepsis.has.intervention,nonshock.has.intervention,shock.has.intervention)]
