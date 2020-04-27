rm(list=ls())
source("src/R/eicu/load_reference_data_2.r")
source("src/R/eicu/functions/eval_carry_forward.R")

library(kernlab)
library(tidyr)
library(ggplot2)

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

num.cores = detectCores()
cluster = makeCluster(num.cores)
clusterExport(cluster,c("eval.early.prediction.timestamps.combined.rf.comorbidities","eval.table.with.sofa","eval.table.with.sofa.comorbidities","eval.carry.forward","eval.sum.in.past","eval.max.in.past","eval.interval","means","model","clinical.data","sampling.rate.data","lengths","has.sepsis","has.shock","preshock.sample","shock.onsets","comorbidities"))
clusterEvalQ(cluster,library(glmnet))
tic("Early prediction, parallel (all, with timestamps)")
shock.predictions = parLapply(cluster, 1:sum(has.shock), function(x) eval.early.prediction.timestamps.combined.rf.comorbidities(model,clinical.data[lengths>0][has.shock][[x]],sampling.rate.data[lengths>0][has.shock][[x]],comorbidities[,has.shock][,x],means))
toc()
stopCluster(cluster)


shock.patients = sapply(clinical.data[lengths>0][has.shock], function(x) x$subject.id)
shock.patient.indices = sapply(shock.patients, function(x) which(patient.result$patientunitstayid==x))
shock.hospital.discharge = patient.result$hospitaldischargeoffset[shock.patient.indices]
shock.icu.discharge = patient.result$unitdischargeoffset[shock.patient.indices]


# Generate time series
num.cores = detectCores()
cluster = makeCluster(num.cores)
clusterExport(cluster, c("eval.carry.forward","shock.predictions","shock.onsets"))
tic("Creating timeseries from risk scores")
risk.score.ts = parLapply(cluster, 1:sum(has.shock), function(x) {
    eval.carry.forward(shock.onsets[x]+(0:12)*60, shock.predictions[[x]]$timestamps, shock.predictions[[x]]$predictions)
})
toc()
stopCluster(cluster)



mortality.labels = has.mortality[has.shock]

x = risk.score.ts.table

tic("spectral clustering")
sc2 = specc(x,4)
toc()

saveRDS(sc2,file="post.shock.rds")

sc2.confusion = array(NA, dim=c(4,2))
sc2.confusion[1,1] = sum(sc2@.Data == 1 & mortality.labels)
sc2.confusion[1,2] = sum(sc2@.Data == 1 & !mortality.labels)
sc2.confusion[2,1] = sum(sc2@.Data == 2 & mortality.labels)
sc2.confusion[2,2] = sum(sc2@.Data == 2 & !mortality.labels)
sc2.confusion[3,1] = sum(sc2@.Data == 3 & mortality.labels)
sc2.confusion[3,2] = sum(sc2@.Data == 3 & !mortality.labels)
sc2.confusion[4,1] = sum(sc2@.Data == 4 & mortality.labels)
sc2.confusion[4,2] = sum(sc2@.Data == 4 & !mortality.labels)

tp.fraction = apply(sc2.confusion,MARGIN=1,function(x) x[1]/sum(x))
cluster.indices = order(tp.fraction,decreasing=T)

s <- function(x1, x2, alpha=1) {
    exp(- alpha * norm(as.matrix(x1-x2), type="F"))
}

make.similarity <- function(my.data, similarity) {
    N <- nrow(my.data)
    S <- matrix(rep(NA,N^2), ncol=N)
    for(i in 1:N) {
        for(j in 1:N) {
            S[i,j] <- similarity(my.data[i,], my.data[j,])
        }
    }
    S
}

S <- make.similarity(x, s)

make.affinity <- function(S, n.neighboors=2) {
    N <- length(S[,1])
    
    if (n.neighboors >= N) {  # fully connected
        A <- S
    } else {
        A <- matrix(rep(0,N^2), ncol=N)
        for(i in 1:N) { # for each line
            # only connect to those points with larger similarity 
            best.similarities <- sort(S[i,], decreasing=TRUE)[1:n.neighboors]
            for (s in best.similarities) {
                j <- which(S[i,] == s)
                A[i,j] <- S[i,j]
                A[j,i] <- S[i,j] # to make an undirected graph, ie, the matrix becomes symmetric
            }
        }
    }
    A  
}

A <- make.affinity(S, 3)  # use 3 neighboors (includes self)
D <- diag(apply(A, 1, sum))
U <- D - A
evL <- eigen(U, symmetric=TRUE)
plot(1:10, rev(evL$values)[1:10], log="y")
abline(v=4.5, col="red", lty=2)

# Visualize
windows = 0:12
risk.score.centroids = do.call(rbind,lapply(1:4, function(y) {
    colMeans(x[sc2@.Data==y,],na.rm=T)
}))
risk.score.sds = do.call(rbind,lapply(1:4, function(y) {
    colSds(x[sc2@.Data==y,],na.rm=T)
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

ggplot(trajectory.confidence.bounds,aes(group=cluster,x=t,y=center,color=cluster))+geom_line(size=2)+geom_ribbon(aes(ymin=lower,ymax=upper,alpha=0.3,fill=cluster))+
    labs(fill="Cluster", color="Cluster")+guides(alpha=FALSE)+xlab("Time (hrs relative to shock onset)")+ylab("Risk Score")

# Visualize ext
# ext.windows = -12:12
ext.windows = 0:48

num.cores = detectCores()
cluster = makeCluster(num.cores)
clusterExport(cluster, c("eval.carry.forward","shock.predictions","shock.onsets","ext.windows"))
tic("Creating timeseries from risk scores")
risk.score.ts.ext = parLapply(cluster, 1:sum(has.shock), function(x) {
    eval.carry.forward(shock.onsets[x]+(ext.windows)*60, shock.predictions[[x]]$timestamps, shock.predictions[[x]]$predictions)
})
toc()
stopCluster(cluster)

risk.score.ts.table = do.call(rbind,risk.score.ts)
risk.score.ts.ext.table = do.call(rbind,risk.score.ts.ext)

risk.score.centroids = do.call(rbind,lapply(1:4, function(y) {
    colMeans(risk.score.ts.ext.table[sc2@.Data==y,],na.rm=T)
}))
risk.score.sds = do.call(rbind,lapply(1:4, function(y) {
    colSds(risk.score.ts.ext.table[sc2@.Data==y,],na.rm=T)
}))

trajectory.confidence.bounds = data.frame(t=rep(NA,length(ext.windows)*4),cluster=rep(NA,length(ext.windows)*4),center=rep(NA,length(ext.windows)*4),lower=rep(NA,length(ext.windows)*4),upper=rep(NA,length(ext.windows)*4))
for (i in 1:4) {
    indices = (i-1)*length(ext.windows)+(1:length(ext.windows))
    trajectory.confidence.bounds$t[indices] = ext.windows
    trajectory.confidence.bounds$cluster[indices] = which(cluster.indices==i)
    trajectory.confidence.bounds$center[indices] = risk.score.centroids[i,]
    trajectory.confidence.bounds$lower[indices] = risk.score.centroids[i,] - risk.score.sds[i,]*1
    trajectory.confidence.bounds$upper[indices] = risk.score.centroids[i,] + risk.score.sds[i,]*1
}
trajectory.confidence.bounds$cluster=as.factor(trajectory.confidence.bounds$cluster)

ggplot(trajectory.confidence.bounds,aes(group=cluster,x=t,y=center,color=cluster))+geom_line(size=2)+geom_ribbon(aes(ymin=lower,ymax=upper,alpha=0.3,fill=cluster))+
    labs(fill="Cluster", color="Cluster")+guides(alpha=FALSE)+xlab("Time (hrs relative to shock onset)")+ylab("Risk Score")

