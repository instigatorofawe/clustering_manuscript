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
#shock.detection.times = sapply(shock.predictions[has.detection.shock], function(x) x$timestamps[min(which(x$predictions>=threshold))])

shock.detection.times = sapply(shock.predictions[has.detection.shock], function(x) x$timestamps[which.min(x$predictions>=threshold)])
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

windows = seq(from=0,to=12,by=1)

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

s <- function(x1, x2, alpha=1) {
    exp(- alpha * norm(as.matrix(x1-x2), type="F"))
}

make.similarity <- function(my.data, similarity) {
    N <- nrow(my.data)
    S <- matrix(rep(NA,N^2), ncol=N)
    diag(S) = 1
    for(i in 1:(N-1)) {
        for (j in (i+1):N) {
            S[i,j] = similarity(my.data[i,],my.data[j,])
            S[j,i] = S[i,j]
        }
    }
    S
}

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
tic()
S <- make.similarity(trajectories, s)
toc()
A <- make.affinity(S, 3)  # use 3 neighboors (includes self)
D <- diag(apply(A, 1, sum))
U <- D - A
evL <- eigen(U, symmetric=TRUE)
plot(1:10, rev(evL$values)[1:10], log="y")
abline(v=4.5, col="red", lty=2)

tic()
sc = specc(trajectories,centers=4)
toc()
tic()
sc.5 = specc(trajectories,centers=5)
toc()
tic()
sc.3 = specc(trajectories,centers=3)
toc()

save(sc,sc.3,sc.5,file="spec.clustering.detection.rdata")
#tic()
#ts = diss(trajectories,"DTWARP")
#toc()

#hc = hclust(ts)
#dtwc = cutree(hc,k=4)
#confusion = array(NA,dim=c(4,3))
#for (i in 1:4) {
#    for (j in 1:3) {
#        confusion[i,j] = sum(dtwc==i&labels==j)
#    }
#}

#centroids = as.data.frame(cbind(0:12,sapply(1:4,function(x) colMeans(trajectories[dtwc==x,]))))
#colnames(centroids) = c("t",1,2,3,4)
#centroids.processed = gather(data=centroids,"centroid","value",-t)
#ggplot(centroids.processed,aes(x=t,y=value,group=centroid,color=centroid))+geom_line(size=1)


# Visualize
centroids = as.data.frame(cbind(windows,t(sc@centers)))
colnames(centroids) = c("t",1,2,3,4)
centroids.processed = gather(data=centroids,"centroid","value",-t)
ggplot(centroids.processed,aes(x=t,y=value,group=centroid,color=centroid))+geom_line(size=1)


confusion = array(NA,dim=c(4,3))
for (i in 1:4) {
    for (j in 1:3) {
        confusion[i,j] = sum(sc@.Data==i&labels==j)
    }
}


#windows = 0:12

# Compute means, confidence bounds for each cluster
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
    trajectory.confidence.bounds$cluster[indices] = i
    trajectory.confidence.bounds$center[indices] = sc@centers[i,]
    trajectory.confidence.bounds$lower[indices] = sc@centers[i,] - risk.score.sds[i,]*1
    trajectory.confidence.bounds$upper[indices] = sc@centers[i,] + risk.score.sds[i,]*1
}
trajectory.confidence.bounds$cluster=as.factor(trajectory.confidence.bounds$cluster)

# Construct confidence bound plot
ggplot(trajectory.confidence.bounds,aes(group=cluster,x=t,y=center,color=cluster))+geom_line(size=2)+geom_ribbon(aes(ymin=lower,ymax=upper,alpha=0.3,fill=cluster))


# K-nearest neighbors classifier
# Arbitrarily choose k = 5
trajectories.train = rbind(nonsepsis.post.table[nonsepsis.sample[has.detection.nonsepsis],],
    nonshock.post.table[nonshock.sample[has.detection.nonshock],],
    shock.post.table[preshock.sample[has.detection.shock],])

trajectories.test = rbind(nonsepsis.post.table[!nonsepsis.sample[has.detection.nonsepsis],],
    nonshock.post.table[!nonshock.sample[has.detection.nonshock],],
    shock.post.table[!preshock.sample[has.detection.shock],])

labels.train = c(rep(1,sum(nonsepsis.sample[has.detection.nonsepsis])),
    rep(2,sum(nonshock.sample[has.detection.nonshock])),
    rep(3,sum(preshock.sample[has.detection.shock])))

labels.test = c(rep(1,sum(!nonsepsis.sample[has.detection.nonsepsis])),
    rep(2,sum(!nonshock.sample[has.detection.nonshock])),
    rep(3,sum(!preshock.sample[has.detection.shock])))

tic()
sc.train = specc(trajectories.train, centers=4)
toc()

fpr = array(NA,dim=c(13,4))
for (i in 1:13) {
    knn.prediction = knn(trajectories.train[,1:i,drop=F],trajectories.test[,1:i,drop=F],sc.train@.Data,5)
    for (j in 1:4) {
        fpr[i,j] = sum(knn.prediction==j&labels.test!=3)/sum(knn.prediction==j)
    }
}
data = data.frame(t=0:12,fpr)
colnames(data) = c("t","1","2","3","4")
reshaped.data = gather(data,"cluster","value",-t)
ggplot(reshaped.data,aes(x=t,y=value,group=cluster,color=cluster))+geom_line()+xlab("Hours of post-detection data")+ylab("Proportion FP")

# npv = rep(NA,13)
# size.knn = rep(NA,13)

# for (i in 1:13) {
#     knn.prediction = knn(trajectories.train[,1:i,drop=F],trajectories.test[,1:i,drop=F],sc.train@.Data,5)
#     npv[i] = sum(knn.prediction==3&labels.test!=3)/sum(knn.prediction==3)
#     size.knn[i] = sum(knn.prediction==3)
# }


# npv.rf = rep(NA,13)
# size.rf = rep(NA,13)

# train.data = data.frame(trajectories.train,y=as.factor(sc.train@.Data))
# test.data = data.frame(trajectories.test,y=as.factor(labels.test))

# for (i in 1:13) {
   
#     rf = ranger(y~.,train.data[,c(1:i,14)])
#     z = predict(rf,test.data[,c(1:i,14)])
#     npv.rf[i] = sum(z$predictions==3&labels.test!=3)/sum(z$predictions==3)
#     size.rf[i] = sum(z$predictions==3)
# }


fpr = array(NA,dim=c(13,4))

train.data = data.frame(trajectories.train,y=as.factor(sc.train@.Data))
test.data = data.frame(trajectories.test,y=as.factor(labels.test))

for (i in 1:13) {
   
    rf = ranger(y~.,train.data[,c(1:i,14)])
    z = predict(rf,test.data[,c(1:i,14)])
    for (j in 1:4) {
        fpr[i,j] = sum(z$predictions==j&labels.test!=3)/sum(z$predictions==j)

    }   
}
data = data.frame(t=0:12,fpr)
colnames(data) = c("t","1","2","3","4")
reshaped.data = gather(data,"cluster","value",-t)
ggplot(reshaped.data,aes(x=t,y=value,group=cluster,color=cluster))+geom_line()+xlab("Hours of post-detection data")+ylab("Proportion FP")




plot.data = data.frame(t=rep(0:12,2),npv=c(npv,npv.rf),group=c(rep("knn",13),rep("rf",13)))
ggplot(plot.data,aes(x=t,y=npv,group=group,color=group))+geom_line()+xlab("Hours of data post-detection")+ylab("Proportion false positive of removed set")

plot.data = data.frame(t=rep(0:12,2),size=c(size.knn,size.rf),group=c(rep("knn",13),rep("rf",13)))
ggplot(plot.data,aes(x=t,y=size,group=group,color=group))+geom_line()+xlab("Hours of data post-detection")+ylab("Size of removed set")

qplot(0:12,npv) + xlab("Hours of data post-detection") + ylab("Proportion false positive of removed set")


# Build trajectories: post-icu admission
nonsepsis.post.icu.trajectories = lapply(nonsepsis.predictions, function(x) {
    eval.timestamps = (0:12)*60
    return(eval.carry.forward(eval.timestamps, x$timestamps, x$predictions))
})

nonshock.post.icu.trajectories = lapply(nonshock.predictions, function(x) {
    eval.timestamps = (0:12)*60
    return(eval.carry.forward(eval.timestamps, x$timestamps, x$predictions))
})

shock.post.icu.trajectories = lapply(shock.predictions, function(x) {
    eval.timestamps = (0:12)*60
    return(eval.carry.forward(eval.timestamps, x$timestamps, x$predictions))
})

nonsepsis.post.icu.table = do.call(rbind,nonsepsis.post.icu.trajectories)
nonshock.post.icu.table = do.call(rbind,nonshock.post.icu.trajectories)
shock.post.icu.table = do.call(rbind,shock.post.icu.trajectories)

trajectories.icu = rbind(nonsepsis.post.icu.table,nonshock.post.icu.table,shock.post.icu.table)
means = colMeans(trajectories.icu,na.rm=T)
for (i in 1:dim(trajectories.icu)[2]) {
    trajectories.icu[is.na(trajectories.icu[,i]),i] = means[i]
}

S <- make.similarity(trajectories.icu, s)
A <- make.affinity(S, 3)  # use 3 neighboors (includes self)
D <- diag(apply(A, 1, sum))
U <- D - A
evL <- eigen(U, symmetric=TRUE)
plot(1:10, rev(evL$values)[1:10], log="y")
abline(v=4.5, col="red", lty=2)


# Analysis of by ewts

ewts.confusion = rep(NA,4)

ewts.confusion[1] = median(shock.ewts[sc@.Data[labels==3]==1&shock.ewts<0])
ewts.confusion[2] = median(shock.ewts[sc@.Data[labels==3]==2&shock.ewts<0])
ewts.confusion[3] = median(shock.ewts[sc@.Data[labels==3]==3&shock.ewts<0])
ewts.confusion[4] = median(shock.ewts[sc@.Data[labels==3]==4&shock.ewts<0])