rm(list=ls())
source("src/R/eicu/load_reference_data_2.R")
load("data/eicu/post_detection_clustering.RData")

library(kernlab)
library(ggplot2)
library(tidyr)
library(tictoc)
library(class)
library(matrixStats)

training = c(nonsepsis.sample[has.detection.nonsepsis],nonshock.sample[has.detection.nonshock],preshock.sample[has.detection.shock])
test.clusters = sc@.Data[!training]

confusion = array(NA,dim=c(4,3))
confusion.train = array(NA,dim=c(4,3))
for (i in 1:4) {
    for (j in 1:3) {
        confusion[i,j] = sum(sc@.Data==i&labels==j)
        confusion.train[i,j] = sum(sc.train@.Data==i&labels.train==j)
    }
}

cluster.indices = order(confusion[,3]/rowSums(confusion),decreasing=T)
cluster.indices.train = order(confusion.train[,3]/rowSums(confusion.train),decreasing=T)

acc = rep(NA,13)
acc.exclusive = rep(NA,13)

acc.tertiary = rep(NA,13)
sens.tertiary = rep(NA,13)
spec.tertiary = rep(NA,13)

for (i in 1:13) {
    knn.prediction = knn(trajectories.train[,1:i,drop=F],trajectories.test[,1:i,drop=F],sc.train@.Data,5)
    prediction.indices = sapply(knn.prediction, function(x) which(cluster.indices.train==x))
    translated.clusters = cluster.indices[prediction.indices]
    acc[i] = sum(translated.clusters==test.clusters)/length(test.clusters)
    acc.exclusive[i] = sum(translated.clusters[test.clusters==cluster.indices[1]|test.clusters==cluster.indices[4]]==test.clusters[test.clusters==cluster.indices[1]|test.clusters==cluster.indices[4]])/sum(test.clusters==cluster.indices[1]|test.clusters==cluster.indices[4])

    acc.tertiary[i] = (sum(translated.clusters[test.clusters==cluster.indices[1]]==cluster.indices[1])+sum(translated.clusters[test.clusters!=cluster.indices[1]]!=cluster.indices[1]))/length(translated.clusters)
    sens.tertiary[i] = sum(translated.clusters[test.clusters==cluster.indices[1]]==cluster.indices[1])/sum(test.clusters==cluster.indices[1])
    spec.tertiary[i] = sum(translated.clusters[test.clusters!=cluster.indices[1]]!=cluster.indices[1])/sum(test.clusters!=cluster.indices[1])
}

#qplot(0:12,acc,geom="line")+xlab("Data (hrs)")+ylab("Accuracy")+ylim(c(0,1))

#qplot(0:12,acc.exclusive,geom="line")+xlab("Data (hrs)")+ylab("Accuracy")+ylim(c(0,1))

qplot(0:12,acc.tertiary,geom="line")+xlab("Data (hrs)")+ylab("Accuracy")+ylim(0.5,1)+
    scale_x_continuous(breaks=c(0:12))+
    theme(axis.text.x = element_text(size=12),
          axis.title.x = element_text(size=14),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=14),
          legend.text = element_text(size=12),
          legend.title = element_text(size=14))

#qplot(0:12,sens.tertiary,geom="line")+xlab("Data (hrs)")+ylab("Sensitivity")+ylim(0.4,1)
#qplot(0:12,spec.tertiary,geom="line")+xlab("Data (hrs)")+ylab("Specificity")+ylim(0.5,1)