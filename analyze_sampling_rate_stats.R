rm(list=ls())
library(matrixStats)

clinical.data = readRDS("processed/clinical_data_icd9_sofa_vent.rds")
source("functions/generate_sampling_rate_table.R")

tic("Generating sampling rate data")
sampling.rate.data = lapply(clinical.data, generate.sampling.rate.table)
toc()

sampling.rate = t(sapply(sampling.rate.data, function(x) {
    sapply(2:28, function(y) {
        if (!is.null(x[[y]])) {
            return(mean(x[[y]]$values,na.rm=T))
        } else {
            return(NA)
        }
    })
}))

names = names(clinical.data[[1]])
sampling.rate.mean = colMeans(sampling.rate,na.rm=T)
sampling.rate.median = colMedians(sampling.rate,na.rm=T)
names(sampling.rate.mean) = names[-1]
names(sampling.rate.median) = names[-1]