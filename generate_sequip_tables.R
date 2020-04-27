rm(list=ls())

load("data/sequip/early_prediction_data_sofa.v2.alt.rdata")
sequip.names = colnames(nonshock.data[[1]])

load("data/eicu/reference_dataset_rx_combined_2.rdata")
eicu.names = colnames(nonsepsis.data[[1]])

indices = sapply(sequip.names, function(x) which(eicu.names==x))
indices$lact = which(eicu.names=="lactate")
indices = do.call(c,indices)

nonsepsis.train = lapply(nonsepsis.data, function(x) x[,indices])
nonshock.train = lapply(nonshock.data, function(x) x[,indices])
preshock.train = lapply(preshock.data, function(x) x[,indices])

saveRDS(list(nonsepsis.train=nonsepsis.train,nonshock.train=nonshock.train,preshock.train=preshock.train),
        file="data/eicu/sequip.training.rds")