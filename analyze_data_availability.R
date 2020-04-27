rm(list=ls())
clinical.data = readRDS("data/eicu/clinical_data_icd9_sofa_vent.rds")

names = names(clinical.data[[1]])

availability = t(sapply(clinical.data, function(x) {
    sapply(2:28, function(y) {
        !is.null(x[[y]])
    })
}))

availability.fraction = apply(availability, 2, function(x) sum(x)/length(x))
names(availability.fraction) = names[-1]