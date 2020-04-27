rm(list=ls())
intake.data = readRDS("intake.data.rds")
# For each patient

library(pracma)
library(tictoc)

unique.patients = sort(unique(intake.data$patientunitstayid))
fluid.data = vector(mode="list",length=length(unique.patients))

interval = 1000

for (i in 1:length(unique.patients)) {
    if (i %% interval == 0) {
        fprintf("%d of %d...\n",i,length(unique.patients))
    }

    current.data = intake.data[intake.data$patientunitstayid==unique.patients[i],]
    offsets = sort(unique(current.data$intakeoutputoffset))
    values = sapply(offsets, function(x) max(current.data$intaketotal[current.data$intakeoutputoffset==x]))

    fluid.data[[i]] = list(subject.id = unique.patients[i], timestamps = offsets, values = values)
    
}

saveRDS(fluid.data,file="fluid.data.rds")