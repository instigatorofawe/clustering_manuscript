rm(list=ls())

library(pracma)

load("processed/diagnosis_query.rdata")
has.infection.icd9 = readRDS("processed/has_infection_icd9.rds")
clinical.data = readRDS("processed/clinical_data_icd9_sofa_vent.rds")
patient.data = readRDS("processed/patient_data.rds")
sofa.scores = readRDS("processed/sofa_scores.rds")
query.subjects = subjects[has.infection.icd9]

# Analyze unit type
unit.types = unique(patient.data$unittype)
unit.counts = sapply(unit.types, function(x) sum(patient.data$unittype==x))
unit.percentages = unit.counts/sum(unit.counts)

lengths = sapply(sofa.scores, function(x) length(x$timestamps))
sepsis.labels = sapply(sofa.scores[lengths>0], function(x) rowSums(x[2:7])>=2)
has.sepsis = sapply(sepsis.labels, any)

sepsis.timestamps = mapply(function(x,y) x$timestamps[y], sofa.scores[lengths>0][has.sepsis],sepsis.labels[has.sepsis])

# Determine shock onsets
shock.labels = mapply(function(x,y) x&y$lactate&y$vasopressors, sepsis.labels, sofa.scores[lengths>0])
has.shock = sapply(shock.labels, function(x) any(x,na.rm=T)) #& has.dx[lengths>0]

clinical.icustays = sapply(clinical.data[lengths>0], function(x) x$subject.id)
sepsis.icustays = clinical.icustays[has.sepsis]
shock.icustays = clinical.icustays[has.shock]

patients.has.shock = is.element(patient.data$patientunitstayid, shock.icustays)
patients.has.sepsis = is.element(patient.data$patientunitstayid, sepsis.icustays)

# Categorize patients
# Number of patients by category
unique.subjects = unique(patient.data$uniquepid)
shock.subjects = unique(patient.data$uniquepid[patients.has.shock])
sepsis.subjects = unique(patient.data$uniquepid[patients.has.sepsis])
sepsis.subjects = sepsis.subjects[!is.element(sepsis.subjects, shock.subjects)]
nonsepsis.subjects = unique.subjects[!is.element(unique.subjects,sepsis.subjects)&!is.element(unique.subjects,shock.subjects)]

shock.mortality = sum(patient.data$unitdischargestatus[patients.has.shock]=="Expired")/sum(patients.has.shock)
sepsis.mortality = sum(patient.data$unitdischargestatus[patients.has.sepsis&!patients.has.shock]=="Expired")/sum(patients.has.sepsis&!patients.has.shock)
nonsepsis.mortality = sum(patient.data$unitdischargestatus[!patients.has.sepsis&!patients.has.shock]=="Expired")/sum(!patients.has.sepsis&!patients.has.shock)

shock.gender = sum(patient.data$gender[patients.has.shock]=="Male")/sum(patients.has.shock)
sepsis.gender = sum(patient.data$gender[patients.has.sepsis&!patients.has.shock]=="Male")/sum(patients.has.sepsis&!patients.has.shock)
nonsepsis.gender = sum(patient.data$gender[!patients.has.sepsis&!patients.has.shock]=="Male")/sum(!patients.has.sepsis&!patients.has.shock)

# Manual conversion of > 89 to 89
patient.age.numeric = sapply(patient.data$age, function(x) {
    if (x == "") {
        return(NA)
    } else if (x == "> 89") {
        return(89)
    } else {
        return(as.numeric(x))
    }
})
shock.ages = mean(patient.age.numeric[patients.has.shock],na.rm=T)
sepsis.ages = mean(patient.age.numeric[patients.has.sepsis&!patients.has.shock],na.rm=T)
nonsepsis.ages = mean(patient.age.numeric[!patients.has.sepsis&!patients.has.sepsis],na.rm=T)

shock.ages.sd = sd(patient.age.numeric[patients.has.shock],na.rm=T)
sepsis.ages.sd = sd(patient.age.numeric[patients.has.sepsis&!patients.has.shock],na.rm=T)
nonsepsis.ages.sd = sd(patient.age.numeric[!patients.has.sepsis&!patients.has.sepsis],na.rm=T)

# ICU Stay length
shock.lengths = median(patient.data$unitdischargeoffset[patients.has.shock],na.rm=T)
sepsis.lengths = median(patient.data$unitdischargeoffset[patients.has.sepsis&!patients.has.shock],na.rm=T)
nonsepsis.lengths = median(patient.data$unitdischargeoffset[!patients.has.sepsis&!patients.has.shock],na.rm=T)