rm(list=ls())
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

patients.has.shock = is.element(subjects, shock.icustays)
patients.has.sepsis = is.element(subjects, sepsis.icustays)

comorbidities = array(F,dim=c(length(subjects),17))

# 1. Myocardial infarction (1)
comorbidities[,1] = sapply(codes, function(x) any((x>=410&x<411)|(x>=412&x<413)))

# 2. Congestive heart failure (1)
comorbidities[,2] = sapply(codes, function(x) any((x>=398.91&x<398.92)|(x>=402.01&x<402.02)|(x>=402.11&x<402.12)|
    (x>=402.91&x<402.92)|(x>=404.01&x<404.01)|(x>=404.03&x<404.04)|(x>=404.11&x<404.12)|(x>=404.13&x<404.14)|
    (x>=404.91&x<404.92)|(x>=404.93&x<404.94)|(x>=425.4&x<426)|(x>=428&x<429)))

# 3. Peripheral vascular disease (1)
comorbidities[,3] = sapply(codes, function(x) any((x>=93&x<94)|(x>=437.3&x<437.4)|(x>=440&x<441)|(x>=441&x<442)|
    (x>=443.1&x<444)|(x>=47.1&x<47.2)|(x>=557.1&x<557.2)|(x>=557.9&x<558)))

# 4. Cerebrovascular disease (1)
comorbidities[,4] = sapply(codes, function(x) any((x>=362.34&x<362.35)|(x>=430&x<439)))

# 5. Dementia (1)
comorbidities[,5] = sapply(codes, function(x) any((x>=290&x<291)|(x>=294.1&x<294.2)|(x>=331.2&x<331.3)))

# 6. Chronic pulmonary disease (1)
comorbidities[,6] = sapply(codes, function(x) any((x>=416.8&x<416.9)|(x>=416.9&x<417)|(x>=490&x<506)|
    (x>=506.4&x<506.5)|(x>=508.1&x<508.2)|(x>=508.8&x<508.9)))

# 7. Rheumatic disease (1)
comorbidities[,7] = sapply(codes, function(x) any((x>=446.5&x<446.6)|(x>=710&x<710.5)|(x>=714&x<714.3)|
    (x>=714.8&x<714.9)|(x>=725&x<726)))

# 8. Peptic ulcer disease (1)
comorbidities[,8] = sapply(codes, function(x) any((x>=531&x<535)))

# 9. Mild liver disease (1)
comorbidities[,9] = sapply(codes, function(x) any((x>=70.22&x<70.23)|(x>=70.23&x<70.24)|(x>=70.32&x<70.33)|
    (x>=70.33&x<70.34)|(x>=70.44&x<70.45)|(x>=70.54&x<70.55)|(x>=70.6&x<70.7)|(x>=70.9&x<71)|(x>=570&x<571)|
    (x>=571&x<572)|(x>=573.3&x<573.4)|(x>=573.4&x<573.5)|(x>=573.8&x<573.9)|(x>=573.9&x<574)))

# 10. Diabetes without chronic complication (1)
comorbidities[,10] = sapply(codes, function(x) any((x>=250&x<250.4)|(x>=250.8&x<250.9)|(x>=250.9&x<260)))

# 11. Diabetes with chronic complication (2)
comorbidities[,11] = sapply(codes, function(x) any(x>=250.4&x<250.8))

# 12. Hemiplegia or paraplegia (2)
comorbidities[,12] = sapply(codes, function(x) any((x>=334.1&x<334.2)|(x>=342&x<343)|(x>=343&x<344)|(x>=344&x<344.7)|
    (x>=344.9&x<345)))

# 13. Renal disease (2)
comorbidities[,13] = sapply(codes, function(x) any((x>=403.01&x<403.02)|(x>=403.11&x<403.12)|(x>=403.91&x<403.92)|
    (x>=404.02&x<404.03)|(x>=404.03&x<404.04)|(x>=404.12&x<404.14)|(x>=404.92&x<404.94)|(x>=582&x<583.8)|
    (x>=585&x<587)|(x>=588&x<589)))

# 14. Any malignancy (2)
comorbidities[,14] = sapply(codes, function(x) any((x>=140&x<173)|(x>=174&x<195.9)|(x>=200&x<209)|(x>=238.6&x<238.7)))

# 15. Moderate or severe liver disease (3)
comorbidities[,15] = sapply(codes, function(x) any((x>=456&x<456.3)|(x>=572.2&x<572.9)))

# 16. Metastatic solid tumor (6)
comorbidities[,16] = sapply(codes, function(x) any((x>=196&x<200)))
# 17. AIDS/HIV (6)
comorbidities[,17] = sapply(codes, function(x) any((x>=42&x<45)))


# Compute aggregate score
score = comorbidities[,1] + comorbidities[,2] + comorbidities[,3] + comorbidities[,4] + comorbidities[,5] +
    comorbidities[,6] + comorbidities[,7] + comorbidities[,8] + comorbidities[,9] + comorbidities[,10] +
    comorbidities[,11] * 2 + comorbidities[,12] * 2 + comorbidities[,13] * 2 + comorbidities[,14] * 2 +
    comorbidities[,15] * 3 + comorbidities[,16] * 6 + comorbidities[,17] * 6

