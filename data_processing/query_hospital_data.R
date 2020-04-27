rm(list=ls())

library(RPostgreSQL)
library(pracma)
library(tictoc)

user = "postgres"
password = "postgres"
db = "eicu"

tic("Hospital query")
query = "select * from hospital"
connection = dbConnect(PostgreSQL(), user=user, password=password, dbname=db)
hospital.result = dbGetQuery(connection, query)
dbDisconnect(connection)
toc()

saveRDS(hospital.result, file="processed/hospital_data.rds")

# Analyze hospital id distribution

icd.subjects = readRDS("processed/icd9_subjects.rds")
patient.result = readRDS("processed/patient_data.rds")
patient.result = patient.result[is.element(patient.result$patientunitstayid,icd.subjects),]


hospital.ids = unique(patient.result$hospitalid)
counts = sapply(hospital.ids, function(x) sum(patient.result$hospitalid == x))
fractions = counts/sum(counts)

patient.counts = sapply(hospital.ids, function(x) length(unique(patient.result$uniquepid[patient.result$hospitalid==x])))

