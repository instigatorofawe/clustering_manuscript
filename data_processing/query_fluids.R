rm(list=ls())

library(RPostgreSQL)
library(tictoc)

user = "postgres"
password = "postgres"
db = "eicu"

query = "select * from intakeoutput where intaketotal > 0"
tic("Fluid intake query")
connection = dbConnect(PostgreSQL(), user=user, password=password, dbname=db)
intake.data = dbGetQuery(connection, query)
dbDisconnect(connection)
toc()

saveRDS(intake.data,file="intake.data.rds")