rm(list=ls())

library(RPostgreSQL)
library(tictoc)

user = "postgres"
password = "postgres"
db = "eicu"

query = "select * from pasthistory"

tic("Past history query")
connection = dbConnect(PostgreSQL(), user=user, password=password, dbname=db)
past.history = dbGetQuery(connection, query)
dbDisconnect(connection)
toc()

saveRDS(past.history,file="past.history.rds")