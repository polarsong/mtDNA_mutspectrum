rm(list = ls(all=TRUE))
library(ggfortify)
library(ggplot2)
library(ggpubr)
library(ggbiplot)
library(phytools)
library(nlme)
library(geiger)
library(ggtree)
library(stringr)
library(dplyr)
library(RSQLite)
rangeMapper <- "rangeMapper_ecog.00929.sqlite"
sqlite.driver <- dbDriver("SQLite")
db <- dbConnect(sqlite.driver,
                dbname = rangeMapper)

dbListTables(db) #Выводит список таблиц в SQL базе
BIO_Life_History <- dbReadTable(db, "BIO_Life_History")
write.csv2(BIO_Life_History, file="BIO_Life_History.csv")
if (file.exists("rangeMapper_ecog.00929.sqlite")) file.remove("rangeMapper_ecog.00929.sqlite")
