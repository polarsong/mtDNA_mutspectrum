rm(list = ls(all=TRUE))
library(ggplot2)
library(ggpubr)
library(dplyr)
library(readr)
table1 = read_tsv('expected_freqs.tsv')
table2 = table1[table1$Label == 'syn',]
