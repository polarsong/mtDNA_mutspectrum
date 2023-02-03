rm(list=ls(all=TRUE))

df_mtdna = read.csv('../../Head/2Scripts/Birds_dataset_paper.csv')
mut_data = read.table("C:/Users/User/Desktop/Birds mutspec results from Bogdan/mutations.tsv", header = TRUE)
df1 = read.table("C:/Users/User/Desktop/Birds mutspec results from Bogdan/expected_mutations.tsv", header = TRUE)
df2 = read.table("C:/Users/User/Desktop/Birds mutspec results from Bogdan/mutspec192.tsv", header = TRUE)
