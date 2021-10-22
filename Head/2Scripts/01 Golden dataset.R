rm(list = ls(all.names = TRUE))
gc() 
df = read.table("../../Body/3Results/Birds codus + neutral.csv", header = TRUE, sep = ',') #opening dirty dataset
df1 = df[df$wrong.amino..==0,]
df1$temp = 1
df2 = aggregate(df1$temp,by = list(df1$Species.name), FUN = sum)
df3 = df2[df2$x != 13,] #checking wrong birds
df4 = df2[df2$x == 13,] #checking right birds
unique(df$Gene.name)#looking for wrong gene names
which(df == "[cytb]", arr.ind=TRUE) #looking for positions
df[11264, "Gene.name"] = "[CYTB]" 
df[11277, "Gene.name"] = "[CYTB]" 
df[11290, "Gene.name"] = "[CYTB]" 
df[11303, "Gene.name"] = "[CYTB]" 
df[11316, "Gene.name"] = "[CYTB]" 
which(df == "[COXI]", arr.ind=TRUE) #looking for positions
df[10969, "Gene.name"] = "[COX1]"
which(df == "[COXII]", arr.ind=TRUE) #looking for positions
df[10970, "Gene.name"] = "[COX2]"
which(df == "[COXIII]", arr.ind=TRUE) #looking for positions
df[10973, "Gene.name"] = "[COX3]"

#NA deleting
nal = df[is.na(df$Gene.name),] #24 birds
which(is.na(df), arr.ind=TRUE)
for (i in unique(nal$Species.name)) #removing NA birds
{
df = df[!(df$Species.name == i),]  
}
unique(df$Gene.name) #checking
write.csv(df, file = "../../Body/3Results/Golden birds dataset.csv", sep = ',', row.names = FALSE)

