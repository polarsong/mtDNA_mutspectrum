################################
################################

rm(list=ls(all=TRUE))

############ Syn mut
SynNuc = read.table("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t')
if (file.exists("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")){file.remove("../../Body/3Results/AllGenesCodonUsageNoOverlap.txt")}

names(SynNuc)
############ AnAge
AA = read.table("../../Body/1Raw/anage_data.txt", header = TRUE, sep = '\t')
AA$Species = paste(AA$Genus,AA$Species,sep = '_')

Alldata =  merge(SynNuc, AA, by = 'Species')
names(Alldata)
Mammalia = Alldata[Alldata$Class.y == "Mammalia",]
table(Mammalia$Order)


vec_of_Marsupials_orders = c("Dasyuromorphia", 'Didelphimorphia', "Diprotodontia", "Microbiotheria", "Notoryctemorphia", "Paucituberculata", "Peramelemorphia")

for (i in 1:nrow(Mammalia)){
  if (Mammalia$Order[i] %in% vec_of_Marsupials_orders){
    Mammalia$Subdivision[i] = "Marsupials"
  }else{
    Mammalia$Subdivision[i] = "nonMarsupials" 
    }  
}
for (i in 1:nrow(Mammalia)){
  if (Mammalia$Subdivision[i] == "nonMarsupials"){
    if (Mammalia$Order[i] == "Monotremata"){
      Mammalia$Subdivision[i] = "Monotremata" 
    }else{
      Mammalia$Subdivision[i] = "Placentals" 
  }
  }  
}


names(Mammalia)
table(Mammalia$Subdivision)

Mammalia= Mammalia[!is.na(Mammalia$Maximum.longevity..yrs.),]

library(ggpubr)

pdf("../../Body/4Figures/WholeGenomeAnalyses.NoOverlap.Placentals_Monotremata_Marsupials_A_G_Longevity.pdf", height = 20, width = 40)


ggboxplot(Mammalia, "Subdivision", "NeutralT",
          fill = "Subdivision", palette =c("#2b8cbe", "#a6bddb", "#ece7f2"))

ggboxplot(Mammalia, "Subdivision", "NeutralC", fill = "Subdivision", palette =c("#2ca25f", "#99d8c9", "#e5f5f9"))

ggboxplot(Mammalia, "Subdivision", "Maximum.longevity..yrs.", fill = "Subdivision", palette =c("#e34a33", "#fdbb84", "#fee8c8"))

dev.off()

