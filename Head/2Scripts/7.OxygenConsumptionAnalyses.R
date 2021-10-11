rm(list=ls(all=T))

library(dplyr)
library(ggplot2)

MetabolicData = read.csv("../../Body/2Derived/FishBaseMetabolicData.csv", sep = ";")

MetabolicData$OxyCon.mg.kg.h.= gsub(",", ".", MetabolicData$OxyCon.mg.kg.h.)
MetabolicData$OnyCon.at20C.= gsub(",", ".", MetabolicData$OnyCon.at20C.)

MetabolicData$OxyCon.mg.kg.h.= as.numeric(MetabolicData$OxyCon.mg.kg.h.)
MetabolicData$OnyCon.at20C.= as.numeric(MetabolicData$OnyCon.at20C.)

length(unique(MetabolicData$Spece)) #206 species

### NEXT USE DPLYR, EVERYONE IT'S FUTURE

MD_standard  = MetabolicData %>% filter(Activity == 'standard' & Applied_stress == 'none_specified') ## take only routine and nonspecified
MD_standard = MD_standard %>% 
  group_by(Spece) %>% 
  mutate(median.oxyconc = median(OxyCon.mg.kg.h.)) %>% 
  select(Spece,median.oxyconc) %>% 
  distinct()

nrow(MD_standard) ### 69 Species for standard

muttemp = read.table('../../Body/2Derived/Supplementary_file_2.txt', header =TRUE)

MD_mut_standard = merge(muttemp, MD_standard, by.y ='Spece', by.x ='Species') # only 10 sps left ;(

MD_mut_standard = MD_mut_standard %>%  select(-c(14,15,16))

MD_mut_standard$AGTC = MD_mut_standard$T_C / MD_mut_standard$A_G
#draw plot cor with T_C and A_G

cor.test(MD_mut_standard$median.oxyconc , MD_mut_standard$A_G, method = 'spearman') ### p = 0.49

pdf('../../Body/4Figures/Supplementary_figure_7.pdf')
ggplot(MD_mut_standard, aes(x = median.oxyconc, y = A_G, label = 1:10))+
  geom_point()+
  theme_bw()+
  geom_text(hjust=-0.40, vjust=-0.40)+
  labs(x = 'Median Oxygen Consumption', y = 'Th to Ch')+
  geom_smooth(method ="lm", color = 'red',size = 0.6, se =F)

cor.test(MD_mut_standard$median.oxyconc , MD_mut_standard$T_C, method = 'spearman') ### p  = 0.068 still bad

ggplot(MD_mut_standard, aes(x = median.oxyconc, y = T_C, label = 1:10))+
  geom_point()+
  theme_bw()+
  labs(x = 'Median Oxygen Consumption', y = 'Ah to Gh')+
  geom_smooth(method ="lm", color = 'red',size = 0.6, se =F)+
  geom_text(hjust=1.3, vjust=-0.20)
  
cor.test(MD_mut_standard$median.oxyconc , MD_mut_standard$AGTC, method = 'spearman') ### p  = 0.14 still bad

ggplot(MD_mut_standard, aes(x = median.oxyconc, y = AGTC, label= 1:10))+
  geom_point()+
  theme_bw()+
  labs(x = 'Median Oxygen Consumption', y = 'A G / T C')+
  geom_smooth(method ="lm", color = 'red',size = 0.6, se =F)+
  geom_text(hjust=-0.4, vjust=0.50)

dev.off()

