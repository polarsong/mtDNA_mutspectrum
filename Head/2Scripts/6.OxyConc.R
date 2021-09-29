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

MD  = MetabolicData %>% filter(Activity == 'standard' & Applied_stress == 'none_specified') ## take only routine and nonspecified
MD = MD %>% 
  group_by(Spece) %>% 
  mutate(median.oxyconc = median(OxyCon.mg.kg.h.)) %>% 
  select(Spece,median.oxyconc) %>% 
  distinct()

nrow(MD) ### 311 Species

muttemp = read.table('../../Body/2Derived/Supplementary_table_2.txt', header =TRUE)

MD_mut = merge(muttemp, MD, by.y ='Spece', by.x ='Species') # only 20 sps left ;(

MD_mut = MD_mut %>%  select(-c(14,15,16))

#draw plot cor with T_C and A_G

cor.test(MD_mut$median.oxyconc , MD_mut$A_G, method = 'spearman') ### p = 0.49

ggplot(MD_mut, aes(x = median.oxyconc, y = A_G))+
  geom_point()+
  theme_bw()+
  labs(x = 'Median Oxygen Consumption', y = 'T to C')+
  geom_smooth(method ="lm", color = 'red',size = 0.6, se =F)

cor.test(MD_mut$median.oxyconc , MD_mut$T_C, method = 'spearman') ### p  = 0.068 still bad

ggplot(MD_mut, aes(x = median.oxyconc, y = T_C))+
  geom_point()+
  theme_bw()+
  labs(x = 'Median Oxygen Consumption', y = 'A to G')+
  geom_smooth(method ="lm", color = 'red',size = 0.6, se =F)



