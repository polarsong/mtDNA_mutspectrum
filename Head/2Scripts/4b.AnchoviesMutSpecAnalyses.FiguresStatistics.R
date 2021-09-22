rm(list=ls(all=TRUE))

library(ggplot2)

fishA = read.table('../../Body/2Derived/Anchovies_Cytb_Clade_a.csv', sep = ';')
names(fishA) = c('codon_number','sub','freq','anc','dec','ancA','decA','typeSub')

fishB = read.table('../../Body/2Derived/Anchovies_Cytb_Clade_b.csv', sep = ';')
names(fishB) = c('codon_number','sub','freq','anc','dec','ancA','decA','typeSub')

### delete nonSyn Subs
fishA = fishA[fishA$typeSub == 'Syn',]

fishB = fishB[fishB$typeSub == 'Syn',]

###make counting

dfA = data.frame(table(fishA$sub))

dfB = data.frame(table(fishB$sub))

dfA$percent = prop.table(dfA$Freq) * 100

dfB$percent = prop.table(dfB$Freq) * 100

dfA$class = 'A'
dfB$class = 'B'

allinone = rbind(dfA,dfB)

### GRAPH FOR C>T on heavy chain
ggplot(data = allinone[allinone$Var1 == 'G>A',], aes(x = class, y = percent)) +
  geom_bar(position="dodge", stat="identity", fill = '#055088')+
  theme_void()+
  geom_text(x = 'A', y = 10, label = '15.6%', size = 5, col = 'white')+
  geom_text(x = 'A', y = 8, label = '(36/230)', size = 5, col = 'white')+
  geom_text(x = 'B', y = 10, label = '15.3%', size = 5, col = 'white')+
  geom_text(x = 'B', y = 8 , label = '(18/118)', size = 5, col = 'white')+
  geom_text(x = 'A', y = 17 , label = 'clade A', size = 5)+
  geom_text(x = 'B', y = 16.5 , label = 'clade B', size = 5)+
  ylim(0,50)

### GRAPH FOR G>A on heavy chain
ggplot(data = allinone[allinone$Var1 == 'C>T',], aes(x = class, y = percent)) +
  geom_bar(position="dodge", stat="identity", fill = '#9c3d37')+
  theme_void()+
  geom_text(x = 'A', y = 10, label = '19.6%', size = 5, col = 'white')+
  geom_text(x = 'A', y = 8, label = '(45/230)', size = 5, col = 'white')+
  geom_text(x = 'B', y = 9, label = '16.1%', size = 5, col = 'white')+
  geom_text(x = 'B', y = 7 , label = '(19/118)', size = 5, col = 'white')+
  geom_text(x = 'A', y = 21 , label = 'clade A',size = 5)+
  geom_text(x = 'B', y = 17.5 , label = 'clade B', size = 5)+
  ylim(0,50)


### GRAPH FOR A>G on heavy chain
ggplot(data = allinone[allinone$Var1 == 'T>C',], aes(x = class, y = percent)) +
  geom_bar(position="dodge", stat="identity", fill = '#73514f')+
  theme_void()+
  geom_text(x = 'A', y = 12, label = '23.5%', size = 5, col = 'white')+
  geom_text(x = 'A', y = 10, label = '(54/230)', size = 5, col = 'white')+
  geom_text(x = 'B', y = 9, label = '14.4%', size = 5, col = 'white')+
  geom_text(x = 'B', y = 7 , label = '(17/118)', size = 5, col = 'white')+
  geom_text(x = 'A', y = 25 , label = 'clade A', size = 5)+
  geom_text(x = 'B', y = 16 , label = 'clade B', size = 5)+
  ylim(0,50)

### GRAPH FOR T>C on heavy chain
ggplot(data = allinone[allinone$Var1 == 'A>G',], aes(x = class, y = percent)) +
  geom_bar(position="dodge", stat="identity", fill = '#036a5b')+
  theme_void()+
  geom_text(x = 'A', y = 15, label = '27.4%', size = 5, col = 'white')+
  geom_text(x = 'A', y = 13, label = '(63/230)', size = 5, col = 'white')+
  geom_text(x = 'B', y = 20, label = '44.1%', size = 5, col = 'white')+
  geom_text(x = 'B', y = 18 , label = '(52/118)', size = 5, col = 'white')+
  geom_text(x = 'A', y = 29 , label = 'clade A', size = 5)+
  geom_text(x = 'B', y = 46 , label = 'clade B', size = 5)+
  ylim(0,50)


## create df for A>G and T>C heavy chain to create mosaic plot

#install.packages("devtools")
#devtools::install_github("haleyjeppson/ggmosaic")
library(ggmosaic)

mos = allinone[allinone$Var1 == 'T>C' | allinone$Var1 == 'A>G', ]

#### TRANSLATE FROM LIGHT TO HEAVY CHAIN TO DRAW GRAPHIC

mos$Var1 = as.character(mos$Var1)

mos[1,1] = 'T>C'
mos[2,1] = 'A>G'
mos[3,1] = 'T>C'
mos[4,1] = 'A>G'


ggplot(data = mos)+
  geom_mosaic(aes(x = product(class,Var1), weight = percent), fill = c('#73514f','#73514f','#036a5b','#036a5b'))+
  theme_mosaic()+
  geom_text(x = 0.17, y = 0.82, label = '14.4%', col = 'white')+
  geom_text(x = 0.17, y = 0.79, label = '(17/118)', col = 'white')+
  geom_text(x = 0.17, y = 0.32, label = '23.5%', col = 'white')+
  geom_text(x = 0.17, y = 0.29, label = '(54/230)', col = 'white')+
  geom_text(x = 0.67, y = 0.72, label = '44.1%', col = 'white')+
  geom_text(x = 0.67, y = 0.69, label = '(52/118)', col = 'white')+
  geom_text(x = 0.67, y = 0.20, label = '27.4%', col = 'white')+
  geom_text(x = 0.67, y = 0.17, label = '(63/230)', col = 'white')+
  labs(x = ' ', y = ' ')

### calculate transversions

### for clade A
list_of_trans = c('A>C','C>A', 'A>T','T>A','G>T','T>G','C>G','G>C')
transA = 0
for ( i in 1:nrow(dfA))
  {
    if (dfA$Var1[i] %in% list_of_trans) {transA = transA + dfA$Freq[i]}
    
  }

transA # 32 transv subs 

percentA = transA/sum(dfA$Freq)*100
percentA = round(percentA, digits = 1) ## round off 

### for clade B

transB = 0
for ( i in 1:nrow(dfB))
{
  if (dfB$Var1[i] %in% list_of_trans) {transB = transB + dfB$Freq[i]}
  
}

transB # 12 transv subs 

percentB = transB/sum(dfB$Freq)*100
percentB = round(percentB, digits = 1) ## round off


### draw graph for transversions

transversion = data.frame(c(percentA, percentB), c('clade A', 'clade B'))
names(transversion) = c('Percentage', 'clade')

ggplot(data = transversion, aes(x = clade, y = Percentage)) +
  geom_bar(position="dodge", stat="identity", fill = '#055088')+
  theme_void()+
  geom_text(x = 'clade A', y = 7, label = '13.9%', col = 'white')+
  geom_text(x = 'clade A', y = 6.5, label = '(32/230)', col = 'white')+
  geom_text(x = 'clade B', y = 5, label = '10.2%', col = 'white')+
  geom_text(x = 'clade B', y = 4.5, label = '(12/118)', col = 'white')+
  geom_text(x = 'clade A', y = 14.5 , label = 'clade A', size = 5)+
  geom_text(x = 'clade B', y = 10.7 , label = 'clade B', size = 5)+
  ylim(0,20)

### make fishers tests

#### Th>Ch btw cladeA and cladeB 
fisher.test(cbind(c(dfA$Freq[dfA$Var1 == 'A>G'], sum(dfA$Freq[dfA$Var1 != 'A>G'])), 
                  c(dfB$Freq[dfB$Var1 == 'A>G'], sum(dfB$Freq[dfB$Var1 != 'A>G'])))) ### p = 0.0025, odds=0.48

#### Ah>Gh btw cladeA and cladeB 
fisher.test(cbind(c(dfA$Freq[dfA$Var1 == 'T>C'], sum(dfA$Freq[dfA$Var1 != 'T>C'])), 
                  c(dfB$Freq[dfB$Var1 == 'T>C'], sum(dfB$Freq[dfB$Var1 != 'T>C'])))) ### p = 0.049 ; odds= 1.81

#### Ch>Th btw cladeA and cladeB 
fisher.test(cbind(c(dfA$Freq[dfA$Var1 == 'G>A'], sum(dfA$Freq[dfA$Var1 != 'G>A'])), 
                  c(dfB$Freq[dfB$Var1 == 'G>A'], sum(dfB$Freq[dfB$Var1 != 'G>A'])))) ### p = 1 ; odds = 1.03

#### Gh>Ah btw cladeA and cladeB 
fisher.test(cbind(c(dfA$Freq[dfA$Var1 == 'C>T'], sum(dfA$Freq[dfA$Var1 != 'C>T'])), 
                  c(dfB$Freq[dfB$Var1 == 'C>T'], sum(dfB$Freq[dfB$Var1 != 'C>T'])))) ### p = 0.4678 ; odds = 1.27

#### transversions btw cladeA and cladeB 
fisher.test(cbind(c(transA, sum(dfA$Freq)-transA), 
                  c(transB, sum(dfB$Freq)-transB))) ### p = 0.395 ; odds = 0.142
