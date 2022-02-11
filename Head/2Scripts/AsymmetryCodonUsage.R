rm(list=ls(all=T))

library(dplyr)
library(ggplot2)

df = read.table('../../Body/3Results/AllGenesCodonUsageNoOverlap.txt',
                header = TRUE, sep = '\t')
df = df[,-c(3,4,5,6,7,8)]
gen_len = read.table("../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt",
                     header = TRUE, sep = '\t')
gen_len$Species = gsub(' ','_',gen_len$Scientific_name)
gen_len = gen_len[,c(11,13)]

df_gen = merge(df,gen_len, by = 'Species')
df_gen = df_gen[df_gen$Gene != 'ND6',] ## delete ND6

### take sps that have all 12 needed genes

for (sp in unique(df_gen$Species))
{
  sp_genes = df_gen[df_gen$Species == sp, ]
  if (nrow(sp_genes) != 12) {df_gen = df_gen[df_gen$Species != sp, ]}
}

nrow(df_gen) ## all have 12 genes, we have 648 mammals

## make this big vector to order for normalization
vec_all = c('TTC','TTT','TCC','TCT','TAC','TAT','TGC','TGT',
            'TTA','TTG','TCA','TCG','TAA','TAG','TGA','TGG',
            'CTC','CTT','CCC','CCT','CAC','CAT','CGC','CGT',
            'CTA','CTG','CCA','CCG','CAA','CAG','CGA','CGG',
            'ATC','ATT','ACC','ACT','AAC','AAT','AGC','AGT',
            'ATA','ATG','ACA','ACG','AAA','AAG','AGA','AGG',
            'GTC','GTT','GCC','GCT','GAC','GAT','GGC','GGT',
            'GTA','GTG','GCA','GCG','GAA','GAG','GGA','GGG')

needed_codons = c('TTC','TCC','TAC','TGC',
              'TTA','TCA','TAA','TGA',
              'CTC','CCC','CAC','CGC',
              'CTA','CCA','CAA','CGA',
              'ATC','ACC','AAC','AGC',
              'ATA','ACA','AAA','AGA',
              'GTC','GCC','GAC','GGC',
              'GTA','GCA','GAA','GGA')


df_codons = df_gen %>% select(Species, Gene, GenerationLength_d,all_of(vec_all))

sp_sum_gen = data.frame(unique(df_codons$Species))

for ( codon in vec_all)
{
  
  sum_of_codon = aggregate(df_codons[ ,codon], by = list(df_codons$Species), FUN = 'sum')[2]
  sp_sum_gen = cbind(sp_sum_gen, sum_of_codon)

}

names(sp_sum_gen) = c('Species', vec_all)

# sp_sum_gen = merge(sp_sum_gen, gen_len)

codon_norm = data.frame()

for (i in 1:nrow(sp_sum_gen))
{
  org_gen = sp_sum_gen[i,]
  org_gen = as.vector(org_gen)
  df_out= data.frame(sp_sum_gen[i,]$Species) 
  for (codon in seq(from = 2, to = 65, by = 2))
  {
    if (org_gen[1,codon] == 0) {df_out = cbind(df_out, 0)}
    else {
      norm_cod = org_gen[1,codon] / (org_gen[1,codon+1] + org_gen[1,codon])
      df_out = cbind(df_out, norm_cod)
    }
  }
  names(df_out) = c('Species', needed_codons)
  codon_norm = rbind(codon_norm,df_out)
}


### decline stop codons
codon_norm = codon_norm %>% select(-c('TAA','AGA'))


## merge with Generation Leght

codon_norm = merge(codon_norm, gen_len)

### calculate median for each organism 
### for triplets ending on C or A

final = data.frame()
for (org in 1:nrow(codon_norm))
{
  sp_r = codon_norm[org,]

  vec_of_c = sp_r %>% select(TTC, TCC, TAC, TGC, CTC, CCC, CAC, CGC,
                               ATC, ACC, AAC, AGC, GTC, GCC, GAC, GGC)
  vec_of_a = sp_r %>% select(TTA, TCA, TGA, CTA, CCA, CAA, CGA,
                               ATA, ACA, AAA, GTA, GCA, GAA, GGA)
    
  med_c = median(as.numeric(vec_of_c), na.rm = TRUE)
  med_a = median(as.numeric(vec_of_a), na.rm = TRUE)
  sp_out = data.frame(sp_r$Species, med_c, med_a, sp_r$GenerationLength_d) 
  final = rbind(final,sp_out)
}

  
names(final) = c('Species', 'med_c', 'med_a', 'GL')

cor.test(final$med_c,final$GL, method = 'spearman')
cor.test(final$med_a,final$GL, method = 'spearman')

a = lm(final$median.vec_of_c. ~ final$g_sp_r.GenerationLength_d) ; summary(a)
summary(lm(final$median.vec_of_a. ~ final$g_sp_r.GenerationLength_d))

df_med_c = final[,-3]
df_med_c$type = 'XXC'
names(df_med_c) = c('Species', 'value', 'GL', 'type')

df_med_a = final[,-2]
df_med_a$type = 'XXA'
names(df_med_a) = c('Species', 'value', 'GL', 'type')

draw_hist = rbind(df_med_a, df_med_c)

ggplot(data = draw_hist, aes(x = value, fill = type))+
  geom_density(alpha = .8)+
  labs(x = 'Value of median', y = 'Frequency')+
  theme_classic()+
  guides(fill=guide_legend(title="Codon"))

wilcox.test(final$med_c, final$med_a, paired = TRUE)

summary(lm(final$med_c~log2(final$GL)))

summary(lm(final$med_a~log2(final$GL)))


final$medcdivmeda = final$med_c / final$med_a

### PGLS
library(ape)
library(geiger)
library(caper)

tree <- read.tree("/Users/diliushchenko/mito/mtDNA_mutspectrum/Body/1Raw/mtalign.aln.treefile.rooted")

Mam = final
row.names(Mam) = Mam$Species

tree_w = treedata(tree, Mam, sort=T, warnings=T)$phy

data<-as.data.frame(treedata(tree_w, Mam, sort=T, warnings=T)$data)

data$Species = as.character(data$Species)

data$med_c = as.numeric(as.character(data$med_c))
data$med_a = as.numeric(as.character(data$med_a))

data$GL = as.numeric(as.character(data$GL))

data$medcdivmeda = as.numeric(data$medcdivmeda)

### to discuss during meeting change position in PGLS
### WHY log2 doesn't work with PGLS med_a

MutComp = comparative.data(tree_w, data, Species, vcv=TRUE)



### SIGNIFICANT 
model = pgls(scale(med_c) ~ log2(GL), MutComp, lambda="ML")
summary(model)


### NOT SIGNIFICANT BUT TO PAPER
model = pgls(scale(med_a) ~ log2(GL), MutComp, lambda="ML")
summary(model)

#### dev c on a

model = pgls(scale(medcdivmeda) ~ log2(GL), MutComp, lambda="ML")
summary(model)

