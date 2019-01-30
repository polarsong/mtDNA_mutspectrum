###################################
###### 
###################################

#### 
rm(list=ls(all=TRUE))

############ read cancer fold changes
unzip("../../Body/1Raw/DEG_GlycOxph_ALL.zip", exdir = "../../Body/1Raw/")
Tr = read.table("../../Body/1Raw/DEG_GlycOxph_ALL.csv", header = TRUE, sep = ';')
if (file.exists("../../Body/1Raw/DEG_GlycOxph_ALL.csv")) file.remove("../../Body/1Raw/DEG_GlycOxph_ALL.csv")

# list of genes involved in Glycolysis (KEGG) from Suppl table 2 (NatureCommunication: Pan-cancer analysis of transcriptional metabolic dysregulation using The Cancer Genome Atlas)
GlycVec = c('ACSS1','ACSS2','ADH1A','ADH1B','ADH1C','ADH4','ADH5','ADH6','ADH7','ADPGK','AKR1A1','ALDH1A3','ALDH1B1','ALDH2','ALDH3A1','ALDH3A2','ALDH3B1','ALDH3B2','ALDH7A1','ALDH9A1','ALDOA','ALDOB','ALDOC','BPGM','DLAT','DLD','ENO1','ENO2','ENO3','FBP1','FBP2','G6PC','G6PC2','GALM','GAPDH','GAPDHS','GCK','GPI','HK1','HK2','HK3','HKDC1','LDHA','LDHAL6A','LDHAL6B','LDHB','LDHC','PANK1','PCK1','PCK2','PDHA1','PDHA2','PDHB','PFKFB1','PFKFB2','PFKFB3','PFKFB4','PFKL','PFKM','PFKP','PGAM1','PGAM2','PGAM4','PGK1','PGK2','PGM1','PGM2','PKLR','PKM','SLC2A2','TPI1')

Tr = Tr[Tr$Gene %in% GlycVec,] 

# find genes with minimal p values -> ADH and ALD
AGG = aggregate(-log10(Tr$adj.P.Val), by = list(Tr$Gene), FUN = median); names(AGG)=c('Gene','MedianLogP'); AGG = AGG[order(AGG$MedianLogP),]

# lets start from well known markers of glycolysis: GAPDH, PKM,  ENO, PFK, PGM, HK, PDC
# GAPDH: Glyceraldehyde 3-phosphate dehydrogenase. Is an enzyme of ~37kDa that catalyzes the sixth step of glycolysis and thus serves to break down glucose for energy and carbon molecules. In addition to this long established metabolic function, GAPDH has recently been implicated in several non-metabolic processes, including transcription activation, initiation of apoptosis,[5] ER to Golgi vesicle shuttling, and fast axonal, or axoplasmic transport.[6] In sperm, a testis-specific isoenzyme GAPDHS is expressed. 
# PKM: Pyruvate Kinase M1/2. This gene encodes a protein involved in glycolysis. The encoded protein is a pyruvate kinase that catalyzes the transfer of a phosphoryl group from phosphoenolpyruvate to ADP, generating ATP and pyruvate. This protein has been shown to interact with thyroid hormone and may mediate cellular metabolic effects induced by thyroid hormones. 
# ENO: enolase (there are three enolase isoenzymes in mammals). They are glycolytic enzymes.
# PFK: Phosphofructokinase (M - muscle; P = platelet, L = liver; PFKL','PFKM','PFKP'); Three phosphofructokinase isozymes exist in humans: muscle, liver and platelet. These isozymes function as subunits of the mammalian tetramer phosphofructokinase, which catalyzes the phosphorylation of fructose-6-phosphate to fructose-1,6-bisphosphate
# PGM: (1,2) Phosphoglucomutase. There are several PGM isozymes, which are encoded by different genes and catalyze the transfer of phosphate between the 1 and 6 positions of glucose
# HK(1,2,3): Hexokinases phosphorylate glucose to produce glucose-6-phosphate, the first step in most glucose metabolism pathways
# PFKFB(1,2,3,4): 6-Phosphofructo-2-Kinase/Fructose-2,6-Biphosphatase. These proteins belong to a family of bifunctional proteins that are involved in both the synthesis and degradation of fructose-2,6-bisphosphate, a regulatory molecule that controls glycolysis in eukaryotes

# http://www.ijcep.com/files/ijcep0063886.pdf

## aggregate Score or logFC. Score works better!

GAPDH = Tr[Tr$Gene == 'GAPDH',]; 
GAPDH = aggregate(GAPDH$Score, by = list(GAPDH$Cancer), FUN = mean); names(GAPDH)=c('CancerType','GAPDH')

ENO   = Tr[Tr$Gene %in% c('ENO1','ENO2','ENO3'),] 
ENO   = aggregate(ENO$Score, by = list(ENO$Cancer), FUN = mean); names(ENO)=c('CancerType','ENO')

PFK  = Tr[Tr$Gene %in% c('PFKL','PFKM','PFKP'),] 
PFK   = aggregate(PFK$Score, by = list(PFK$Cancer), FUN = mean); names(PFK)=c('CancerType','PFK')

PGM  = Tr[Tr$Gene %in% c('PGM1','PGM2'),]
PGM   = aggregate(PGM$Score, by = list(PGM$Cancer), FUN = mean); names(PGM)=c('CancerType','PGM')

HK    = Tr[Tr$Gene %in% c('HK1','HK2','HK3'),]
HK   = aggregate(HK$Score, by = list(HK$Cancer), FUN = mean); names(HK)=c('CancerType','HK')

ADH  = Tr[Tr$Gene %in% c('ADH1A','ADH1B','ADH1C','ADH4','ADH5','ADH6','ADH7'),]
ADH  = aggregate(ADH$Score, by = list(ADH$Cancer), FUN = mean); names(ADH)=c('CancerType','ADH')

ALD = Tr[Tr$Gene %in% c('ALDH1A3','ALDH1B1','ALDH2','ALDH3A1','ALDH3A2','ALDH3B1','ALDH3B2','ALDH7A1','ALDH9A1','ALDOA','ALDOB','ALDOC'),]
ALD  = aggregate(ALD$Score, by = list(ALD$Cancer), FUN = mean); names(ALD)=c('CancerType','ALD')

Glyc = merge(GAPDH,ENO, all = TRUE)
Glyc = merge(Glyc,PFK, all = TRUE)
Glyc = merge(Glyc,PGM, all = TRUE)
Glyc = merge(Glyc,HK, all = TRUE)
Glyc = merge(Glyc,ADH, all = TRUE)
Glyc = merge(Glyc,ALD, all = TRUE)

cor.test(Glyc$GAPDH,Glyc$ENO, method = 'spearman') # positive
cor.test(Glyc$GAPDH,Glyc$PFK, method = 'spearman') # positive

cor.test(Glyc$GAPDH,Glyc$PGM, method = 'spearman') # nothing
cor.test(Glyc$GAPDH,Glyc$HK, method = 'spearman')  # nothing

cor.test(Glyc$GAPDH,Glyc$ADH, method = 'spearman')  # negative
cor.test(Glyc$GAPDH,Glyc$ALD, method = 'spearman')  # nothing

# PCA!?

write.table(Glyc, file = '../../Body/2Derived/Cancer.DeriveLevelOfGlycolysisForEachCancerType.txt')







