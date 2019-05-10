rm(list=ls(all=TRUE))

### non independence of mutations within the same patient / tissue ??? => nested logistic (linear) regression (glm) with classes | should we????

Som = read.table("../../Body/2Derived/HealthySomaticHumanMutations.GTEx.PatientSpecificDerive.txt", header = TRUE, sep = '\t')

### derive MutSpek as a set of Dummy variables:

Som$G_A = 0; Som[Som$Substitution == 'G_A',]$G_A = 1;
Som$T_C = 0; Som[Som$Substitution == 'T_C',]$T_C = 1;
Som$A_G = 0; Som[Som$Substitution == 'A_G',]$A_G = 1;
Som$C_T = 0; Som[Som$Substitution == 'C_T',]$C_T = 1;
Som$Ts = 0;  Som[Som$Substitution == 'G_A' | Som$Substitution == 'T_C' | Som$Substitution == 'A_G' | Som$Substitution == 'C_T',]$Ts = 1;
Som$Tv = 1; Som$Tv = Som$Tv - Som$Ts;  

## many G_A (transitions) in short living cells, 
## few C_T in short living cells.
a<-glm(Som$G_A ~ Som$TurnOverRate, family = binomial); summary(a) # -3.321e-05  6.832e-06  -4.861 1.17e-06 ***
a<-glm(Som$A_G ~ Som$TurnOverRate, family = binomial); summary(a) # 3.438e-06  6.420e-06   0.536    0.592 
a<-glm(Som$T_C ~ Som$TurnOverRate, family = binomial); summary(a) # -1.539e-05  6.202e-06  -2.482   0.0131 *
a<-glm(Som$C_T ~ Som$TurnOverRate, family = binomial); summary(a) # 2.476e-05  6.476e-06   3.824 0.000131 ***
a<-glm(Som$Ts ~ Som$TurnOverRate, family = binomial); summary(a)  # -3.321e-05  6.236e-06  -5.325 1.01e-07 ***

## The level of asymmetry is changing: in fast  dividing mutations are more asymmetrical (high G_A/C_T)
a<-glm(Som[Som$Substitution == 'G_A'|Som$Substitution == 'C_T',]$G_A ~ Som[Som$Substitution == 'G_A'|Som$Substitution == 'C_T',]$TurnOverRate, family = binomial); summary(a) # -4.598e-05  8.497e-06  -5.411 6.26e-08 ***
a<-glm(Som[Som$Substitution == 'T_C'|Som$Substitution == 'A_G',]$T_C ~ Som[Som$Substitution == 'T_C'|Som$Substitution == 'A_G',]$TurnOverRate, family = binomial); summary(a) # -1.391e-05  7.812e-06  -1.780   0.0751 .

### add many potential counfounders: AF (fraction of heteroplsmy), PercentMito (relative expression level of Mitochondrial genes), Age, Sex, Position

#### recode age intervals to Dummy variables:
table(Som$age)
#20-29 30-39 40-49 50-59 60-69 70-79 
#195   156   439   865   880    30
Som$dummy20 <- as.numeric(Som$age == '20-29')
Som$dummy30 <- as.numeric(Som$age == '30-39')
Som$dummy40 <- as.numeric(Som$age == '40-49')
Som$dummy50 <- as.numeric(Som$age == '50-59')
Som$dummy60 <- as.numeric(Som$age == '60-69') 

#### sex - recode sex:
table(Som$sex)
Som$DummyMale <-as.numeric(Som$sex == '1')

a<-glm(Som$G_A ~ Som$TurnOverRate + Som$AF + Som$mitoReads + Som$percentMito + Som$Position + Som$dummy20 + Som$dummy30 + Som$dummy40 + Som$dummy50 + Som$dummy60 + Som$DummyMale, family = binomial); summary(a)
a<-glm(Som$T_C ~ Som$TurnOverRate + Som$AF + Som$mitoReads + Som$percentMito + Som$Position + Som$dummy20 + Som$dummy30 + Som$dummy40 + Som$dummy50 + Som$dummy60 + Som$DummyMale, family = binomial); summary(a)
a<-glm(Som$C_T ~ Som$TurnOverRate + Som$AF + Som$mitoReads + Som$percentMito + Som$Position + Som$dummy20 + Som$dummy30 + Som$dummy40 + Som$dummy50 + Som$dummy60 + Som$DummyMale, family = binomial); summary(a)
a<-glm(Som$A_G ~ Som$TurnOverRate + Som$AF + Som$mitoReads + Som$percentMito + Som$Position + Som$dummy20 + Som$dummy30 + Som$dummy40 + Som$dummy50 + Som$dummy60 + Som$DummyMale, family = binomial); summary(a)
a<-glm(Som$Ts ~ Som$TurnOverRate + Som$AF + Som$mitoReads + Som$percentMito + Som$Position + Som$dummy20 + Som$dummy30 + Som$dummy40 + Som$dummy50 + Som$dummy60 + Som$DummyMale, family = binomial); summary(a)

##### opposite - turnover as a function of substitutions:

a<-lm(Som$TurnOverRate ~ Som$G_A + Som$T_C + Som$C_T + Som$A_G + Som$AF + Som$mitoReads + Som$percentMito+ Som$dummy20 + Som$dummy30 + Som$dummy40 + Som$dummy50 + Som$dummy60 + Som$DummyMale); summary(a)
a<-glm(Som$TurnOverRate ~ G_A + T_C + C_T + A_G + AF + mitoReads + percentMito, data = Som); summary(a)
