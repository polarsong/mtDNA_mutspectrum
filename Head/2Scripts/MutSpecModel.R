#### 0: one mutation per unit of time, but unit of time is different for three categories of substitutions:

#C_T => AbsoluteTime (in days)
#A_G => NormoxicTime (in the number of days, necessy to have oxygen enough for one A>G mutation: high number for hypoxic, low number for normoxic)
#TheRest => ReplicationTime (in days, needed for replication)

# NormoxicTime miht be ~ 1/ReplicationTime

ReplicationTime = 10
NormoxicTime = 5

EndOfTime = 100 # years
C_T = 0
A_G = 0
TheRest = 0
for (i in (1:EndOfTime))
{
  C_T = C_T + 1
  if (i/NormoxicTime == (round(i/NormoxicTime)))   {A_G = A_G + 1}
  if (i/ReplicationTime == (round(i/ReplicationTime)))   {TheRest = TheRest + 1}
  
  OneLine = data.frame(C_T,A_G,TheRest)
  if (i == 1) {Final = OneLine}
  if (i >  1) {Final = rbind(Final,OneLine)}
}

Final$days = row.names(Final)
Final$Fr.C_T = Final$C_T/(Final$C_T + Final$A_G + Final$TheRest)
Final$Fr.A_G = Final$A_G/(Final$C_T + Final$A_G + Final$TheRest)
Final$Fr.TheRest = Final$TheRest/(Final$C_T + Final$A_G + Final$TheRest)

plot(Final$days,Final$Fr.C_T, ylim = c(0,1), xlim = c(1,EndOfTime), col = 'green', ylab = 'freq'); par(new=TRUE) # 
plot(Final$days,Final$Fr.A_G, ylim = c(0,1), xlim = c(1,EndOfTime), col = 'red', ylab = 'freq'); par(new=TRUE)
plot(Final$days,Final$Fr.TheRest, ylim = c(0,1), xlim = c(1,EndOfTime), col = 'grey', ylab = 'freq');

#### 1: different tissues have different ReplicationTime and NormoxicTime  

ReplicationTimeVec = c(1,10,20,30,40,50,60,70,80,90,100)
NormoxicTimeVec = c(1,10,20,30,40,50,60,70,80,90,100)
for (count1 in 1:length(ReplicationTimeVec))
  { # count1 = 1
  for (count2 in 1:length(NormoxicTimeVec))
    { # count2 = 1
    ReplicationTime = ReplicationTimeVec[count1]
    NormoxicTime = NormoxicTimeVec[count2]
    EndOfTime = 100 # years
    C_T = 0
    A_G = 0
    TheRest = 0
    for (i in (1:EndOfTime))
    {
      C_T = C_T + 1
      if (i/NormoxicTime == (round(i/NormoxicTime)))   {A_G = A_G + 1}
      if (i/ReplicationTime == (round(i/ReplicationTime)))   {TheRest = TheRest + 1}
  
      TempOneLine = data.frame(ReplicationTime,NormoxicTime,C_T,A_G,TheRest)
      if (i == 1) {TempFinal = TempOneLine}
      if (i >  1) {TempFinal = rbind(TempFinal,TempOneLine)}
    }

  TempFinal$days = row.names(TempFinal)
  TempFinal$Fr.C_T = TempFinal$C_T/(TempFinal$C_T + TempFinal$A_G + TempFinal$TheRest)
  TempFinal$Fr.A_G = TempFinal$A_G/(TempFinal$C_T + TempFinal$A_G + TempFinal$TheRest)
  TempFinal$Fr.TheRest = TempFinal$TheRest/(TempFinal$C_T + TempFinal$A_G + TempFinal$TheRest)
  TempFinal$TotalMut = TempFinal$C_T + TempFinal$A_G + TempFinal$TheRest
  
  if (count1 ==1 & count2 == 1) {Final = TempFinal}
  if (count1 > 1 | count2 > 1)  {Final = rbind(Final,TempFinal)}
  }
}

ToDraw = Final[Final$days == 100,]
plot(ToDraw$ReplicationTime,ToDraw$NormoxicTime, cex = c(ToDraw$Fr.A_G)*5)
plot(ToDraw$ReplicationTime,ToDraw$NormoxicTime, cex = c(ToDraw$Fr.C_T)*5)
plot(ToDraw$ReplicationTime,ToDraw$NormoxicTime, cex = c(ToDraw$Fr.TheRest)*5)
plot(ToDraw$ReplicationTime,ToDraw$NormoxicTime, cex = c(ToDraw$TotalMut)/20)


ToDraw2 = ToDraw[ToDraw$ReplicationTime == 1,]
par(mfrow=c(1,2))
plot(ToDraw2$NormoxicTime,ToDraw2$Fr.A_G, col = 'green', ylim = c(0,1), pch = 16, ylab = 'fractions'); par(new = TRUE);
plot(ToDraw2$NormoxicTime,ToDraw2$Fr.C_T, col = 'red', ylim = c(0,1), ylab = 'fractions'); par(new = TRUE);
plot(ToDraw2$NormoxicTime,ToDraw2$Fr.TheRest, col = 'black', ylim = c(0,1), ylab = 'fractions');
plot(ToDraw2$NormoxicTime,ToDraw2$TotalMut)

ToDraw2 = ToDraw[ToDraw$NormoxicTime == 1,]
par(mfrow=c(1,2))
plot(ToDraw2$ReplicationTime,ToDraw2$Fr.A_G, col = 'green', ylim = c(0,1), pch = 16, ylab = 'fractions'); par(new = TRUE);
plot(ToDraw2$ReplicationTime,ToDraw2$Fr.C_T, col = 'red', ylim = c(0,1), ylab = 'fractions'); par(new = TRUE);
plot(ToDraw2$ReplicationTime,ToDraw2$Fr.TheRest, col = 'black', ylim = c(0,1), ylab = 'fractions');
plot(ToDraw2$ReplicationTime,ToDraw2$TotalMut)


### 2: ReplicationTime and NormoxicTime can change during agind / cancer development


##### OOOOOLLLLLLLDDDDDDDD

#### 1: simple model 

CoeffOfTime =  0.001 # 3/36000 => so that 3 mutations in mtDNA per lifetime
CoeffOfNormoxia = 0.0005
CellLongevity = 50
ReplMutagenesis = 0.0001
EndOfTime = 36000 # days ~ human lifespan 360
C_T = 0
A_G = 0
TheRest = 0
for (i in (1:EndOfTime))
{
  C_T = C_T + CoeffOfTime
  A_G = A_G + CoeffOfNormoxia
  if (i/CellLongevity == (round(i/CellLongevity)))
  {TheRest = TheRest + ReplMutagenesis}
  OneLine = data.frame(C_T,A_G,TheRest)
  if (i == 1) {Final = OneLine}
  if (i >  1) {Final = rbind(Final,OneLine)}
}

Final$days = row.names(Final)
Final$Fr.C_T = Final$C_T/(Final$C_T + Final$A_G + Final$TheRest)
Final$Fr.A_G = Final$A_G/(Final$C_T + Final$A_G + Final$TheRest)
Final$Fr.TheRest = Final$TheRest/(Final$C_T + Final$A_G + Final$TheRest)

plot(Final$days,Final$Fr.C_T, ylim = c(0,1), xlim = c(1,EndOfTime), col = 'green', ylab = 'freq'); par(new=TRUE) # 
plot(Final$days,Final$Fr.A_G, ylim = c(0,1), xlim = c(1,EndOfTime), col = 'red', ylab = 'freq'); par(new=TRUE)
plot(Final$days,Final$Fr.TheRest, ylim = c(0,1), xlim = c(1,EndOfTime), col = 'grey', ylab = 'freq');


#### 2: with age cells start to divide slower and they become more normoxic

CoeffOfTime =  0.001 # 3/36000 => so that 3 mutations in mtDNA per lifetime
CoeffOfNormoxia = 0.0005
CellLongevity = 50
ReplMutagenesis = 0.0001
EndOfTime = 36000 # days ~ human lifespan 360
C_T = 0
A_G = 0
TheRest = 0
for (i in (1:EndOfTime))
{
  C_T = C_T + CoeffOfTime
  A_G = A_G + CoeffOfNormoxia + i*0.0001*CoeffOfNormoxia # cells are becoming more normoxic 
  if (i/CellLongevity == (round(i/CellLongevity)))
  {TheRest = TheRest + ReplMutagenesis - ReplMutagenesis*i*0.0001}
  OneLine = data.frame(C_T,A_G,TheRest)
  if (i == 1) {Final = OneLine}
  if (i >  1) {Final = rbind(Final,OneLine)}
}

Final$days = row.names(Final)
Final$Fr.C_T = Final$C_T/(Final$C_T + Final$A_G + Final$TheRest)
Final$Fr.A_G = Final$A_G/(Final$C_T + Final$A_G + Final$TheRest)
Final$Fr.TheRest = Final$TheRest/(Final$C_T + Final$A_G + Final$TheRest)
Final$Total = Final$C_T + Final$A_G + Final$TheRest

plot(Final$days,Final$Fr.C_T, ylim = c(0,1), xlim = c(1,EndOfTime), col = 'green', ylab = 'freq'); par(new=TRUE) # 
plot(Final$days,Final$Fr.A_G, ylim = c(0,1), xlim = c(1,EndOfTime), col = 'red', ylab = 'freq'); par(new=TRUE)
plot(Final$days,Final$Fr.TheRest, ylim = c(0,1), xlim = c(1,EndOfTime), col = 'grey', ylab = 'freq');

Ymin = 0; Ymax=max(Final$Total);
plot(Final$days,Final$Total, xlim = c(1,EndOfTime), ylim=c(Ymin,Ymax), col = 'black', ylab = 'total'); par(new = TRUE)
plot(Final$days,Final$C_T, xlim = c(1,EndOfTime), ylim=c(Ymin,Ymax), col = 'green', ylab = 'total'); par(new = TRUE)
plot(Final$days,Final$A_G, xlim = c(1,EndOfTime), ylim=c(Ymin,Ymax), col = 'red', ylab = 'total'); par(new = TRUE)
plot(Final$days,Final$TheRest, xlim = c(1,EndOfTime), ylim=c(Ymin,Ymax), col = 'grey', ylab = 'total');


# if CoeffOfNormoxia is negatively proportional to CellLongevity - do I see that fraction of A>G is higher in slow dividing tissues?
# add absolute number of mutations:

