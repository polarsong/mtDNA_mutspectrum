rm(list=ls(all=TRUE))

Time = seq(1:100)
GaIntercept = 9
TcIntercept = 0
GaSlope = 0.01
TcSlope = 0.05

pdf("../../Body/4Figures/MutSpecModelA.R01.pdf")
par(mfrow=c(2,2))
MutSpec = data.frame(Time)
MutSpec$G_A = GaIntercept + MutSpec$Time*GaSlope
MutSpec$T_C = TcIntercept + MutSpec$Time*TcSlope
plot(MutSpec$Time,MutSpec$G_A, ylim =c(0,max(c(MutSpec$G_A,MutSpec$T_C))), col='red', ylab = '', xlab = ''); par(new=TRUE)
plot(MutSpec$Time,MutSpec$T_C, ylim =c(0,max(c(MutSpec$G_A,MutSpec$T_C))), col='green', ylab = 'substitution probability', xlab = 'time being single stranded');
MutSpec$FrGa = MutSpec$G_A/(MutSpec$G_A + MutSpec$T_C)
MutSpec$FrTc = MutSpec$T_C/(MutSpec$T_C + MutSpec$G_A)
MutSpec$FrTcGa = MutSpec$FrTc/(MutSpec$FrGa) 
plot(MutSpec$Time,MutSpec$FrGa, ylim =c(0,max(c(MutSpec$FrGa,MutSpec$FrTc))), col='red', ylab = '', xlab = ''); par(new=TRUE)
plot(MutSpec$Time,MutSpec$FrTc, ylim =c(0,max(c(MutSpec$FrGa,MutSpec$FrTc))), col='green', ylab = 'fraction of the substitution', xlab = 'time being single stranded');
MutSpec$FrTcGa = MutSpec$FrTc/(MutSpec$FrGa) 
plot(MutSpec$Time,MutSpec$FrTcGa, col='grey', ylab = '', xlab = '');
MutSpec$FrTcGa = MutSpec$FrTc/(MutSpec$FrGa + MutSpec$FrTc) 
plot(MutSpec$Time,MutSpec$FrTcGa, col='grey', ylab = '', xlab = '');

dev.off()