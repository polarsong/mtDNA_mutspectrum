BodyMassKg = seq(1,1000,1) # kilogramms

##### table 1, Lindstedt 1981, Body size, physiological time, and longevity of homeothermic animals
LifeSpanInCaptivityMinuts = 6.10*10^6*BodyMassKg^0.20 # minutes versus Kg
plot(log(BodyMassKg),log(LifeSpanInCaptivityMinuts))
plot(BodyMassKg,LifeSpanInCaptivityMinuts/525600) # FromMinutesToYear 60*24*365 = 525600

##### back to BodyMass from LifeSpanInCaptivity
BodyMassKgPred = ((LifeSpanInCaptivityMinuts/525600)/(6.10*10^6))^(1/0.2)
plot(BodyMassKg,BodyMassKgPred) # FromMinutesToYear 60*24*365 = 525600

LifeSpan = BodyMass^(1/4); # if we add any linear coefficient => it doesn't changes
plot(log(BodyMass),log(GenTime))

help(exp)

#### eq 8, Lindstedt 1981, Body size, physiological time, and longevity of homeothermic animals
TimeOfMat = 2.93*10^5*BodyMass^0.18 # minutes !! mistake!! should be 10^6!!!
plot(log(BodyMass),log(TimeOfMat))
plot(BodyMass,TimeOfMat/525600) # FromMinutesToYear 60*24*365 = 525600

TimeOfMat = 2.93*10^6*BodyMass^0.18 # minutes !! mistake!! should be 10^6!!!
plot(BodyMass,TimeOfMat/525600) # FromMinutesToYear 60*24*365 = 525600
