# Write .pin file that lists initial values for each estimated parameter

setwd("~/Documents/MyDocuments/GradSchoolUW/Dissertation/BaMMUpdates/TwoSourceR/EmilysStreams/Pioneer")
setwd("~/Documents/MyDocuments/GradSchoolUW/Dissertation/BaMMUpdates/2013Streams/Moose")

library(R2admb)

data=read.csv("Beaver.csv", sep=",", header=TRUE)
head(data)

#### Format for writing .pin file
# Initial values
Pmax=312
IHalfSat=2800
alphaPI=5.4
Rref=120
Eb=0.32
Ep=0.65
alphaProd=0.1
betaProd=0.7
k20=0.24
alphaR=1
delO18=-14.14
O2init=12.4
delO18init=21
sigmaO2=0.114
sigmadelO18=0.288

## Write .pin file for BaMM
write_pin('BaMM', list(Pmax=log(Pmax),IHalfSat=log(IHalfSat),alphaPI=log(alphaPI),
                       Rref=log(Rref),Eb=log(Eb),Ep=log(Ep),alphaProd=log(alphaProd),
                       betaProd=log(betaProd),k20=log(k20),alphaR=alphaR,delO18=delO18,O2=log(O2init),
                       initO18=log(delO18init),sigmaO2=log(sigmaO2),sigmaO18=log(sigmadelO18)))



