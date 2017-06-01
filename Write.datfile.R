# Attempting to run ADMB 

setwd("~/Documents/MyDocuments/GradSchoolUW/Dissertation/BaMMUpdates/TwoSourceR/EmilysStreams/Pioneer")
setwd("~/Documents/MyDocuments/GradSchoolUW/Dissertation/BaMMUpdates/2013Streams/Moose")

library(R2admb)

data=read.csv("Beaver.csv", sep=",", header=TRUE)
head(data)

## Write .dat file for BaMM
altitude=1211
aspect=0.001
DST=0
Latitude=44.6791
Longitude=-114.4217
Slope=0.001
SolarConst=1367
TimeZone=-8
Transmiss=0.8
Area=1
Depth=0.12
Salinity=0.00
alphaP=1
alphaGK=0.9972
SMeasFlag=1
LightSatFlag=2
SScale=0.27626576
RKstep=2

# Nobs
n_irr=length(data$TempC)
n_tC=length(data$TempC)
n_O2conc=length(data$TempC)
n_d18O=2

# O2, temp, light, time and date data
O2conc=data$O2Conc

DOY=data$JulianDate
TOD=data$DecimalTime
Interval=data$TimeSinceStart
irr=data$PAR
tC=data$TempC
delta_18O_O2=c(16,16)
DOY_18O=c(head(DOY, n=1), tail(DOY, n=1))
TOD_18O=c(head(TOD, n=1), tail(TOD, n=1))
Int_18O=c(head(Interval, n=1), tail(Interval, n=1))
endofdata=-999

write_dat("BaMM", list(altitude=altitude, aspect=aspect, DST=DST, Latitude=Latitude,
                       Longitude=Longitude,Slope=Slope, SolarConst=SolarConst,
                       TimeZone=TimeZone,Transmiss=Transmiss, Area=Area, Depth=Depth,
                       Salinity=Salinity, alphaP=alphaP, alphaGK=alphaGK, SMeasFlag=SMeasFlag,
                       LightSatFlag=LightSatFlag,SScale=SScale, RKstep=RKstep,
                       n_irr=n_irr, n_tC=n_tC, n_O2conc=n_O2conc, n_d18O=n_d18O, irr=irr, DOY, TOD, Interval,
                       temp=tC, DOY, TOD, Interval, O2data=O2conc,DOY,TOD,Interval,
                       O18data=delta_18O_O2,DOY_18O,TOD_18O,Int_18O, endofdata=endofdata))

