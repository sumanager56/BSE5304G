# Solution to HW 2
#if (!require("pacman")) install.packages("pacman")
pacman::p_load(rnoaa,EcoHydRology,lattice)
# For a station near Blacksburg we look at
# https://maps.waterdata.usgs.gov/mapper/index.html
# or https://waterdata.usgs.gov/nwis/rt
# and like gage:

# Finding weather station with names:
# waterdata.usgs.gov
myflowgage_id="04218518"
myflowgage=get_usgs_gage(myflowgage_id,
              begin_date="2016-01-01",end_date="2022-02-01")
# Now use the functions meteo_distance and ghcnd_stations to start 
# looking for weather stations nearby
# 
stns=meteo_distance(
  station_data=ghcnd_stations(),
  lat=myflowgage$declat,
  long=myflowgage$declon,
  units = "deg",
  radius = 30,
  limit = NULL
)
# We are looking for stations with elements that have PRCP, TMAX and TMIN 
# and current data (i.e. Year 2021). 
WXStn=stns[stns$element=="TMAX"&stns$last_year>=2021,]$id[1]
WXData=meteo_pull_monitors(
  monitors=WXStn,
  keep_flags = FALSE,
  date_min = "2016-01-01",
  date_max = NULL,
  var = c("TMAX","TMIN","PRCP")
)
summary(WXData)  
plot(WXData$date,WXData$prcp)
# Creat an aligned modeldata data frame to build our model in
modeldata=merge(WXData,myflowgage$flowdata,by.x="date",by.y="mdate")
summary(modeldata)  #
# Convert Units and Area normalize flow to match (depth)
# flow(m^3/day) / (area(km^2) * 1000^2m/km^2) * 1000mm/m = flow(mm/day)
modeldata$Qmm = modeldata$flow/myflowgage$area/10^3
# It is good practice to use similar object names to the 
# values in the documentation of the model (P,Q,MaxTemp,MinTemp)
modeldata$MaxTemp=modeldata$tmax/10 # Converting to C
modeldata$MinTemp=modeldata$tmin/10 # Converting to C
modeldata$P=modeldata$prcp/10 # Converting to mm
# View(modeldata)  
# Compare your precipitation to the flow out of your basin
mean(modeldata$Qmm)
mean(modeldata$P)
modeldata$P[is.na(modeldata$P)]=0
modeldata$MinTemp[is.na(modeldata$MinTemp)]=0
modeldata$MaxTemp[is.na(modeldata$MaxTemp)]=modeldata$MinTemp[is.na(modeldata$MaxTemp)] +1
summary(modeldata)

TMWB=modeldata
  
summary(TMWB)

# Building our Soil Wetting and Drying Functions
#
soilwetting<-function(AWprev,dP_func,AWC_func){
  AW_func<-AWprev+dP_func
  excess_func<-0.0
  c(AW_func,excess_func)
} 

soildrying<-function(AWprev,dP_func,AWC_func){
  AW_func=AWprev*exp(dP_func/AWC_func)
  excess_func<-0.0
  c(AW_func,excess_func)
}
# soil_wetting_above_capacity function
soil_wetting_above_capacity<-function(AWprev,dP_func,AWC_func){
  AW_func<-AWC_func
  excess_func<-AWprev+dP_func-AWC_func
  c(AW_func,excess_func)
}

#From here onwards, we will estimate PET differently.
TMWB$Albedo=.23
TMWB$Albedo[TMWB$SNO>0]=.95
?PET_fromTemp
attach(TMWB)
PET=PET_fromTemp(Jday=(1+as.POSIXlt(date)$yday),Tmax_C = MaxTemp,Tmin_C = MinTemp,albedo=Albedo,lat_radians = myflowgage$declat*pi/180) * 1000
TMWB$PET=PET
plot(date,PET)
detach(TMWB)
rm(list=c("PET"))

myflowgage$FldCap=.45
myflowgage$WiltPt=.15
myflowgage$Z=1000
#TMWB$AWC=(myflowgage$FldCap-myflowgage$WiltPt)*myflowgage$Z # 
TMWB$AWC=350
TMWB$dP = 0 # Initializing Net Precipitation
TMWB$ET = 0 # Initializing ET
TMWB$AW = 0 # Initializing AW
TMWB$Excess = 0 # Initializing Excess

#Some initializations based on SWAT model 
TMWB$AvgTemp=(TMWB$MaxTemp+TMWB$MinTemp)/2
SFTmp = 1  # referred to as SFTMP in SWAT input (Table 1)
bmlt6 = 7   # referred to as SMFMX in SWAT input (Table 1)
bmlt12 = 1.4  # referred to as SMFMN in SWAT input adjusted for season
Tmlt = SFTmp  # Assumed to be same as SnowFall Temperature
Tlag = 0.5  # referred to as TIMP in SWAT input (Table 1)
TMWB$bmlt = (bmlt6 + bmlt12)/2 + (bmlt6 - bmlt12)/2 *  sin(2*pi/365*(julian(TMWB$date,origin = as.Date("2000-01-01"))-81))
# Initialize SNO, Tsno as well as the first values of each
TMWB$SNO = 0  # Snow Depth (mm)
TMWB$Tsno = 0  # Snow Temp (C)
TMWB$SNOmlt = 0  # Snow Melt (mm)

attach(TMWB)
for (t in 2:length(date)){
  Tsno[t]= Tsno[t-1] * (1.0-Tlag) +  AvgTemp[t] * Tlag
  if(AvgTemp[t] < SFTmp){
    SNO[t]= SNO[t-1] + P[t]
  }  else {
    SNOmlt[t]= bmlt[t] * SNO[t-1] * ((Tsno[t]+MaxTemp[t])/2 - Tmlt) 
    SNOmlt[t]= min(SNOmlt[t],SNO[t-1])
    SNO[t]= SNO[t-1] -SNOmlt[t]
  }
  print(t)
}
lines(date,SNO,col="red")
detach(TMWB)
TMWB$Tsno=Tsno
TMWB$SNO=SNO
TMWB$SNOmlt=SNOmlt
rm(list=c("SNO", "SNOmlt", "Tsno"))


# Loop to calculate AW and Excess
attach(TMWB)
for (t in 2:length(AW)){
  # This is where Net Precipitation is now calculated
  # Do you remember what Net Precip is? Refer to week 2 notes
  ET[t] = min (AW[t-1],PET[t])
  ET[t] = (AW[t-1]/AWC[t-1])*PET[t] # New Model
  if(AvgTemp[t] >= SFTmp){
    dP[t] = P[t] - ET[t] + SNOmlt[t] 
  }  else {
    dP[t] = ET[t]
  }
  # From here onward, everything is the same as Week2’s lab
  if (dP[t]<=0) {
    values<-soildrying(AW[t-1],dP[t],AWC[t])
  } else if((dP[t]>0) & (AW[t-1]+dP[t])<=AWC[t]) {
    values<-soilwetting(AW[t-1],dP[t],AWC[t])
  } else {
    values<-soil_wetting_above_capacity(AW[t-1],dP[t],AWC[t])
  }
  AW[t]<-values[1]
  Excess[t]<-values[2]
  print(t)
}
TMWB$AW=AW
TMWB$Excess=Excess
TMWB$dP=dP
rm(list=c("AW","dP","ET", "Excess"))
detach(TMWB) # IMPORTANT TO DETACH


TMWB$Qpred=NA
TMWB$Qpred[1]=0
TMWB$S=NA
TMWB$S[1]=0
attach(TMWB)
fcres=.1   #changing from 0.1 to 0.5
for (t in 2:length(date)){
  S[t]=S[t-1]+Excess[t]     
  Qpred[t]=fcres*S[t]
  S[t]=S[t]-Qpred[t]
}
TMWB$S=S
TMWB$Qpred=Qpred # UPDATE vector BEFORE DETACHING

#Make a plot that has Qmm, P,and Qpred over time
plot(date,P,col="black")
lines(date,Qmm,type = "l",col="black")
lines(date,Qpred,col="blue")
detach(TMWB) # IMPORTANT TO DETACH
rm(list=c("Qpred","S"))

#Calculated NSE for calibrated model
NSE=function(Yobs,Ysim){
  return(1-sum((Yobs-Ysim)^2, na.rm=TRUE)/sum((Yobs-mean(Yobs, na.rm=TRUE))^2, na.rm=TRUE))
}
# Yobs and Ysim are vectors with the same length and order
# Yobs is observed streamflow
# Ysim is modeled

#Calling function to calculate  NSE
NSE(TMWB$Qmm,TMWB$Qpred)

