# Solution to HW 2
if (!require("pacman")) install.packages("pacman")
pacman::p_load(rnoaa,EcoHydRology,lattice)
# For a station near Blacksburg we look at
# https://maps.waterdata.usgs.gov/mapper/index.html
# or https://waterdata.usgs.gov/nwis/rt
# and like gage:

# Finding weather station with names:
# waterdata.usgs.gov
myflowgage_id="03173000"
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
WXStn=stns[stns$element=="TMAX"&stns$last_year>=2021,]$id[2]
WXData=meteo_pull_monitors(
  monitors=WXStn,
  keep_flags = FALSE,
  date_min = "2016-01-01",
  date_max = NULL,
  var = c("TMAX","TMIN","PRCP")
)
summary(WXData)  #
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

# Add some some soil parameters to be associated with the area 
# above the flow gage.

myflowgage$FldCap=.45
myflowgage$WiltPt=.15
myflowgage$Z=1000
TMWB$AWC=(myflowgage$FldCap-myflowgage$WiltPt)*myflowgage$Z # 


TMWB$PET = mean(TMWB$P,na.rm=T)-mean(TMWB$Qmm,na.rm=T)  # in mm/day
TMWB$ET = TMWB$PET # in mm/day
TMWB$dP = TMWB$P - TMWB$PET
TMWB$AW=NA  #Assigns all values in column with “NA” (Not available)
TMWB$AW[1]=250
TMWB$Excess=NA
TMWB$Excess[1]=0
head(TMWB)
attach(TMWB)
for (t in 2:length(date)){
  if (dP[t]< 0) {  
    values<-soildrying(AW[t-1],dP[t],AWC[t])
  } else if (AW[t-1]+dP[t]>AWC[t]) {
    values<-soil_wetting_above_capacity(AW[t-1],dP[t],AWC[t])
  } else {
    values<-soilwetting (AW[t-1],dP[t],AWC[t])
  }
  AW[t]<-values[1]
  Excess[t]<-values[2]
  
}
detach(TMWB)
TMWB$AW <-AW
TMWB$Excess<-Excess
rm(list=c("AW","Excess"))
TMWB$Qpred=NA
TMWB$Qpred[1]=0
TMWB$S=NA
TMWB$S[1]=0

attach(TMWB)
fcres=.3   # reservoir coefficient
for (t in 2:length(date)){
  S[t]=S[t-1]+Excess[t]     
  Qpred[t]=fcres*S[t]
  S[t]=S[t]-Qpred[t]
}
detach(TMWB) # IMPORTANT TO DETACH
TMWB$S=S
TMWB$Qpred=Qpred # UPDATE vector BEFORE DETACHING
rm(list=c("S","Qpred"))
summary(TMWB)
# View(TMWB)
dev.off()
####################################
###let's start to plot
#####################################
#pcp and Qmm and Qpred
maxRange <- 1.1*(max(TMWB$P,na.rm = T) + max(TMWB$Qmm,na.rm = T))
p1<- ggplot() +
  # Use geom_tile to create the inverted hyetograph. geom_tile has a bug that displays a warning message for height and width, you can ignore it.
  geom_tile(data = TMWB, aes(x=date,y = -1*(P/2-maxRange), # y = the center point of each bar
                               height = P,
                               width = 1),
            fill = "black",
            color = "black") +
  # Plot your discharge data
  geom_line(data=TMWB,aes(x=date, y = Qmm, colour ="Qmm"), size=1) +
  geom_line(data=TMWB,aes(x=date, y = Qpred, colour= "Qpred"), size=1) +
  scale_colour_manual("", 
                      breaks = c("Qmm", "Qpred"),
                      values = c("red", "blue")) +
  # Create a second axis with sec_axis() and format the labels to display the original precipitation units.
  scale_y_continuous(name = "Discharge (mm/day)",
                     sec.axis = sec_axis(trans = ~-1*(.-maxRange),
                                         name = "Precipitation (mm/day)"))+
  scale_x_continuous(name = NULL,labels = NULL)+
  ggtitle(myflowgage$gagename)

#ET and Excess
p2 <- ggplot(TMWB, aes(x=date)) +
  geom_line(aes(y=Excess, colour="Excess"), size=1)+
  geom_line(aes(y=ET, colour="ET"), size=1) +
  scale_colour_manual("", 
                      breaks = c("ET", "Excess"),
                      values = c("red", "blue")) +
  scale_y_continuous(name = "Depth (mm/day)",) +
  scale_x_continuous(name = NULL,labels = NULL) 
#FOR AW
p3 <- ggplot(TMWB, aes(x=date)) +
  geom_line(aes(y=AW,colour="AW"), size=1) +
  scale_colour_manual("", 
                      breaks = c("AW"),
                      values = c("black")) +
  scale_y_continuous(
    # Features of the first axis
    name = "AW (mm)",
    
  )

p1 + p2 + p3+ plot_layout(ncol = 1, widths = c(2,2,1))



SFTmp = 1  # referred to as SFTMP in SWAT input (Table 1)
bmlt6 = 4.5   # referred to as SMFMX in SWAT input (Table 1)
bmlt12 = 0.0  # referred to as SMFMN in SWAT input adjusted for season
Tmlt = SFTmp  # Assumed to be same as SnowFall Temperature
Tlag = 1  # referred to as TIMP in SWAT input (Table 1)
TMWB$AvgTemp=(TMWB$MaxTemp+TMWB$MinTemp)/2
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
plot(date,SNO,type="l")
detach(TMWB)
TMWB$Tsno=Tsno
TMWB$SNO=SNO
TMWB$SNOmlt=SNOmlt
rm(list=c("SNO", "SNOmlt", "Tsno"))


bmlt12 = 0.0
bmlt6 = 3
SFTmp = 7
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

# notice that there is an Energy Balance based Snow Accumulation 
# and Melt model in the EcoHydRology package.

?SnowMelt
attach(TMWB)
SNO_Energy=SnowMelt(date, P, MaxTemp-3, MinTemp-3, myflowgage$declat, 
                      slope = 0,
                      aspect = 0, tempHt = 1, windHt = 2, groundAlbedo = 0.25,
                      SurfEmissiv = 0.95, windSp = 2, forest = 0, startingSnowDepth_m = 0,
                      startingSnowDensity_kg_m3=450)
# How do we know what units slope and aspect are in? (view function)
lines(date,SNO_Energy$SnowWaterEq_mm)
detach(TMWB)

#edit(SnowMelt)

TMWB$Albedo=.23
TMWB$Albedo[TMWB$SNO>0]=.95
?PET_fromTemp
attach(TMWB)
PET=PET_fromTemp(Jday=(1+as.POSIXlt(date)$yday),Tmax_C = MaxTemp,Tmin_C = MinTemp,albedo=Albedo,lat_radians = myflowgage$declat*pi/180) * 1000
TMWB$PET=PET
plot(date,PET)
detach(TMWB)
rm(list=c("PET"))

TMWB$AWC=???? #Fld Cap = .45, Wilt Pt = .15, z=1000mm
myflowgage$FldCap=.45
myflowgage$WiltPt=.15
myflowgage$Z=1000
TMWB$AWC=(myflowgage$FldCap-myflowgage$WiltPt)*myflowgage$Z # 
TMWB$dP = 0 # Initializing Net Precipitation
TMWB$ET = 0 # Initializing ET
TMWB$AW = 0 # Initializing AW
TMWB$Excess = 0 # Initializing Excess


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
fcres=.3
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
