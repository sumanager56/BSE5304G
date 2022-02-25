options(repos ="http://cran.us.r-project.org")  # required to get latest libs
# Installing the packages we will play with today
if (!require("pacman")) install.packages("pacman")
pacman::p_load(elevatr,soilDB,rgdal,raster,ggplot2,patchwork)
# From previous weeks
pacman::p_load(EcoHydRology,rnoaa,curl,httr)
# 3 Functions to calculate SWE and excess when soil is drying, 
#   wetting, and wetting above capacity
#
#browseURL("https://github.com/vtdrfuka/BSE5304_2022/tree/main/functions")
#browseURL("https://github.com/vtdrfuka/BSE5304_2022/blob/main/functions/TMWBmodel.R")
#browseURL("https://raw.githubusercontent.com/vtdrfuka/BSE5304_2022/main/functions/NSE.R")
source("https://raw.githubusercontent.com/vtdrfuka/BSE5304_2022/main/functions/TMWBmodel.R")
#CN model here
CNmodel<-function(CNmodeldf, CNavg = 75,IaFrac = 0.05,fnc_slope=0, 
                  fnc_aspect=0,func_DAWC=.3,func_z=1000,fnc_fcres=.3) {
  CNmodeldf=modeldata
  #CNavg=75;
  #IaFrac=0.1; 
  #fnc_slope=0; fnc_aspect=0;func_DAWC=.3;
  #func_z=1000;fnc_fcres=.3
  # Energy Balance based Snow Accumulation 
  # and Melt model from the EcoHydRology package.
  attach(CNmodeldf)
  SNO_Energy=SnowMelt(date, P, MaxTemp-3, MinTemp-3, myflowgage$declat, 
                      slope = fnc_slope, aspect = fnc_aspect, tempHt = 1, 
                      windHt = 2, groundAlbedo = 0.25,SurfEmissiv = 0.95, windSp = 2, 
                      forest = 0, startingSnowDepth_m = 0,startingSnowDensity_kg_m3=450)
  # We will update the -3 in the above to be a lapse rate adjustment
  detach(CNmodeldf)
  CNmodeldf$SNO=SNO_Energy$SnowWaterEq_mm
  CNmodeldf$SNOmlt=SNO_Energy$SnowMelt_mm
  CNmodeldf$SnowfallWatEq_mm=SNO_Energy$SnowfallWatEq_mm
  CNmodeldf$SnowMelt_mm=SNO_Energy$SnowMelt_mm
  attach(CNmodeldf)
  CNmodeldf$Albedo=.23
  CNmodeldf$Albedo[CNmodeldf$SNO>0]=.95
  PET=PET_fromTemp(Jday=(1+as.POSIXlt(date)$yday),
                   Tmax_C = MaxTemp,Tmin_C = MinTemp,
                   lat_radians = myflowgage$declat*pi/180) * 1000
  CNmodeldf$PET=PET
  detach(CNmodeldf)
  rm(list="PET")
  
  CNmodeldf$AWC=func_DAWC*func_z
  # Oh, this we want to vary some of these around our watershed!
  CNmodeldf$dP = 0 # Initializing Net Precipitation
  CNmodeldf$ET = 0 # Initializing ET
  CNmodeldf$AW = 0 # Initializing AW
  CNmodeldf$Excess = 0 # Initializing Excess
  CNmodeldf$S =0 # Initializing S
  CNmodeldf$Qpred=0 # Initializing Qpred
  attach(CNmodeldf)
  SSCNavg=(1000/CNavg-10)*25.4
  SSCN=SoilStorage(S_avg=SSCNavg, field_capacity=func_DAWC*.9,
                   soil_water_content=0.1*func_DAWC, porosity=func_DAWC)
  Ia_init=IaFrac*SSCN   
  CNmodeldf$CNavg = CNavg
  CNmodeldf$SSCNavg = SSCNavg
  CNmodeldf$SSCN = SSCN
  detach(CNmodeldf)
  rm(list=c("CNavg", "SSCN", "SSCNavg"))
  CNmodeldf$Ia = Ia_init
  attach(CNmodeldf)
  # Those processes that are dependant on prior days conditions, we run as a 
  # loop through each of the days.
  for (t in 2:length(AW)){
    ET[t] = AW[t-1]/AWC[t-1]*PET[t]
    # Calculating Net Precipitation which adds in slope above's Excess
    dP[t] = SNO_Energy$Rain_mm[t] - ET[t] + 
      SNO_Energy$SnowMelt_mm[t]    # CN Solution
    # Is the soil saturated, and thus can't take more dP? 
    if (AW[t-1] + dP[t]>=AWC[t]){
      Excess[t]=AW[t-1] + dP[t] -AWC[t]
      AW[t]=AWC[t]
      # Otherwise, if dP is less than the initial abstraction? 
      # https://en.wikipedia.org/wiki/Runoff_curve_number#Definition
    } else if (dP[t]<=Ia[t]) {
      Excess[t]=0.0
      AW[t]=AW[t-1] + dP[t]
    } else {
      Excess[t]=(dP[t]-Ia[t])^2/(dP[t]-Ia[t]+SSCN[t])
      AW[t]=AW[t-1] + dP[t] -Excess[t]
    }
    S[t]=S[t-1]+Excess[t]
    Qpred[t]=fnc_fcres*S[t]
    S[t]=S[t]-Qpred[t]
  }
  CNmodeldf$ET=ET
  CNmodeldf$dP=dP
  CNmodeldf$AW=AW
  CNmodeldf$Excess=Excess
  CNmodeldf$S=S
  CNmodeldf$Qpred=Qpred # UPDATE vector BEFORE DETACHING
  rm(list=c("AW", "dP", "ET", "Excess", "Qpred", "S"))
  detach(CNmodeldf)
  return(CNmodeldf)
}
# Download a soils dataset for your basin based on the WebSoilSurvey method 
# and replace this url with your own
#url="https://websoilsurvey.sc.egov.usda.gov/DSD/Download/AOI/toxbqa5gcd1swjxqyryf50pz/wss_aoi_2022-02-24_10-23-42.zip"
#download.file(url,"mysoil.zip")
#unzip("mysoil.zip")
# https://cran.r-project.org/web/packages/elevatr/elevatr.pdf
# https://cran.r-project.org/web/packages/raster/raster.pdf
# https://cran.r-project.org/web/packages/soilDB/soilDB.pdf
# https://cran.r-project.org/web/packages/rgdal/rgdal.pdf
# use the function to get data from USGS 02049500, 
#BLACKWATER RIVER NEAR FRANKLIN, VA
myflowgage_id="02047000"
myflowgage=get_usgs_gage(myflowgage_id,begin_date = "2015-01-01",end_date = "2022-03-01")
# Note that flow returned is in m3/day, but we want mm/day for the basin
myflowgage$flowdata$Qmm = myflowgage$flowdata$flow/myflowgage$area/10^3

stns=meteo_distance(
  station_data=ghcnd_stations(),
  lat=myflowgage$declat,
  long=myflowgage$declon,
  units = "deg",
  radius = 30,
  limit = NULL
)
# We are looking for stations with elements that have PRCP, TMAX and TMIN 
# and current data (i.e. Year 2022). 3
WXStn=stns[stns$element=="TMAX"&stns$last_year>=2021,]$id[1]
WXData=meteo_pull_monitors(
  monitors=WXStn,
  keep_flags = FALSE,
  date_min = "2016-01-01",
  date_max = NULL,
  var = c("TMAX","TMIN","PRCP") 
)
summary(WXData)  #

# Create an aligned modeldata data frame to build our model in
modeldata=merge(WXData,myflowgage$flowdata,by.x="date",by.y="mdate")
summary(modeldata)  #
modeldata$MaxTemp=modeldata$tmax/10 # Converting to C
modeldata$MinTemp=modeldata$tmin/10 # Converting to C
modeldata$P=modeldata$prcp/10 # Converting to mm
# View(modeldata)  
# Compare your precipitation to the flow out of your basin
modeldata$P[is.na(modeldata$P)]=0
modeldata$MinTemp[is.na(modeldata$MinTemp)]=0
modeldata$MaxTemp[is.na(modeldata$MaxTemp)]=
  modeldata$MinTemp[is.na(modeldata$MaxTemp)] +1
modeldata$MaxTemp[modeldata$MaxTemp<=modeldata$MinTemp]=
  modeldata$MinTemp[modeldata$MaxTemp<=modeldata$MinTemp]+1
modeldata$AvgTemp=(modeldata$MaxTemp+modeldata$MinTemp)/2.0

summary(modeldata)
modeldata[is.na(modeldata)]=0 # A Quick BUT sloppy removal of NAs
TMWB=modeldata
# Last weeks homework example
source("https://raw.githubusercontent.com/vtdrfuka/BSE5304_2022/main/functions/TMWBmodel.R")
# Calibrating the parameters one at a time
for (fcres in seq(.1,.5,.1)){
  TMWBnew=TMWBmodel(TMWB=TMWB,fcres=fcres)
  print(paste(fcres,NSE(TMWBnew$Qmm,TMWBnew$Qpred)))
}
# fcres=.1 provides the best NSE
for (SFTmp in seq(-5,20)){
  TMWBnew=TMWBmodel(TMWB=TMWB,fcres=.1,SFTmp = SFTmp)
  print(paste(SFTmp,NSE(TMWBnew$Qmm,TMWBnew$Qpred)))
}
# SFTMp=10 provides the best NSE
for(AWCval in seq(50,350,50)){
  TMWBnew=TMWBmodel(TMWB=TMWB,fcres=.1,SFTmp = 10,Tlag = .5,AWCval = AWCval)
  print(paste(AWCval,NSE(TMWBnew$Qmm,TMWBnew$Qpred)))
}
# AWCval=350 provides the best NSE
# Best result for "LICK RUN ABOVE PATTON AVENUE AT ROANOKE, VA" NSE = .34 . 
TMWBnew=TMWBmodel(TMWB=TMWB,fcres=.1,SFTmp = 10,Tlag = .5,AWCval = 350)
print(paste(AWCval,NSE(TMWBnew$Qmm,TMWBnew$Qpred)))


# Model Performance 
plot(TMWBnew$date,TMWBnew$Qpred,type="l")
NSeff(TMWBnew$Qmm,TMWBnew$Qpred)


#CN model calling
CNmodeldf <- CNmodel(CNmodeldf=modeldata)
plot(CNmodeldf$date,CNmodeldf$Qpred,type="l")
NSeff(CNmodeldf$Qmm,CNmodeldf$Qpred)

#Trying to calibrate CNmodeldf
for (CNavg in seq(40,95,5)){
  FinalmodelCNnew=CNmodel(CNmodeldf=modeldata,CNavg=CNavg)
  print(paste(CNavg,NSE(FinalmodelCNnew$Qmm,FinalmodelCNnew$Qpred)))
}

#Best NSE=-4.36 obtained using CNavg=85
for (IaFrac in seq(0.05,0.5,0.05)){
  FinalmodelCNnew=CNmodel(CNmodeldf=modeldata,CNavg=85,IaFrac=IaFrac)
  print(paste(IaFrac,NSE(FinalmodelCNnew$Qmm,FinalmodelCNnew$Qpred)))
}
#Best NSE obtained using IaFrac=0.05
for (fnc_fcres in seq(0.1,0.5,0.05)){
  FinalmodelCNnew=CNmodel(CNmodeldf=modeldata,CNavg=85,IaFrac=0.05,fnc_fcres=fnc_fcres)
  print(paste(fnc_fcres,NSE(FinalmodelCNnew$Qmm,FinalmodelCNnew$Qpred)))
}
#Best NSE (-0.72) obtained using fnc_fcres=0.1 


#plotting using ggplot
colors <- c("Predicted flow" = "magenta", "Observed flow" = "cyan")
Date <- FinalmodelCNnew$date
Qpred <- FinalmodelCNnew$Qpred
Qobs <- FinalmodelCNnew$Qmm
#
p1<- ggplot() +
  # Plot your discharge data
  geom_line(aes(x=Date, y = Qpred), color="magenta", size=0.5)+
  geom_line(aes(x=Date, y = Qobs), color="cyan",size=0.5) +
  labs(x = "",
       y = "Discharge, mm") +
  ggtitle("CN model - NOTTOWAY RIVER NEAR SEBRELL, VA")+
  scale_color_identity(name = "Legend",
                       breaks = c("magenta", "cyan"),
                       labels = c("Predicted", "Observed"),
                       guide = "legend")+
  theme(
    plot.title = element_text(size=12, face="bold"),
    axis.title.y = element_text(size=11, face="bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size=11, face="bold"),
    axis.line.x = element_line(color="black",size=0.3),
    axis.line.y = element_line(color="black",size=0.3),
    panel.border = element_rect(color="black",fill=NA,size=0.3)
  )

#Plotting obs vs predicted for TMWB model

Date <- TMWBnew$date
Qobs <- TMWBnew$Qmm
Qpred_TMWB <- TMWBnew$Qpred
p2 <- ggplot() +
  # Plot your discharge data
  geom_line(aes(x=Date, y = Qpred_TMWB, color="magenta"),size=0.5)+
  geom_line(aes(x=Date, y = Qobs, color="cyan"),size=0.5) +
  labs(x = "Date",
       y = "Discharge, mm") +
  ggtitle("TMWB model")+
  scale_color_identity(name = "Legend",
                       breaks = c("magenta", "cyan"),
                       labels = c("Predicted", "Observed"),
                       guide = "legend")+
  theme(
    plot.title = element_text(size=12, face="bold"),
    axis.title.y = element_text(size=11, face="bold"),
    axis.title.x = element_text(size=11, face="bold"),
    axis.text.x = element_text(size=11, face="bold"),
    #axis.ticks.x = element_blank(),
    axis.text.y = element_text(size=11, face="bold"),
    axis.line.x = element_line(color="black",size=0.3),
    axis.line.y = element_line(color="black",size=0.3),
    panel.border = element_rect(color="black",fill=NA,size=0.3)
  )

#dev.off()
p1 + p2 +  plot_layout(nrow = 2) 

#Plotting Snowmelt vs snowwater equivalent for CN model
colors <- c("Snowmelt" = "magenta", "Snowwater_equivalent" = "cyan")
Date <- FinalmodelCNnew$date
SnowMelt <- FinalmodelCNnew$SnowMelt_mm
SNW_eq <- FinalmodelCNnew$SnowfallWatEq_mm
#
p1<- ggplot() +
  # Plot your discharge data
  geom_line(aes(x=Date, y = SnowMelt), color="magenta", size=0.5)+
  geom_line(aes(x=Date, y = SNW_eq), color="cyan",size=0.5) +
  labs(x = "",
       y = "Depth, mm") +
  ggtitle("CN model - SnowMelt vs Snowwater equivalent")+
  scale_color_identity(name = "Legend",
                       breaks = c("magenta", "cyan"),
                       labels = c("SnowMelt", "SNW_eq"),
                       guide = "legend")+
  theme(
    plot.title = element_text(size=12, face="bold"),
    axis.title.y = element_text(size=11, face="bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size=11, face="bold"),
    axis.line.x = element_line(color="black",size=0.3),
    axis.line.y = element_line(color="black",size=0.3),
    panel.border = element_rect(color="black",fill=NA,size=0.3)
  )

#Plotting Snowmelt vs snowwater equivalent for TMWB model
Date <- TMWBnew$date
SnowMeltTMWB <- TMWBnew$SNOmlt
SNW_eqTMWB <- TMWBnew$SNO

p2 <- ggplot() +
  # Plot your discharge data
  geom_line(aes(x=Date, y = SnowMeltTMWB, color="magenta"),size=0.5)+
  geom_line(aes(x=Date, y = SNW_eqTMWB, color="cyan"),size=0.5) +
  labs(x = "Date",
       y = "Depth, mm") +
  ggtitle("TMWB model - SnowMelt vs Snowwater equivalent")+
  scale_color_identity(name = "Legend",
                       breaks = c("magenta", "cyan"),
                       labels = c("SnowMelt", "SNW_eq"),
                       guide = "legend")+
  theme(
    plot.title = element_text(size=12, face="bold"),
    axis.title.y = element_text(size=11, face="bold"),
    axis.title.x = element_text(size=11, face="bold"),
    axis.text.x = element_text(size=11, face="bold"),
    #axis.ticks.x = element_blank(),
    axis.text.y = element_text(size=11, face="bold"),
    axis.line.x = element_line(color="black",size=0.3),
    axis.line.y = element_line(color="black",size=0.3),
    panel.border = element_rect(color="black",fill=NA,size=0.3)
  )

#dev.off()
p1 + p2 +  plot_layout(nrow = 2)
