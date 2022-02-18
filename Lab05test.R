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

# Download a soils dataset for your basin based on the WebSoilSurvey method 
# and replace this url with your own
url="https://websoilsurvey.sc.egov.usda.gov/DSD/Download/AOI/bky5q4i2wkdjmvm3t0d1gniq/wss_aoi_2022-02-11_09-31-24.zip"
download.file(url,"mysoil.zip")
unzip("mysoil.zip")
# https://cran.r-project.org/web/packages/elevatr/elevatr.pdf
# https://cran.r-project.org/web/packages/raster/raster.pdf
# https://cran.r-project.org/web/packages/soilDB/soilDB.pdf
# https://cran.r-project.org/web/packages/rgdal/rgdal.pdf
# use the function to get data from USGS 0205551460 
#LICK RUN ABOVE PATTON AVENUE AT ROANOKE, VA
myflowgage_id="0205551460"
myflowgage=get_usgs_gage(myflowgage_id,begin_date = "2015-01-01",end_date = "2022-03-01")
# Note that flow returned is in m3/day, but we want mm/day for the basin
myflowgage$flowdata$Qmm = myflowgage$flowdata$flow/myflowgage$area/10^3
# the soil data using the soilDB package
# What is this returning? Why do we care?
# This needs to be completed based on your download
mysoil=readOGR("wss_aoi_2022-02-11_09-31-24/spatial/soilmu_a_aoi.shp")    
# Explore the mysoil dataset which is returned
mybbox=c(mysoil@bbox)
# First associate mukey with cokey from component
mysoil$mukey=mysoil$MUKEY  # or rename the column
mukey_statement = format_SQL_in_statement(unique(mysoil$mukey))
print(mukey_statement)
q_mu2co = paste("SELECT mukey,cokey FROM component WHERE mukey IN ", mukey_statement, sep="")
print(q_mu2co)
mu2co = SDA_query(q_mu2co)
# Second associate cokey with ksat_r,awc_r,hzdepb_r from chorizon
cokey_statement = format_SQL_in_statement(unique(mu2co$cokey))
q_co2ch = paste("SELECT cokey,ksat_r,awc_r,hzdepb_r  FROM chorizon WHERE cokey IN ", cokey_statement, sep="")
print(q_co2ch)
co2ch = SDA_query(q_co2ch)
# Last, bring them back together, and aggregate based on max values
# of ksat_r,awc_r, and hzdepb_r
mu2ch=merge(mu2co,co2ch)
summary(mu2ch)
mu2chmax=aggregate(mu2ch,list(mu2ch$mukey),max)

proj4_ll = "+proj=longlat"
proj4string(mysoil) = proj4_ll
mydem=get_elev_raster(locations=mysoil, 
                        z = 11, prj =proj4string(mysoil) ,
                        src ="aws",clip="bbox",expand = 0.001)

summary(terrain(mydem, opt='slope',unit = "degrees"))
# What is this 'slope'? Use the man page for the terrain() function to answer
plot(terrain(mydem, opt='TPI',unit = "degrees"))
# What is this 'TPI'? 
summary(terrain(mydem, opt='TRI',unit = "degrees"))
plot(terrain(mydem, opt='TRI',unit = "degrees"))
# What is this 'TRI'? 

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
# fcres=.3 is the highest NSE
for (SFTmp in seq(-5,20)){
  TMWBnew=TMWBmodel(TMWB=TMWB,fcres=.3,SFTmp = SFTmp)
  print(paste(SFTmp,NSE(TMWBnew$Qmm,TMWBnew$Qpred)))
}
for(AWCval in seq(50,350,50)){
  TMWBnew=TMWBmodel(TMWB=TMWB,fcres=.3,SFTmp = 11,Tlag = .5,AWCval = AWCval)
  print(paste(AWCval,NSE(TMWBnew$Qmm,TMWBnew$Qpred)))
}
# Best result for "LICK RUN ABOVE PATTON AVENUE AT ROANOKE, VA" NSE = .34 . 
TMWBnew=TMWBmodel(TMWB=TMWB,fcres=.3,SFTmp = 9,Tlag = .5,AWCval = 100)
print(paste(AWCval,NSE(TMWBnew$Qmm,TMWBnew$Qpred)))

# For a simple spatial model, use TMWB, to initialize 3 
# slope components for Top, Mid, and Bottom of the hillside. 

#Initializing the model and making a column to put data
summary(modeldata)
TopSlope=modeldata
MidSlope=modeldata
BotSlope=modeldata

# Running the model with the three HRU dataframes
# Low slope but highest ksat
#############################################
##GETTING ARGUMENT VALUES BASED ON EACH MODEL
#############################################
summary(mu2chmax)
##topslope has the lowest Z, lowest slope as botslope, midslope has the #highest slope, we assume that awc in % is 
#equal for all three models
awcpercent=mean(mu2chmax$awc_r,na.rm=TRUE)
### for TopSlope model lowest depth
TopslopeZ=min(mu2chmax$hzdepb_r)*10 #in mm
MidslopeZ=mean(mu2chmax$hzdepb_r)*10
BotslopeZ=max(mu2chmax$hzdepb_r)*10
#####################AWCVAL is mm of awc= depth*awc(%)
AWCvalTop=TopslopeZ*awcpercent
AWCvalMid=MidslopeZ*awcpercent
AWCvalBot=BotslopeZ*awcpercent
####################################calculation of slope from terrain #function,note that the unit is in degree
summary(terrain(mydem, opt='slope',unit = "degrees"))
SlopeTop=1.837328 #degree
SlopeBot=1.837328 #degree
SlopeMid=40.600004 #degree
# Running the model with the three HRU dataframes
# Low slope but highest ksat
# These 3 function calls are what you will vary for the Lab 04 homework 
TopSlope = TMWBmodel(TMWB = TopSlope,SFTmp = 1, 
                       AWCval = AWCvalTop,
                       Tlag = .5,fcres=.3,Slope = atan(SlopeTop/100))
MidSlope$P=TopSlope$Excess+MidSlope$P
# Higher slope, medium ksat, fcres=0.5 
MidSlope = TMWBmodel(TMWB = MidSlope,SFTmp = 1, 
                       AWCval = AWCvalMid,
                       Tlag = .5,fcres=0.5,Slope = atan(SlopeMid/100))
# Low Slope and lowest ksat, $fcres=0.2
BotSlope$P=MidSlope$Excess+BotSlope$P
BotSlope = TMWBmodel(TMWB = BotSlope,SFTmp = 1, 
                       AWCval = AWCvalBot,
                       Tlag = .5,fcres=0.2,Slope = atan(SlopeBot/100))
##############AW Plots
p1=ggplot() +
  geom_line(data=BotSlope,aes(x=date, y = AW,colour="BotSlope")) +
  geom_line(data=MidSlope,aes(x=date, y = AW,colour="MidSlope")) +
  geom_line(data=TopSlope,aes(x=date, y = AW,colour="TopSlope")) +
  labs(x = 'Date', y = 'AW (mm)')+
  scale_colour_manual("", 
                      breaks = c("BotSlope", "MidSlope", "TopSlope"),
                      values = c("black", "blue","red"))+
  theme(text = element_text(size = 15))+
  ggtitle("(a)")
##############Excess
p2=ggplot() +
  geom_line(data=BotSlope,aes(x=date, y = Excess,colour="BotSlope")) +
  geom_line(data=MidSlope,aes(x=date, y = Excess,colour="MidSlope")) +
  geom_line(data=TopSlope,aes(x=date, y = Excess,colour="TopSlope")) +
  labs(x = 'Date', y = 'Excess (mm)')+
  scale_colour_manual("", 
                      breaks = c("BotSlope", "MidSlope", "TopSlope"),
                      values = c("black", "blue","red"))+
  theme(text = element_text(size = 15))+
  ggtitle("(b)")
p1 + p2 + plot_layout(ncol = 1, widths = c(1, 1))

# These 3 function calls are what you will vary for the Lab 04 homework 

# Model Performance 
plot(BotSlope$date,BotSlope$Qpred,type="l")
NSeff(BotSlope$Qmm,BotSlope$Qpred)

# Publishing - for when you want to create your own CRAN Package and become
# CRAN FAMOUS!
package.skeleton("BSEHydroModels",list=c("soil_wetting_above_capacity",
                                         "soilwetting","soildrying","TMWBmodel","NSE"))

url="https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHD/HU8/HighResolution/Shape/NHD_H_03010101_HU8_Shape.zip"
download.file(url,"NHD_H_03010101_HU8_Shape.zip")
unzip("NHD_H_03010101_HU8_Shape.zip",exdir="03010101")
streams=readOGR("03010101/Shape/NHDFlowline.dbf")
mystream=subset(streams,gnis_id=="01478950")
mybbox=c(mystream@bbox)
mysoil = mapunit_geom_by_ll_bbox(mybbox)
mukey_statement = format_SQL_in_statement(unique(mysoil$mukey))
q_mu2co = paste("SELECT mukey,cokey FROM component WHERE mukey IN ", mukey_statement, sep="")
mu2co = SDA_query(q_mu2co)
cokey_statement = format_SQL_in_statement(unique(mu2co$cokey))
# THE ONLY DIFFERENCE BETWEEN LAST WEEKS LAB AND THE HW3 IS:
# We modify these following lines from our Lab4 exercise to access 
# the "resdept" from the "corestrictions" table as described in the PDF.
# So, from last week we modify this line:
q_co2ch = paste("SELECT cokey,ksat_r,awc_r,hzdepb_r  FROM chorizon WHERE cokey IN ", cokey_statement, sep="")
# To HW3 needs (here we are using "cr" for corestrictions where before 
# we used "ch" for "chorizon"... naming should reflect the data (though 
# many would use the full name so they can remember when read their code)
q_co2cr = paste("SELECT cokey,resdept_r FROM corestrictions WHERE cokey IN ", cokey_statement, sep="")
co2cr = SDA_query(q_co2cr)
View(co2cr)

# Last, bring them back together, and aggregate based on max values
# of ksat_r,awc_r, and hzdepb_r
mu2cr=merge(mu2co,co2cr)

mu2crmax=aggregate(mu2cr,list(mu2cr$mukey),max)
depth2restrict<-merge(mysoil,mu2crmax,by="mukey")
# And then a simple plot will do...
plot(mystream,col="red",lwd=4)
plot(depth2restrict,col=topo.colors(10),add=T)
plot(mystream,col="red",lwd=4,add=T)
# Though you want to explore nicer plots
pols1 <- list("sp.lines", as(mystream, 'SpatialLines'), col = "red", lwd = 4)
spplot(depth2restrict, sp.layout=list(pols1), zcol="resdept_r", xlab="Longitude", ylab="Latitude", main="Depth to Restrictive Layer", colorkey=T, col.regions=colorRampPalette(c("grey87","royalblue4"))(100), checkEmptyRC=T, add=T,xlim=mystream@bbox[1,],ylim=mystream@bbox[2,])
# or other colormap
spplot(depth2restrict, sp.layout=list(pols1), zcol="resdept_r", xlab="Longitude", ylab="Latitude", main="Depth to Restrictive Layer", colorkey=T, col.regions=topo.colors(100), checkEmptyRC=T, add=T,xlim=mystream@bbox[1,],ylim=mystream@bbox[2,])

####CN
#CNModel
#
CNmodel<-function(CNmodeldf, CNavg = 75,IaFrac = 0.05,fnc_slope=0, 
                  fnc_aspect=0,func_DAWC=.3,func_z=1000,fnc_fcres=.3) {
  CNmodeldf=modeldata
  CNavg=75;IaFrac=0.05; fnc_slope=0; fnc_aspect=0;func_DAWC=.3;
  func_z=1000;fnc_fcres=.3
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
#
# Like before, initializing the 3 hillslope classes
#
TopSlopeCN=modeldata
MidSlopeCN=modeldata
BotSlopeCN=modeldata
# Call the new CNmodel() function with Top,Mid,BotSlope HRU objects,
# passing the Qpred into the lower HRUs HillslopeAboveExcess (as area scaled flow)
TopSlopeCN=CNmodel(TopSlopeCN, CNavg = 60)
TopSlopeCN = CNmodel(CNmodeldf = TopSlopeCN, CNavg = 60,fnc_slope=0,
                     fnc_aspect=0,func_DAWC=.3,
                     func_z=500,fnc_fcres=.3)
MidSlope$P=TopSlope$Excess+MidSlope$P
# Higher slope, medium ksat, fcres=0.5 
MidSlopeCN = CNmodel(CNmodeldf = MidSlopeCN, CNavg = 60,fnc_slope=0, 
                     fnc_aspect=0,func_DAWC=.3,
                     func_z=750,fnc_fcres=.5)
# Low Slope and lowest ksat, $fcres=0.2
BotSlope$P=MidSlope$Excess+BotSlope$P
BotSlopeCN = CNmodel(CNmodeldf = BotSlopeCN, CNavg = 60,fnc_slope=0, 
                     fnc_aspect=0,func_DAWC=.3,
                     func_z=1000,fnc_fcres=.2)


#AW Plots HW1
plot(BotSlopeCN$date,BotSlopeCN$AW,type="l",col=1,xlab="Date",ylab="AW (mm)")
lines(MidSlopeCN$date,MidSlopeCN$AW,type="l",col=2)
lines(TopSlopeCN$date,TopSlopeCN$AW,type="l",col=3)
# Excess Plots HW1
plot(BotSlopeCN$date,BotSlopeCN$Excess,type="l",col=1,xlab="Date",ylab="Excess (mm)")
lines(MidSlopeCN$date,MidSlopeCN$Excess,type="l",col=2)
lines(TopSlopeCN$date,TopSlopeCN$Excess,type="l",col=3)

# PET and ET HW2
plot(BotSlopeCN$date,BotSlopeCN$PET,type="l",col=1,xlab="Date",ylab="(P)ET (mm)")
lines(BotSlopeCN$date,BotSlopeCN$ET,type="l",col=2)
lines(MidSlopeCN$date,MidSlopeCN$ET,type="l",col=3)
lines(TopSlopeCN$date,TopSlopeCN$ET,type="l",col=4)
# or as cumulative summations
plot(TopSlopeCN$date,cumsum(BotSlopeCN$PET),type="l",
     xlab="Date",ylab="(P)ET")
lines(TopSlopeCN$date,cumsum(TopSlopeCN$ET),col="red")
lines(MidSlopeCN$date,cumsum(MidSlopeCN$ET),col="green")
lines(BotSlopeCN$date,cumsum(BotSlopeCN$ET),col="blue")


# Cumulative Summary of QPred is very informative
plot(BotSlopeCN$date,cumsum(BotSlopeCN$Qpred),type="l",
     xlab="Date",ylab="Flow Q Cumulative Summary (mm)")
lines(MidSlopeCN$date,cumsum(MidSlopeCN$Qpred),col="red")
lines(TopSlopeCN$date,cumsum(TopSlopeCN$Qpred),col="green")

# Model Performance 
plot(BotSlopeCN$date,BotSlopeCN$Qpred,type="l")
NSeff(BotSlopeCN$Qmm,BotSlopeCN$Qpred)

# finish building all the hillslope HRUs….
