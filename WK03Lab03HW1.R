#if (!require("pacman")) install.packages("pacman")
pacman::p_load(rnoaa,EcoHydRology,lattice)
library(patchwork)
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

#HOmework 1
?SnowMelt
attach(TMWB)
SNO_Energy_0=SnowMelt(date, P, MaxTemp-3, MinTemp-3, myflowgage$declat, 
                      slope = 0,
                      aspect = 0, tempHt = 1, windHt = 2, groundAlbedo = 0.25,
                      SurfEmissiv = 0.95, windSp = 2, forest = 0, startingSnowDepth_m = 0,
                      startingSnowDensity_kg_m3=450)

SNO_Energy_10N=SnowMelt(date, P, MaxTemp-3, MinTemp-3, myflowgage$declat, 
                        slope = atan(10/100),
                        aspect = 0, tempHt = 1, windHt = 2, groundAlbedo = 0.25,
                        SurfEmissiv = 0.95, windSp = 2, forest = 0, startingSnowDepth_m = 0,
                        startingSnowDensity_kg_m3=450)

SNO_Energy_10S=SnowMelt(date, P, MaxTemp-3, MinTemp-3, myflowgage$declat, 
                        slope = atan(10/100),
                        aspect = 180*pi/180, tempHt = 1, windHt = 2, groundAlbedo = 0.25,
                        SurfEmissiv = 0.95, windSp = 2, forest = 0, startingSnowDepth_m = 0,
                        startingSnowDensity_kg_m3=450)

SNO_Energy_45NW=SnowMelt(date, P, MaxTemp-3, MinTemp-3, myflowgage$declat, 
                         slope = atan(45/100),
                         aspect = 315*pi/180, tempHt = 1, windHt = 2, groundAlbedo = 0.25,
                         SurfEmissiv = 0.95, windSp = 2, forest = 0, startingSnowDepth_m = 0,
                         startingSnowDensity_kg_m3=450)

SNO_Energy_45SW=SnowMelt(date, P, MaxTemp-3, MinTemp-3, myflowgage$declat, 
                         slope = atan(45/100),
                         aspect = 225*pi/180, tempHt = 1, windHt = 2, groundAlbedo = 0.25,
                         SurfEmissiv = 0.95, windSp = 2, forest = 0, startingSnowDepth_m = 0,
                         startingSnowDensity_kg_m3=450)
detach(TMWB)

#let's save only whats required
Date=SNO_Energy_0$Date
SWE_mm_0=SNO_Energy_0$SnowWaterEq_mm
SWE_mm_10N=SNO_Energy_10N$SnowWaterEq_mm
SWE_mm_10S=SNO_Energy_10S$SnowWaterEq_mm
SWE_mm_45NW=SNO_Energy_45NW$SnowWaterEq_mm
SWE_mm_45SW=SNO_Energy_45SW$SnowWaterEq_mm

Slope_Aspect <- c("flat slope","10% North","10% South",
                 "45% Northwest","45% Southwest")
mean_SWE <- c(mean(SWE_mm_0),mean(SWE_mm_10N),
              mean(SWE_mm_10S),mean(SWE_mm_45NW),mean(SWE_mm_45SW))

dataframe <- data.frame(Slope_Aspect,mean_SWE)

write.csv(dataframe,"Mean_SWEs.csv") 
#Let's create a dataframe consisting these values
SWEmm <- data.frame(Date,SWE_mm_0,SWE_mm_10N,SWE_mm_10S,SWE_mm_45NW,SWE_mm_45SW)

p1<- ggplot() +
  # Plot your discharge data
  geom_line(data=SWEmm,aes(x=Date, y = SWE_mm_0), color="magenta",size=0.5)+
  xlab("") +
  ylab("")+
  ggtitle("Flat slope")+
  theme(
    plot.title = element_text(size=12, face="bold"),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size=11, face="bold"),
    axis.line.x = element_line(color="black",size=0.3),
    axis.line.y = element_line(color="black",size=0.3),
    panel.border = element_rect(color="black",fill=NA,size=0.3)
  )

p2<- ggplot() +
  geom_line(data=SWEmm,aes(x=Date, y = SWE_mm_10N), color="magenta",size=0.5)+
  xlab("")+
  ylab("")+
  ggtitle("10% slope North facing")+
  theme(
    plot.title = element_text(size=12, face="bold"),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size=11, face="bold"),
    axis.line.x = element_line(color="black",size=0.3),
    axis.line.y = element_line(color="black",size=0.3),
    panel.border = element_rect(color="black",fill=NA,size=0.3)
  )
p3<- ggplot() +
  geom_line(data=SWEmm,aes(x=Date, y = SWE_mm_10S), color="magenta",size=0.5)+
  xlab("")+
  ylab("Snow water eq, mm")+
  ggtitle("10% slope South facing")+
  theme(
    plot.title = element_text(size=12, face="bold"),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size=11, color="black",face="bold"),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size=12, face="bold"),
    axis.line.x = element_line(color="black",size=0.3),
    axis.line.y = element_line(color="black",size=0.3),
    panel.border = element_rect(color="black",fill=NA,size=0.3)
  )

p4<- ggplot() +
  geom_line(data=SWEmm,aes(x=Date, y = SWE_mm_45NW),color="magenta", size=0.5)+
  xlab("")+
  ylab("")+
  ggtitle("45% slope Northwest facing")+
  theme(
    plot.title = element_text(size=12, face="bold"),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size=11, face="bold"),
    axis.line.x = element_line(color="black",size=0.3),
    axis.line.y = element_line(color="black",size=0.3),
    panel.border = element_rect(color="black",fill=NA,size=0.3)
  )


p5<- ggplot() +
  geom_line(data=SWEmm,aes(x=Date, y = SWE_mm_45SW),color="magenta", size=0.5)+
  xlab("Date")+
  ylab("")+
  ggtitle("45% slope Southwest facing")+
  theme(
    plot.title = element_text(size=12, face="bold"),
    axis.title.x = element_text(size=12,color="black",face="bold"),
    axis.text.x = element_text(size=12,face="bold"),
    axis.text.y = element_text(size=11, face="bold"),
    axis.line.x = element_line(color="black",size=0.3),
    axis.line.y = element_line(color="black",size=0.3),
    panel.border = element_rect(color="black",fill=NA,size=0.3)
  )

p1 + p2 + p3 + p4 + p5 + plot_layout(nrow = 5) 

