# Get rid of plots and Views from Lab08
#
dir.create("~/Week09Th")
setwd("~/Week09Th")
url="https://raw.githubusercontent.com/vtdrfuka/BSE5304/main/Lab08Sol.R"
download.file(url,"Lab08Sol.R")
file.edit("Lab08Sol.R")   # Run to get the TIC0[1:5] data objects
#
TIC01$tici=1
TIC02$tici=2
TIC03$tici=3
TIC04$tici=4
TIC05$tici=5
TIC01$area=myflowgage$area/nTIclass
TIC02$area=myflowgage$area/nTIclass
TIC03$area=myflowgage$area/nTIclass
TIC04$area=myflowgage$area/nTIclass
TIC05$area=myflowgage$area/nTIclass


P_Model_Top=function(TIC,tau=9.3,fertrate=5.4*10^(-4),TempBias=0,
                     kF=.015,dt=1,nTIclass=5){
#TIC=TIC05;tau=9.3;fertrate=5.4*10^(-4);TempBias=0;kF=.015;dt=1;nTIclass=5

# 
# Build a dataframe for your DP Load model in TI Class 05,(you will need to do 
#this 5 times for TIC 1-4 where you will need the date, Qpred, and Average 
# Temperature
DPTI=data.frame(date=TIC$date,
                  Rt=TIC$Qpred/1000,        #m^3/m^2/day from mm/day
                  Tavg=TempBias+TIC$AvgTemp)
# Take it to volumetric water content rather # than volumetric soil content, 
# should be in m of water
#tau=9.3  # days
#dt=1     # days time step
#kF=.015  # Table 2
# Initialize MF and DF
DPTI$MF=0
DPTI$DF=0
# Spread your P Fertilizer on ~May 1, ~August 1, and ~October 1
DPTI$MF[(format(DPTI$date,"%j") %in% c(121,213,274))]=fertrate
# Remember what we say about attaching! 
attach(DPTI)
#
# Loop to solve MF and DF
for (i in 2:length(date)){
  if(MF[i]<=MF[i-1]){
    MF[i]=MF[i-1]*exp(-dt/tau)-DF[i-1]
  }
  DF[i]=MF[i]*(kF*MF[i]*Rt[i]/(1+kF*MF[i]*Rt[i]))
}
DPTI$MF=MF
DPTI$DF=DF
detach(DPTI)
rm(list=c("MF","DF")) # Clean up the environment 
#dev.off() #reset graphics device
#plot(DPTI$date,DPTI$MF)
#plot(DPTI$date,DPTI$DF)

# Calculate your Export Coef from Easton et al 2007 Figure 2 using TI Class 5
# and note that the bold gives the TI Class. This figure gives a range of 
# 0 - 520 micrograms/litre 
# For TIC=5
muTS_TI=(((520-0)/nTIclass)*TIC$tici+0)/1000000 # use VSAsol$TIClass table to 
# micrograms/litre * 10^3 litre/m^3 / 10^9micrograms/kg
# assign TIC to calc (remember TIC05 is in location 1 in the VSAsol$TIClass 
# table
# Setting range of Soil P values (MS) using the range given on page 7 of 
# Easton et al. 2007 (3.7-18.5mg/kg), then 
# Moore 1993… assume soil density is 2000kg/m^3 of soil
MS_TI=(((18.5-3.7)/nTIclass)*TIC$tici+3.7)*2000/10^6  # range from Easton et 
#      mg/kg * 2000kg/m^3 / 10^6 mg/kg    has to be in kg/m^3
# al. 2007, pg 7. Moore suggests linear relationship
# We will take care of all of TIClass 05 now as will so 
# it makes sense when you repeat for TI Classes 1-4
# You will use muTS_TI01 thru muTS_TI04 to finish this lab
#
QS= 3.0 # A guess using the middle of the range 1-5
TR=20   # reference Temperature from Table 2.
DPTI$muS= muTS_TI*QS^((DPTI$Tavg-TR)/10)  # Eq. 5
DPTI$DS=(DPTI$muS*MS_TI*DPTI$Rt)          # Eq. 4
#plot(DPTI$date,DPTI$DS)
return(DPTI)

}

# Completed model with default parameters
DPTI01=P_Model_Top(TIC01)
DPTI02=P_Model_Top(TIC02)
DPTI03=P_Model_Top(TIC03)
DPTI04=P_Model_Top(TIC04)
DPTI05=P_Model_Top(TIC05)
# Initialize a dataframe for the total P losses, and 
# build the P Loss in Base Flow
DPLT=data.frame(date=TIC05$date,
                Rt=TIC05$Qpred/1000.0,
                Tavg=(TIC05$AvgTemp+5))
DPLT$B=min(TIC05$Qmm)*myflowgage$area*1000*1000/1000 # m^3/day
muTB=2.1*10^(-5) # Easton Table 2
QB=2.2           # Easton Table 2
TB=17            # Easton Table 2
DPLT$muB=muTB*QB^((DPLT$Tavg-TB)/10)  # Easton eq. 10
DPLT$LB=DPLT$muB*DPLT$B     # Easton eq. 9
# Sum the total Stream DP Loss (Baseflow + Fert + Soil/Plant)
DPLT$LT=DPLT$LB +
  TIC05$area*1000*1000*(DPTI05$DF + DPTI05$DS)+
  TIC04$area*1000*1000*(DPTI04$DF + DPTI04$DS)+
  TIC03$area*1000*1000*(DPTI03$DF + DPTI03$DS)+
  TIC02$area*1000*1000*(DPTI02$DF + DPTI02$DS)+
  TIC01$area*1000*1000*(DPTI01$DF + DPTI01$DS)

plot(DPLT$date,DPLT$LT, type="l",col=2)
mean(DPLT$LT)*365   # kg/year into stream 
# [1] 165.5578

# HW1 Set fertrate to 0, and remove baseflow
DPTI01=P_Model_Top(TIC01,fertrate = 0)
DPTI02=P_Model_Top(TIC02,fertrate = 0)
DPTI03=P_Model_Top(TIC03,fertrate = 0)
DPTI04=P_Model_Top(TIC04,fertrate = 0)
DPTI05=P_Model_Top(TIC05,fertrate = 0)
DPLT=data.frame(date=TIC05$date,
                Rt=TIC05$Qpred/1000.0,
                Tavg=(TIC05$AvgTemp+5))
DPLT$B=min(TIC05$Qmm)*myflowgage$area*1000*1000/1000 # m^3/day
muTB=2.1*10^(-5) # Easton Table 2
QB=2.2           # Easton Table 2
TB=17            # Easton Table 2
DPLT$muB=muTB*QB^((DPLT$Tavg-TB)/10)  # Easton eq. 10
DPLT$LB=DPLT$muB*DPLT$B     # Easton eq. 9
# Just remove DPLT$LB from the calculation
DPLT$LT=      # DPLT$LB +
  TIC05$area*1000*1000*(DPTI05$DF + DPTI05$DS)+
  TIC04$area*1000*1000*(DPTI04$DF + DPTI04$DS)+
  TIC03$area*1000*1000*(DPTI03$DF + DPTI03$DS)+
  TIC02$area*1000*1000*(DPTI02$DF + DPTI02$DS)+
  TIC01$area*1000*1000*(DPTI01$DF + DPTI01$DS)

mean(DPLT$LT)*365   # kg/year into stream if no fert added
# [1] 91.1213

(myflowgage$area*1000*1000)*(.2)*(2000) * # kg of soil in top 20cm at 2000kg/m^3
    mean(c(3.7,18.5))/1000/1000 *.1 / # 10 % of kg P / kg soil
    mean(DPLT$LT)/365 # daily mean kg/day into years
# [1] 63.24278

#HW2) add 5deg
DPTI01=P_Model_Top(TIC01,TempBias=5)
DPTI02=P_Model_Top(TIC02,TempBias=5)
DPTI03=P_Model_Top(TIC03,TempBias=5)
DPTI04=P_Model_Top(TIC04,TempBias=5)
DPTI05=P_Model_Top(TIC05,TempBias=5)
DPLT=data.frame(date=TIC05$date,
                Rt=TIC05$Qpred/1000.0,
                Tavg=(TIC05$AvgTemp+5))
DPLT$B=min(TIC05$Qmm)*myflowgage$area*1000*1000/1000 # m^3/day
muTB=2.1*10^(-5) # Easton Table 2
QB=2.2           # Easton Table 2
TB=17            # Easton Table 2
DPLT$muB=muTB*QB^((DPLT$Tavg-TB)/10)  # Easton eq. 10
DPLT$LB=DPLT$muB*DPLT$B     # Easton eq. 9
DPLT$LT=DPLT$LB +
  TIC05$area*1000*1000*(DPTI05$DF + DPTI05$DS)+
  TIC04$area*1000*1000*(DPTI04$DF + DPTI04$DS)+
  TIC03$area*1000*1000*(DPTI03$DF + DPTI03$DS)+
  TIC02$area*1000*1000*(DPTI02$DF + DPTI02$DS)+
  TIC01$area*1000*1000*(DPTI01$DF + DPTI01$DS)

plot(DPLT$date,DPLT$LT, type="l",ylim=c(1,10))
mean(DPLT$LT)*365  # 
#[1] 232.2632
# First model run with no bias gave [1] 165.5578


#HW3
DPTI01=P_Model_Top(TIC01,tau = 2)
DPTI02=P_Model_Top(TIC02,tau = 2)
DPTI03=P_Model_Top(TIC03,tau = 2)
DPTI04=P_Model_Top(TIC04,tau = 2)
DPTI05=P_Model_Top(TIC05,tau = 2)
DPLT=data.frame(date=TIC05$date,
                Rt=TIC05$Qpred/1000.0,
                Tavg=(TIC05$AvgTemp+5))
DPLT$B=min(TIC05$Qmm)*myflowgage$area*1000*1000/1000 # m^3/day
muTB=2.1*10^(-5) # Easton Table 2
QB=2.2           # Easton Table 2
TB=17            # Easton Table 2
DPLT$muB=muTB*QB^((DPLT$Tavg-TB)/10)  # Easton eq. 10
DPLT$LB=DPLT$muB*DPLT$B     # Easton eq. 9
DPLT$LT=DPLT$LB +
  TIC05$area*1000*1000*(DPTI05$DF + DPTI05$DS)+
  TIC04$area*1000*1000*(DPTI04$DF + DPTI04$DS)+
  TIC03$area*1000*1000*(DPTI03$DF + DPTI03$DS)+
  TIC02$area*1000*1000*(DPTI02$DF + DPTI02$DS)+
  TIC01$area*1000*1000*(DPTI01$DF + DPTI01$DS)
lines(DPTI05$date,DPTI05$DF, type="l",col=1)
mean(DPLT$LT)*365 # Does very little difference to overall P budget

# then skip to the 
DPTI01=P_Model_Top(TIC01,tau = 9.3)
DPTI02=P_Model_Top(TIC02,tau = 9.3)
DPTI03=P_Model_Top(TIC03,tau = 9.3)
DPTI04=P_Model_Top(TIC04,tau = 9.3)
DPTI05=P_Model_Top(TIC05,tau = 9.3)

DPLT=data.frame(date=TIC05$date,
                Rt=TIC05$Qpred/1000.0,
                Tavg=(TIC05$AvgTemp+5))
DPLT$B=min(TIC05$Qmm)*myflowgage$area*1000*1000/1000 # m^3/day
muTB=2.1*10^(-5) # Easton Table 2
QB=2.2           # Easton Table 2
TB=17            # Easton Table 2
DPLT$muB=muTB*QB^((DPLT$Tavg-TB)/10)  # Easton eq. 10
DPLT$LB=DPLT$muB*DPLT$B     # Easton eq. 9
DPLT$LT=DPLT$LB +
  TIC05$area*1000*1000*(DPTI05$DF + DPTI05$DS)+
  TIC04$area*1000*1000*(DPTI04$DF + DPTI04$DS)+
  TIC03$area*1000*1000*(DPTI03$DF + DPTI03$DS)+
  TIC02$area*1000*1000*(DPTI02$DF + DPTI02$DS)+
  TIC01$area*1000*1000*(DPTI01$DF + DPTI01$DS)
lines(DPTI05$date,DPTI05$DF, type="l",col=1)
mean(DPLT$LT)*365 # Does very little difference to overall P budget

# Basically, there is a lot of P in the soil already
(myflowgage$area*1000*1000)*(.2)*(2000) * # kg of soil in top 20cm at 2000kg/m^3
  mean(c(3.7,18.5))/1000/1000  # of kg P / kg soil

