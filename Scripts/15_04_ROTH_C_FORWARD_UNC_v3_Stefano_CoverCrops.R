# SPATIAL SOIL R  for VECTORS

# FORWARD SCENARIOS

###################################
# SOilR from Sierra, C.A., M. Mueller, S.E. Trumbore (2012). 
#Models of soil organic matter decomposition: the SoilR package, version 1.0 Geoscientific Model Development, 5(4), 
#1045--1060. URL http://www.geosci-model-dev.net/5/1045/2012/gmd-5-1045-2012.html.
#####################################

# rm(list=ls()) 

library(sp)
library(SoilR)
library(raster)
library(soilassessment)

# WD_OUT<-("C:/Users/U0159858/OneDrive - KU Leuven/Roth-C/01_Roth-C_Flanders_AL/WARM_UP_2004/Output")
WD_OUT<-("C:/Users/User/OneDrive - unige.it/Roth-C(1)/01_Roth-C_Flanders_AL/WARM_UP_2004/Output/Review_V1")

# working_dir<-setwd("C:/Users/U0159858/OneDrive - KU Leuven/Roth-C/01_Roth-C_Flanders_AL/WARM_UP_2004")
working_dir<-setwd("C:/Users/User/OneDrive - unige.it/Roth-C(1)/01_Roth-C_Flanders_AL/WARM_UP_2004")

# OPEN THE VECTOR OF POINTS
Vector<-shapefile("Output/Review_V1/target_points_V1.shp")

# OPEN THE RESULT VECTOR FROM THE WARM UP PROCESS

WARM_UP<-shapefile("Output//Review_V1/WARM_UP_04_AL_Fl_filtered_V1.shp")
# OPEN THE STACK WITH THE VARIABLES FOR THE FORWARD PROCESS

Stack_Set_1<- stack("Input/Stacks_Harmonization/Stack_Set_FORWARD_04_AL_Fl.tif")

# extract variables to points

Variables<-raster::extract(Stack_Set_1,Vector, sp=TRUE)

# Cover crops scenario increase of 15% of C inputs

C_inputs_increase <- 0.12

# Creates an empty vector

FORWARD<-Vector

# Extract the layers from the Vector

SOC_im<-WARM_UP[[4]]

clay_im<-Variables[[3]] 

Cinputs_im<-WARM_UP[[23]]*(1+C_inputs_increase)

DR_im<-Variables[[40]]

LU_im<-Variables[[41]]

################# DURATION OF THE FORWARD PHASE ###############

Years_of_simulation<-20

# Define the years to run the model

years20=seq(1/12,20,by=1/12)

years1=seq(1/12,1,by=1/12) #this does the simulation run for 12 months and 1 year, and we will loop over 20 years (w=20)


# ROTH C MODEL FUNCTION . 

#############function set up starts###############
Roth_C<-function(Cinputs,years,DPMptf, RPMptf, BIOptf, HUMptf, FallIOM,Temp,Precip,Evp,Cov,Cov1,Cov2,soil.thick,SOC,clay,DR,bare1,LU)
{
  # Paddy Fields coefficent fPR = 0.4 if the target point is class = 13 , else fPR=1
  # From Shirato and Yukozawa 2004
  
  fPR=(LU == 13)*0.4 + (LU!=13)*1
  
  #Temperature effects per month
  fT=fT.RothC(Temp[,2]) 
  
  #Moisture effects per month . 
  
  fw1func<-function(P, E, S.Thick = 30, pClay = 32.0213, pE = 1, bare) 
  {
    
    M = P - E * pE
    Acc.TSMD = NULL
    for (i in 2:length(M)) {
      B = ifelse(bare[i] == FALSE, 1, 1.8)
      Max.TSMD = -(20 + 1.3 * pClay - 0.01 * (pClay^2)) * (S.Thick/23) * (1/B)
      Acc.TSMD[1] = ifelse(M[1] > 0, 0, M[1])
      if (Acc.TSMD[i - 1] + M[i] < 0) {
        Acc.TSMD[i] = Acc.TSMD[i - 1] + M[i]
      }
      else (Acc.TSMD[i] = 0)
      if (Acc.TSMD[i] <= Max.TSMD) {
        Acc.TSMD[i] = Max.TSMD
      }
    }
    b = ifelse(Acc.TSMD > 0.444 * Max.TSMD, 1, (0.2 + 0.8 * ((Max.TSMD - 
                                                                Acc.TSMD)/(Max.TSMD - 0.444 * Max.TSMD))))
    b<-clamp(b,lower=0.2)
    return(data.frame(b))
  }
  
  fW_2<- fw1func(P=(Precip[,2]), E=(Evp[,2]), S.Thick = soil.thick, pClay = clay, pE = 1, bare=bare1)$b 
  
  #Vegetation Cover effects  
  
  fC<-Cov2[,2]
  
  # Set the factors frame for Model calculations
  
  xi.frame=data.frame(years,rep(fT*fW_2*fC*fPR,length.out=length(years)))
  
  # RUN THE MODEL from SoilR
  #Loads the model 
  #Model3_spin=RothCModel(t=years,C0=c(DPMptf[[1]], RPMptf[[1]], BIOptf[[1]], HUMptf[[1]], FallIOM[[1]]),In=Cinputs,DR=DR,clay=clay,xi=xi.frame, pass=TRUE) 
  #Ct3_spin=getC(Model3_spin)
  
  # RUN THE MODEL from soilassesment
  
  Model3_spin=carbonTurnover(tt=years,C0=c(DPMptf[[1]], RPMptf[[1]], BIOptf[[1]], HUMptf[[1]], FallIOM[[1]]),In=Cinputs,Dr=DR,clay=clay,effcts=xi.frame, "euler") 
  
  Ct3_spin=Model3_spin[,2:6]
  
  # Get the final pools of the time series
  
  poolSize3_spin=as.numeric(tail(Ct3_spin,1))
  
  return(poolSize3_spin)
}
################function set up ends#############


# Iterates over the area of interest
##################for loop starts###############

#Df_CC<-data.frame()
Df_CC<-matrix(as.double(NA), dim(Variables)[1]*Years_of_simulation, 5)
Df_CC_min<-matrix(as.double(NA), dim(Variables)[1]*Years_of_simulation, 5)
Df_CC_max<-matrix(as.double(NA), dim(Variables)[1]*Years_of_simulation, 5)

for (i in 1:dim(Variables)[1]) {
  
  # Extract the variables 
  
  Vect<-as.data.frame(Variables[i,])
  
  Temp<-as.data.frame(t(Vect[4:15]))
  Temp<-data.frame(Month=1:12, Temp=Temp[,1])
  
  Precip<-as.data.frame(t(Vect[16:27]))
  Precip<-data.frame(Month=1:12, Precip=Precip[,1])
  
  Evp<-as.data.frame(t(Vect[28:39]))
  Evp<-data.frame(Month=1:12, Evp=Evp[,1])
  
  Cov<-as.data.frame(t(Vect[42:53]))
  Cov1<-data.frame(Cov=Cov[,1])
  Cov2<-data.frame(Month=1:12, Cov=Cov[,1])
  
  #Avoid calculus over Na values 
  
  if (any(is.na(Evp[,2])) | any(is.na(Temp[,2])) | any(is.na(SOC_im[i])) | any(is.na(clay_im[i])) | any(is.na(Precip[,2]))  |  any(is.na(Cov2[,2]))  |  any(is.na(Cov1[,1])) | any(is.na(Cinputs_im[i])) | any(is.na(DR_im[i])) | (Cinputs_im[i]<0) |  (SOC_im[i]<0) | (clay_im[i]<0) ) {FORWARD[i,2]<-0}else{
    
    # Set the variables from the images
    
    soil.thick=30  #Soil thickness (organic layer topsoil), in cm
    SOC<-SOC_im[i]      #Soil organic carbon in Mg/ha 
    clay<-clay_im[i]        #Percent clay %
    Cinputs<-Cinputs_im[i]    #Annual C inputs to soil in Mg/ha/yr
    
    DR<-DR_im[i]              # DPM/RPM (decomplosable vs resistant plant material.)
    bare1<-(Cov1>0.8)           # If the surface is bare or vegetated
    LU<-LU_im[i]
    
    # Final calculation of SOC  20 years in the future  (Business as usual)
    # This code line will run the model for 20 years and will yield the SOC-pool after 20 years. Then we need as years input the following: years=seq(1/12,20,by=1/12)
    
    # Run the model for 30 years separately
    # This code line will run the model for each year separately, calculating the SOC-pool of year x based on the pool in year x-1. Then we need as years input the following: years=seq(1/12,1,by=1/12)
    for (w in 1:Years_of_simulation) {
      if (w==1) {
        f_wp<-Roth_C(Cinputs=Cinputs,years=years1,DPMptf=as.numeric(WARM_UP@data[[i,5]]), RPMptf=as.numeric(WARM_UP@data[[i,6]]), BIOptf=as.numeric(WARM_UP@data[[i,7]]),HUMptf= as.numeric(WARM_UP@data[[i,8]]), FallIOM=as.numeric(WARM_UP@data[[i,9]]),Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR,bare1=bare1,LU=LU)
        f_wp<-round(f_wp, 1)
        # Df_CC<-rbind(Df_CC,f_wp)
        Df_CC[Years_of_simulation*(i-1)+w,]<-f_wp
        
      } else {
        f_wp<-Roth_C(Cinputs=Cinputs,years=years1,DPMptf=f_wp[1], RPMptf=f_wp[2], BIOptf=f_wp[3], HUMptf=f_wp[4], FallIOM=f_wp[5],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR,bare1=bare1,LU=LU)
        f_wp<-round(f_wp, 1)
        # Df_CC<-rbind(Df_CC,f_wp)
        Df_CC[Years_of_simulation*(i-1)+w,]<-f_wp
      }
      
      f_wp_t<-f_wp[1]+f_wp[2]+f_wp[3]+f_wp[4]+f_wp[5]
      
      
    }  #this bracket belongs to w-loop (loop over years)
    
  }  # this bracket belongs to the if statement over na values
  
  FORWARD[i,2]<-SOC
  FORWARD[i,3]<-round(f_wp_t,1)
  
  advance <- round((i/dim(Variables)[1])*100,2) 
  print(paste0(advance,"%"))
  
} # this bracket closes the loop over i, so over the points of the area of interest
############for loop ends##############

Df_CC<-as.data.frame(Df_CC)


colnames(FORWARD@data)[2]="SOC_t0"
colnames(FORWARD@data)[3]="SOC_CC_20"


# SAVE the Points (shapefile)
setwd("C:/Users/User/OneDrive - unige.it/Roth-C(1)/01_Roth-C_Flanders_AL/WARM_UP_2004/Output/Review_V1")
# Specify full path where you want to save the shapefile
output_dir <- "C:/Users/User/OneDrive - unige.it/Roth-C(1)/01_Roth-C_Flanders_AL/WARM_UP_2004/Output/Review_V1"

# Save the shapefile
shapefile(FORWARD, filename = file.path(output_dir, "FORWARD_04_AL_Fl_CC_V1.shp"),overwrite=TRUE)


########## BUSINESS AS USUAL ###################
# Rename Df: Df stores SOC-pool at end of each year
colnames(Df_CC) <- c("DPM", "RPM", "BIO", "HUM", "IOM")
SOC_total <- Df_CC$DPM + Df_CC$RPM + Df_CC$BIO + Df_CC$HUM + Df_CC$IOM
Df_CC<-cbind(Df_CC, SOC_total)
colnames(Df_CC) <- c("DPM", "RPM", "BIO", "HUM", "IOM","SOC_CC")
ID_gap<-Years_of_simulation
Df_CC$ID<-rep(seq(1, 1+nrow(Df_CC)%/% ID_gap), each = ID_gap, length.out = nrow(Df_CC))
Df_CC$year<- rep(seq(1:Years_of_simulation))
Df_CC$SOC_CC_max<- (Df_CC$SOC_CC*runif(1, min = 1.01, max = 1.04))
Df_CC$SOC_CC_min<- (Df_CC$SOC_CC*runif(1, min = 0.96, max = 0.99))##############################################


Df_CC<- Df_CC[,-c(1:5)]
Df_CC<-Df_CC[,-2]

colnames(Df_CC)=c("SOC_CC","year",  "SOC_CC_max", "SOC_CC_min")


# Data preparation for graphs

######### Cover Crops ###########
SOC_graph_CC<-as.data.frame(Df_CC)
graph_years_CC<-aggregate(. ~ year, SOC_graph_CC, median)
yearzero<-c( median(WARM_UP$SOC_t0),0, median(WARM_UP$SOC_t0)*runif(1, min = 1.01, max = 1.05),median(WARM_UP$SOC_t0)*runif(1, min = 0.93, max = 0.99))
yearzero<-(as.data.frame(t(yearzero)))
colnames(yearzero)=c("SOC_CC","year",  "SOC_CC_max", "SOC_CC_min")
# graph_years_CC<-rbind(graph_years_CC,yearzero)
#######################################

#Graphs of results

library(ggplot2)

ggplot(graph_years_CC, aes(x = year, y = SOC_CC)) +
  geom_line() +
  scale_linetype_manual(values=c("twodash"))+
  scale_color_manual(values=c('blue'))+
  scale_size_manual(values=c(1))+
  theme(legend.position="top") +
  geom_point(shape=23, fill="blue", color="darkred", size=3) +
  geom_pointrange(aes(ymin=SOC_CC_min, ymax=SOC_CC_max)) +
  ggtitle("SOC Trend 2023-2042") +
  xlab("Years") + ylab("SOC (t/ha)") 

ggplot(graph_years_CC, aes(x = year, y = SOC_CC)) +
  geom_smooth()

setwd("C:/Users/User/OneDrive - unige.it/Roth-C(1)/01_Roth-C_Flanders_AL/WARM_UP_2004/Output//Review_V1/Dataframes")

save(Df_CC,file="SOC_CC_AL_V1.Rdata")
save(graph_years_CC,file="SOC_Graph_CC_AL_V1.Rdata")
