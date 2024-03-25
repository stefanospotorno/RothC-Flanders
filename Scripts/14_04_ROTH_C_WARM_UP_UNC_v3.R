# SPATIAL SOIL R  for VECTORS

#   ROTH C phase 3: WARM UP

###################################
# SOilR from Sierra, C.A., M. Mueller, S.E. Trumbore (2012). 
#Models of soil organic matter decomposition: the SoilR package, version 1.0 Geoscientific Model Development, 5(4), 
#1045--1060. URL http://www.geosci-model-dev.net/5/1045/2012/gmd-5-1045-2012.html.
#####################################

rm(list=ls()) 

library(sp)
library(SoilR)
library(raster)
library(soilassessment)

WD_FOLDER <- "C:/Users/User/..."
setwd(WD_FOLDER)

#Open empty vector

Vector<-shapefile(".../AgriculturalPointsFromR.shp")

# Vector_sub<-Vector[0:100,]

#Open Warm Up Stack

Stack_Set_warmup<- stack(".../Stack_Set_WARM_UP_04_AL_Fl.tif")

plot(Stack_Set_warmup)# Open Result from SPIN UP PROCESS. A vector with 5 columns , one for each pool

Spin_up<-shapefile(".../SPIN_UP_04_AL_Fl_V1.shp")
Spin_up<-as.data.frame(Spin_up)

# Open Precipitation , temperature, and EVapotranspiration file 20 years x 12 = 240 layers x 3

# CRU LAYERS
TEMP<-stack("Input/Climate_BEL/Temp_Stack_05-22_Fl.tif")

PREC<-stack("Input/Climate_BEL/Prec_Stack_05-22_Fl.tif")

PET<-stack("Input/Climate_BEL/PET_Stack_05-22_Fl.tif")

# Set CRU LAYERS

WD_AOI<-("C:/Users/User/.../VL_Flanders")

# Open the shapefile of the region/country
# setwd(WD_AOI)
# AOI<-readOGR("Flanders_31370.shp")
# AOI_WGS84<-readOGR("Flanders_WGS84.shp")
# 
# WD_LU<-("C:/Users/User/OneDrive - unige.it/Roth-C(1)/01_Roth-C_Flanders_AL/WARM_UP_2004/Input/LandUse")
# setwd(WD_LU)
# LU_AOI_WGS84<-raster("ESA_Land_Cover_13clases_FAO_Fl_WGS84.tif")
# LU_AOI<-projectRaster(LU_AOI_WGS84,crs="+proj=lcc +lat_0=90 +lon_0=4.36748666666667 +lat_1=51.1666672333333 +lat_2=49.8333339 +x_0=150000.013 +y_0=5400088.438 +ellps=intl +units=m +no_defs", method = 'ngb')
# 
# TEMP_AOI<-crop(TEMP,AOI_WGS84, snap="out")
# TEMP_AOI<-resample(TEMP_AOI,LU_AOI_WGS84)
# TEMP_AOI<-crop(TEMP_AOI,AOI_WGS84)
# TEMP_AOI<-mask(TEMP_AOI,AOI_WGS84)
# TEMP_AOI<-projectRaster(TEMP_AOI,crs="+proj=lcc +lat_0=90 +lon_0=4.36748666666667 +lat_1=51.1666672333333 +lat_2=49.8333339 +x_0=150000.013 +y_0=5400088.438 +ellps=intl +units=m +no_defs")
# TEMP_AOI<-stack(TEMP_AOI)
# 
# PREC_AOI<-crop(PREC,AOI_WGS84, snap="out")
# PREC_AOI<-resample(PREC_AOI,LU_AOI_WGS84)
# PREC_AOI<-crop(PREC_AOI,AOI_WGS84)
# PREC_AOI<-mask(PREC_AOI,AOI_WGS84)
# PREC_AOI<-projectRaster(PREC_AOI,crs="+proj=lcc +lat_0=90 +lon_0=4.36748666666667 +lat_1=51.1666672333333 +lat_2=49.8333339 +x_0=150000.013 +y_0=5400088.438 +ellps=intl +units=m +no_defs")
# PREC_AOI<-stack(PREC_AOI)
# 
# PET_AOI<-crop(PET,AOI_WGS84, snap="out")
# PET_AOI<-resample(PET_AOI,LU_AOI_WGS84)
# PET_AOI<-crop(PET_AOI,AOI_WGS84)
# PET_AOI<-mask(PET_AOI,AOI_WGS84)
# PET_AOI<-projectRaster(PET_AOI,crs="+proj=lcc +lat_0=90 +lon_0=4.36748666666667 +lat_1=51.1666672333333 +lat_2=49.8333339 +x_0=150000.013 +y_0=5400088.438 +ellps=intl +units=m +no_defs")
# PET_AOI<-stack(PET_AOI)

# setwd("C:/Users/User/OneDrive - unige.it/Roth-C(1)/01_Roth-C_Flanders_AL/WARM_UP_2004/Input/Stacks_Harmonization")
# writeRaster(TEMP_AOI,filename="TEMP_AOI_WARM_UP_04_AL_Fl",format='GTiff',overwrite=TRUE)
# writeRaster(PREC_AOI,filename="PREC_AOI_WARM_UP_04_AL_Fl",format='GTiff',overwrite=TRUE)
# writeRaster(PET_AOI,filename="PET_AOI_WARM_UP_04_AL_Fl",format='GTiff',overwrite=TRUE)
# # 
setwd("C:/Users/User/OneDrive - unige.it/Roth-C(1)/01_Roth-C_Flanders_AL/WARM_UP_2004/Input/Stacks_Harmonization")
TEMP_AOI<-stack("TEMP_AOI_WARM_UP_04_AL_Fl.tif")
PREC_AOI<-stack("PREC_AOI_WARM_UP_04_AL_Fl.tif")
PET_AOI<-stack("PET_AOI_WARM_UP_04_AL_Fl.tif")

#Also check that the number of layers of DR and Land use are the same as the years of climate data. 
setwd(WD_FOLDER)

#Open Mean NPP MIAMI 1981 - 2000

NPP<-raster("Input/NPP/NPP_MIAMI_MEAN_85-04_31370_Fl.tif")

NPP_MEAN_MIN<-raster("Input/NPP/NPP_MIAMI_MEAN_85-04_Fl_MIN_31370.tif")

NPP_MEAN_MAX<-raster("Input/NPP/NPP_MIAMI_MEAN_85-04_Fl_MAX_31370.tif")

#Open LU layer (year 2000).

LU_AOI<-raster("Input/LandUse/ESA_Land_Cover_00_13clases_FAO_Fl_300m.tif")

#Apply NPP coeficientes
NPP<-((LU_AOI==2 | LU_AOI==12 | LU_AOI==13)*NPP*0.53+ (LU_AOI==4)*NPP*0.88 + (LU_AOI==3 | LU_AOI==5 | LU_AOI==6 | LU_AOI==8)*NPP*0.72)
NPP_MEAN_MIN<-((LU_AOI==2 | LU_AOI==12 | LU_AOI==13)*NPP_MEAN_MIN*0.53+ (LU_AOI==4)*NPP_MEAN_MIN*0.88 + (LU_AOI==3 | LU_AOI==5 | LU_AOI==6 | LU_AOI==8)*NPP_MEAN_MIN*0.72)
NPP_MEAN_MAX<-((LU_AOI==2 | LU_AOI==12 | LU_AOI==13)*NPP_MEAN_MAX*0.53+ (LU_AOI==4)*NPP_MEAN_MAX*0.88 + (LU_AOI==3 | LU_AOI==5 | LU_AOI==6 | LU_AOI==8)*NPP_MEAN_MAX*0.72)


# Extract variables to points

Vector_points<-raster::extract(Stack_Set_warmup,Vector,sp=TRUE)
Vector_points<-raster::extract(TEMP,Vector_points,sp=TRUE)
Vector_points<-raster::extract(PREC,Vector_points,sp=TRUE)
Vector_points<-raster::extract(PET,Vector_points,sp=TRUE)
Vector_points<-raster::extract(NPP,Vector_points,sp=TRUE)
Vector_points<-raster::extract(NPP_MEAN_MIN,Vector_points,sp=TRUE)
Vector_points<-raster::extract(NPP_MEAN_MAX,Vector_points,sp=TRUE)

WARM_UP<-Vector


#WARM_UP<-readOGR("WARM_UP_Fl.shp")

# Warm Up number of years simulation: 22

yearsSimulation<-dim(TEMP)[3]/12

clim_layers<-yearsSimulation*12

nppBand<-nlayers(Stack_Set_warmup)+clim_layers*3+2

firstClimLayer<-nlayers(Stack_Set_warmup)+2

nppBand_min<-nppBand+1

nppBand_max<-nppBand+2

nDR_beg<-(16+yearsSimulation)
nDR_end<-nDR_beg+(yearsSimulation-1)

# Extract the layers from the Vector

SOC_im<-Vector_points[[2]] 

clay_im<-Vector_points[[3]] 

LU_im<-Vector_points[[16]]

NPP_im<-Vector_points[[nppBand]]

NPP_im_MIN<-Vector_points[[nppBand_min]]

NPP_im_MAX<-Vector_points[[nppBand_max]]

# Define Years

years=seq(1/12,1,by=1/12)


# ROTH C MODEL FUNCTION . 

###########function set up starts################
Roth_C<-function(Cinputs,years,DPMptf, RPMptf, BIOptf, HUMptf, FallIOM,Temp,Precip,Evp,Cov,Cov1,Cov2,soil.thick,SOC,clay,DR,bare1,LU)

{

# Paddy fields coefficent fPR = 0.4 if the target point is class = 13 , else fPR=1
# From Shirato and Yukozawa 2004

fPR=(LU == 13)*0.4 + (LU!=13)*1

#Temperature effects per month
fT=fT.RothC(Temp[,2]) 

#Moisture effects per month . If evapotranspiration is used pE=1

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

#Vegetation Cover effects  C1: No till Agriculture, C2: Conventional Agriculture, C3: Grasslands and Forests, C4 bareland and Urban

fC<-Cov2[,2]

# Set the factors frame for Model calculations

xi.frame=data.frame(years,rep(fT*fW_2*fC*fPR,length.out=length(years)))

# RUN THE MODEL from SoilR
#Loads the model if pass=TRUE calculates the model even if invalid. 
#Model3_spin=RothCModel(t=years,C0=c(DPMptf, RPMptf, BIOptf, HUMptf, FallIOM),In=Cinputs,DR=DR,clay=clay,xi=xi.frame, pass=TRUE) 
#Calculates stocks for each pool per month
#Ct3_spin=getC(Model3_spin)

# RUN THE MODEL from soilassesment

Model3_spin=carbonTurnover(tt=years,C0=c(DPMptf, RPMptf, BIOptf, HUMptf, FallIOM),In=Cinputs,Dr=DR,clay=clay,effcts=xi.frame, "euler") 

Ct3_spin=Model3_spin[,2:6]

# Get the final pools of the time series

poolSize3_spin=as.numeric(tail(Ct3_spin,1))

return(poolSize3_spin)

}
##############funtion set up ends##########

# Iterates over the area of interest and over 22 years 

Cinputs<-c()
Cinputs_min<-c()
Cinputs_max<-c()
NPP_M_MIN<-c()
NPP_M_MAX<-c()
NPP_M<-c()

############for loop starts################
for (i in 1:(length(Vector_points))) {

gt<-firstClimLayer
gp<-gt+clim_layers
gevp<-gp+clim_layers

for (w in 1:(dim(TEMP)[3]/12)) {

# print(c("year:",w))
# Extract the variables 

Vect<-as.data.frame(Vector_points[i,])

Temp<-as.data.frame(t(Vect[gt:(gt+11)]))
Temp<-data.frame(Month=1:12, Temp=Temp[,1])

Precip<-as.data.frame(t(Vect[gp:(gp+11)]))
Precip<-data.frame(Month=1:12, Precip=Precip[,1])

Evp<-as.data.frame(t(Vect[gevp:(gevp+11)]))
Evp<-data.frame(Month=1:12, Evp=Evp[,1])
	
Cov<-as.data.frame(t(Vect[4:15]))
Cov1<-data.frame(Cov=Cov[,1])
Cov2<-data.frame(Month=1:12, Cov=Cov[,1])

DR_im<-as.data.frame(t(Vect[nDR_beg:nDR_end])) # DR one per year according to LU
DR_im<-data.frame(DR_im=DR_im[,1])

gt<-gt+12
gp<-gp+12
gevp<-gevp+12

#Avoid calculus over Na values 

if (any(is.na(Evp[,2])) | any(is.na(Temp[,2])) | any(is.na(SOC_im[i])) | any(is.na(clay_im[i])) | any(is.na(Spin_up[i,3]))  | any(is.na(NPP_im[i])) | any(is.na(Precip[,2]))  |  any(is.na(Cov2[,2]))  |  any(is.na(Cov1[,1]))  | any(is.na(DR_im[,1]))    |  (SOC_im[i]<0) | (clay_im[i]<0) | (Spin_up[i,3]<=0) ) {WARM_UP[i,2]<-0}else{

# Get the variables from the vector

soil.thick=30  #Soil thickness (organic layer topsoil), in cm
SOC<-SOC_im[i]      #Soil organic carbon in Mg/ha 
clay<-clay_im[i]        #Percent clay %

DR<-DR_im[w,1]              # DPM/RPM (decomplosable vs resistant plant material.)
bare1<-(Cov1>0.8)           # If the surface is bare or vegetated
NPP_81_00<-NPP_im[i]
NPP_81_00_MIN<-NPP_im_MIN[i]
NPP_81_00_MAX<-NPP_im_MAX[i]

# PHASE 2  : WARM UP .  years (w)

# Cinputs 
T<-mean(Temp[,2])
P<-sum(Precip[,2])
NPP_M[w]<-NPPmodel(P,T,"miami")*(1/100)*0.5
NPP_M[w]<-(LU_im[i]==2 | LU_im[i]==12 | LU_im[i]==13)*NPP_M[w]*0.53+ (LU_im[i]==4)*NPP_M[w]*0.88 + (LU_im[i]==3 | LU_im[i]==5 | LU_im[i]==6 | LU_im[i]==8)*NPP_M[w]*0.72

if (w==1) {Cinputs[w]<-(Spin_up[i,3]/NPP_81_00)*NPP_M[w]} else {Cinputs[w]<-(Cinputs[[w-1]]/ NPP_M[w-1]) * NPP_M[w]} 

# Cinputs MIN

Tmin<-mean(Temp[,2]*1.02)
Pmin<-sum(Precip[,2]*0.95)
NPP_M_MIN[w]<-NPPmodel(Pmin,Tmin,"miami")*(1/100)*0.5
NPP_M_MIN[w]<-(LU_im[i]==2 | LU_im[i]==12 | LU_im[i]==13)*NPP_M_MIN[w]*0.53+ (LU_im[i]==4)*NPP_M_MIN[w]*0.88 + (LU_im[i]==3 | LU_im[i]==5 | LU_im[i]==6 | LU_im[i]==8)*NPP_M_MIN[w]*0.72

if (w==1) {Cinputs_min[w]<-(Spin_up[i,10]/NPP_81_00)*NPP_M_MIN[w]} else {Cinputs_min[w]<-(Cinputs_min[[w-1]]/ NPP_M_MIN[w-1]) * NPP_M_MIN[w]} 

# Cinputs MAX

Tmax<-mean(Temp[,2]*0.98)
Pmax<-sum(Precip[,2]*1.05)
NPP_M_MAX[w]<-NPPmodel(Pmax,Tmax,"miami")*(1/100)*0.5
NPP_M_MAX[w]<-(LU_im[i]==2 | LU_im[i]==12 | LU_im[i]==13)*NPP_M_MAX[w]*0.53+ (LU_im[i]==4)*NPP_M_MAX[w]*0.88 + (LU_im[i]==3 | LU_im[i]==5 | LU_im[i]==6 | LU_im[i]==8)*NPP_M_MAX[w]*0.72

if (w==1) {Cinputs_max[w]<-(Spin_up[i,11]/NPP_81_00)*NPP_M_MAX[w]} else {Cinputs_max[w]<-(Cinputs_max[[w-1]]/ NPP_M_MAX[w-1]) * NPP_M_MAX[w]} 

# Run the model for 2001-2022

if (w==1) {
f_wp<-Roth_C(Cinputs=Cinputs[1],years=years,DPMptf=Spin_up[i,5], RPMptf=Spin_up[i,6], BIOptf=Spin_up[i,7], HUMptf=Spin_up[i,8], FallIOM=Spin_up[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR,bare1=bare1,LU=LU_im[i])
} else {
f_wp<-Roth_C(Cinputs=Cinputs[w],years=years,DPMptf=f_wp[1], RPMptf=f_wp[2], BIOptf=f_wp[3], HUMptf=f_wp[4], FallIOM=f_wp[5],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR,bare1=bare1,LU=LU_im[i])
}

f_wp_t<-f_wp[1]+f_wp[2]+f_wp[3]+f_wp[4]+f_wp[5]

# Run the model for minimum values

if (w==1) {
f_wp_min<-Roth_C(Cinputs=Cinputs_min[1],years=years,DPMptf=Spin_up[i,13], RPMptf=Spin_up[i,14], BIOptf=Spin_up[i,15], HUMptf=Spin_up[i,16], FallIOM=Spin_up[i,17],Temp=Temp*1.02,Precip=Precip*0.95,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*0.8,clay=clay*0.9,DR=DR,bare1=bare1,LU=LU_im[i])
} else {
f_wp_min<-Roth_C(Cinputs=Cinputs_min[w],years=years,DPMptf=f_wp_min[1], RPMptf=f_wp_min[2], BIOptf=f_wp_min[3], HUMptf=f_wp_min[4], FallIOM=f_wp_min[5],Temp=Temp*1.02,Precip=Precip*0.95,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*0.8,clay=clay*0.9,DR=DR,bare1=bare1,LU=LU_im[i])
}

f_wp_t_min<-f_wp_min[1]+f_wp_min[2]+f_wp_min[3]+f_wp_min[4]+f_wp_min[5]

# Run the model for maximum values

if (w==1) {
f_wp_max<-Roth_C(Cinputs=Cinputs_max[1],years=years,DPMptf=Spin_up[i,19], RPMptf=Spin_up[i,20], BIOptf=Spin_up[i,21], HUMptf=Spin_up[i,22], FallIOM=Spin_up[i,23],Temp=Temp*0.98,Precip=Precip*1.05,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*1.2,clay=clay*1.1,DR=DR,bare1=bare1,LU=LU_im[i])
} else {
f_wp_max<-Roth_C(Cinputs=Cinputs_max[w],years=years,DPMptf=f_wp_max[1], RPMptf=f_wp_max[2], BIOptf=f_wp_max[3], HUMptf=f_wp_max[4], FallIOM=f_wp_max[5],Temp=Temp*0.98,Precip=Precip*1.05,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*1.2,clay=clay*1.1,DR=DR,bare1=bare1,LU=LU_im[i])
}

f_wp_t_max<-f_wp_max[1]+f_wp_max[2]+f_wp_max[3]+f_wp_max[4]+f_wp_max[5]

# print(w)
#print(c(i,SOC,Spin_up[i,3],NPP_81_00,Cinputs[w],f_wp_t))
# print(c(NPP_M[w],Cinputs[w]))
}
}
if (is.na(mean(Cinputs))){ CinputFORWARD<-NA} else { 

CinputFORWARD<-mean(Cinputs)

CinputFORWARD_min<-mean(Cinputs_min)

CinputFORWARD_max<-mean(Cinputs_max)

WARM_UP[i,2]<-SOC
WARM_UP[i,3]<-Cinputs[18]
WARM_UP[i,4]<-f_wp_t
WARM_UP[i,5]<-f_wp[1]
WARM_UP[i,6]<-f_wp[2]
WARM_UP[i,7]<-f_wp[3]
WARM_UP[i,8]<-f_wp[4]
WARM_UP[i,9]<-f_wp[5]
WARM_UP[i,10]<-CinputFORWARD
WARM_UP[i,11]<-f_wp_t_min
WARM_UP[i,12]<-f_wp_min[1]
WARM_UP[i,13]<-f_wp_min[2]
WARM_UP[i,14]<-f_wp_min[3]
WARM_UP[i,15]<-f_wp_min[4]
WARM_UP[i,16]<-f_wp_min[5]
WARM_UP[i,17]<-f_wp_t_max
WARM_UP[i,18]<-f_wp_max[1]
WARM_UP[i,19]<-f_wp_max[2]
WARM_UP[i,20]<-f_wp_max[3]
WARM_UP[i,21]<-f_wp_max[4]
WARM_UP[i,22]<-f_wp_max[5]
WARM_UP[i,23]<-CinputFORWARD_min
WARM_UP[i,24]<-CinputFORWARD_max

Cinputs<-c()
Cinputs_min<-c()
Cinputs_max<-c()
}
print(i)
}

################for loop ends#############

colnames(WARM_UP@data)[2]="SOC_FAO"
colnames(WARM_UP@data)[3]="Cin_t0"
colnames(WARM_UP@data)[4]="SOC_t0"
colnames(WARM_UP@data)[5]="DPM_w_up"
colnames(WARM_UP@data)[6]="RPM_w_up"
colnames(WARM_UP@data)[7]="BIO_w_up"
colnames(WARM_UP@data)[8]="HUM_w_up"
colnames(WARM_UP@data)[9]="IOM_w_up"
colnames(WARM_UP@data)[10]="Cin_mean"
colnames(WARM_UP@data)[11]="SOC_t0min"
colnames(WARM_UP@data)[12]="DPM_w_min"
colnames(WARM_UP@data)[13]="RPM_w_min"
colnames(WARM_UP@data)[14]="BIO_w_min"
colnames(WARM_UP@data)[15]="HUM_w_min"
colnames(WARM_UP@data)[16]="IOM_w_min"
colnames(WARM_UP@data)[17]="SOC_t0max"
colnames(WARM_UP@data)[18]="DPM_w_max"
colnames(WARM_UP@data)[19]="RPM_w_max"
colnames(WARM_UP@data)[20]="BIO_w_max"
colnames(WARM_UP@data)[21]="HUM_w_max"
colnames(WARM_UP@data)[22]="IOM_w_max"
colnames(WARM_UP@data)[23]="Cin_min"
colnames(WARM_UP@data)[24]="Cin_max"

setwd("C:/Users/User/.../Output")

# Specify full path where you want to save the shapefile
output_dir <- "C:/Users/User/OneDrive - unige.it/Roth-C(1)/01_Roth-C_Flanders_AL/WARM_UP_2004/Output/Review_V1"

# Save the shapefile
shapefile(WARM_UP, filename = file.path(output_dir, "WARM_UP_04_AL_Fl_V1.shp"),overwrite=TRUE)

WARM_UP <- shapefile("WARM_UP_04_AL_Fl_V1.shp")

WARM_UP_filtered <- WARM_UP[!is.na(WARM_UP@data$Cin_t0),]
WARM_UP_filtered <- WARM_UP_filtered[WARM_UP_filtered@data$Cin_t0 != Inf,]

target_points <- WARM_UP_filtered
target_points <- target_points[,-c(2:24)]

shapefile(WARM_UP_filtered, filename = file.path(output_dir, "WARM_UP_04_AL_Fl_filtered_V1.shp"),overwrite=TRUE)
shapefile(target_points, filename = file.path(output_dir, "target_points_V1.shp"),overwrite=TRUE)
