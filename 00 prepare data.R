#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prepare data for L.boscii spatio-temporal INLA model #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Modified 04/06/2022 #
#~~~~~~~~~~~~~~~~~~~~~~~
# Francisco Izquierdo  #
#~~~~~~~~~~~~~~~~~~~~~~~

## Press Ctrl + Shift + O to see the document outline

## In this script we arrange the different survey L.boscii files
## Once we have all in the same dataset we can run the model

# start here! ------------------------------------------------------------------

rm(list=ls()) ## Clean environment

# read data -------------------------------------------------------------------

dat_a<-openxlsx::read.xlsx("./input/datasets/boscii_datos_modelo_cococha.xlsx")
head(dat_a)
summary(dat_a)
## 1993-2020 
## Data without 0's, hauls from 1 to 139
## abundance variable = numero
## unique ID column for haul and year = ID, composed by year and haul

dat_b<-openxlsx::read.xlsx("./input/datasets/LancesN93_N20.xlsx")
head(dat_b)
summary(dat_b)
## 1993-2020
## Data without 0s
## NO abundance var
## NO presence var
## NO ID column

# prepare  -------------------------------------------------------------------------

# To obtain a dataset with abundance and zeros, we need to combine dat_a + dat_b
# In addition, we first have to create dat_b ID column from lance and year
head(dat_a$ID)
head(dat_b)
dat_b$ID=paste0(dat_b$camp,"_",dat_b$lance,sep="")
head(dat_b$ID)

# join -------------------------------------------------------------------------

## We join both datasets in order to get hauls with zeros
head(dat_a)
head(dat_b)

library(conflicted)
library(dplyr)
conflict_prefer("select", "dplyr")

data<- dat_a %>% 
  full_join(dat_b, by = c("lance","lat","long","prof","temp","sali","year","ID")) %>%
  select(ID,year,lance,lat,long,prof,temp,sali,year,numero) # join all hauls

summary(data$numero) # NAs 
data$numero[is.na(data$numero)] = 0 # substitute NAs by 0s
summary(data$numero)  

# set up -------------------------------------------------------------------------

## Set up variables for INLA model
data$year<-as.factor(data$year) # year as variable ID
data$year_id<-as.numeric(data$year) # year as variable ID
data$presence<-ifelse(data$numero>0, 1, 0) # Create presence absence
head(data)

# sediment -------------------------------------------------------------------------

## extract values of variable sediment at sampling locations
## first we need to re-project from UTM to WGS84

library(raster)
sedi<-raster("./input/covariates/sediment/AllPredictionsSedimentWith2019.tif")
extent(sedi)
crs(sedi)
res(sedi)
plot(sedi) # the raster is projected in UTM

## Project Raster
wgs84 <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
sedi_wgs84 <- projectRaster(sedi, crs = wgs84)
extent(sedi_wgs84)
crs(sedi_wgs84)
res(sedi_wgs84)
plot(sedi_wgs84)

## we need to change the extent
ext<-extent(-10,-1,41.5,44.5)
sedi_wgs84 <- crop(sedi_wgs84, ext) # CROP para "cortar" un extent del mapa 
plot(sedi_wgs84)
points(data[,5:4], pch=16, cex=0.6)

## Write the RasterLayer to disk (See datatype documentation for other formats)
writeRaster(sedi_wgs84, filename="./input/covariates/sediment/sediment WGS84.tif", datatype='INT1U', overwrite=TRUE)

## extract sediment values at coordinates
data$sedi=extract(sedi_wgs84, data[,5:4])

# outliers ---------------------------------------------------------------------

# simple check to ensure that are no variables outside their normal range

summary(data$sali) # 0 value incorrect, 25.65 incorrect
data$sali[data$sali < 34] <- NA # Replace 0 with NA

summary(data$temp) # 0 value incorrect
data$temp[data$temp == 0] <- NA # Replace 0 with NA

summary(data$prof) # seems fine

summary(data$sedi) # seems fine

# save -------------------------------------------------------------------------

# with NAs
summary(data) # note that there are NAs in envars
write.table(data,"./input/datasets/boscii data INLA 1993-2020 (with NAs).txt")

# without NAs
data<-data %>%
  na.omit()
write.table(data,"./input/datasets/boscii data INLA 1993-2020.txt")
