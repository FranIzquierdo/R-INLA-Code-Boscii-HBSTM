#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# L.boscii spatio-temporal INLA model #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Modified 04/06/2022 #
#~~~~~~~~~~~~~~~~~~~~~~~
# Francisco Izquierdo  #
#~~~~~~~~~~~~~~~~~~~~~~~

## Press Ctrl + Shift + O to see the document outline

## In this script we standardize the L.boscii survey index  
## A Bayesian Hierarchical Spatio-temporal Hurdle model
## Binomial-Gamma (presence/ausence, abundance) in WGS84 lonlat
## Spatio-temporal progressive structure not shared (Paradinas et al., 2020)
##  "RW2" effect for covariates. Inla.group n=14 knots.
## Year smoothed "RW2" shared effect.
## Predicted output trend with Inla.posterior.sample


# Clean environment
rm(list=ls()) 
library(INLA)
inla.setOption(scale.model.default = TRUE) # set scale.model=TRUE, see scale tutorial

# read data 
data<-read.table(file="./input/datasets/boscii data INLA 1993-2020.txt", dec=".", header=TRUE)

# as. factors
data$year<-as.factor(data$year)

# create dir
dir_final<-paste(getwd(),"/output/03 final model", sep="")
dir.create(dir_final)

# model setup ------------------------------------------------------------------------

## mesh ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Different mesh sizes were tried for a single year
## in order to get a representative resolution in all time steps

## As we have so many years, we choose a moderated resolution mesh
length(unique(data$year))

#data=subset(data,data$year==2018)
coords<-as.matrix(data[,5:4])
bound=inla.nonconvex.hull(as.matrix(data[,5:4]), 
                          convex=-0.025, 
                          eps=0.05, 
                          resolution=40)
mesh <- inla.mesh.2d(loc=coords, 
                     boundary=bound, 
                     max.edge=c(0.05, 0.65), 
                     offset=c(0.2,0.5), 
                     cutoff=0.14, 
                     min.angle = 0.1)
mesh$n

# save jpeg
jpeg(paste(dir_final, "/mesh for all data n=370.jpeg",sep=""), quality = 300, height = 800, width = 1000)
plot(mesh,asp=1, main="")
points(coords, pch=16, cex=0.6)
dev.off()

## spde ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## P(Range < 0.5 degrees (66.6 km)) = 0.001  and P(sigma > 0.6) = 0.05

## Then, probability of spatial range being smaller  than 66.6 km is small
## and probability of standard deviation being higher than 0.6 is small

spde <- inla.spde2.pcmatern(mesh, 
                            prior.range=c(0.5, 0.05), 
                            prior.sigma=c(0.6, 0.05))

## priors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

hyper.prec <- list(theta = list(prior="pc.prec", param = c(0.5, 0.05)))# informative good for bathy

bathy.prec <- list(theta = list(prior="pc.prec", param = c(1, 0.05)))# informative good for bathy

## Amat ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

est.temp <- inla.spde.make.A(mesh, loc=matrix(c(data$long,data$lat),ncol=2),
                             group=data$year_id,
                             n.group=length(unique(data$year_id)))

## idx ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mesh.index.bin<- inla.spde.make.index("i.bin", 
                                      n.spde=spde$n.spde, 
                                      n.group=max(data$year_id))
mesh.index.con<- inla.spde.make.index("i.con", 
                                      n.spde=spde$n.spde, 
                                      n.group=max(data$year_id))

## stack ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

est.bin<-inla.stack(data=list(y=cbind(data$presence,NA)),
                    A=list(est.temp, 1),
                    effects=list(mesh.index.bin,
                                 list(bin.b0=1,
                                      prof.bin=inla.group(data$prof, n=12),
                                      year.bin=data$year_id)),
                    tag='est.bin')

est.con<-inla.stack(data=list(y=cbind(NA,ifelse(data$presence>0,data$numero,NA))),
                    A=list(est.temp, 1),
                    effects=list(mesh.index.con,
                                 list(con.b0=1,
                                      prof.con=inla.group(data$prof, n=12),
                                      year.con=data$year_id)),
                    tag='est.con')

est.stack<-inla.stack(est.bin,est.con)

## link ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

link=c(rep(1,length(est.stack$data$data[,1])/2), rep(2,length(est.stack$data$data[,1])/2))

# run --------------------------------------------------------------------------

## we add year shared rw2 effect in order to derive a shared temporal trend

## shared

final <- y ~ -1 + 
  f(bin.b0,model="linear", prec.linear=.1,mean.linear=-.5) + con.b0 +
  f(prof.bin, model="rw2", hyper=bathy.prec) + 
  f(prof.con, model="rw2", hyper=bathy.prec) +
  f(year.bin, model="rw2") + 
  f(year.con, copy="year.bin", fixed=F) +
  f(i.bin,model=spde,group = i.bin.group,control.group = list(model="ar1")) +
  f(i.con, copy="i.bin", fixed=F)

modfinal <- inla(final,
            family=c('binomial',"gamma"),
            data=inla.stack.data(est.stack), 
            control.compute=list(dic=TRUE,cpo=TRUE,waic=T,config=TRUE),# important config=TRUE
            control.predictor=list(A=inla.stack.A(est.stack), 
                                   compute=TRUE,link=link),
            verbose=TRUE,control.inla = list(strategy = "gaussian"), 
            num.threads = 1)

saveRDS(modfinal, file=paste(dir_final,"/final - prof ns - progressive s - year s .rds", sep=""))


## With verbose=TRUE we must check the internal computation values in order
## to know if the model is valid:
## First, theta values (hyperparameters i.e. spatial effect) must be close to 1
## 0.9 or 1.1 would be the limits.
## Second, the Eigenvectors of the hessian, which defines how many time has been
## INLA looking for the mode of the posterior distribution
## Similar to the convergence in Markov chains
## The values should not be higher than 3000 (e.g., 0, 300, 500, 3000)

# output plots  ---------------------------------------------------------------

shared<-readRDS(file=paste(dir_sel_st,"/r6 - prof & temp ns - progressive s.rds", sep=""))

## Plot all model effects
pdf(paste0(dir_final,"/INLA final model output plots.pdf"), # File name
    width = 8, height = 7, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk",    # Color model (cmyk is required for most publications)
    paper = "A4")          # Paper size

plot(shared)
dev.off() 

# RW2 --------------------------------------------------------------------------------------

library(ggplot2)

col_pre<-"#418FEC"
col_abu<-"#EF8C1D"

suabinm <- shared$summary.random$prof.bin$mean
suabin2 <- shared$summary.random$prof.bin$`0.025quant`
suabin9 <-shared$summary.random$prof.bin$`0.975quant`
suabinID<-shared$summary.random$prof.bin$ID
suabin<-data.frame(suabinm, suabin2,suabin9,suabinID)

m1<-ggplot(data = suabin, aes(x = suabinID, y = suabinm))+
  geom_ribbon(aes(x = suabinID, ymin = suabin2, ymax = suabin9), alpha = 0.6,fill=col_pre)+
  geom_line(aes(x = suabinID, y = suabinm), color="black", size=0.8)+
  geom_hline(yintercept = 0, colour="grey30", linetype="dashed")+
  ggtitle("")+
  xlab("Depth (m)")+ ylim(-15,6)+
  ylab("Depth effect Bernoulli")+
  theme_bw()+theme(axis.text.x= element_text(size=11), axis.text.y= element_text(size=11)) +
  scale_x_continuous(breaks=seq(0,800,100))+ theme(plot.title = element_text(hjust=0.5))

suaconm <- shared$summary.random$prof.con$mean
suacon2 <- shared$summary.random$prof.con$`0.025quant`
suacon9 <-shared$summary.random$prof.con$`0.975quant`
suaconID<-shared$summary.random$prof.con$ID
suacon<-as.data.frame(suaconm, suacon2,suacon9,suaconID)


m2<-ggplot(data = suacon, aes(x = suaconID, y = suaconm))+
  ggtitle("")+
  geom_ribbon(aes(x = suaconID, ymin = suacon2, ymax = suacon9), alpha = 0.6, fill=col_abu)+
  geom_line(aes(x = suaconID, y = suaconm), color="black", size=0.8)+
  geom_hline(yintercept = 0, colour="grey30", linetype="dashed")+
  xlab("Depth (m)")+ ylim(-15,6)+
  ylab("Depth effect Gamma")+
  theme_bw()+theme(axis.text.x= element_text(size=11), axis.text.y= element_text(size=11)) +
  scale_x_continuous(breaks=seq(0,800,100))+ theme(plot.title = element_text(hjust=0.5))

m1
ggsave(paste0(dir_final,"/Bathy RW2 effects Bernoulli.jpg"),dpi=300, height = 4, width = 6)

m2
ggsave(paste0(dir_final,"/Bathy RW2 effects Gamma.jpg"),dpi=300, height = 4, width = 6)

library(patchwork)
m1+m2+plot_layout(ncol=2)
ggsave(paste0(dir_final,"/Bathy RW2 effects.jpg"),dpi=300, height = 3, width = 9)

# spatial effect  -------------------------------------------------------------------------------

## Only the mean spatial effect is plotted
## If desired, select sd instead of mean

# Polygons ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(sp)
library(maptools)
library(rgdal)
library(rgeos)
library(raster)

shape_path <- "./input/oceandatafiles/50m_physical_separated_zip/"
ocean_shapefile <- paste(shape_path, "ne_50m_ocean/ne_50m_ocean.shp", sep="")
layer <- ogrListLayers(ocean_shapefile)
ogrInfo(ocean_shapefile, layer=layer)
ocean_poly <- readOGR(ocean_shapefile, layer=layer)# read the shape file
ocean <- ocean_poly
bbx <- readWKT("POLYGON((-9.65 44.2, -1.1 44.2, -1.1 41.83, -9.65 41.83, -9.65 44.2))") 
proj4string(bbx)<-proj4string(ocean)#Le damos CRS de pen_rec (WGS84)
ocean.cut <- gIntersection(ocean, bbx)
pen.cut<-gDifference(bbx, ocean.cut)

## Prepare map ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(ggplot2)
library(mapdata)
library(marmap)
library("ggspatial")
library("rnaturalearth")
library("rnaturalearthdata")

### Country MAP
world <- ne_countries(scale = "medium", returnclass = "sf")

### Main MAP
b = getNOAA.bathy(lon1 = -9.8, lon2 = -1, lat1 = 44.2, lat2 = 41.2, 
                  resolution = 1)# get bathymetry data
bf = fortify.bathy(b) # convert bathymetry to data frame 
reg = map_data("world2Hires")#names(reg)#table(reg$region)
reg = subset(reg, region %in% c('Spain', 'Portugal', 'France'))
# set map limits
lons = c(-9.8, -1.5)
lats = c(41.6, 44.3)

# create the breaks- and label vectors
ewbrks <- seq(-9,-1,2)
nsbrks <- seq(42,46,1)
ewlbls <- unlist(lapply(ewbrks, function(x) ifelse(x < 0, paste(-x, "ºW"), ifelse(x > 0, paste(x, "ºE"),x))))
nslbls <- unlist(lapply(nsbrks, function(x) ifelse(x < 0, paste(x, "ºS"), ifelse(x > 0, paste(x, "ºN"),x))))

## Projector MATRIX~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(splancs)
box<-bbox(ocean.cut)
Xrange<-c(box[1,1], box[1,2])
Yrange<-c(box[2,1], box[2,2])
nxy=c(120,120)

prj <- inla.mesh.projector(mesh, xlim=range(Xrange),
                           ylim=range(Yrange), dims=nxy) 

m<-mesh$n
k<-which.max(unique(data$year_id))

# bathymetry mask -------------------------------------------------------------------

# we do a mask with bathymetry around the maximum sampled depth

max(data$prof)
bathy<-raster("./input/covariates/bathymetry/bathy.asc")
mask_bathy <- mask(x = bathy , mask = (bathy[[1]]<(950)), maskvalue = 0)
plot(mask_bathy)
writeRaster(mask_bathy, filename="./input/covariates/bathymetry/mask bathy.tif", datatype='INT1U', overwrite=TRUE)

## bathymetry mask
bath<-raster("./input/covariates/bathymetry/mask bathy.tif") # 0-1200 m

# binomial ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

m.prj <- lapply(1:k, function(j) {
  r <- inla.mesh.project(prj,
                         shared$summary.ran$i.bin$mean[1:m + (j - 1) * m])
  return(r) 
})
library(fields)
library(RColorBrewer)
colors <- brewer.pal(9, "YlGnBu")
zlm <- range(unlist(m.prj), na.rm = TRUE)
dev.off()

for (j in 1:length(unique(data$year_id))) {
  
  library(raster)
  sp.mean.raster1<-raster(list(x=prj$x,
                               y = prj$y,
                               z = m.prj[[j]]))
  res<-resample(bath,sp.mean.raster1)
  aa<-mask(sp.mean.raster1,res)
  aa <- rasterToPoints(aa)
  aa<-data.frame(aa)
  colnames(aa)<-c("Longitude","Latitude","MAP")
  #Now make the map
  pr<- ggplot() +
    geom_tile(data=aa, aes(y=Latitude, x=Longitude, fill=MAP)) +
    scale_fill_viridis_c(option = "H",limits = zlm, oob = scales::squish)+
    geom_polygon(data = reg, aes(x = long, y = lat, group = group), 
                 fill= "#D3D3D3", color = "#AEAEAE") +
    coord_map(xlim = lons, ylim = lats)+
    # formatting
    scale_x_continuous(breaks = ewbrks, labels = ewlbls, expand = c(0, 0)) +
    scale_y_continuous(breaks = nsbrks, labels = nslbls, expand = c(0, 0)) +
    ylab(" ")+xlab(" ")+
    theme_bw() + ggtitle(unique(data$year)[j])+ labs(fill = "Presence mean spatial effect") +
    theme(legend.position = "none",axis.text.x= element_text(size=15), axis.text.y= element_text(size=15),
          plot.title = element_text(size=17,hjust = 0.5),
          panel.grid.major = element_line(colour="darkgrey",size=0.4))
  
  ggsave(paste0(dir_final,"/",unique(data$year)[j]," pr mean bin.jpg",sep=""),dpi=300, width = 8, height = 4)
  
}

# gamma ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

m.prj <- lapply(1:k, function(j) {
  r <- inla.mesh.project(prj,
                         shared$summary.ran$i.con$mean[1:m + (j - 1) * m])
  return(r) 
})
library(fields)
library(RColorBrewer)
colors <- brewer.pal(9, "YlGnBu")
zlm <- range(unlist(m.prj), na.rm = TRUE)
dev.off()

for (j in 1:length(unique(data$year_id))) {
  
  library(raster)
  sp.mean.raster1<-raster(list(x=prj$x,
                               y = prj$y,
                               z = m.prj[[j]]))
  res<-resample(bath,sp.mean.raster1)
  aa<-mask(sp.mean.raster1,res)
  aa <- rasterToPoints(aa)
  aa<-data.frame(aa)
  colnames(aa)<-c("Longitude","Latitude","MAP")
  #Now make the map
  pr<- ggplot() +
    geom_tile(data=aa, aes(y=Latitude, x=Longitude, fill=MAP)) +
    scale_fill_viridis_c(option = "H",limits = zlm, oob = scales::squish)+
    geom_polygon(data = reg, aes(x = long, y = lat, group = group), 
                 fill= "#D3D3D3", color = "#AEAEAE") +
    coord_map(xlim = lons, ylim = lats)+
    # formatting
    scale_x_continuous(breaks = ewbrks, labels = ewlbls, expand = c(0, 0)) +
    scale_y_continuous(breaks = nsbrks, labels = nslbls, expand = c(0, 0)) +
    ylab(" ")+xlab(" ")+
    theme_bw() + ggtitle(unique(data$year)[j])+ labs(fill = "Abundance mean spatial effect") +
    theme(legend.position = "none",axis.text.x= element_text(size=15), axis.text.y= element_text(size=15),
          plot.title = element_text(size=17,hjust = 0.5),
          panel.grid.major = element_line(colour="darkgrey",size=0.4))
  
  
  ggsave(paste0(dir_final,"/",unique(data$year)[j]," pr mean con.jpg",sep=""),dpi=300, width = 8, height = 4)
  
}

# prediction ips -------------------------------------------------------------------------

## Inla Posterior Sample (takes long time) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

shared<-readRDS(file=paste(dir_final,"/final - prof ns - progressive s - year s .rds", sep=""))

n=200
#samples = inla.posterior.sample(n=n,shared)
#saveRDS(samples, paste0(dir_final,"/final model IPS.rds"))

## Load posterior samples ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

samples<-readRDS(paste0(dir_final,"/final model IPS.rds"))
contents=shared$misc$configs$contents
contents

##  INTERCEPT
idx_bin.b0=which(contents$tag=="bin.b0")
idx_bin.b0 = contents$start[idx_bin.b0]:contents$start[idx_bin.b0]-1 + (1:contents$length[idx_bin.b0])
idx_bin.b0
bin.b0 = lapply(samples, function(x) x$latent[idx_bin.b0])
bin.b0 = matrix(unlist(bin.b0), ncol = length(idx_bin.b0), byrow = T)
bin.b0

idx_con.b0=which(contents$tag=="con.b0")
idx_con.b0 = contents$start[idx_con.b0]:contents$start[idx_con.b0]-1 + (1:contents$length[idx_con.b0])
idx_con.b0
con.b0 = lapply(samples, function(x) x$latent[idx_con.b0])
con.b0 = matrix(unlist(con.b0), ncol = length(idx_con.b0), byrow = T)
con.b0


## YEAR iid
idx_year.bin=which(contents$tag=="year.bin")
idx_year.bin = contents$start[idx_year.bin]:contents$start[idx_year.bin]-1 + (1:contents$length[idx_year.bin])
idx_year.bin
year.bin = lapply(samples, function(x) x$latent[idx_year.bin])
year.bin = matrix(unlist(year.bin), ncol = length(idx_year.bin), byrow = T)
year.bin

idx_year.con=which(contents$tag=="year.con")
idx_year.con = contents$start[idx_year.con]:contents$start[idx_year.con]-1 + (1:contents$length[idx_year.con])
idx_year.con
year.con = lapply(samples, function(x) x$latent[idx_year.con])
year.con = matrix(unlist(year.con), ncol = length(idx_year.con), byrow = T)
year.con

### Build lineal predictor ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

col_abu<-"#EF8C1D"

dim(bin.b0)   # binomial intercept
dim(year.con) # gamma trend effect
l<-200

int=c()
year<-matrix( nrow = l, ncol = length(unique(data$year))) # n x 19
pred<-matrix( nrow = l, ncol = length(unique(data$year))) # n x 19

for (i in 1:(l)){
  
  for (j in 1:length(unique(data$year))){
    int[i]=con.b0[i,1]
    year[i,j]=year.con[i,j]
    pred[i,j]= int[i] + year.con[i,j]
  }
}

pred # intercept.bin (1) + year.con (20) = 20 values x 1000 samples (1000 x 19)
a_fun=function(x){exp(x)} ## exponent: transform function for gamma
pred_trans<-a_fun(pred) # abundance kg (original scale)

jpeg(paste0(dir_final,"/predicted trend gamma boxplot.jpeg"), quality = 300, width = 1300, height = 800, res=200)
colnames(pred_trans)<-(seq(from=min(as.numeric(as.character(data$year))), to=max(as.numeric(as.character(data$year))), by=1))
boxplot(pred_trans, col=col_abu,main="exp (int + year.gamma), 500 samples", alpha=0.7, ylim=c(0,15))
dev.off()

write.table(pred_trans, file=paste0(dir_final,"/predicted trend gamma all samples.txt"), append = FALSE, sep = "\t", dec = ".",row.names = TRUE, col.names = TRUE)

## ColMeans~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

predict<-colMeans(pred_trans) # INDEX: Mean of 1000 values each level of year
predict<-cbind(predict, seq(from=min(as.numeric(as.character(data$year))), to=max(as.numeric(as.character(data$year))), by=1))
colnames(predict)<-c("pred","year")
predict<-as.data.frame(predict)
summary(predict)
write.table(predict, file=paste0(dir_final,"/predicted trend gamma all samples.txt"), append = FALSE, sep = "\t", dec = ".",row.names = TRUE, col.names = TRUE)

ggplot(data=predict, aes(x=year, y=pred)) +
  geom_line(color=col_abu) + 
  geom_point(color=col_abu)+
  scale_x_discrete(breaks=seq(min(as.numeric(as.character(data$year))),max(as.numeric(as.character(data$year))),1), limits=predict$year)+
  ylab("Predicted abundance index") + xlab("Year") + theme_light() + ylim(0,4.5)+
  theme(axis.text.x = element_text(angle = 90))
  
ggsave(paste0(dir_final,"Predicted abundance trend by year.jpg"))

# prediction total -------------------------------------------------------------

ips<-readRDS(paste0(dir_final,"/final model IPS.rds"))
contents=shared$misc$configs$contents
contents

## we take samples of the posterior distribution for all model components
## we combine the components to build the linear predictor and des-transform

t_fun=function(x){exp(x)} ## exponent for gamma

psam <- sapply(ips, function(x){
  ## we define function inside of sapply ()
  
  # x<-ips[[1]] # ips object
  library(dplyr)
  
  ## intercept
  intercept <- x$latent %>% rownames(.) %>% stringr::str_detect("^con.b0") %>% x$latent[.,]
  
  ## rw2 year effect
  rw_y<-x$latent %>% rownames(.) %>% stringr::str_detect("^year.con") %>% x$latent[.,]
  
  ## rw2 bathy effect
  rw_b<-x$latent %>% rownames(.) %>% stringr::str_detect("^prof.bin") %>% x$latent[.,]
  
  
  ## spatio-temporal effect
  st<-x$latent %>% rownames(.) %>% stringr::str_detect("^i.con") %>% x$latent[.,]
  
  ## structure besag in spatio-temporal matrix
  stmat<-matrix(st, nrow=mesh$n, 
                 ncol=length(unique(data$year_id)) )
  
  ## df predictor
  predictor<-matrix(0, nrow=mesh$n, 
                    ncol=length(unique(data$year_id)) )
  
  ##  loop through each year + iid for each cell
  for (s in 1:mesh$n){
    
    for (t in 1:length(unique(data$year_id))){ 
      
      ## build linear predictor
      predictor[s,t] = intercept + rw_b[t] + rw_y[t] + stmat[s,t]    
      
    }   
  }
  
  ## loop to pass from matrix to vector
  pred=0
  for (i in 1:length(unique(data$year_id))){
    pred <- c(pred, predictor[,i])
    
  }
  ## remove first 0 row
  pred=pred[-1]
  ## apply des-transform function
  t_fun(pred)
  
})

## in each column we have each of the 100 samples
## in each row we loop each value of s and t
psam

## calculate median and quantiles of the 100 samples
q.sam_al_a <- apply(psam, 1, quantile,
                    c(.025, 0.05, 0.5, 0.95, .975), na.rm =TRUE)

## transpose matrix and pass to data.frame
df_plot <- data.frame(t(q.sam_al_a)) # quantile (cols) and space and time in rows
df_plot$year<-sort(rep(unique(data$year), (mesh$n)))
df_plot$meshX<-(rep((mesh$loc[,1]),length(unique(data$year))))
df_plot$meshY<-(rep((mesh$loc[,2]),length(unique(data$year))))

# map --------------------------------------------------------------------------

library(ggplot2)
library(mapdata)
library(marmap)
library("ggspatial")
library("rnaturalearth")
library("rnaturalearthdata")

### Country MAP
world <- ne_countries(scale = "medium", returnclass = "sf")

### Main MAP
b = getNOAA.bathy(lon1 = -9.8, lon2 = -1, lat1 = 44.2, lat2 = 41.2, 
                  resolution = 1)# get bathymetry data
bf = fortify.bathy(b) # convert bathymetry to data frame 
reg = map_data("world2Hires")#names(reg)#table(reg$region)
reg = subset(reg, region %in% c('Spain', 'Portugal', 'France'))
# set map limits
lons = c(-9.8, -1.5)
lats = c(41.6, 44.3)

# create the breaks- and label vectors
ewbrks <- seq(-9,-1,2)
nsbrks <- seq(42,46,1)
ewlbls <- unlist(lapply(ewbrks, function(x) ifelse(x < 0, paste(-x, "ºW"), ifelse(x > 0, paste(x, "ºE"),x))))
nslbls <- unlist(lapply(nsbrks, function(x) ifelse(x < 0, paste(x, "ºS"), ifelse(x > 0, paste(x, "ºN"),x))))

## Projector MATRIX~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(splancs)
box<-bbox(ocean.cut)
Xrange<-c(box[1,1], box[1,2])
Yrange<-c(box[2,1], box[2,2])
nxy=c(120,120)

prj <- inla.mesh.projector(mesh, xlim=range(Xrange),
                           ylim=range(Yrange), dims=nxy) 

m<-mesh$n
k<-which.max(unique(data$year_id))


# gamma ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

m.prj <- lapply(1:k, function(j) {
  r <- inla.mesh.project(prj,
                         shared$summary.ran$i.con$mean[1:m + (j - 1) * m])
  return(r) 
})
library(fields)
library(RColorBrewer)
colors <- brewer.pal(9, "YlGnBu")
zlm <- range(unlist(m.prj), na.rm = TRUE)
dev.off()

for (j in 1:length(unique(data$year_id))) {
  
  library(raster)
  sp.mean.raster1<-raster(list(x=prj$x,
                               y = prj$y,
                               z = m.prj[[j]]))
  res<-resample(bath,sp.mean.raster1)
  aa<-mask(sp.mean.raster1,res)
  aa <- rasterToPoints(aa)
  aa<-data.frame(aa)
  colnames(aa)<-c("Longitude","Latitude","MAP")
  #Now make the map
  pr<- ggplot() +
    geom_tile(data=aa, aes(y=Latitude, x=Longitude, fill=MAP)) +
    scale_fill_viridis_c(option = "H",limits = zlm, oob = scales::squish)+
    geom_polygon(data = reg, aes(x = long, y = lat, group = group), 
                 fill= "#D3D3D3", color = "#AEAEAE") +
    coord_map(xlim = lons, ylim = lats)+
    # formatting
    scale_x_continuous(breaks = ewbrks, labels = ewlbls, expand = c(0, 0)) +
    scale_y_continuous(breaks = nsbrks, labels = nslbls, expand = c(0, 0)) +
    ylab(" ")+xlab(" ")+
    theme_bw() + ggtitle(unique(data$year)[j])+ labs(fill = "Abundance mean spatial effect") +
    theme(legend.position = "none",axis.text.x= element_text(size=15), axis.text.y= element_text(size=15),
          plot.title = element_text(size=17,hjust = 0.5),
          panel.grid.major = element_line(colour="darkgrey",size=0.4))
  
  
  ggsave(paste0(dir_final,"/",unique(data$year)[j]," pr mean con.jpg",sep=""),dpi=300, width = 8, height = 4)
  
}
