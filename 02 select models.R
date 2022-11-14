#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# L.boscii spatio-temporal INLA model selection #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Modified 05/06/2022 #
#~~~~~~~~~~~~~~~~~~~~~~~
# Francisco Izquierdo  #
#~~~~~~~~~~~~~~~~~~~~~~~

## Press Ctrl + Shift + O to see the document outline

## Selection of models
## In this script we standardize the L.boscii survey index  
## Bayesian Hierarchical Spatio-temporal Hurdle models
## Binomial-Gamma (presence/ausence, abundance) in WGS84 lonlat
## Spatio-temporal progressive structure not shared (Paradinas et al., 2020)
##  "RW2" effect for covariates. Inla.group n=14 knots.
## Year smoothed "RW2" shared effect.
## Predicted output trend with Inla.posterior.sample

# start here! ------------------------------------------------------------------

# Clean environment
rm(list=ls()) 
library(INLA)
inla.setOption(scale.model.default = TRUE) # set scale.model=TRUE, see scale tutorial

# read data 
data<-read.table(file="./input/datasets/boscii data INLA 1993-2020.txt", dec=".", header=TRUE)

# as. factors
data$year<-as.factor(data$year)

# create dir
dir_sel<-paste(getwd(),"/output/02 select models", sep="")
dir_sel_l<-paste(getwd(),"/output/02 select models/1) envar linear", sep="")
dir_sel_rw<-paste(getwd(),"/output/02 select models/2) envar rw2", sep="")
dir_sel_st<-paste(getwd(),"/output/02 select models/3) spatio-temporal", sep="")

dir.create(dir_sel)
dir.create(dir_sel_l)
dir.create(dir_sel_rw)
dir.create(dir_sel_st)

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
plot(mesh,asp=1, main="")
points(coords, pch=16, cex=0.6)

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
                                      sali.bin=inla.group(data$sali, n=12),
                                      temp.bin=inla.group(data$temp, n=12),
                                      sedi.bin=inla.group(data$sedi, n=12),
                                      year.bin=data$year_id)),
                    tag='est.bin')

est.con<-inla.stack(data=list(y=cbind(NA,ifelse(data$presence>0,data$numero,NA))),
                    A=list(est.temp, 1),
                    effects=list(mesh.index.con,
                                 list(con.b0=1,
                                      prof.con=inla.group(data$prof, n=12),
                                      sali.con=inla.group(data$sali, n=12),
                                      temp.con=inla.group(data$temp, n=12),
                                      sedi.con=inla.group(data$sedi, n=12),
                                      year.con=data$year_id)),
                    tag='est.con')

est.stack<-inla.stack(est.bin,est.con)

## link ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

link=c(rep(1,length(est.stack$data$data[,1])/2), rep(2,length(est.stack$data$data[,1])/2))

# select envars --------------------------------------------------------------------------------------------------------

# 1) linear ---------------------------------------------------------------------------------------------------------------

## linear environmental variable selection, not shared effects

# l 0 - int ------------------------------------------------------------------------------------------------------

l0 <- y ~ -1 + 
  f(bin.b0,model="linear", prec.linear=.1,mean.linear=-.5) + con.b0 

mod <- inla(l0,
            family=c('binomial',"gamma"),
            data=inla.stack.data(est.stack), 
            control.compute=list(dic=TRUE,cpo=TRUE,waic=T,config=TRUE),# important config=TRUE
            control.predictor=list(A=inla.stack.A(est.stack), 
                                   compute=TRUE,link=link),
            verbose=TRUE,control.inla = list(strategy = "gaussian"), 
            num.threads = 1)

saveRDS(mod, file=paste(dir_sel_l,"/l0 - intercept.rds", sep=""))

# l 1 - prof ------------------------------------------------------------------------------------------------------

## not shared 

l1 <- y ~ -1 + 
  f(bin.b0,model="linear", prec.linear=.1,mean.linear=-.5) + con.b0 +
  f(prof.bin, model="linear") + 
  f(prof.con, model="linear") 

mod <- inla(l1,
            family=c('binomial',"gamma"),
            data=inla.stack.data(est.stack), 
            control.compute=list(dic=TRUE,cpo=TRUE,waic=T,config=TRUE),# important config=TRUE
            control.predictor=list(A=inla.stack.A(est.stack), 
                                   compute=TRUE,link=link),
            verbose=TRUE,control.inla = list(strategy = "gaussian"), 
            num.threads = 1)

saveRDS(mod, file=paste(dir_sel_l,"/l1 - prof ns.rds", sep=""))

# l 2 - sali ------------------------------------------------------------------------------------------------------

## not shared

l2 <- y ~ -1 + 
  f(bin.b0,model="linear", prec.linear=.1,mean.linear=-.5) + con.b0 +
  f(sali.bin, model="linear") + 
  f(sali.con, model="linear") 

mod <- inla(l2,
            family=c('binomial',"gamma"),
            data=inla.stack.data(est.stack), 
            control.compute=list(dic=TRUE,cpo=TRUE,waic=T,config=TRUE),# important config=TRUE
            control.predictor=list(A=inla.stack.A(est.stack), 
                                   compute=TRUE,link=link),
            verbose=TRUE,control.inla = list(strategy = "gaussian"), 
            num.threads = 1)

saveRDS(mod, file=paste(dir_sel_l,"/l2 - sali ns.rds", sep=""))

# l 3 - temp ------------------------------------------------------------------------------------------------------

## not shared

l3 <- y ~ -1 + 
  f(bin.b0,model="linear", prec.linear=.1,mean.linear=-.5) + con.b0 +
  f(temp.bin, model="linear") + 
  f(temp.con, model="linear") 

mod <- inla(l3,
            family=c('binomial',"gamma"),
            data=inla.stack.data(est.stack), 
            control.compute=list(dic=TRUE,cpo=TRUE,waic=T,config=TRUE),# important config=TRUE
            control.predictor=list(A=inla.stack.A(est.stack), 
                                   compute=TRUE,link=link),
            verbose=TRUE,control.inla = list(strategy = "gaussian"), 
            num.threads = 1)

saveRDS(mod, file=paste(dir_sel_l,"/l3 - temp ns.rds", sep=""))

# l 4 - sedi ------------------------------------------------------------------------------------------------------

## not shared

l4 <- y ~ -1 + 
  f(bin.b0,model="linear", prec.linear=.1,mean.linear=-.5) + con.b0 +
  f(sedi.bin, model="linear") + 
  f(sedi.con, model="linear") 

mod <- inla(l4,
            family=c('binomial',"gamma"),
            data=inla.stack.data(est.stack), 
            control.compute=list(dic=TRUE,cpo=TRUE,waic=T,config=TRUE),# important config=TRUE
            control.predictor=list(A=inla.stack.A(est.stack), 
                                   compute=TRUE,link=link),
            verbose=TRUE,control.inla = list(strategy = "gaussian"), 
            num.threads = 1)

saveRDS(mod, file=paste(dir_sel_l,"/l4 - sedi ns.rds", sep=""))

# compare --------------------------------------------------------------------------------------------------------

mod0<- readRDS(paste(dir_sel_l,"./l0 - intercept.rds",sep=""))
mod1<- readRDS(paste(dir_sel_l,"./l1 - prof ns.rds",sep=""))
mod2<- readRDS(paste(dir_sel_l,"./l2 - sali ns.rds",sep=""))
mod3<- readRDS(paste(dir_sel_l,"./l3 - temp ns.rds",sep=""))
mod4<- readRDS(paste(dir_sel_l,"./l4 - sedi ns.rds",sep=""))


SEL_ENVAR<-data.frame(
  model=c("intercept","prof ns","sali ns","temp ns","sedi ns"),
  dic=c(mod0$dic$dic,mod1$dic$dic,mod2$dic$dic,mod3$dic$dic,mod4$dic$dic),
  waic=c(mod0$waic$waic,mod1$waic$waic,mod2$waic$waic,mod3$waic$waic,mod4$waic$waic),
  lcpo=c(-mean(log(mod0$cpo$cpo),na.rm=T),-mean(log(mod1$cpo$cpo),na.rm=T),-mean(log(mod2$cpo$cpo),na.rm=T),-mean(log(mod3$cpo$cpo),na.rm=T),-mean(log(mod4$cpo$cpo),na.rm=T)),
  failure=c(sum((mod0$cpo$failure>0)*1,na.rm=T),sum((mod1$cpo$failure>0)*1,na.rm=T),sum((mod2$cpo$failure>0)*1,na.rm=T),sum((mod3$cpo$failure>0)*1,na.rm=T),sum((mod4$cpo$failure>0)*1,na.rm=T)),
  time=c(mod0$cpu.used[4],mod1$cpu.used[4],mod2$cpu.used[4],mod3$cpu.used[4],mod4$cpu.used[4]))

SEL_ENVAR

write.table(SEL_ENVAR, paste(dir_sel_l, "model selection envars linear.txt", sep=""))

# 2) RW2 ---------------------------------------------------------------------------------------------------------------

## smoothed rw2 environmental variable selection

# run 0 - int ------------------------------------------------------------------------------------------------------

r0 <- y ~ -1 + 
  f(bin.b0,model="linear", prec.linear=.1,mean.linear=-.5) + con.b0 

mod <- inla(r0,
            family=c('binomial',"gamma"),
            data=inla.stack.data(est.stack), 
            control.compute=list(dic=TRUE,cpo=TRUE,waic=T,config=TRUE),# important config=TRUE
            control.predictor=list(A=inla.stack.A(est.stack), 
                                   compute=TRUE,link=link),
            verbose=TRUE,control.inla = list(strategy = "gaussian"), 
            num.threads = 1)

saveRDS(mod, file=paste(dir_sel_rw,"/r0 - intercept.rds", sep=""))

# run 1 - prof ------------------------------------------------------------------------------------------------------

## shared

r1 <- y ~ -1 + 
  f(bin.b0,model="linear", prec.linear=.1,mean.linear=-.5) + con.b0 +
  f(prof.bin, model="rw2", hyper=bathy.prec) + 
  f(prof.con, copy="prof.bin", fixed=F) 

mod <- inla(r1,
                  family=c('binomial',"gamma"),
                  data=inla.stack.data(est.stack), 
                  control.compute=list(dic=TRUE,cpo=TRUE,waic=T,config=TRUE),# important config=TRUE
                  control.predictor=list(A=inla.stack.A(est.stack), 
                                         compute=TRUE,link=link),
                  verbose=TRUE,control.inla = list(strategy = "gaussian"), 
                  num.threads = 1)

saveRDS(mod, file=paste(dir_sel_rw,"/r1 - prof s.rds", sep=""))

## not shared

r1 <- y ~ -1 + 
  f(bin.b0,model="linear", prec.linear=.1,mean.linear=-.5) + con.b0 +
  f(prof.bin, model="rw2", hyper=hyper.prec) + 
  f(prof.con, model="rw2", hyper=hyper.prec) 

mod <- inla(r1,
                  family=c('binomial',"gamma"),
                  data=inla.stack.data(est.stack), 
                  control.compute=list(dic=TRUE,cpo=TRUE,waic=T,config=TRUE),# important config=TRUE
                  control.predictor=list(A=inla.stack.A(est.stack), 
                                         compute=TRUE,link=link),
                  verbose=TRUE,control.inla = list(strategy = "gaussian"), 
                  num.threads = 1)

saveRDS(mod, file=paste(dir_sel_rw,"/r1 - prof ns.rds", sep=""))

# run 2 - sali ------------------------------------------------------------------------------------------------------

## shared
r2 <- y ~ -1 + 
  f(bin.b0,model="linear", prec.linear=.1,mean.linear=-.5) + con.b0 +
  f(sali.bin, model="rw2") + 
  f(sali.con, copy="sali.bin", fixed=F) 

mod <- inla(r2,
                  family=c('binomial',"gamma"),
                  data=inla.stack.data(est.stack), 
                  control.compute=list(dic=TRUE,cpo=TRUE,waic=T,config=TRUE),# important config=TRUE
                  control.predictor=list(A=inla.stack.A(est.stack), 
                                         compute=TRUE,link=link),
                  verbose=TRUE,control.inla = list(strategy = "gaussian"), 
                  num.threads = 1)

saveRDS(mod, file=paste(dir_sel_rw,"/r2 - sali s.rds", sep=""))

## not shared

r2 <- y ~ -1 + 
  f(bin.b0,model="linear", prec.linear=.1,mean.linear=-.5) + con.b0 +
  f(sali.bin, model="rw2", hyper=hyper.prec) + 
  f(sali.con, model="rw2", hyper=hyper.prec) 

mod <- inla(r2,
                  family=c('binomial',"gamma"),
                  data=inla.stack.data(est.stack), 
                  control.compute=list(dic=TRUE,cpo=TRUE,waic=T,config=TRUE),# important config=TRUE
                  control.predictor=list(A=inla.stack.A(est.stack), 
                                         compute=TRUE,link=link),
                  verbose=TRUE,control.inla = list(strategy = "gaussian"), 
                  num.threads = 1)

saveRDS(mod, file=paste(dir_sel_rw,"/r2 - sali ns.rds", sep=""))

# run 3 - temp ------------------------------------------------------------------------------------------------------

## shared
r3 <- y ~ -1 + 
  f(bin.b0,model="linear", prec.linear=.1,mean.linear=-.5) + con.b0 +
  f(temp.bin, model="rw2", hyper=bathy.prec) + 
  f(temp.con, copy="temp.bin", fixed=F) 

mod <- inla(r3,
                  family=c('binomial',"gamma"),
                  data=inla.stack.data(est.stack), 
                  control.compute=list(dic=TRUE,cpo=TRUE,waic=T,config=TRUE),# important config=TRUE
                  control.predictor=list(A=inla.stack.A(est.stack), 
                                         compute=TRUE,link=link),
                  verbose=TRUE,control.inla = list(strategy = "gaussian"), 
                  num.threads = 1)

saveRDS(mod, file=paste(dir_sel_rw,"/r3 - temp s.rds", sep=""))

## not shared

r3 <- y ~ -1 + 
  f(bin.b0,model="linear", prec.linear=.1,mean.linear=-.5) + con.b0 +
  f(temp.bin, model="rw2", hyper=hyper.prec) + 
  f(temp.con, model="rw2", hyper=hyper.prec) 

mod <- inla(r3,
                  family=c('binomial',"gamma"),
                  data=inla.stack.data(est.stack), 
                  control.compute=list(dic=TRUE,cpo=TRUE,waic=T,config=TRUE),# important config=TRUE
                  control.predictor=list(A=inla.stack.A(est.stack), 
                                         compute=TRUE,link=link),
                  verbose=TRUE,control.inla = list(strategy = "gaussian"), 
                  num.threads = 1)

saveRDS(mod, file=paste(dir_sel_rw,"/r3 - temp ns.rds", sep=""))

# run 4 - sedi ------------------------------------------------------------------------------------------------------

## shared
r4 <- y ~ -1 + 
  f(bin.b0,model="linear", prec.linear=.1,mean.linear=-.5) + con.b0 +
  f(sedi.bin, model="rw2") + 
  f(sedi.con, copy="sedi.bin", fixed=F) 

mod <- inla(r4,
                  family=c('binomial',"gamma"),
                  data=inla.stack.data(est.stack), 
                  control.compute=list(dic=TRUE,cpo=TRUE,waic=T,config=TRUE),# important config=TRUE
                  control.predictor=list(A=inla.stack.A(est.stack), 
                                         compute=TRUE,link=link),
                  verbose=TRUE,control.inla = list(strategy = "gaussian"), 
                  num.threads = 1)

saveRDS(mod, file=paste(dir_sel_rw,"/r4 - sedi s.rds", sep=""))

## not shared

r4 <- y ~ -1 + 
  f(bin.b0,model="linear", prec.linear=.1,mean.linear=-.5) + con.b0 +
  f(sedi.bin, model="rw2", hyper=hyper.prec) + 
  f(sedi.con, model="rw2", hyper=hyper.prec) 

mod <- inla(r4,
                  family=c('binomial',"gamma"),
                  data=inla.stack.data(est.stack), 
                  control.compute=list(dic=TRUE,cpo=TRUE,waic=T,config=TRUE),# important config=TRUE
                  control.predictor=list(A=inla.stack.A(est.stack), 
                                         compute=TRUE,link=link),
                  verbose=TRUE,control.inla = list(strategy = "gaussian"), 
                  num.threads = 1)

saveRDS(mod, file=paste(dir_sel_rw,"/r4 - sedi ns.rds", sep=""))

# compare --------------------------------------------------------------------------------------------------------

mod0<- readRDS(paste(dir_sel_rw,"./r0 - intercept.rds",sep=""))
mod1<- readRDS(paste(dir_sel_rw,"./r1 - prof s.rds",sep=""))
mod2<- readRDS(paste(dir_sel_rw,"./r1 - prof ns.rds",sep=""))
mod3<- readRDS(paste(dir_sel_rw,"./r2 - sali s.rds",sep=""))
mod4<- readRDS(paste(dir_sel_rw,"./r2 - sali ns.rds",sep=""))
mod5<- readRDS(paste(dir_sel_rw,"./r3 - temp s.rds",sep=""))
mod6<- readRDS(paste(dir_sel_rw,"./r3 - temp ns.rds",sep=""))
mod7<- readRDS(paste(dir_sel_rw,"./r4 - sedi s.rds",sep=""))
mod8<- readRDS(paste(dir_sel_rw,"./r4 - sedi ns.rds",sep=""))


SEL_ENVAR<-data.frame(
  model=c("intercept","prof s","prof ns","sali s","sali ns","temp s","temp ns","sedi s","sedi ns"),
  dic=c(mod0$dic$dic,mod1$dic$dic,mod2$dic$dic,mod3$dic$dic,mod4$dic$dic,mod5$dic$dic,mod6$dic$dic,mod7$dic$dic,mod8$dic$dic),
  waic=c(mod0$waic$waic,mod1$waic$waic,mod2$waic$waic,mod3$waic$waic,mod4$waic$waic,mod5$waic$waic,mod6$waic$waic,mod7$waic$waic,mod8$waic$waic),
  lcpo=c(-mean(log(mod0$cpo$cpo),na.rm=T),-mean(log(mod1$cpo$cpo),na.rm=T),-mean(log(mod2$cpo$cpo),na.rm=T),-mean(log(mod3$cpo$cpo),na.rm=T),-mean(log(mod4$cpo$cpo),na.rm=T),-mean(log(mod5$cpo$cpo),na.rm=T),-mean(log(mod6$cpo$cpo),na.rm=T),-mean(log(mod7$cpo$cpo),na.rm=T),-mean(log(mod8$cpo$cpo),na.rm=T)),
  failure=c(sum((mod0$cpo$failure>0)*1,na.rm=T),sum((mod1$cpo$failure>0)*1,na.rm=T),sum((mod2$cpo$failure>0)*1,na.rm=T),sum((mod3$cpo$failure>0)*1,na.rm=T),sum((mod4$cpo$failure>0)*1,na.rm=T),sum((mod5$cpo$failure>0)*1,na.rm=T),sum((mod6$cpo$failure>0)*1,na.rm=T),sum((mod7$cpo$failure>0)*1,na.rm=T),sum((mod8$cpo$failure>0)*1,na.rm=T)),
  time=c(mod0$cpu.used[4],mod1$cpu.used[4],mod2$cpu.used[4],mod3$cpu.used[4],mod4$cpu.used[4],mod5$cpu.used[4],mod6$cpu.used[4],mod7$cpu.used[4],mod8$cpu.used[4]))

SEL_ENVAR

## not shared effects give lower DIC values for all variables unless sali
## best variables = prof > temp > sedi > sali

write.table(SEL_ENVAR, paste(dir_sel_rw, "model selection envars RW2", sep=""))

# 2.1) Combine rw2 --------------------------------------------------------------------------------------------------------


# run 5 - prof + temp ns -----------------------------------------------------------------------------------------------

## we combine the 2 most relevant variables
## not shared

r5 <- y ~ -1 + 
  f(bin.b0,model="linear", prec.linear=.1,mean.linear=-.5) + con.b0 +
  f(prof.bin, model="rw2", hyper=bathy.prec) + 
  f(prof.con, model="rw2", hyper=bathy.prec) +
  f(temp.bin, model="rw2", hyper=bathy.prec) + 
  f(temp.con, model="rw2", hyper=bathy.prec) 

mod <- inla(r5,
            family=c('binomial',"gamma"),
            data=inla.stack.data(est.stack), 
            control.compute=list(dic=TRUE,cpo=TRUE,waic=T,config=TRUE),# important config=TRUE
            control.predictor=list(A=inla.stack.A(est.stack), 
                                   compute=TRUE,link=link),
            verbose=TRUE,control.inla = list(strategy = "gaussian"), 
            num.threads = 1)

saveRDS(mod, file=paste(dir_sel_rw,"/r5 - prof & temp ns.rds", sep=""))

shared<-readRDS(file=paste(dir_sel_rw,"/r4 - sedi ns.rds", sep=""))

pdf(paste0(dir_sel_rw,"/INLA model output plots - r5 prof s & temp s.pdf"), # File name
    width = 8, height = 7, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk",    # Color model (cmyk is required for most publications)
    paper = "A4")          # Paper size

plot(shared)
dev.off() 

# run 5 - prof + temp s -----------------------------------------------------------------------------------------------

## we combine the 2 most relevant variables
## not shared

r5 <- y ~ -1 + 
  f(bin.b0,model="linear", prec.linear=.1,mean.linear=-.5) + con.b0 +
  f(prof.bin, model="rw2", hyper=bathy.prec) + 
  f(prof.con, copy="prof.bin", fixed=F) +
  f(temp.bin, model="rw2", hyper=bathy.prec) + 
  f(temp.con, copy="temp.bin", fixed=F) 

mod <- inla(r5,
            family=c('binomial',"gamma"),
            data=inla.stack.data(est.stack), 
            control.compute=list(dic=TRUE,cpo=TRUE,waic=T,config=TRUE),# important config=TRUE
            control.predictor=list(A=inla.stack.A(est.stack), 
                                   compute=TRUE,link=link),
            verbose=TRUE,control.inla = list(strategy = "gaussian"), 
            num.threads = 1)

saveRDS(mod, file=paste(dir_sel_rw,"/r5 - prof & temp s.rds", sep=""))

shared<-readRDS(file=paste(dir_sel_rw,"/r9 - prof s & temp s - progressive s.rds", sep=""))

pdf(paste0(dir_sel_st,"/INLA model output plots - r5 prof s & temp s.pdf"), # File name
    width = 8, height = 7, # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk",    # Color model (cmyk is required for most publications)
    paper = "A4")          # Paper size

plot(`r9 - prof ns - progressive s`)
dev.off() 


# run 5 - prof + sedi ns -----------------------------------------------------------------------------------------------

r5 <- y ~ -1 + 
  f(bin.b0,model="linear", prec.linear=.1,mean.linear=-.5) + con.b0 +
  f(prof.bin, model="rw2", hyper=bathy.prec) + 
  f(prof.con, copy="prof.bin", fixed=F) +
  f(sedi.bin, model="rw2", hyper=bathy.prec) + 
  f(sedi.con, model="rw2", hyper=bathy.prec) 

mod <- inla(r5,
            family=c('binomial',"gamma"),
            data=inla.stack.data(est.stack), 
            control.compute=list(dic=TRUE,cpo=TRUE,waic=T,config=TRUE),# important config=TRUE
            control.predictor=list(A=inla.stack.A(est.stack), 
                                   compute=TRUE,link=link),
            verbose=TRUE,control.inla = list(strategy = "gaussian"), 
            num.threads = 1)

saveRDS(mod, file=paste(dir_sel_rw,"/r5 - prof & sedi ns.rds", sep=""))

r5 <- y ~ -1 + 
  f(bin.b0,model="linear", prec.linear=.1,mean.linear=-.5) + con.b0 +
  f(prof.bin, model="rw2", hyper=bathy.prec) + 
  f(prof.con, copy="prof.bin", fixed=F) +
  f(sedi.bin, model="rw2", hyper=bathy.prec) + 
  f(sedi.con, copy="sedi.bin", fixed=F) 

mod <- inla(r5,
            family=c('binomial',"gamma"),
            data=inla.stack.data(est.stack), 
            control.compute=list(dic=TRUE,cpo=TRUE,waic=T,config=TRUE),# important config=TRUE
            control.predictor=list(A=inla.stack.A(est.stack), 
                                   compute=TRUE,link=link),
            verbose=TRUE,control.inla = list(strategy = "gaussian"), 
            num.threads = 1)

saveRDS(mod, file=paste(dir_sel_rw,"/r5 - prof & sedi s.rds", sep=""))


# compare  ----------------------------------------------------------------------------------------------------

mod1<- readRDS(paste(dir_sel,"/3) spatio-temporal/r9 - prof ns & temp (l) s - progressive s.rds",sep=""))
mod2<- readRDS(paste(dir_sel,"/3) spatio-temporal/r9 - prof ns & temp (l) ns - progressive s.rds",sep=""))
mod3<- readRDS(paste(dir_sel,"/3) spatio-temporal/r9 - prof ns - progressive s - year s .rds",sep=""))

mod0<- readRDS(paste(dir_sel,"/3) spatio-temporal/r9 - prof ns - progressive s.rds",sep=""))
mod1 <- readRDS("D:/FRAN/D/R (all)/R code boscii st/output/02 select models/3) spatio-temporal/r6 - prof & temp ns - progressive s.rds")
mod2 <- readRDS("D:/FRAN/D/R (all)/R code boscii st/output/02 select models/3) spatio-temporal/r9 - prof s & sedi s - progressive s.rds")
mod3 <- readRDS("D:/FRAN/D/R (all)/R code boscii st/output/02 select models/3) spatio-temporal/r9 - prof s & temp s - progressive s.rds")

SEL_ENVAR<-data.frame(
  model=c("prof ns - progressive s","prof ns temp ns - progressive s","prof s & sedi s - progressive s", "prof s & temp s - progressive s"),
  dic=c(mod0$dic$dic,mod1$dic$dic,mod2$dic$dic,mod3$dic$dic),
  waic=c(mod0$waic$waic,mod1$waic$waic,mod2$waic$waic,mod3$waic$waic),
  lcpo=c(-mean(log(mod0$cpo$cpo),na.rm=T),-mean(log(mod1$cpo$cpo),na.rm=T),-mean(log(mod2$cpo$cpo),na.rm=T),-mean(log(mod3$cpo$cpo),na.rm=T)),
  failure=c(sum((mod0$cpo$failure>0)*1,na.rm=T),sum((mod1$cpo$failure>0)*1,na.rm=T),sum((mod2$cpo$failure>0)*1,na.rm=T),sum((mod3$cpo$failure>0)*1,na.rm=T)),
  time=c(mod0$cpu.used[4],mod1$cpu.used[4],mod2$cpu.used[4],mod3$cpu.used[4]))

SEL_ENVAR

write.table(SEL_ENVAR, paste(dir_sel_st, "/extra progressive rw2 combined models.txt", sep=""))


# 3) ST --------------------------------------------------------------------------------------------------------

## Try different spatio-temporal structures from Paradinas et al., 2017 
## Progressive, persistent and opportunistic

# run 6 - progressive ----------------------------------------------------------------------------------------------------

## Amat 
est.temp <- inla.spde.make.A(mesh, loc=matrix(c(data$long,data$lat),ncol=2),
                             group=data$year_id,
                             n.group=length(unique(data$year_id)))

## idx 
mesh.index.bin<- inla.spde.make.index("i.bin", 
                                      n.spde=spde$n.spde, 
                                      n.group=max(data$year_id))
mesh.index.con<- inla.spde.make.index("i.con", 
                                      n.spde=spde$n.spde, 
                                      n.group=max(data$year_id))

## stack 

est.bin<-inla.stack(data=list(y=cbind(data$presence,NA)),
                    A=list(est.temp, 1),
                    effects=list(mesh.index.bin,
                                 list(bin.b0=1,
                                      prof.bin=inla.group(data$prof, n=12),
                                      temp.bin=inla.group(data$temp, n=12),
                                      year.bin=data$year_id)),
                    tag='est.bin')

est.con<-inla.stack(data=list(y=cbind(NA,ifelse(data$presence>0,data$numero,NA))),
                    A=list(est.temp, 1),
                    effects=list(mesh.index.con,
                                 list(con.b0=1,
                                      prof.con=inla.group(data$prof, n=12),
                                      temp.con=inla.group(data$temp, n=12),
                                      year.con=data$year_id)),
                    tag='est.con')

est.stack<-inla.stack(est.bin,est.con)

## link
link=c(rep(1,length(est.stack$data$data[,1])/2), rep(2,length(est.stack$data$data[,1])/2))


## shared

r6 <- y ~ -1 + 
  f(bin.b0,model="linear", prec.linear=.1,mean.linear=-.5) + con.b0 +
  f(prof.bin, model="rw2", hyper=bathy.prec) + 
  f(prof.con, model="rw2", hyper=bathy.prec) +
  f(temp.bin, model="rw2", hyper=bathy.prec) + 
  f(temp.con, model="rw2", hyper=bathy.prec) +
  f(i.bin,model=spde,group = i.bin.group,control.group = list(model="ar1")) +
  f(i.con, copy="i.bin", fixed=F)

mod <- inla(r6,
            family=c('binomial',"gamma"),
            data=inla.stack.data(est.stack), 
            control.compute=list(dic=TRUE,cpo=TRUE,waic=T,config=TRUE),# important config=TRUE
            control.predictor=list(A=inla.stack.A(est.stack), 
                                   compute=TRUE,link=link),
            verbose=TRUE,control.inla = list(strategy = "gaussian"), 
            num.threads = 1)

saveRDS(mod, file=paste(dir_sel,"/r6 - prof & temp ns - progressive s.rds", sep=""))

## not shared

r6b <- y ~ -1 + 
  f(bin.b0,model="linear", prec.linear=.1,mean.linear=-.5) + con.b0 +
  f(prof.bin, model="rw2", hyper=bathy.prec) + 
  f(prof.con, model="rw2", hyper=bathy.prec) +
  f(temp.bin, model="rw2", hyper=bathy.prec) + 
  f(temp.con, model="rw2", hyper=bathy.prec) +
  f(i.bin,model=spde,group = i.bin.group,control.group = list(model="ar1")) +
  f(i.con,model=spde,group = i.bin.group,control.group = list(model="ar1")) 
  
modb <- inla(r6b,
            family=c('binomial',"gamma"),
            data=inla.stack.data(est.stack), 
            control.compute=list(dic=TRUE,cpo=TRUE,waic=T,config=TRUE),# important config=TRUE
            control.predictor=list(A=inla.stack.A(est.stack), 
                                   compute=TRUE,link=link),
            verbose=TRUE,control.inla = list(strategy = "gaussian"), 
            num.threads = 1)

saveRDS(modb, file=paste(dir_sel,"/r6 - prof & temp ns - progressive ns.rds", sep=""))

# run 7 - opportunistic ----------------------------------------------------------------------------------------------------

# This model includes the TEMPORAL component associated to the SPATIAL component as a replica

# A. Matrix
bin.temp <- inla.spde.make.A(mesh, loc=cbind(data$long, data$lat),
                             repl=data$year_id,
                             n.repl=length(unique(data$year_id)))

con.temp <- inla.spde.make.A(mesh, loc=cbind(data$long, data$lat),
                             repl=data$year_id,
                             n.repl=length(unique(data$year_id)))

# Index
mesh.index.bin<- inla.spde.make.index("i.bin", n.spde=spde$n.spde, n.repl=length(unique(data$year_id)))
mesh.index.con<- inla.spde.make.index("i.con", n.spde=spde$n.spde, n.repl=length(unique(data$year_id)))

# Stack
bin.stack.repl<-inla.stack(data=list(y=cbind(data$presence,NA)),
                           A=list(bin.temp, 1),
                           effects=list(mesh.index.bin,
                                        list(bin.b0=1,
                                             prof.bin=inla.group(data$prof, n=12),
                                             temp.bin=inla.group(data$temp, n=12))),
                           tag='est_bin')

con.stack.repl<-inla.stack(data=list(y=cbind(NA,ifelse(data$numero>0,data$numero,NA))),
                           A=list(con.temp, 1),
                           effects=list(mesh.index.con,
                                        list(con.b0=1,
                                             prof.con=inla.group(data$prof, n=12),
                                             temp.con=inla.group(data$temp, n=12))),                           
                           tag='est_con')

repl.stack=inla.stack(bin.stack.repl,con.stack.repl)

# Link
link=c(rep(1,length(repl.stack$data$data[,1])/2),rep(2,length(repl.stack$data$data[,1])/2))

## shared

r7 <- y ~ -1 + 
  f(bin.b0,model="linear", prec.linear=.1,mean.linear=-.5) + con.b0 +
  f(prof.bin, model="rw2", hyper=bathy.prec) + 
  f(prof.con, model="rw2", hyper=bathy.prec) +
  f(temp.bin, model="rw2", hyper=bathy.prec) + 
  f(temp.con, model="rw2", hyper=bathy.prec) +
  f(i.bin, model=spde, replicate=i.bin.group) +
  f(i.con, copy="i.bin", fixed=F)

mod <- inla(r7,
            family=c('binomial',"gamma"),
            data=inla.stack.data(repl.stack), 
            control.compute=list(dic=TRUE,cpo=TRUE,waic=T,config=TRUE),# important config=TRUE
            control.predictor=list(A=inla.stack.A(repl.stack), 
                                   compute=TRUE,link=link),
            verbose=TRUE,control.inla = list(strategy = "gaussian"), 
            num.threads = 1)

saveRDS(mod, file=paste(dir_sel,"/r7 - prof & temp ns - year s - opportunistic s.rds", sep=""))

## not shared

r7 <- y ~ -1 + 
  f(bin.b0,model="linear", prec.linear=.1,mean.linear=-.5) + con.b0 +
  f(prof.bin, model="rw2", hyper=bathy.prec) + 
  f(prof.con, model="rw2", hyper=bathy.prec) +
  f(temp.bin, model="rw2", hyper=bathy.prec) + 
  f(temp.con, model="rw2", hyper=bathy.prec) +
  f(i.bin,model=spde, replicate=i.bin.group) +
  f(i.con,model=spde, replicate=i.con.group) 
  
  mod <- inla(r7,
              family=c('binomial',"gamma"),
              data=inla.stack.data(repl.stack), 
              control.compute=list(dic=TRUE,cpo=TRUE,waic=T,config=TRUE),# important config=TRUE
              control.predictor=list(A=inla.stack.A(repl.stack), 
                                     compute=TRUE,link=link),
              verbose=TRUE,control.inla = list(strategy = "gaussian"), 
              num.threads = 1)

saveRDS(mod, file=paste(dir_sel,"/r7 - prof & temp ns - year s - opportunistic ns.rds", sep=""))

# this model showed estimation problems finding the hessian

# run 8 - persistent ----------------------------------------------------------------------------------------------------

# This model includes the TEMPORAL component as identical independent random effects (iid)

# A. Matrix
A.bin <- inla.spde.make.A(mesh, loc=cbind(data$long, data$lat))
A.con <- inla.spde.make.A(mesh, loc=cbind(data$long, data$lat))

# Stack
est.bin<-inla.stack(data=list(y=cbind(data$presence,NA)),
                    A=list(A.bin, 1),
                    effects=list(i.bin=1:spde$n.spde,
                                 list(bin.b0=1,
                                      prof.bin=inla.group(data$prof, n=12),
                                      temp.bin=inla.group(data$temp, n=12),
                                      year.bin=data$year_id)),
                    tag='est.bin')

est.con<-inla.stack(data=list(y=cbind(NA,ifelse(data$numero>0,data$numero,NA))),
                    A=list(A.con, 1),
                    effects=list(i.con=1:spde$n.spde,
                                 list(con.b0=1,
                                      prof.con=inla.group(data$prof, n=12),
                                      temp.con=inla.group(data$temp, n=12),
                                      year.con=ifelse(data$numero>0,data$year_id,NA))),
                    tag='est.con')

pers.stack=inla.stack(est.bin,est.con)

# Link
link=c(rep(1,length(pers.stack$data$data[,1])/2),rep(2,length(pers.stack$data$data[,1])/2))

## shared

r8 <- y ~ -1 + 
  f(bin.b0,model="linear", prec.linear=.1,mean.linear=-.5) + con.b0 +
  f(prof.bin, model="rw2", hyper=bathy.prec) + 
  f(prof.con, model="rw2", hyper=bathy.prec) +
  f(temp.bin, model="rw2", hyper=bathy.prec) + 
  f(temp.con, model="rw2", hyper=bathy.prec) +
  f(year.bin, model="iid", hyper=hyper.prec) + 
  f(year.con, copy="year.bin", fixed=F) +
  f(i.bin,model=spde) +
  f(i.con, copy="i.bin", fixed=F)

mod <- inla(r8,
            family=c('binomial',"gamma"),
            data=inla.stack.data(pers.stack), 
            control.compute=list(dic=TRUE,cpo=TRUE,waic=T,config=TRUE),# important config=TRUE
            control.predictor=list(A=inla.stack.A(pers.stack), 
                                   compute=TRUE,link=link),
            verbose=TRUE,control.inla = list(strategy = "gaussian"), 
            num.threads = 1)

saveRDS(mod, file=paste(dir_sel,"/r8 - prof & temp ns - year s - persistent s.rds", sep=""))

## not shared

r8 <- y ~ -1 + 
  f(bin.b0,model="linear", prec.linear=.1,mean.linear=-.5) + con.b0 +
  f(prof.bin, model="rw2", hyper=bathy.prec) + 
  f(prof.con, model="rw2", hyper=bathy.prec) +
  f(temp.bin, model="rw2", hyper=bathy.prec) + 
  f(temp.con, model="rw2", hyper=bathy.prec) +
  f(year.bin, model="iid", hyper=hyper.prec) + 
  f(year.con, model="iid", hyper=hyper.prec) + 
  f(i.bin,model=spde) +
  f(i.con,model=spde) +
  
  mod <- inla(r8,
              family=c('binomial',"gamma"),
              data=inla.stack.data(est.stack), 
              control.compute=list(dic=TRUE,cpo=TRUE,waic=T,config=TRUE),# important config=TRUE
              control.predictor=list(A=inla.stack.A(est.stack), 
                                     compute=TRUE,link=link),
              verbose=TRUE,control.inla = list(strategy = "gaussian"), 
              num.threads = 1)

saveRDS(mod, file=paste(dir_sel,"/r8 - prof & temp ns - year s - persistent ns.rds", sep=""))


# compare  ----------------------------------------------------------------------------------------------------

mod0<- readRDS(paste(dir_sel,"/3) spatio-temporal/r5 - prof & temp ns.rds",sep=""))
mod1<- readRDS(paste(dir_sel,"/3) spatio-temporal/r6 - prof & temp ns - progressive s.rds",sep=""))
mod2<- readRDS(paste(dir_sel,"/3) spatio-temporal/r6 - prof & temp ns - progressive ns.rds",sep=""))
mod3<- readRDS(paste(dir_sel,"/3) spatio-temporal/r7 - prof & temp ns - opportunistic s.rds",sep=""))
mod4<- readRDS(paste(dir_sel,"/3) spatio-temporal/r7 - prof & temp ns - opportunistic ns.rds",sep=""))
mod5<- readRDS(paste(dir_sel,"/3) spatio-temporal/r8 - prof & temp ns - persistent s.rds",sep=""))
mod6<- readRDS(paste(dir_sel,"/3) spatio-temporal/r8 - prof & temp ns - persistent s.rds",sep="")) # repetir!
mod7<- readRDS(paste(dir_sel,"/3) spatio-temporal/r9 - prof & temp ns - progressive s - year s .rds" ,sep="")) 


SEL_ENVAR<-data.frame(
  model=c("p&t ns","p&t ns prog s","p&t ns prog ns","p&t ns oppo s","p&t ns oppo ns","p&t ns pers s","p&t ns pers ns","p&t ns prog s year s"),
  dic=c(mod0$dic$dic,mod1$dic$dic,mod2$dic$dic,mod3$dic$dic,mod4$dic$dic,mod5$dic$dic,mod6$dic$dic,mod7$dic$dic),
  waic=c(mod0$waic$waic,mod1$waic$waic,mod2$waic$waic,mod3$waic$waic,mod4$waic$waic,mod5$waic$waic,mod6$waic$waic,mod7$waic$waic),
  lcpo=c(-mean(log(mod0$cpo$cpo),na.rm=T),-mean(log(mod1$cpo$cpo),na.rm=T),-mean(log(mod2$cpo$cpo),na.rm=T),-mean(log(mod3$cpo$cpo),na.rm=T),-mean(log(mod4$cpo$cpo),na.rm=T),-mean(log(mod5$cpo$cpo),na.rm=T),-mean(log(mod6$cpo$cpo),na.rm=T),-mean(log(mod7$cpo$cpo),na.rm=T)),
  failure=c(sum((mod0$cpo$failure>0)*1,na.rm=T),sum((mod1$cpo$failure>0)*1,na.rm=T),sum((mod2$cpo$failure>0)*1,na.rm=T),sum((mod3$cpo$failure>0)*1,na.rm=T),sum((mod4$cpo$failure>0)*1,na.rm=T),sum((mod5$cpo$failure>0)*1,na.rm=T),sum((mod6$cpo$failure>0)*1,na.rm=T),sum((mod7$cpo$failure>0)*1,na.rm=T)),
  time=c(mod0$cpu.used[4],mod1$cpu.used[4],mod2$cpu.used[4],mod3$cpu.used[4],mod4$cpu.used[4],mod5$cpu.used[4],mod6$cpu.used[4],mod7$cpu.used[4]))

SEL_ENVAR

write.table(SEL_ENVAR, paste(dir_sel_st, "model selection ST.txt", sep=""))


# run 9 - check ----------------------------------------------------------------------------------------------------------------------

## when including st effect, some rw2 variable effects change their shapes
## variable temp with rw2 shows a plane effect but is an important variable to reduce WAIC
## then we include temp as linear effect and compare it with the model without linear effect


## Amat 
est.temp <- inla.spde.make.A(mesh, loc=matrix(c(data$long,data$lat),ncol=2),
                             group=data$year_id,
                             n.group=length(unique(data$year_id)))

## idx 
mesh.index.bin<- inla.spde.make.index("i.bin", 
                                      n.spde=spde$n.spde, 
                                      n.group=max(data$year_id))
mesh.index.con<- inla.spde.make.index("i.con", 
                                      n.spde=spde$n.spde, 
                                      n.group=max(data$year_id))

## stack 

est.bin<-inla.stack(data=list(y=cbind(data$presence,NA)),
                    A=list(est.temp, 1),
                    effects=list(mesh.index.bin,
                                 list(bin.b0=1,
                                      prof.bin=inla.group(data$prof, n=12),
                                      temp.bin=data$temp,
                                      year.bin=data$year_id)),
                    tag='est.bin')

est.con<-inla.stack(data=list(y=cbind(NA,ifelse(data$presence>0,data$numero,NA))),
                    A=list(est.temp, 1),
                    effects=list(mesh.index.con,
                                 list(con.b0=1,
                                      prof.con=inla.group(data$prof, n=12),
                                      temp.con=data$temp,
                                      year.con=data$year_id)),
                    tag='est.con')

est.stack<-inla.stack(est.bin,est.con)

## link
link=c(rep(1,length(est.stack$data$data[,1])/2), rep(2,length(est.stack$data$data[,1])/2))


## temp linear ns

r9 <- y ~ -1 + 
  f(bin.b0,model="linear", prec.linear=.1,mean.linear=-.5) + con.b0 +
  f(prof.bin, model="rw2", hyper=bathy.prec) + 
  f(prof.con, model="rw2", hyper=bathy.prec) +
  f(temp.bin, model="linear") + 
  f(temp.con, model="linear") +
  f(i.bin,model=spde,group = i.bin.group,control.group = list(model="ar1")) +
  f(i.con, copy="i.bin", fixed=F)

mod <- inla(r9,
            family=c('binomial',"gamma"),
            data=inla.stack.data(est.stack), 
            control.compute=list(dic=TRUE,cpo=TRUE,waic=T,config=TRUE),# important config=TRUE
            control.predictor=list(A=inla.stack.A(est.stack), 
                                   compute=TRUE,link=link),
            verbose=TRUE,control.inla = list(strategy = "gaussian"), 
            num.threads = 1)

saveRDS(mod, file=paste(dir_sel_st,"/r9 - prof ns & temp (l) ns - progressive s.rds", sep=""))


## temp linear s

r9 <- y ~ -1 + 
  f(bin.b0,model="linear", prec.linear=.1,mean.linear=-.5) + con.b0 +
  f(prof.bin, model="rw2", hyper=bathy.prec) + 
  f(prof.con, model="rw2", hyper=bathy.prec) +
  f(temp.bin, model="linear") + 
  f(temp.con, copy="temp.bin") +
  f(i.bin,model=spde,group = i.bin.group,control.group = list(model="ar1")) +
  f(i.con, copy="i.bin", fixed=F)

mod <- inla(r9,
            family=c('binomial',"gamma"),
            data=inla.stack.data(est.stack), 
            control.compute=list(dic=TRUE,cpo=TRUE,waic=T,config=TRUE),# important config=TRUE
            control.predictor=list(A=inla.stack.A(est.stack), 
                                   compute=TRUE,link=link),
            verbose=TRUE,control.inla = list(strategy = "gaussian"), 
            num.threads = 1)

saveRDS(mod, file=paste(dir_sel_st,"/r9 - prof ns & temp (l) s - progressive s.rds", sep=""))

## prof s temp s
 
## without temp

r9 <- y ~ -1 + 
  f(bin.b0,model="linear", prec.linear=.1,mean.linear=-.5) + con.b0 +
  f(prof.bin, model="rw2", hyper=bathy.prec) + 
  f(prof.con, model="rw2", hyper=bathy.prec) +
  f(i.bin,model=spde,group = i.bin.group,control.group = list(model="ar1")) +
  f(i.con, copy="i.bin", fixed=F)

mod <- inla(r9,
            family=c('binomial',"gamma"),
            data=inla.stack.data(est.stack), 
            control.compute=list(dic=TRUE,cpo=TRUE,waic=T,config=TRUE),# important config=TRUE
            control.predictor=list(A=inla.stack.A(est.stack), 
                                   compute=TRUE,link=link),
            verbose=TRUE,control.inla = list(strategy = "gaussian"), 
            num.threads = 1)

saveRDS(mod, file=paste(dir_sel_st,"/r9 - prof ns - progressive s.rds", sep=""))


# add year -----------------------------------------------------------------------------------------------------------------------------

# run 10 - year s ----------------------------------------------------------------------------------------------------------------------

## when including st effect, some rw2 variable effects change their shapes
## variable temp with rw2 shows a plane effect but is an important variable to reduce WAIC
## then we include temp as linear effect and compare it with the model without linear effect


## Amat 
est.temp <- inla.spde.make.A(mesh, loc=matrix(c(data$long,data$lat),ncol=2),
                             group=data$year_id,
                             n.group=length(unique(data$year_id)))

## idx 
mesh.index.bin<- inla.spde.make.index("i.bin", 
                                      n.spde=spde$n.spde, 
                                      n.group=max(data$year_id))
mesh.index.con<- inla.spde.make.index("i.con", 
                                      n.spde=spde$n.spde, 
                                      n.group=max(data$year_id))

## stack 

est.bin<-inla.stack(data=list(y=cbind(data$presence,NA)),
                    A=list(est.temp, 1),
                    effects=list(mesh.index.bin,
                                 list(bin.b0=1,
                                      prof.bin=inla.group(data$prof, n=12),
                                      temp.bin=data$temp,
                                      year.bin=data$year_id)),
                    tag='est.bin')

est.con<-inla.stack(data=list(y=cbind(NA,ifelse(data$presence>0,data$numero,NA))),
                    A=list(est.temp, 1),
                    effects=list(mesh.index.con,
                                 list(con.b0=1,
                                      prof.con=inla.group(data$prof, n=12),
                                      temp.con=data$temp,
                                      year.con=data$year_id)),
                    tag='est.con')

est.stack<-inla.stack(est.bin,est.con)

## link
link=c(rep(1,length(est.stack$data$data[,1])/2), rep(2,length(est.stack$data$data[,1])/2))

## we add year shared rw2 effect in order to derive a shared temporal trend

## shared

r10 <- y ~ -1 + 
  f(bin.b0,model="linear", prec.linear=.1,mean.linear=-.5) + con.b0 +
  f(prof.bin, model="rw2", hyper=bathy.prec) + 
  f(prof.con, model="rw2", hyper=bathy.prec) +
  f(temp.bin, model="linear") + 
  f(temp.con, model="linear") +
  f(year.bin, model="rw2") + 
  f(year.con, copy="year.bin", fixed=F) +
  f(i.bin,model=spde,group = i.bin.group,control.group = list(model="ar1")) +
  f(i.con, copy="i.bin", fixed=F)

mod <- inla(r10,
            family=c('binomial',"gamma"),
            data=inla.stack.data(est.stack), 
            control.compute=list(dic=TRUE,cpo=TRUE,waic=T,config=TRUE),# important config=TRUE
            control.predictor=list(A=inla.stack.A(est.stack), 
                                   compute=TRUE,link=link),
            verbose=TRUE,control.inla = list(strategy = "gaussian"), 
            num.threads = 1)

saveRDS(mod, file=paste(dir_sel_st,"/r10 - prof & temp (l) ns - progressive s - year s .rds", sep=""))


# compare  ----------------------------------------------------------------------------------------------------

mod0<- readRDS(paste(dir_sel,"/3) spatio-temporal/r9 - prof ns - progressive s.rds",sep=""))
mod1<- readRDS(paste(dir_sel,"/3) spatio-temporal/r9 - prof ns & temp (l) s - progressive s.rds",sep=""))
mod2<- readRDS(paste(dir_sel,"/3) spatio-temporal/r9 - prof ns & temp (l) ns - progressive s.rds",sep=""))
mod3<- readRDS(paste(dir_sel,"/3) spatio-temporal/r9 - prof ns - progressive s - year s .rds",sep=""))


SEL_ENVAR<-data.frame(
  model=c("prof ns", "prof ns temp(l) s","prof ns temp(l) ns","prof ns temp(l) ns & year"),
  dic=c(mod0$dic$dic,mod1$dic$dic,mod2$dic$dic,mod3$dic$dic),
  waic=c(mod0$waic$waic,mod1$waic$waic,mod2$waic$waic,mod3$waic$waic),
  lcpo=c(-mean(log(mod0$cpo$cpo),na.rm=T),-mean(log(mod1$cpo$cpo),na.rm=T),-mean(log(mod2$cpo$cpo),na.rm=T),-mean(log(mod3$cpo$cpo),na.rm=T)),
  failure=c(sum((mod0$cpo$failure>0)*1,na.rm=T),sum((mod1$cpo$failure>0)*1,na.rm=T),sum((mod2$cpo$failure>0)*1,na.rm=T),sum((mod3$cpo$failure>0)*1,na.rm=T)),
  time=c(mod0$cpu.used[4],mod1$cpu.used[4],mod2$cpu.used[4],mod3$cpu.used[4]))

SEL_ENVAR

write.table(SEL_ENVAR, paste(dir_sel_st, "model selection final add year.txt", sep=""))

# alternative  ----------------------------------------------------------------------------------------------------
