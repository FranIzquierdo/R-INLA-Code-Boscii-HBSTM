#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# L.boscii spatio-temporal INLA exploratory analysis #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Modified 14/06/2022 #
#~~~~~~~~~~~~~~~~~~~~~~~
# Francisco Izquierdo  #
#~~~~~~~~~~~~~~~~~~~~~~~

## Press Ctrl + Shift + O to see the document outline

# start here! ------------------------------------------------------------------

# clean environment
rm(list=ls()) 

# create dir
dir_ea<-paste(getwd(),"/output/01 exploratoty analysis", sep="")
dir.create(dir_ea)

# read data 
data<-read.table(file="./input/datasets/boscii data INLA 1993-2020.txt", dec=".", header=TRUE)
data_na<-read.table(file="./input/datasets/boscii data INLA 1993-2020 (with NAs).txt", dec=".", header=TRUE)

# set variable type
data$ID<-as.factor(data$ID)
data$year<-as.factor(data$year)
data$year_id<-as.factor(data$year_id)
data$presence<-as.factor(data$presence)

data_na$ID<-as.factor(data_na$ID)
data_na$year<-as.factor(data_na$year)
data_na$year_id<-as.factor(data_na$year_id)
data_na$presence<-as.factor(data_na$presence)

# EA ---------------------------------------------------------------------------

# NAs --------------------------------------------------------------------------

visdat::vis_dat(data_na)
ggsave(paste0(dir_ea, "/variable type with NAs.jpeg"),dpi=300)

visdat::vis_miss(data_na) # missing values per var
ggsave(paste0(dir_ea, "/variable NAs (missing).jpeg"),dpi=300)

# stats table -----------------------------------------------------------------

psych::describe(data)
write.csv(psych::describe(data),paste(dir_ea, "/describe data stats.csv", sep=""))

# Discrete vars -----------------------------------------------------------------

table(data$year) 
table(data$lance)
table(data$presence)

library(inspectdf)
data %>% inspect_cat %>% show_plot
ggsave(paste0(dir_ea, "/discrete variables.jpeg"),dpi=300)

library(ggplot2)
ggplot(data = data, aes(x = year, fill=presence))  +  
  geom_bar(alpha=0.6) + 
  scale_fill_manual(values=c("0"="grey48", "1"="#418FEC"))+
  theme(axis.text.x = element_text(angle = 90)) + 
  ggtitle("Hauls by year")  + 
  ylab("Nº hauls")+xlab("Year") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90))
ggsave(paste0(dir_ea,"/discrete variables - hauls vs presence by year.jpg",sep=""), width=8, height=8, dpi=300)

# Continuous vars ---------------------------------------------------------

library(inspectdf)
data %>% inspect_num %>% show_plot
ggsave(paste0(dir_ea, "/continuous variables.jpeg"),dpi=300)

# Histograms
p0 <- qplot(data$lance, geom="histogram") + xlab("Nº hauls") + theme_light() 
p1 <- qplot(data$numero, geom="histogram") + xlab("Nº fish by 30 min trawl") + theme_light() 
p2 <- qplot(data$temp, geom="histogram")+ xlab("Temperature (ºC)") + theme_light() 
p3 <- qplot(data$sali, geom="histogram")+ xlab("Salinity (psu)") + theme_light() 
p4 <- qplot(data$prof, geom="histogram")+ xlab("Depth (m)") + theme_light() 
p5 <- qplot(data$sedi, geom="histogram")+ xlab("Sediment (?)") + theme_light() 

library(gridExtra)
grid.arrange(p0,p1, p2, p3, p4,p5, nrow = 2)
ggsave(paste0(dir_ea,"/continuous vars histogram.jpg",sep=""), width=8, height=8, dpi=300)

# Correlation ------------------------------------------------------------------

# Corr >0.7 not into the model
sdmdata<-data[,c(9,6,7,8,12)]

# Corrtest 
library(Hmisc)
library(corrplot)
matrix<-rcorr(as.matrix(sdmdata[,-1], type = "pearson"))
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}


jpeg(paste0(dir_ea,"/corrplot.jpeg"), width=1200, height=900, quality = 300, pointsize = 28)
corrplot(matrix$r, type="lower", tl.col = "black",method="number",
         p.mat = matrix$P, sig.level = 0.05)
dev.off()

# performanceanalytics
jpeg(paste0(dir_ea,"/corrplot analytics.jpeg"), width=1200, height=900, quality = 300, pointsize = 28)
library(PerformanceAnalytics)
chart.Correlation(sdmdata,
                  method="pearson",
                  histogram=TRUE,
                  pch=16)
dev.off()

# GGpairs ----------------------------------------------------------------------

# pairs
l<-length(unique(data$year))
jpeg(paste0(dir_ea,"/pairs colour by year.jpeg"), width=1200, height=900, quality = 300, pointsize = 28)
pairs(sdmdata, smooth=TRUE, density=TRUE, method="pearson",
      hist.col="cyan", col=hcl.colors(l, "Temps"[data$year]),
      bg=hcl.colors(l, "Temps"[data$year]))
dev.off()

# ggpairs
library(GGally)
my_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point() + 
    geom_smooth(method=loess, fill="red", color="red", ...) +
    geom_smooth(method=lm, fill="blue", color="blue", ...)
  p
}
g = ggpairs(sdmdata, lower = list(continuous = my_fn))
g
ggsave(paste(dir_ea,"/pairs gg and correlation.jpg",sep=""), width=8, height=8, dpi=300)

# Collinearity -----------------------------------------------------------------

#Support function for corvif
myvif <- function(mod) {
  v <- vcov(mod)
  assign <- attributes(model.matrix(mod))$assign
  if (names(coefficients(mod)[1]) == "(Intercept)") {
    v <- v[-1, -1]
    assign <- assign[-1]
  } else warning("No intercept: vifs may not be sensible.")
  terms <- labels(terms(mod))
  n.terms <- length(terms)
  if (n.terms < 2) stop("The model contains fewer than 2 terms")
  if (length(assign) > dim(v)[1] ) {
    diag(tmp_cor)<-0
    if (any(tmp_cor==1.0)){
      return("Sample size is too small, 100% collinearity is present")
    } else {
      return("Sample size is too small")
    }
  }
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/2Df)")
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] <- length(subs)
  }
  if (all(result[, 2] == 1)) {
    result <- data.frame(GVIF=result[, 1])
  } else {
    result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  }
  invisible(result)
}

corvif <- function(dataz) {
  dataz <- as.data.frame(dataz)
  
  #vif part
  form    <- formula(paste("fooy ~ ",paste(strsplit(names(dataz)," "),collapse=" + ")))
  dataz   <- data.frame(fooy=1 + rnorm(nrow(dataz)) ,dataz)
  lm_mod  <- lm(form,dataz)
  
  cat("\n\nVariance inflation factors\n\n")
  print(myvif(lm_mod))
}

# Multicollinearity
corvif(sdmdata[,-1]) # GVIF > 3 not into the model

# MAP  -------------------------------------------------------------------------

col_pre<-"#418FEC"
col_abu<-"#EF8C1D"

## MAP sampling points
library(ggplot2)
library(mapdata)
library(marmap)
library("ggspatial")
library("rnaturalearth")
library("rnaturalearthdata")

### Mini MAP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

world <- ne_countries(scale = "medium", returnclass = "sf")

ply = data.frame(
  lon = c(-9.8,-9.8,-9.8, -1, -1, -1),
  lat = c(44.3,  43, 42, 42, 43, 44.3))# build a polygon of study area

ggplot(data = world) +
  geom_sf(fill="#D3D3D3", colour="#AEAEAE") + 
  coord_sf(xlim = c(-19, 19), ylim = c(32, 54), expand = FALSE, datum=NA) +
  geom_polygon(data = ply, aes(x = lon, y = lat), color = "black", 
               fill=col_pre, alpha = 0.2, size=0.1) +
  xlab(" ") + ylab(" ") + theme_bw()+  theme(panel.spacing = unit(0.1, "cm"))

ggsave(paste(dir_ea,"/study area minimap.jpg",sep=""), width=8, height=8, dpi=300)

### Main MAP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

b = getNOAA.bathy(lon1 = -9.8, lon2 = -1, lat1 = 44.2, lat2 = 41.2, 
                  resolution = 1)# get bathymetry data
bf = fortify.bathy(b) # convert bathymetry to data frame 
reg = map_data("world2Hires")#names(reg)#table(reg$region)
reg = subset(reg, region %in% c('Spain', 'Portugal', 'France'))
# set map limits
lons = c(-9.8, -1)
lats = c(41.5, 44.5)

# create the breaks- and label vectors
ewbrks <- seq(-8,-1,2)
nsbrks <- seq(42,46,1)
ewlbls <- unlist(lapply(ewbrks, function(x) ifelse(x < 0, paste(-x, "ºW"), ifelse(x > 0, paste(x, "ºE"),x))))
nslbls <- unlist(lapply(nsbrks, function(x) ifelse(x < 0, paste(x, "ºS"), ifelse(x > 0, paste(x, "ºN"),x))))

ggplot()+
  
  # add 100m contour
  geom_contour(data = bf, 
               aes(x=x, y=y, z=z),
               breaks=c(-200),
               size=c(0.3),
               colour="grey")+ 
 
  # add 250m contour
  geom_contour(data = bf, 
               aes(x=x, y=y, z=z),
               breaks=c(-800),
               size=c(0.6),
               colour="grey")+ 
  # add coastline
  geom_polygon(data = reg, aes(x = long, y = lat, group = group), 
               fill= "#D3D3D3", color = "#AEAEAE") +
  
  # Ciudades
  #annotate("text", x=-8.2, y=43.1, label= "La Coruña", size=3, fontface= "bold") + 
  #annotate("text", x=-6, y=43.3, label= "Gijón", size=3, fontface= "bold") + 
  #annotate("text", x=-3.9, y=43.2, label= "Santander", size=3 ,fontface= "bold") + 
  #annotate("text", x=-2.55, y=43, label= "Bilbao", size=3,fontface= "bold") + 
  #annotate("text", x=-8.3, y=42.3, label= "Vigo", size=3 ,fontface= "bold") + 
  
  # add points
  geom_point(data = data, aes(x = long, y = lat, color=factor(presence) ),
             stroke = .5, size = 0.8, 
             alpha = 1, shape = 19)+
  scale_color_manual(values=c("0"="black", "1"="#418FEC"))+
  
  # configure projection and plot domain
  coord_map(xlim = lons, ylim = lats) +
  
  # formatting
  scale_x_continuous(breaks = ewbrks, labels = ewlbls, expand = c(0, 0)) +
  scale_y_continuous(breaks = nsbrks, labels = nslbls, expand = c(0, 0)) +
  ylab(" ")+xlab(" ")+ 
  theme_bw() + 
  theme(legend.position = "bottom" ,axis.text.x= element_text(size=10), axis.text.y= element_text(size=10))+
  labs(subtitle = "All sampling points 1993-2020", col="Presence")

ggsave(paste(dir_ea,"/study area main.jpg",sep=""), width=8, height=8, dpi=300)

## Presence/absence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ewbrks <- seq(-8,-1,3)
nsbrks <- seq(42,46,2)
ewlbls <- unlist(lapply(ewbrks, function(x) ifelse(x < 0, paste(-x, "ºW"), ifelse(x > 0, paste(x, "ºE"),x))))
nslbls <- unlist(lapply(nsbrks, function(x) ifelse(x < 0, paste(x, "ºS"), ifelse(x > 0, paste(x, "ºN"),x))))

ggplot()+
  
   # add coastline
  geom_polygon(data = reg, aes(x = long, y = lat, group = group), 
               fill= "#D3D3D3", color = "#AEAEAE") +
  

  # add points
  geom_point(data = data, aes(x = long, y = lat, color=factor(presence) ),
             stroke = .5, size = 0.5, 
             alpha = 1, shape = 19)+
  scale_color_manual(values=c("0"="black", "1"="#418FEC"))+
  # configure projection and plot domain
  coord_map(xlim = lons, ylim = lats) +
  # formatting
  scale_x_continuous(breaks = ewbrks, labels = ewlbls, expand = c(0, 0)) +
  scale_y_continuous(breaks = nsbrks, labels = nslbls, expand = c(0, 0)) +
  ylab(" ")+xlab(" ")+ facet_wrap(~year,dir="v")+
  theme_bw() + theme(legend.position = "bottom",axis.text.x= element_text(size=8), axis.text.y= element_text(size=8),
                     strip.text = element_text(size = 6))+
  labs(subtitle = "", col="Presence")+
  theme(strip.background = element_rect(fill = "white"))

ggsave(paste(dir_ea,"/map presence by year.jpg",sep=""), width=8, height=8, dpi=300)

## Abundance ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Filter desired years
#library(dplyr)
#data<- data %>% filter( year==2020 | year==2021)

ggplot()+
  
  # add coastline
  geom_polygon(data = reg, aes(x = long, y = lat, group = group), 
               fill= "#D3D3D3", color = "#AEAEAE") +
  
  
  
   # add points
  geom_point(data = data, aes(x = long, y = lat, color=factor(presence), size=numero ),
             stroke = .5, 
             alpha = 0.5, shape = 19)+
  
  scale_color_manual(values = c("grey", "#EF8C1D", "#FC4E07")) +
  scale_size(range = c(0, 10)) +
  #scale_size(range = c(0, 188), name="Nº fish by 30 min trawl")+
  #scale_color_gradient(low="blue", high="red")+

  #scale_color_manual(values=c("0"="black", "1"="orange"))+
  
  # configure projection and plot domain
  coord_map(xlim = lons, ylim = lats) +
  
  # formatting
  scale_x_continuous(breaks = ewbrks, labels = ewlbls, expand = c(0, 0)) +
  scale_y_continuous(breaks = nsbrks, labels = nslbls, expand = c(0, 0)) +
  ylab(" ")+xlab(" ")+ facet_wrap(~year,dir="v")+
  theme_bw() + theme(plot.caption = element_text(hjust = 0),legend.position="none",strip.text = element_text(size = 6),axis.text.x= element_text(size=8), axis.text.y= element_text(size=8))+
  labs(caption="Abundance (number of fish) - bubble size range (1, 188)")+
  theme(strip.background = element_rect(fill = "white"))

ggsave(paste(dir_ea,"/map abundance by year.jpg",sep=""), width=8, height=8, dpi=300)

