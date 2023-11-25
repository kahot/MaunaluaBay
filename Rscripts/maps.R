library(maptools)
library(ggplot2)
library(sp)
#library(gpclib)
library(base)
library(rgeos)
gpclibPermit()
#theme_set(theme_bw())

library(maps)
library(mapdata)
library(scales) 
library(devtools)
library(dplyr)
library(stringr)

library(ggspatial)

setwd("~/Dropbox/R")
#########################################################


samps <- read.csv("/Users/kahotisthammer/Dropbox/R/FieldSamples.csv")  
coordi<-samps[,c(2:3)]
mbay<-samps[c(1:2),c(2:3)]

#include Kewalo
samps2 <- read.csv("/Users/kahotisthammer/Dropbox/R/FieldSamples2.csv")  
coordi2<-samps2[,c(2:3)]

#Oahu map
hawaii_map<-map_data("worldHires",region=c("Hawaii:Oahu"))
ggplot(hawaii_map, aes(x=long, y=lat, group=group))+geom_polygon(fill="#F4F1D5",colour="black")+
        coord_equal() + ylab("Latitude") + xlab("Longitude")+ 
        theme_bw()+theme(axis.text=element_text(size=12), axis.title=element_text(size=12))+
        theme(panel.background = element_rect(fill = '#D4F5F5', color = 'purple'))+
        annotation_scale(location = "tr", width_hint = 0.5) +
        annotation_north_arrow(location = "tr", which_north = "true", 
                               pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                               style = north_arrow_fancy_orienteering) 


# Mbay only
world.shp <- readShapeSpatial("~/Dropbox/R/TM_WORLD/TM_WORLD_BORDERS_SIMPL-0.3.shp")
hawaiicoast <- readShapeSpatial("~/Dropbox/R/hawaii_coastline/hawaii_coastline.shp")

sites <- read.csv("Data/mbay_sites.csv")  
coordi<-sites[,c(2:3)]



#
mbay<-map(hawaiicoast, xlim=c(-157.74,-157.69),ylim=c(21.255,21.3),col="gray20")
map.axes()
points(sites$lon, sites$lat, pch=16, col="blue", cex=1)

mbay<-map(hawaiicoast,fill=TRUE, col="gray90")
map.axes()
points(sites$lon, sites$lat, pch=16, col="blue", cex=1)

mbay<-map("worldHires", "Hawaii:Oahu", xlim=c(-157.76,-157.69),ylim=c(21.26,21.31), fill=TRUE, col="gray90")
points(samps2$lon, samps2$lat, pch=16, col="black", cex=0.9)
map.axes()

map('state')



ggsave("Output/Oahu_map.png", height = 4, width = 5, units="in", dpi=300)


hawaii_map<-map_data("worldHires",region=c("Hawaii:Oahu","Hawaii:Maui","Hawaii:Lanai","Hawaii:Molokai"))
ggplot(hawaii_map, aes(x=long, y=lat, group=group))+geom_polygon(fill="gray",colour="black")+
        coord_equal() +geom_point(data=coordi, aes(group=NULL), colour="black", shape=16,cex=2.5)+ ylab("Latitude") + xlab("Longitude") 

hawaii_map<-map_data("worldHires",region=c("Hawaii:Oahu"))
ggplot(hawaii_map, aes(x=long, y=lat, group=group))+geom_polygon(fill="gray",colour="black")+
        coord_equal() +geom_point(data=mbay, aes(group=NULL), colour="black", shape=16,cex=2.5)+ ylab("Latitude") + xlab("Longitude")+ theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))

hawaii_map<-map_data("worldHires",region=c("Hawaii:Oahu","Hawaii:Maui","Hawaii:Lanai","Hawaii:Molokai","Hawaii:Kauai","Hawaii:Hawaii","Hawaii:Niihau"))
ggplot(hawaii_map, aes(x=long, y=lat, group=group))+geom_polygon(fill="gray",colour="black")+
        coord_equal() + ylab("Latitude") + xlab("Longitude")+ theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))


hawaii_map<-map_data("worldHires",region=c("Hawaii:Oahu","Hawaii:Maui","Hawaii:Lanai","Hawaii:Molokai","Hawaii:Kauai"))
ggplot(hawaii_map, aes(x=long, y=lat, group=group))+geom_polygon(fill="gray",colour="black")+
        coord_equal() +geom_point(data=coordi, aes(group=NULL), colour="black", shape=16,cex=2.5)+ ylab("Latitude") + xlab("Longitude") 


#### MAP Wahikuli and La Perouse

samps3 <- read.csv("/Users/kahotisthammer/Dropbox/R/FieldSamples_Maui.csv")  
coordi3<-samps3[,c(2:3)]

hawaii_map<-map_data("worldHires",region=c("Hawaii:Maui"))
ggplot(hawaii_map, aes(x=long, y=lat, group=group))+geom_polygon(fill="gray",colour="black")+
        coord_equal() +geom_point(data=coordi3, aes(group=NULL), colour="black", shape=16,cex=1)+ ylab("Latitude") + xlab("Longitude") 



#####
library(sf)
library(raster)
library(dplyr)
library(spData)
library(spDataLarge)

library(tmap)    # for static and interactive maps
#library(leaflet) # for interactive maps
library(ggplot2)


tm_shape(hawaii)+
        tm_fill() +tm_borders() 
tm_shape(nz)+ tm_fill() 
