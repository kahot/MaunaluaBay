library(ggplot2)
library(dplyr)
library(plotrix)
library(cowplot)
library(reshape2)

#color schemes for site
# Ref, N, 3,2.5,K,T
colors<-c("#0070C0","#843C0B","#C65A12","#F4B183","#F9CBAD","#FBE6D6")

sediments<-read.csv("Data/sediment_analysis_results.csv")
sedi<-melt(sediments, id.vars=c("type", "compound"))
colnames(sedi)[3:4]<-c("Site","Level")
sedi$Level[is.na(sedi$Level)]<-0
sedi$Level<-factor(sedi$Level, levels=c(0,1,2,3))

ggplot(sedi, aes(x=Site, y=compound, fill=Level)) + 
    geom_tile(color="white", size=0.05)+
    scale_fill_manual(values=c("white","#FFFECD","#FFD876","#FE8E3C"),name="",labels=c("Not detected","Suspected","Presence", "Hit"))+
    theme_bw()
    theme(axis.text.x = element_text(angle=45, hjust=1))+
    theme(legend.title = element_text(size=10), axis.title=element_blank())+
    facet_grid(treatment~origin, scales="free",space="free", labeller=plot_labeller)+
    theme(strip.text.y = element_text(angle = 0), strip.background.x = element_rect(color="gray20"))+
    theme(axis.text.y=element_text(size=7, hjust=0))+
    theme(panel.spacing.x=unit(0, "lines"))+
    theme(legend.text.align = 0, legend.text = element_text(size=8))

    
    
