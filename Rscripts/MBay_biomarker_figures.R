library(ggplot2)
library(dplyr)
library(plotrix)
library(cowplot)
library(gridExtra)

#color schemes for site
# Ref, N, 3,2.5,K,T
colors<-c("#0070C0","#843C0B","#C65A12","#F4B183","#F9CBAD","#FBE6D6")

#Read the data
#Maunalua Bay
biomarkers<-read.csv("Data/biomarkers.csv", fill = T) 
#Maui Means and SE
#Maui Means and SE
maui<-read.csv("Data/maui_means.csv")


# Create figures for all biomarkers

markers<-colnames(biomarkers)[3:13]
# sort the order by appearance in the manuscript text:
markers<-markers[c(3,4,1,6,7,9,11,8,5,2,10)]


plots<-list()
for (i in 1:11){
    marker<-markers[i]
    df<-biomarkers[,c("Site", marker)]
    colnames(df)[2]<-"amount"
    
    #add maui for now till obtain the data
    df<-rbind(df, c("Maui",0))
    df$amount<-as.numeric(df$amount)
    
    means<-df %>% 
        group_by(Site) %>% # our group
        summarise( # summarise operation by group
            mean = mean(amount, na.rm = T),
            std = sd(amount, na.rm = T),
            SE=std.error(amount, na.rm = T))
    
    #Add maui mean & SE
    if(i!=11) {
        means$mean[1]<-maui$mean[maui$Biomarker==marker]
        means$SE[1]<-maui$SE[maui$Biomarker==marker]
    }
    means$Site<-factor(means$Site, levels=c("Maui","SiteN","Site3","Site2","SiteT",'SiteK'))
    
    if (i!= 10) unit<-"fmole/mg TSP"
    if (i==10) unit<-expression(paste("AP sites per", 10^5,"bp"))

    if (i==10) marker = "DNA AP sites"
    if (i==11) marker = "Heme oxygenase"
    
    if (i<9){
        plots[[i]]<-ggplot(means, aes(x=Site,y=mean,ymin=mean-SE, ymax=mean+SE, fill=Site))+
            geom_bar(stat = "identity", width=0.8)+
            geom_errorbar(width=.2, color="gray30", linewidth=0.2)+
            xlab("")+ylab(unit)+
            scale_fill_manual(values=colors, guide="none")+
            theme_classic()+
            theme(axis.text.x=element_blank(),axis.title.y = element_text(size=9))+
            ggtitle(marker)+
            theme(plot.title = element_text(size = 10))
    }
    
    if (i>=9){
    plots[[i]]<-ggplot(means, aes(x=Site,y=mean,ymin=mean-SE, ymax=mean+SE, fill=Site))+
        geom_bar(stat = "identity", width=0.8)+
        geom_errorbar(width=.2, color="gray30", linewidth=0.2)+
        xlab("")+ylab(unit)+
        scale_fill_manual(values=colors, guide="none")+
        theme_classic()+
        theme(axis.text.x=element_text(size=9, angle=45, hjust=1), axis.title.y = element_text(size=9))+
        ggtitle(marker)+
        theme(plot.title = element_text(size = 10))+
        scale_x_discrete(labels=c("SiteN"="Site N","Site2"="Site 2.5","SiteT"="Site T","Site3"="Site 3","SiteK"="Site K"))
    }
    #ggsave(paste0("Output/",marker,".png"), width = 5.5, height=4, dpi=300)
}


#png(paste0("Output/biomarkers_barcharts.png"), height = 7, width = 6.5, res=300, units = "in")
#do.call(grid.arrange, c(plots, ncol=3))
#dev.off()

#png(paste0("Output/biomarkers_barcharts_all.png"), height = 7, width = 6.5, res=300, units = "in")
pdf(paste0("Output/biomarkers_barcharts_all.pdf"), height = 7, width = 6.5)

ggdraw()+
    draw_plot(plots[[1]],x=0,   y=0.65,width=0.33,height=0.2)+
    draw_plot(plots[[2]],x=0.33,y=0.65,width=0.33,height=0.2)+
    draw_plot(plots[[3]],x=0.66,y=0.65,width=0.33,height=0.2)+
    draw_plot(plots[[4]],x=0,   y=0.45,width=0.33,height=0.2)+
    draw_plot(plots[[5]],x=0.33,y=0.45,width=0.33,height=0.2)+
    draw_plot(plots[[6]],x=0.66,y=0.45,width=0.33,height=0.2)+
    draw_plot(plots[[7]],x=0,   y=0.25,width=0.33,height=0.2)+
    draw_plot(plots[[8]],x=0.33,y=0.25,width=0.33,height=0.2)+
    draw_plot(plots[[9]],x=0.66,y=0.2,width=0.33,height=0.25)+
    draw_plot(plots[[10]],x=0,  y=0,width=0.33,height=0.25)+
    draw_plot(plots[[11]],x=0.33,y=0,width=0.33,height=0.25)+
draw_plot_label(c("A","B","C","D","E","F","G","H","I","J","K"), 
                c(0,0.33,0.66,0,0.33,0.66,0,0.33,0.66,0,0.33),
                c(0.85,0.85,0.85,0.65,0.65,0.65,0.45,0.45,0.45,0.25,0.25), size = 11)

dev.off()

#stats

pretest<-data.frame(marker=markers)
for (i in 1:11){
    marker<-markers[i]
    df<-biomarkers[,c("Site", marker)]
    colnames(df)[2]<-"amount"
    
    re1<-shapiro.test(df$amount) 
    pretest$shapiro.test[i]<-re1$p.value
    pretest$shapiro_pass[i]<-ifelse(re1$p.value>0.05, "Y", "N")
    re2<-bartlett.test(amount ~ Site, data=df)
    pretest$barlett.test[i]<-re2$p.value
    pretest$barlett_pass[i]<-ifelse(re2$p.value>0.05, "Y", "N")
}


#           marker shapiro.test shapiro_pass barlett.test barlett_pass
#1       Porphyrin 7.825081e-02            Y 4.643598e-01            Y
#2    DNA.AP.sites 8.350118e-02            Y 1.382346e-01            Y
#3       Ubiquitin 4.122826e-01            Y 4.895295e-01            Y
#4           Hsp70 8.111860e-02            Y 3.935800e-01            Y
#5           Hsp60 2.065272e-05            N 2.915311e-06            N
#6             MxR 3.661060e-04            N 1.073331e-05            N
#7             GST 1.929391e-02            N 3.492335e-01            Y
#8            MutY 4.809779e-01            Y 6.508116e-01            Y
#9        Cu.ZnSOD 1.524208e-03            N 2.927210e-03            N
#10 Heme.Oxygenase 1.491656e-01            Y 6.034722e-01            Y
#11            GPx 7.817460e-02            Y 2.082087e-01            Y




# Test correlation with distance

markers<-colnames(biomarkers)[3:13]
cortest<-c(1,1,1,1,2,2,2,2,2,1,1)
correlation<-data.frame(marker=markers)
Markers<-list()
for (i in 1:11){
    marker<-markers[i]
    df<-biomarkers[,c("Site", marker)]
    colnames(df)[2]<-"amount"
    
    means<-df %>% 
        group_by(Site) %>% # our group
        summarise( # summarise operation by group
            mean = mean(amount, na.rm = T),
            std = sd(amount, na.rm = T),
            SE=std.error(amount, na.rm = T))
    means$distance<-c(0.84, 0.52, 2.16, 0.23, 1.56)
    #means$Site<-factor(means$Site, levels=c("Maui","SiteN","Site3","Site2","SiteK",'SiteT'))
    Markers[[i]]<-means
    names(Markers)[i]<-marker
    
    if(cortest[i]==1) test<-'pearson'
    if(cortest[i]==2) test<-"spearman"
    results<-cor.test(means$mean, means$distance, method=test)
    
    correlation$method[i]<-test
    correlation$p.val[i]<-results$p.value
    correlation$coeff[i]<-results$estimate
}
correlation
#           marker   method       p.val       coeff
#1       Porphyrin  pearson 0.388597414 -0.5021825
#2    DNA.AP.sites  pearson 0.735487472 -0.2092857
#3       Ubiquitin  pearson 0.006942922 -0.9676741*
#4           Hsp70  pearson 0.040054028 -0.8952459*
#5           Hsp60 spearman 0.350000000 -0.6000000
#6             MxR spearman 0.083333333 -0.9000000.
#7             GST spearman 0.233333333 -0.7000000
#8            MutY spearman 0.950000000  0.1000000
#9        Cu.ZnSOD spearman 0.133333333 -0.8000000
#10 Heme.Oxygenase  pearson 0.054634439 -0.8708322.
#11            GPx  pearson 0.533818911 -0.3751319


# Create correlation plot for significant marker
corrs<-correlation[correlation$p.val<0.1,]
corrs<-corrs[order(corrs$p.val),]
corplots<-list()
for (i in 1:nrow(corrs)){
    marker<-corrs$marker[i]
    df<-Markers[[marker]]
    df$Site<-factor(df$Site,levels=c("SiteN","Site3","Site2","SiteT",'SiteK'))
    
    ymax<-max(df$mean)*1.05
    ymin<-min(df$mean)*0.8
    pval<-round(corrs$p.val[i], digits = 3)
    r=round(corrs$coeff[i], digits=3)
    if(pval<0.05&pval>0.01) star="*"
    if(pval<0.01) star="**"
    if (pval>0.05) star="."
    if (i<3) xlab=""
    if (i>2) xlab="Distance (km)"
    if (i==3) marker="Heme oxygenase"
    
    
    corplots[[i]]<-ggplot(df, aes(x=distance, y=mean, fill=Site))+
        geom_point(size=3.8, shape=21, color="gray50")+
        theme_classic()+xlab(xlab)+
        ylab('fmole/mg TSP')+
        ylim(ymin,ymax)+
        ggtitle(marker)+theme(plot.title = element_text(size = 10))+
        scale_fill_manual(values=colors[2:6], guide="none")+
        theme(axis.text=element_text(size=9))+
        annotate("text", x=1.8, y=ymax, label=paste0("r = ", r,star))
        #geom_smooth(method="lm", se=FALSE,color="gray", linewidth=0.5, linetype=2)
    
}

# Create legend for the correlation plot
ggplot(df, aes(x=distance, y=mean, fill=Site))+
    geom_point(size=3.8, shape=21, color="gray50")+
    theme_classic()+xlab(xlab)+
    ylab('fmole/mg TSP')+
    ggtitle(marker)+theme(plot.title = element_text(size = 10))+
    scale_fill_manual(values=colors[2:6], labels=c("Site N","Site 3","Site 2.5","Site T",'Site K'))+
    theme(axis.text=element_text(size=9))+
    annotate("text", x=1.8, y=ymax, label=paste0("r = ", pval,star))+
    theme(legend.title=element_blank(),
        legend.position = c(.95, .85),
        legend.justification = c("right", "top")
     )
ggsave("Output/corplot_legend.png", width = 5, height = 4, dpi=300)

#png(paste0("Output/biomarkers_correlation_plots.png"), height = 5, width = 7, res=300, units = "in")
pdf(paste0("Output/biomarkers_correlation_plots2.pdf"), height = 5, width = 7)

ggdraw()+
    draw_plot(corplots[[1]],x=0,y=0.5,width=0.5,height=0.49)+
    draw_plot(corplots[[2]],x=0.5,y=0.5,width=0.5,height=0.49)+
    draw_plot(corplots[[3]],x=0,y=0,width=0.5,height=0.5)+
    draw_plot(corplots[[4]],x=0.5, y=0,width=0.5,height=0.5)+
    draw_plot_label(c("A","B","C","D"), 
                    c(0,0.5,0,0.5),
                    c(0.99,0.99,0.5,0.5), size = 11)

dev.off()

