library(ggplot2)
library(dplyr)
library(plotrix)
library(cowplot)

#color schemes for site
# N, 3,2.5,K,T, Ref
colors<-c("#843C0B","#C65A12","#F4B183","#F9CBAD","#FBE6D6","#0070C0")

# Biomarker data for Maunalua Bay
biomarkers<-read.csv("Data/biomarkers.csv", fill = T) 

#Biomarker data for Maui
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
    means$Site<-factor(means$Site, levels=c("SiteN","Site3","Site2","SiteT",'SiteK',"Maui"))
    
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
}



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
    

#Stats

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


# Apply anova with tukey-posthoc for passed biomarkers

markers1<-markers[c(1,2,3,4,8,10,11)]
markers2<-markers[c(5,6,7,9)]

stats1<-data.frame(marker=markers1)
stats_pairwise1<-list()
for (i in 1:7){
    df<-biomarkers[,c("Site", markers1[i])]
    colnames(df)[2]<-"amount"
    df<-df[!is.na(df$amount),]
    anova = aov(amount~Site,df)
    sum<-summary(anova)[[1]]
    
    stats1$anova_pval[i]<-sum$`Pr(>F)`[1]
    
    post<-TukeyHSD(x=anova,'Site', conf.level=0.95)
    stats_pairwise1[[i]]<-post$Site
    names(stats_pairwise1)[i]<-markers1[i]
}
stats1
#           marker   anova_pval
#1      Porphyrin 5.250592e-01
#2   DNA.AP.sites 2.066118e-02
#3      Ubiquitin 5.994201e-05
#4          Hsp70 4.977158e-09
#5           MutY 5.382348e-08
#6 Heme.Oxygenase 1.229648e-07
#7            GPx 8.211623e-04

stats_pairwise1
#$Porphyrin
#                   diff       lwr      upr     p adj
#Site3-Site2  0.002456668 -5.273169 5.278082 1.0000000
#SiteK-Site2 -2.393913668 -7.398811 2.610984 0.6028272
#SiteN-Site2 -1.420683332 -6.696309 3.854942 0.9210280
#SiteT-Site2 -0.411099000 -5.415997 4.593799 0.9990471
#SiteK-Site3 -2.396370336 -7.401268 2.608527 0.6019343
#SiteN-Site3 -1.423140000 -6.698765 3.852485 0.9205766
#SiteT-Site3 -0.413555667 -5.418453 4.591342 0.9990245
#SiteN-SiteK  0.973230336 -4.031667 5.978128 0.9745113
#SiteT-SiteK  1.982814668 -2.735848 6.701477 0.7072025
#SiteT-SiteN  1.009584333 -3.995313 6.014482 0.9709100
#
#$DNA.AP.sites
#               diff        lwr       upr      p adj
#Site3-Site2  291.73250  -360.9803  944.4453 0.66429132
#SiteK-Site2 -227.11475  -846.3325  392.1030 0.79970309
#SiteN-Site2  -27.57125  -680.2840  625.1415 0.99993362
#SiteT-Site2  424.48958  -171.3529 1020.3321 0.24109129
#SiteK-Site3 -518.84725 -1138.0650  100.3705 0.12671898
#SiteN-Site3 -319.30375  -972.0165  333.4090 0.58785386
#SiteT-Site3  132.75708  -463.0854  728.5996 0.95962523
#SiteN-SiteK  199.54350  -419.6742  818.7612 0.86313594
#SiteT-SiteK  651.60433    92.6545 1210.5542 0.01797114
#SiteT-SiteN  452.06083  -143.7817 1047.9034 0.19202317
#
#$Ubiquitin
#               diff         lwr       upr        p adj
#Site3-Site2   7.916667  -62.880560  78.71389 9.972297e-01
#SiteK-Site2 -77.833333 -141.156298 -14.51037 1.098192e-02
#SiteN-Site2  55.500000   -7.822965 118.82296 1.053152e-01
#SiteT-Site2 -27.500000  -90.822965  35.82296 7.033919e-01
#SiteK-Site3 -85.750000 -156.547227 -14.95277 1.242193e-02
#SiteN-Site3  47.583333  -23.213893 118.38056 3.034152e-01
#SiteT-Site3 -35.416667 -106.213893  35.38056 5.856996e-01
#SiteN-SiteK 133.333333   70.010369 196.65630 2.177098e-05
#SiteT-SiteK  50.333333  -12.989631 113.65630 1.655715e-01
#SiteT-SiteN -83.000000 -146.322965 -19.67704 6.228851e-03
#
#$Hsp70
#               diff        lwr        upr        p adj
#Site3-Site2   4.000000  -5.530845  13.530845 7.283886e-01
#SiteK-Site2  -9.833333 -18.357980  -1.308687 1.837921e-02
#SiteN-Site2  15.500000   6.975353  24.024647 1.652346e-04
#SiteT-Site2 -13.500000 -22.024647  -4.975353 8.905669e-04
#SiteK-Site3 -13.833333 -23.364178  -4.302489 2.299002e-03
#SiteN-Site3  11.500000   1.969155  21.030845 1.281905e-02
#SiteT-Site3 -17.500000 -27.030845  -7.969155 1.454236e-04
#SiteN-SiteK  25.333333  16.808687  33.857980 7.866270e-08
#SiteT-SiteK  -3.666667 -12.191313   4.857980 7.105707e-01
#SiteT-SiteN -29.000000 -37.524647 -20.475353 6.561907e-09
#
#$MutY
#                 diff         lwr         upr        p adj
#Site3-Site2   1.166667  -9.5720721  11.9054054 9.975235e-01
#SiteK-Site2  -9.333333 -18.9383533   0.2716866 5.962695e-02
#SiteN-Site2 -27.000000 -36.6050200 -17.3949800 2.097591e-07
#SiteT-Site2 -17.166667 -26.7716866  -7.5616467 2.062423e-04
#SiteK-Site3 -10.500000 -21.2387388   0.2387388 5.743602e-02
#SiteN-Site3 -28.166667 -38.9054054 -17.4279279 6.860021e-07
#SiteT-Site3 -18.333333 -29.0720721  -7.5945946 3.662021e-04
#SiteN-SiteK -17.666667 -27.2716866  -8.0616467 1.421690e-04
#SiteT-SiteK  -7.833333 -17.4383533   1.7716866 1.481376e-01
#SiteT-SiteN   9.833333   0.2283134  19.4383533 4.302629e-02
#
#$Heme.Oxygenase
#                  diff         lwr         upr        p adj
#Site3-Site2  -63.25000 -151.736923   25.236923 2.485954e-01
#SiteK-Site2 -161.50000 -240.645110  -82.354890 3.426190e-05
#SiteN-Site2   19.00000  -60.145110   98.145110 9.521784e-01
#SiteT-Site2 -174.83333 -253.978443  -95.688224 1.067065e-05
#SiteK-Site3  -98.25000 -186.736923   -9.763077 2.452036e-02
#SiteN-Site3   82.25000   -6.236923  170.736923 7.708954e-02
#SiteT-Site3 -111.58333 -200.070256  -23.096411 8.811185e-03
#SiteN-SiteK  180.50000  101.354890  259.645110 6.552324e-06
#SiteT-SiteK  -13.33333  -92.478443   65.811776 9.867492e-01
#SiteT-SiteN -193.83333 -272.978443 -114.688224 2.123782e-06
#
#$GPx
#                  diff       lwr        upr       p adj
#Site3-Site2   8.833333 -10.89156 28.5582290 0.679717199
#SiteK-Site2 -12.333333 -29.97582  5.3091497 0.267929782
#SiteN-Site2 -16.833333 -34.47582  0.8091497 0.066419102
#SiteT-Site2 -20.333333 -37.97582 -2.6908503 0.018502887
#SiteK-Site3 -21.166667 -40.89156 -1.4417710 0.031316779
#SiteN-Site3 -25.666667 -45.39156 -5.9417710 0.006657782
#SiteT-Site3 -29.166667 -48.89156 -9.4417710 0.001892002
#SiteN-SiteK  -4.500000 -22.14248 13.1424830 0.941044529
#SiteT-SiteK  -8.000000 -25.64248  9.6424830 0.669778391
#SiteT-SiteN  -3.500000 -21.14248 14.1424830 0.975774397


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

