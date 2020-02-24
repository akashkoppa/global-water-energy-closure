# Script to plot cluster metrics in the Budyko space
# Author - Akash Koppa
# Date - 2019-08-05

# clear workspace
rm(list=ls())

# load required libraries
library(ggplot2)
library(raster)
library(stringr)
library(maptools)
library(plyr)
library(RColorBrewer)
library(reshape2)
library(viridis)

# set working directory
setwd("/home/akash/Documents/FifthPaper/Budyko_Analysis_New/")

# read in the shapefile
shp = shapefile("Input/HydroSHEDS_All_Level05/hybas_all_lev05_v1c.shp")

# read in the distance metric results
dist_metric = read.table("Output/Budyko_Distance_HydroSHEDS_Lvl05_Rn-CERES_noSM2RAIN.txt",header=TRUE)

# Plot of ET cluster density vs P cluster density with aridity
cluster_sens = read.table("Output/Budyko_Cluster_Density_Sensitivity_noSM2RAIN.txt",header=TRUE)

# separate the area and AI data
hydrosheds_id = as.character(dist_metric$HydroSHEDS_ID)
area = dist_metric$Area
ai = dist_metric$AI
dist_metric = dist_metric[,-c(1:3)]

# normalize the cluster metrics by the aridity index
dist_metric = dist_metric/ai

# normalize the sensitivity metric by aridity index
cluster_sens = cluster_sens/ai
data_cluster_sens = data.frame(P = cluster_sens$CD_P_Mean,
                               ET = cluster_sens$CD_ET_Mean,
                               Aridity = ai)
#data_cluster_sens_ratio = data.frame(Ratio = data_)

# replace radius metric > 100 as NA
#cluster_metric$Max_Radius[which(cluster_metric$Max_Radius > 10)] = NA
#cluster_metric$Mean_Radius[which(cluster_metric$Mean_Radius > 10)] = NA
#cluster_metric$Min_Radius[which(cluster_metric$Min_Radius > 10)] = NA
#cluster_metric$Sd_Radius[which(cluster_metric$Sd_Radius > 10)] = NA
#cluster_metric$Density_Max_Radius[which(cluster_metric$Density_Max_Radius > 500)] = NA
#cluster_metric$Density_Mean_Radius[which(cluster_metric$Density_Mean_Radius > 500)] = NA

#cluster_metric = apply(cluster_metric, MARGIN = 2, FUN = function(x)(x - min(x, na.rm=TRUE))/diff(range(x,na.rm=TRUE)))
#cluster_metric = data.frame(cluster_metric)

# calculate mean, standard deviation, min, and max
cluster_metric = matrix(data = NA, nrow=nrow(dist_metric),ncol = 4)
cluster_metric = data.frame(cluster_metric)
colnames(cluster_metric) = c("Max_Radius","Mean_Radius","Min_Radius","Sd_Radius")
cluster_metric$Max_Radius  = apply(dist_metric,1,max,na.rm=TRUE)
cluster_metric$Max_Radius[which(cluster_metric$Max_Radius==Inf)] = NA
cluster_metric$Max_Radius[which(cluster_metric$Max_Radius==-Inf)] = NA
cluster_metric$Mean_Radius = apply(dist_metric,1,mean,na.rm=TRUE)
cluster_metric$Min_Radius  = apply(dist_metric,1,min,na.rm=TRUE)
cluster_metric$Min_Radius[which(cluster_metric$Min_Radius==Inf)] = NA
cluster_metric$Min_Radius[which(cluster_metric$Min_Radius==-Inf)] = NA
cluster_metric$Sd_Radius   = apply(dist_metric,1,sd,na.rm=TRUE)

# create different classes for each of the metric (maybe five or six classes each) and then assign colors to them.
# I think that will create a better map than what I have now. 
# The bins are created based on quantiles
# Mean Radius
quantiles_meanrad = quantile(cluster_metric$Mean_Radius, probs=c(0.1,0.25,0.5,0.75,0.9),na.rm=TRUE)
bins_meanrad = rep(NA,nrow(cluster_metric))
bins_meanrad[which(cluster_metric$Mean_Radius <= quantiles_meanrad[1])]                                                     = paste0("0.0","-", round(quantiles_meanrad[1],1))
bins_meanrad[which(cluster_metric$Mean_Radius <= quantiles_meanrad[2] & cluster_metric$Mean_Radius > quantiles_meanrad[1])] = paste0(round(quantiles_meanrad[1],1),"-", round(quantiles_meanrad[2],1))
bins_meanrad[which(cluster_metric$Mean_Radius <= quantiles_meanrad[3] & cluster_metric$Mean_Radius > quantiles_meanrad[2])] = paste0(round(quantiles_meanrad[2],1),"-", round(quantiles_meanrad[3],1))
bins_meanrad[which(cluster_metric$Mean_Radius <= quantiles_meanrad[4] & cluster_metric$Mean_Radius > quantiles_meanrad[3])] = paste0(round(quantiles_meanrad[3],1),"-", round(quantiles_meanrad[4],1))
bins_meanrad[which(cluster_metric$Mean_Radius <= quantiles_meanrad[5] & cluster_metric$Mean_Radius > quantiles_meanrad[4])] = paste0(round(quantiles_meanrad[4],1),"-", round(quantiles_meanrad[5],1))
bins_meanrad[which(cluster_metric$Mean_Radius >  quantiles_meanrad[5])]                                                     = paste0("> ",round(quantiles_meanrad[5],1))

# Standard deviation of Radius
quantiles_stdrad = quantile(cluster_metric$Sd_Radius, probs=c(0.1,0.25,0.5,0.75,0.9),na.rm=TRUE)
bins_stdrad = rep(NA,nrow(cluster_metric))
bins_stdrad[which(cluster_metric$Sd_Radius <= quantiles_stdrad[1])]                                                     = paste0("0.0","-", round(quantiles_stdrad[1],3))
bins_stdrad[which(cluster_metric$Sd_Radius <= quantiles_stdrad[2] & cluster_metric$Sd_Radius > quantiles_stdrad[1])] = paste0(round(quantiles_stdrad[1],2),"-", round(quantiles_stdrad[2],2))
bins_stdrad[which(cluster_metric$Sd_Radius <= quantiles_stdrad[3] & cluster_metric$Sd_Radius > quantiles_stdrad[2])] = paste0(round(quantiles_stdrad[2],2),"-", round(quantiles_stdrad[3],2))
bins_stdrad[which(cluster_metric$Sd_Radius <= quantiles_stdrad[4] & cluster_metric$Sd_Radius > quantiles_stdrad[3])] = paste0(round(quantiles_stdrad[3],2),"-", round(quantiles_stdrad[4],2))
bins_stdrad[which(cluster_metric$Sd_Radius <= quantiles_stdrad[5] & cluster_metric$Sd_Radius > quantiles_stdrad[4])] = paste0(round(quantiles_stdrad[4],2),"-", round(quantiles_stdrad[5],2))
bins_stdrad[which(cluster_metric$Sd_Radius >  quantiles_stdrad[5])]                                                     = paste0("> ",round(quantiles_stdrad[5],2))

# Min Radius
quantiles_minrad = quantile(cluster_metric$Min_Radius, probs=c(0.1,0.25,0.5,0.75,0.9),na.rm=TRUE)
bins_minrad = rep(NA,nrow(cluster_metric))
bins_minrad[which(cluster_metric$Min_Radius <= quantiles_minrad[1])]                                                     = paste0("0.0","-", round(quantiles_minrad[1],2))
bins_minrad[which(cluster_metric$Min_Radius <= quantiles_minrad[2] & cluster_metric$Min_Radius > quantiles_minrad[1])] = paste0(round(quantiles_minrad[1],2),"-", round(quantiles_minrad[2],2))
bins_minrad[which(cluster_metric$Min_Radius <= quantiles_minrad[3] & cluster_metric$Min_Radius > quantiles_minrad[2])] = paste0(round(quantiles_minrad[2],2),"-", round(quantiles_minrad[3],2))
bins_minrad[which(cluster_metric$Min_Radius <= quantiles_minrad[4] & cluster_metric$Min_Radius > quantiles_minrad[3])] = paste0(round(quantiles_minrad[3],2),"-", round(quantiles_minrad[4],2))
bins_minrad[which(cluster_metric$Min_Radius <= quantiles_minrad[5] & cluster_metric$Min_Radius > quantiles_minrad[4])] = paste0(round(quantiles_minrad[4],2),"-", round(quantiles_minrad[5],2))
bins_minrad[which(cluster_metric$Min_Radius >  quantiles_minrad[5])]                                                     = paste0("> ",round(quantiles_minrad[5],2))


# Cluster Sensitivity for Precipitation
quantiles_sens_p = quantile(data_cluster_sens$P, probs=c(0.1,0.25,0.5,0.75,0.9),na.rm=TRUE)
bins_sens_p = rep(NA,nrow(data_cluster_sens))
bins_sens_p[which(data_cluster_sens$P <= quantiles_sens_p[1])]                                                     = paste0("0.0","-", round(quantiles_sens_p[1],2))
bins_sens_p[which(data_cluster_sens$P <= quantiles_sens_p[2] & data_cluster_sens$P > quantiles_sens_p[1])] = paste0(round(quantiles_sens_p[1],2),"-", round(quantiles_sens_p[2],2))
bins_sens_p[which(data_cluster_sens$P <= quantiles_sens_p[3] & data_cluster_sens$P > quantiles_sens_p[2])] = paste0(round(quantiles_sens_p[2],2),"-", round(quantiles_sens_p[3],2))
bins_sens_p[which(data_cluster_sens$P <= quantiles_sens_p[4] & data_cluster_sens$P > quantiles_sens_p[3])] = paste0(round(quantiles_sens_p[3],2),"-", round(quantiles_sens_p[4],2))
bins_sens_p[which(data_cluster_sens$P <= quantiles_sens_p[5] & data_cluster_sens$P > quantiles_sens_p[4])] = paste0(round(quantiles_sens_p[4],2),"-", round(quantiles_sens_p[5],2))
bins_sens_p[which(data_cluster_sens$P >  quantiles_sens_p[5])]                                                     = paste0("> ",round(quantiles_sens_p[5],2))


# Cluster Sensitivity for Evapotranspiration
quantiles_sens_et = quantile(data_cluster_sens$ET, probs=c(0.1,0.25,0.5,0.75,0.9),na.rm=TRUE)
bins_sens_et = rep(NA,nrow(data_cluster_sens))
#bins_sens_et[which(data_cluster_sens$ET <= quantiles_sens_et[1])]                                                     = paste0("0.0","-", round(quantiles_sens_et[1],2))
#bins_sens_et[which(data_cluster_sens$ET <= quantiles_sens_et[2] & data_cluster_sens$ET > quantiles_sens_et[1])] = paste0(round(quantiles_sens_et[1],2),"-", round(quantiles_sens_et[2],2))
#bins_sens_et[which(data_cluster_sens$ET <= quantiles_sens_et[3] & data_cluster_sens$ET > quantiles_sens_et[2])] = paste0(round(quantiles_sens_et[2],2),"-", round(quantiles_sens_et[3],2))
#bins_sens_et[which(data_cluster_sens$ET <= quantiles_sens_et[4] & data_cluster_sens$ET > quantiles_sens_et[3])] = paste0(round(quantiles_sens_et[3],2),"-", round(quantiles_sens_et[4],2))
#bins_sens_et[which(data_cluster_sens$ET <= quantiles_sens_et[5] & data_cluster_sens$ET > quantiles_sens_et[4])] = paste0(round(quantiles_sens_et[4],2),"-", round(quantiles_sens_et[5],2))
#bins_sens_et[which(data_cluster_sens$ET >  quantiles_sens_et[5])]                                                     = paste0("> ",round(quantiles_sens_et[5],2))
# harmonize the legend
bins_sens_et[which(data_cluster_sens$ET <= quantiles_sens_p[1])]                                                     = paste0("0.0","-", round(quantiles_sens_p[1],2))
bins_sens_et[which(data_cluster_sens$ET <= quantiles_sens_p[2] & data_cluster_sens$ET > quantiles_sens_p[1])] = paste0(round(quantiles_sens_p[1],2),"-", round(quantiles_sens_p[2],2))
bins_sens_et[which(data_cluster_sens$ET <= quantiles_sens_p[3] & data_cluster_sens$ET > quantiles_sens_p[2])] = paste0(round(quantiles_sens_p[2],2),"-", round(quantiles_sens_p[3],2))
bins_sens_et[which(data_cluster_sens$ET <= quantiles_sens_p[4] & data_cluster_sens$ET > quantiles_sens_p[3])] = paste0(round(quantiles_sens_p[3],2),"-", round(quantiles_sens_p[4],2))
bins_sens_et[which(data_cluster_sens$ET <= quantiles_sens_p[5] & data_cluster_sens$ET > quantiles_sens_p[4])] = paste0(round(quantiles_sens_p[4],2),"-", round(quantiles_sens_p[5],2))
bins_sens_et[which(data_cluster_sens$ET >  quantiles_sens_p[5])]                                                     = paste0("> ",round(quantiles_sens_p[5],2))

# ratio of mean P to mean ET
quantiles_sens_p2et = quantile(data_cluster_sens$P/data_cluster_sens$ET, probs=c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm=TRUE)
bins_sens_p2et = rep(NA,nrow(data_cluster_sens))
bins_sens_p2et[which(data_cluster_sens$P/data_cluster_sens$ET <= 1)]                                                     = paste0("0.0","-", "1.00")
bins_sens_p2et[which(data_cluster_sens$P/data_cluster_sens$ET <= quantiles_sens_p2et[2] & data_cluster_sens$P/data_cluster_sens$ET > 1)] = paste0("1.00","-", round(quantiles_sens_p2et[2],2))
bins_sens_p2et[which(data_cluster_sens$P/data_cluster_sens$ET <= quantiles_sens_p2et[3] & data_cluster_sens$P/data_cluster_sens$ET > quantiles_sens_p2et[2])] = paste0(round(quantiles_sens_p2et[2],2),"-", round(quantiles_sens_p2et[3],2))
bins_sens_p2et[which(data_cluster_sens$P/data_cluster_sens$ET <= quantiles_sens_p2et[4] & data_cluster_sens$P/data_cluster_sens$ET > quantiles_sens_p2et[3])] = paste0(round(quantiles_sens_p2et[3],2),"-", round(quantiles_sens_p2et[4],2))
bins_sens_p2et[which(data_cluster_sens$P/data_cluster_sens$ET <= quantiles_sens_p2et[5] & data_cluster_sens$P/data_cluster_sens$ET > quantiles_sens_p2et[4])] = paste0(round(quantiles_sens_p2et[4],2),"-", round(quantiles_sens_p2et[5],2))
bins_sens_p2et[which(data_cluster_sens$P/data_cluster_sens$ET >  quantiles_sens_p2et[5])]                                                     = paste0("> ",round(quantiles_sens_p2et[5],2))

flag = 0 
if(flag>1){
# Max Radius
quantiles_maxrad = quantile(cluster_metric$Max_Radius, probs=c(0.1,0.25,0.5,0.75,0.9),na.rm=TRUE)
bins_maxrad = rep(NA,nrow(cluster_metric))
bins_maxrad[which(cluster_metric$Max_Radius <= quantiles_maxrad[1])] = paste0("0.0"," - ", round(quantiles_maxrad[1],1))
bins_maxrad[which(cluster_metric$Max_Radius <= quantiles_maxrad[2] & cluster_metric$Max_Radius > quantiles_maxrad[1])] = paste0(round(quantiles_maxrad[1],1)," - ", round(quantiles_maxrad[2],1))
bins_maxrad[which(cluster_metric$Max_Radius <= quantiles_maxrad[3] & cluster_metric$Max_Radius > quantiles_maxrad[2])] = paste0(round(quantiles_maxrad[2],1)," - ", round(quantiles_maxrad[3],1))
bins_maxrad[which(cluster_metric$Max_Radius <= quantiles_maxrad[4] & cluster_metric$Max_Radius > quantiles_maxrad[3])] = paste0(round(quantiles_maxrad[3],1)," - ", round(quantiles_maxrad[4],1))
bins_maxrad[which(cluster_metric$Max_Radius <= quantiles_maxrad[5] & cluster_metric$Max_Radius > quantiles_maxrad[4])] = paste0(round(quantiles_maxrad[4],1)," - ", round(quantiles_maxrad[5],1))
bins_maxrad[which(cluster_metric$Max_Radius > quantiles_maxrad[5])] = paste0("> ",round(quantiles_maxrad[5],1))


}
# plot the 1) max, 2) mean, 3) min, and 4) standard deviation of cluster radius, 
# 5) Density (max radius), 6) Density (mean radius)

# Prepare the data for ggplot
#shp$Cluster_Radius_Max  = cluster_metric$Max_Radius
#shp$Cluster_Radius_Mean = cluster_metric$Mean_Radius
#shp$Cluster_Radius_Min  = cluster_metric$Min_Radius
#shp$Cluster_Radius_SD   = cluster_metric$Sd_Radius
#shp$Density_Max_Radius  = cluster_metric$Density_Max_Radius
#shp$Density_Mean_Radius = cluster_metric$Density_Mean_Radius

# based on bins
#shp$Cluster_Radius_Max  = bins_maxrad
shp$Cluster_Radius_Mean = bins_meanrad
shp$Cluster_Radius_Min  = bins_minrad
shp$Cluster_Radius_SD   = bins_stdrad
shp$Sensitivity_P       = bins_sens_p
shp$Sensitivity_ET      = bins_sens_et
shp$Sensitivity_P2ET    = bins_sens_p2et

shp@data$id = rownames(shp@data)
shp_df = fortify(shp)
shp_df = join(shp_df,shp@data,by="id")

flag = 0
if (flag>0){
# Sensitivity of Cluster Radius (Mean P)
levels_reqd_sens_p = c(paste0("0.0","-", round(quantiles_sens_p[1],2)),
                       paste0(round(quantiles_sens_p[1],2),"–", round(quantiles_sens_p[2],2)),
                       paste0(round(quantiles_sens_p[2],2),"-", round(quantiles_sens_p[3],2)),
                       paste0(round(quantiles_sens_p[3],2),"-", round(quantiles_sens_p[4],2)),
                       paste0(round(quantiles_sens_p[4],2),"-", round(quantiles_sens_p[5],2)),
                       paste0("> ",round(quantiles_sens_p[5],2)))
levels_reqd_sens_p = c("0.0-0.1","0.1-0.14","0.14-0.23","0.23-0.39","0.39-0.7", "> 0.7")
shp_df$Sensitivity_P = factor(shp_df$Sensitivity_P, levels=rev(levels_reqd_sens_p))

# color for the six bins (may be we need more bins)
colors_reqd = brewer.pal(length(unique(bins_sens_p))-1, "RdBu")

plot_cluster_sens_p = ggplot()+
  geom_polygon(data=shp_df, aes(x=long,y=lat,group=group,fill=Sensitivity_P),colour="black",size=0.1/2)+
  scale_fill_manual(breaks=as.character(levels(shp_df$Sensitivity_P)),values=colors_reqd,na.value="lightgrey",name="Closure Radius (Mean P)")+
  theme(axis.line = element_line(color="black",size=0.1/2))+
  theme(legend.position = "bottom")+
  theme(text = element_text(family="Times", size= 12))+
  theme(legend.text = element_text(size=10))+
  theme(legend.title = element_text(face="bold"))+
  theme(axis.text = element_text(size=12))+
  guides(fill=guide_legend(title="Cluster Radius (Mean P):", face="bold",nrow=1,byrow = TRUE))+
  theme(strip.background = element_blank())+
  theme(panel.background = element_rect(fill=alpha("lightgrey",0.3/2)))+
  theme(panel.grid.major = element_line(colour="black",size=0.1/2),
        panel.grid.minor = element_line(colour="black",size=0.1/2),
        panel.border = element_rect(colour="black",fill=NA))+
  labs(x = "Longitude", y="Latitude")
  #guides(fill=guide_legend(title="Cluster Radius (Mean P)"))
#ggsave(plot_cluster_sens_p,filename = "Graphs_For_Manuscript_Revised/Figure_5/Cluster_Sensitivity_Mean_P_noSM2RAIN_withBoundary.png", height=115, width=190,dpi=300, units="mm")

# Sensitivity of Cluster Radius (Mean P)
#levels_reqd_sens_et = c(paste0("0.0","-", round(quantiles_sens_et[1],2)),
#                      paste0(round(quantiles_sens_et[1],2),"-", round(quantiles_sens_et[2],2)),
#                        paste0(round(quantiles_sens_et[2],2),"-", round(quantiles_sens_et[3],2)),
#                        paste0(round(quantiles_sens_et[3],2),"-", round(quantiles_sens_et[4],2)),
#                        paste0(round(quantiles_sens_et[4],2),"-", round(quantiles_sens_et[5],2)),
#                        paste0("> ",round(quantiles_sens_et[5],2)))
levels_reqd_sens_et = c(paste0("0.0","-", round(quantiles_sens_p[1],2)),
                       paste0(round(quantiles_sens_p[1],2),"–", round(quantiles_sens_p[2],2)),
                       paste0(round(quantiles_sens_p[2],2),"-", round(quantiles_sens_p[3],2)),
                       paste0(round(quantiles_sens_p[3],2),"-", round(quantiles_sens_p[4],2)),
                       paste0(round(quantiles_sens_p[4],2),"-", round(quantiles_sens_p[5],2)),
                       paste0("> ",round(quantiles_sens_p[5],2)))
levels_reqd_sens_et = c("0.0-0.1","0.1-0.14","0.14-0.23","0.23-0.39","0.39-0.7", "> 0.7")
shp_df$Sensitivity_ET = factor(shp_df$Sensitivity_ET, levels=rev(levels_reqd_sens_et))

plot_cluster_sens_et = ggplot()+
  geom_polygon(data=shp_df, aes(x=long,y=lat,group=group,fill=Sensitivity_ET),colour="black",size=0.1/2)+
  scale_fill_manual(breaks=as.character(levels(shp_df$Sensitivity_ET)),values=colors_reqd,na.value="lightgrey",name="Closure Radius (Mean ET)")+
  theme(axis.line = element_line(color="black",size=0.1/2))+
  theme(legend.position = "bottom")+
  theme(text = element_text(family="Times", size= 12))+
  theme(legend.text = element_text(size=10))+
  theme(legend.title = element_text(face="bold"))+
  theme(axis.text = element_text(size=12))+
  guides(fill=guide_legend(title="Cluster Radius (Mean ET):", face="bold",nrow=1,byrow = TRUE))+
  theme(strip.background = element_blank())+
  theme(panel.background = element_rect(fill=alpha("lightgrey",0.3/2)))+
  theme(panel.grid.major = element_line(colour="black",size=0.1/2),
        panel.grid.minor = element_line(colour="black",size=0.1/2),
        panel.border = element_rect(colour="black",fill=NA))+
  labs(x = "Longitude", y="Latitude")
  #guides(fill=guide_legend(title="Cluster Radius (Mean ET)"))
#ggsave(plot_cluster_sens_et,filename = "Graphs_For_Manuscript_Revised/Figure_5/Cluster_Sensitivity_Mean_ET_noSM2RAIN_withBoundary.png", height=115, width=195,dpi=300, units="mm")

# Ratio of P to ET
levels_reqd_sens_p2et = c(paste0("0.0","-", "1.00"),
                        paste0("1.00","–", round(quantiles_sens_p2et[2],2)),
                        paste0(round(quantiles_sens_p2et[2],2),"-", round(quantiles_sens_p2et[3],2)),
                        paste0(round(quantiles_sens_p2et[3],2),"-", round(quantiles_sens_p2et[4],2)),
                        paste0(round(quantiles_sens_p2et[4],2),"-", round(quantiles_sens_p2et[5],2)),
                        paste0("> ",round(quantiles_sens_p2et[5],2)))
levels_reqd_sens_p2et = c("0.0-1.00","1.00-1.84","1.84-2.7","2.7-3.69","3.69-4.73", "> 4.73")
shp_df$Sensitivity_P2ET = factor(shp_df$Sensitivity_P2ET, levels=rev(levels_reqd_sens_p2et))
plot_cluster_sens_p2et = ggplot()+
  geom_polygon(data=shp_df, aes(x=long,y=lat,group=group,fill=Sensitivity_P2ET),colour="black",size=0.1/2)+
  scale_fill_manual(breaks=as.character(levels(shp_df$Sensitivity_P2ET)),values=colors_reqd,na.value="lightgrey",name="Closure Radius (Mean ET)")+
  theme(axis.line = element_line(color="black",size=0.1/2))+
  theme(legend.position = "bottom")+
  theme(text = element_text(family="Times", size= 12))+
  theme(legend.text = element_text(size=10))+
  theme(legend.title = element_text(face="bold"))+
  theme(axis.text = element_text(size=12))+
  guides(fill=guide_legend(title="CR (Mean P)/CR (Mean ET):", face="bold",nrow=1,byrow = TRUE))+
  theme(strip.background = element_blank())+
  theme(panel.background = element_rect(fill=alpha("lightgrey",0.3/2)))+
  theme(panel.grid.major = element_line(colour="black",size=0.1/2),
        panel.grid.minor = element_line(colour="black",size=0.1/2),
        panel.border = element_rect(colour="black",fill=NA))+
  labs(x = "Longitude", y="Latitude")
#guides(fill=guide_legend(title="Cluster Radius (Mean ET)"))
ggsave(plot_cluster_sens_p2et,filename = "Graphs_For_Manuscript_Revised/Figure_5/Cluster_Sensitivity_Mean_P2ET_noSM2RAIN_withBoundary.png", height=115, width=205,dpi=300, units="mm")


}

flag=2
if(flag > 0){
  
  # Mean Radius 
  levels_reqd_meanrad = c(paste0("0.0","-", round(quantiles_meanrad[1],1)),
                          paste0(round(quantiles_meanrad[1],1),"-", round(quantiles_meanrad[2],1)),
                          paste0(round(quantiles_meanrad[2],1),"-", round(quantiles_meanrad[3],1)),
                          paste0(round(quantiles_meanrad[3],1),"-", round(quantiles_meanrad[4],1)),
                          paste0(round(quantiles_meanrad[4],1),"-", round(quantiles_meanrad[5],1)),
                          paste0("> ",round(quantiles_meanrad[5],1)))
  shp_df$Cluster_Radius_Mean = factor(shp_df$Cluster_Radius_Mean, levels=rev(levels_reqd_meanrad))
  
  
  # color for the six bins (may be we need more bins)
  colors_reqd = brewer.pal(length(unique(bins_meanrad))-1, "RdBu")
  
  plot_cluster_mean_radius = ggplot()+
    geom_polygon(data=shp_df, aes(x=long,y=lat,group=group,fill=Cluster_Radius_Mean))+
    geom_polygon(data=shp_df, aes(x=long,y=lat,group=group),fill=NA,colour="black",size=0.1/2)+
    #  scale_fill_viridis()+
    #scale_fill_brewer(type="qual",palette = "Paired")+
    scale_fill_manual(breaks=as.character(levels(shp_df$Cluster_Radius_Mean)),values=colors_reqd,name="Closure Error (Mean)")+
    theme(axis.line = element_line(color="black",size=0.1/2))+
    theme(strip.background = element_blank())+
    theme(panel.background = element_rect(fill=alpha("lightgrey",0.3/2)))+
    theme(panel.grid.major = element_line(colour="black",size=0.1/2),
          panel.grid.minor = element_line(colour="black",size=0.1/2),
          panel.border = element_rect(colour="black",fill=NA))+
    theme(legend.position = "bottom")+
    theme(text = element_text(family="Times", size= 12))+
    theme(legend.text = element_text(size=10))+
    theme(legend.title = element_text(face="bold"))+
    theme(axis.text = element_text(size=12))+
    guides(fill=guide_legend(title="Closure Error (Mean):", face="bold",nrow=1,byrow = TRUE))+
    labs(x = "Longitude", y="Latitude")
  #guides(fill=guide_legend(title="Budyko Error (Mean)"))
  #ggsave(plot_cluster_mean_radius,filename = "Graphs_For_Manuscript_Revised/Figure_3/Dist_Metric_Mean_Normalized_noSM2RAIN_withBoundary.png", height=115, width=195,dpi=300, units="mm")
  
  # Standard Deviation of Radius
  levels_reqd_stdrad = c(paste0("0.0"," - ", round(quantiles_stdrad[1],3)),
                         paste0(round(quantiles_stdrad[1],2),"-", round(quantiles_stdrad[2],2)),
                         paste0(round(quantiles_stdrad[2],2),"-", round(quantiles_stdrad[3],2)),
                         paste0(round(quantiles_stdrad[3],2),"-", round(quantiles_stdrad[4],2)),
                         paste0(round(quantiles_stdrad[4],2),"-", round(quantiles_stdrad[5],2)),
                         paste0("> ",round(quantiles_stdrad[5],2)))
  shp_df$Cluster_Radius_SD = factor(shp_df$Cluster_Radius_SD, levels=rev(levels_reqd_stdrad))
  
  plot_cluster_sd_radius = ggplot()+
    geom_polygon(data=shp_df, aes(x=long,y=lat,group=group,fill=Cluster_Radius_SD),colour="black",size=0.1/2)+
    scale_fill_manual(breaks=as.character(levels(shp_df$Cluster_Radius_SD)),values=colors_reqd,na.value="lightgrey",name="Closure Error (SD)")+
    theme(axis.line = element_line(color="black",size=0.1/2))+
    theme(strip.background = element_blank())+
    theme(panel.background = element_rect(fill=alpha("lightgrey",0.3/2)))+
    theme(panel.grid.major = element_line(colour="black",size=0.1/2),
          panel.grid.minor = element_line(colour="black",size=0.1/2),
          panel.border = element_rect(colour="black",fill=NA))+
    theme(legend.position = "bottom")+
    theme(text = element_text(family="Times", size= 12))+
    theme(legend.text = element_text(size=10))+
    theme(legend.title = element_text(face="bold"))+
    theme(axis.text = element_text(size=12))+
    guides(fill=guide_legend(title="Closure Error (SD):", face="bold",nrow=1,byrow = TRUE))+
    labs(x = "Longitude", y="Latitude")
# ggsave(plot_cluster_sd_radius,filename = "Graphs_For_Manuscript_Revised/Figure_3/Dist_Metric_SD_Normalized_noSM2RAIN_withBoundary.png", height=115, width=190,dpi=300, units = "mm")
  
  # Min Radius
  levels_reqd_minrad = c(paste0("0.0"," - ", round(quantiles_minrad[1],2)),
                         paste0(round(quantiles_minrad[1],2),"-", round(quantiles_minrad[2],2)),
                         paste0(round(quantiles_minrad[2],2),"-", round(quantiles_minrad[3],2)),
                         paste0(round(quantiles_minrad[3],2),"-", round(quantiles_minrad[4],2)),
                         paste0(round(quantiles_minrad[4],2),"-", round(quantiles_minrad[5],2)),
                         paste0("> ",round(quantiles_minrad[5],2)))
  levels_reqd_minrad = c("0.02-0.04", "0.04-0.1","0.1-0.18","0.18-0.39","> 0.39")
  shp_df$Cluster_Radius_Min = factor(shp_df$Cluster_Radius_Min, levels=rev(levels_reqd_minrad))
  
  plot_cluster_min_radius = ggplot()+
    geom_polygon(data=shp_df, aes(x=long,y=lat,group=group,fill=Cluster_Radius_Min),colour="black",size=0.1/2)+
    scale_fill_manual(breaks=as.character(levels(shp_df$Cluster_Radius_Min)),values=colors_reqd[-5],na.value="lightgrey",name="Closure Error (Min)")+
    theme(axis.line = element_line(color="black",size=0.1/2))+
    theme(strip.background = element_blank())+
    theme(panel.background = element_rect(fill=alpha("lightgrey",0.3/2)))+
    theme(panel.grid.major = element_line(colour="black",size=0.1/2),
          panel.grid.minor = element_line(colour="black",size=0.1/2),
          panel.border = element_rect(colour="black",fill=NA))+
    theme(legend.position = "bottom")+
    theme(text = element_text(family="Times", size= 12))+
    theme(legend.text = element_text(size=10))+
    theme(legend.title = element_text(face="bold"))+
    theme(axis.text = element_text(size=12))+
    guides(fill=guide_legend(title="Closure Error (Min):", face="bold",nrow=1,byrow = TRUE))+
    labs(x = "Longitude", y="Latitude")
  ggsave(plot_cluster_min_radius,filename = "Graphs_For_Manuscript_Revised/Supplementary/Dist_Metric_Minimum_Normalized_noSM2RAIN_withBoundary.pdf", height=115, width=190,dpi=300, units="mm")

}


flag=0
if(flag > 0){
  # Max Radius
  levels_reqd_maxrad = c(paste0("0.0"," - ", round(quantiles_maxrad[1],1)),
                         paste0(round(quantiles_maxrad[1],1)," - ", round(quantiles_maxrad[2],1)),
                         paste0(round(quantiles_maxrad[2],1)," - ", round(quantiles_maxrad[3],1)),
                         paste0(round(quantiles_maxrad[3],1)," - ", round(quantiles_maxrad[4],1)),
                         paste0(round(quantiles_maxrad[4],1)," - ", round(quantiles_maxrad[5],1)),
                         paste0("> ",round(quantiles_maxrad[5],1)))
  shp_df$Cluster_Radius_Max = factor(shp_df$Cluster_Radius_Max, levels=levels_reqd_maxrad)
  
  plot_cluster_max_radius = ggplot()+
    geom_polygon(data=shp_df, aes(x=long,y=lat,group=group,fill=factor(Cluster_Radius_Max)))+
    #scale_fill_viridis()+
    scale_fill_brewer(type="qual",palette = "Paired")+
    theme(axis.line = element_line(color="black"))+
    theme(strip.background = element_blank())+
    theme(panel.grid.major = element_line(colour="black"),
          panel.grid.minor = element_line(colour="grey"),
          panel.border = element_rect(colour="black",fill=NA))+
    theme(panel.background = element_rect(fill=alpha("#C8AD7F",0.5)))+
    guides(fill=guide_legend(title="Budyko Error (Maximum)"))
  #ggsave(plot_cluster_max_radius,filename = "Output/Dist_Metric_Maximum_Normalized.png", height=6, width=10,dpi=300)
}


