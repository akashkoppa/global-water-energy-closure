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

# read in the clustering metric results
cluster_metric = read.table("Output/Budyko_Cluster_Metrics_RS_noSM2RAIN.txt",header=TRUE)

# subset the problem catchments in Amazonia
#catch_reqd = c(6050207670, 6050156880, 6050207660, 6050287830, 6050266740, 6050266880, 6050294360, 6050294270, 6050247260, 6050247270)
#index_reqd = which(cluster_metric$HydroSHEDS_ID %in% catch_reqd)

# subset ei and ai data for the Amazon region
#cluster_metric_reqd = cluster_metric[index_reqd,]

# subset a few catchments with higher uncertainty in Amazonia
#catch_reqd_high = c(6050315910, 6050315920, 6050029730)
#index_reqd_high = which(cluster_metric$HydroSHEDS_ID %in% catch_reqd_high)
#cluster_metric_reqd_high = cluster_metric[index_reqd_high,]

# separate the area and AI data
hydrosheds_id = as.character(cluster_metric$HydroSHEDS_ID)
area = cluster_metric$Area
ai = cluster_metric$AI
cluster_metric = cluster_metric[,-c(1:3)]

# normalize the cluster metrics by the aridity index
cluster_metric = cluster_metric/ai

# replace radius metric > 100 as NA
#cluster_metric$Max_Radius[which(cluster_metric$Max_Radius > 10)] = NA
#cluster_metric$Mean_Radius[which(cluster_metric$Mean_Radius > 10)] = NA
#cluster_metric$Min_Radius[which(cluster_metric$Min_Radius > 10)] = NA
#cluster_metric$Sd_Radius[which(cluster_metric$Sd_Radius > 10)] = NA
#cluster_metric$Density_Max_Radius[which(cluster_metric$Density_Max_Radius > 500)] = NA
#cluster_metric$Density_Mean_Radius[which(cluster_metric$Density_Mean_Radius > 500)] = NA

#cluster_metric = apply(cluster_metric, MARGIN = 2, FUN = function(x)(x - min(x, na.rm=TRUE))/diff(range(x,na.rm=TRUE)))
#cluster_metric = data.frame(cluster_metric)

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
bins_stdrad[which(cluster_metric$Sd_Radius <= quantiles_stdrad[1])]                                                     = paste0("0.0","-", round(quantiles_stdrad[1],2))
bins_stdrad[which(cluster_metric$Sd_Radius <= quantiles_stdrad[2] & cluster_metric$Sd_Radius > quantiles_stdrad[1])] = paste0(round(quantiles_stdrad[1],2),"-", round(quantiles_stdrad[2],2))
bins_stdrad[which(cluster_metric$Sd_Radius <= quantiles_stdrad[3] & cluster_metric$Sd_Radius > quantiles_stdrad[2])] = paste0(round(quantiles_stdrad[2],2),"-", round(quantiles_stdrad[3],2))
bins_stdrad[which(cluster_metric$Sd_Radius <= quantiles_stdrad[4] & cluster_metric$Sd_Radius > quantiles_stdrad[3])] = paste0(round(quantiles_stdrad[3],2),"-", round(quantiles_stdrad[4],2))
bins_stdrad[which(cluster_metric$Sd_Radius <= quantiles_stdrad[5] & cluster_metric$Sd_Radius > quantiles_stdrad[4])] = paste0(round(quantiles_stdrad[4],2),"-", round(quantiles_stdrad[5],2))
bins_stdrad[which(cluster_metric$Sd_Radius >  quantiles_stdrad[5])]                                                     = paste0("> ",round(quantiles_stdrad[5],2))

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

# Min Radius
quantiles_minrad = quantile(cluster_metric$Min_Radius, probs=c(0.1,0.25,0.5,0.75,0.9),na.rm=TRUE)
bins_minrad = rep(NA,nrow(cluster_metric))
bins_minrad[which(cluster_metric$Min_Radius <= quantiles_minrad[1])]                                                     = paste0("0.0"," - ", round(quantiles_minrad[1],3))
bins_minrad[which(cluster_metric$Min_Radius <= quantiles_minrad[2] & cluster_metric$Min_Radius > quantiles_minrad[1])] = paste0(round(quantiles_minrad[1],3)," - ", round(quantiles_minrad[2],3))
bins_minrad[which(cluster_metric$Min_Radius <= quantiles_minrad[3] & cluster_metric$Min_Radius > quantiles_minrad[2])] = paste0(round(quantiles_minrad[2],3)," - ", round(quantiles_minrad[3],3))
bins_minrad[which(cluster_metric$Min_Radius <= quantiles_minrad[4] & cluster_metric$Min_Radius > quantiles_minrad[3])] = paste0(round(quantiles_minrad[3],3)," - ", round(quantiles_minrad[4],3))
bins_minrad[which(cluster_metric$Min_Radius <= quantiles_minrad[5] & cluster_metric$Min_Radius > quantiles_minrad[4])] = paste0(round(quantiles_minrad[4],3)," - ", round(quantiles_minrad[5],3))
bins_minrad[which(cluster_metric$Min_Radius >  quantiles_minrad[5])]                                                     = paste0("> ",round(quantiles_minrad[5],3))

# Cluster Density (Based on Max Radius)
quantiles_dmaxrad = quantile(cluster_metric$Density_Max_Radius, probs=c(0.1,0.25,0.5,0.75,0.9),na.rm=TRUE)
bins_dmaxrad = rep(NA,nrow(cluster_metric))
bins_dmaxrad[which(cluster_metric$Density_Max_Radius <= quantiles_dmaxrad[1])]                                                     = paste0("0.0"," - ", round(quantiles_dmaxrad[1],1))
bins_dmaxrad[which(cluster_metric$Density_Max_Radius <= quantiles_dmaxrad[2] & cluster_metric$Density_Max_Radius > quantiles_dmaxrad[1])] = paste0(round(quantiles_dmaxrad[1],1)," - ", round(quantiles_dmaxrad[2],1))
bins_dmaxrad[which(cluster_metric$Density_Max_Radius <= quantiles_dmaxrad[3] & cluster_metric$Density_Max_Radius > quantiles_dmaxrad[2])] = paste0(round(quantiles_dmaxrad[2],1)," - ", round(quantiles_dmaxrad[3],1))
bins_dmaxrad[which(cluster_metric$Density_Max_Radius <= quantiles_dmaxrad[4] & cluster_metric$Density_Max_Radius > quantiles_dmaxrad[3])] = paste0(round(quantiles_dmaxrad[3],1)," - ", round(quantiles_dmaxrad[4],1))
bins_dmaxrad[which(cluster_metric$Density_Max_Radius <= quantiles_dmaxrad[5] & cluster_metric$Density_Max_Radius > quantiles_dmaxrad[4])] = paste0(round(quantiles_dmaxrad[4],1)," - ", round(quantiles_dmaxrad[5],1))
bins_dmaxrad[which(cluster_metric$Density_Max_Radius >  quantiles_dmaxrad[5])]                                                     = paste0("> ",round(quantiles_dmaxrad[5],1))

# Cluster Density (Based on Mean Radius)
quantiles_dmeanrad = quantile(cluster_metric$Density_Mean_Radius, probs=c(0.1,0.25,0.5,0.75,0.9),na.rm=TRUE)
bins_dmeanrad = rep(NA,nrow(cluster_metric))
bins_dmeanrad[which(cluster_metric$Density_Mean_Radius <= quantiles_dmeanrad[1])]                                                     = paste0("0.0"," - ", round(quantiles_dmeanrad[1],1))
bins_dmeanrad[which(cluster_metric$Density_Mean_Radius <= quantiles_dmeanrad[2] & cluster_metric$Density_Mean_Radius > quantiles_dmeanrad[1])] = paste0(round(quantiles_dmeanrad[1],1)," - ", round(quantiles_dmeanrad[2],1))
bins_dmeanrad[which(cluster_metric$Density_Mean_Radius <= quantiles_dmeanrad[3] & cluster_metric$Density_Mean_Radius > quantiles_dmeanrad[2])] = paste0(round(quantiles_dmeanrad[2],1)," - ", round(quantiles_dmeanrad[3],1))
bins_dmeanrad[which(cluster_metric$Density_Mean_Radius <= quantiles_dmeanrad[4] & cluster_metric$Density_Mean_Radius > quantiles_dmeanrad[3])] = paste0(round(quantiles_dmeanrad[3],1)," - ", round(quantiles_dmeanrad[4],1))
bins_dmeanrad[which(cluster_metric$Density_Mean_Radius <= quantiles_dmeanrad[5] & cluster_metric$Density_Mean_Radius > quantiles_dmeanrad[4])] = paste0(round(quantiles_dmeanrad[4],1)," - ", round(quantiles_dmeanrad[5],1))
bins_dmeanrad[which(cluster_metric$Density_Mean_Radius >  quantiles_dmeanrad[5])]                                                    = paste0("> ",round(quantiles_dmeanrad[5],1))
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
#shp$Cluster_Radius_Min  = bins_minrad
shp$Cluster_Radius_SD   = bins_stdrad
#shp$Density_Max_Radius  = bins_dmaxrad
#shp$Density_Mean_Radius = bins_dmeanrad

# try sequential colors once more
#cluster_mean_radius = cluster_metric$Mean_Radius
#cluster_mean_radius[which(cluster_mean_radius > 1)] = NA
#shp$Cluster_Radius_Mean = cluster_mean_radius
shp@data$id = rownames(shp@data)
shp_df = fortify(shp)
shp_df = join(shp_df,shp@data,by="id")

# Mean Radius 
levels_reqd_meanrad = c(paste0("0.0","-", round(quantiles_meanrad[1],1)),
                       paste0(round(quantiles_meanrad[1],1),"-", round(quantiles_meanrad[2],1)),
                       paste0(round(quantiles_meanrad[2],1),"-", round(quantiles_meanrad[3],1)),
                       paste0(round(quantiles_meanrad[3],1),"-", round(quantiles_meanrad[4],1)),
                       paste0(round(quantiles_meanrad[4],1),"-", round(quantiles_meanrad[5],1)),
                       paste0("> ",round(quantiles_meanrad[5],1)))
shp_df$Cluster_Radius_Mean = factor(shp_df$Cluster_Radius_Mean, levels=levels_reqd_meanrad)

# centroid coordinates of the basins
#shp_centroid = as.data.frame(coordinates(shp))
#hybas = shp@data$HYBAS_ID
#shp_centroid = cbind(shp_centroid, hybas)
#colnames(shp_centroid) = c("long","lat","hybas_id")

# reorder the factors 
shp_df$Cluster_Radius_Mean = factor(shp_df$Cluster_Radius_Mean, levels=rev(levels(shp_df$Cluster_Radius_Mean)))

# color for the six bins (may be we need more bins)
#colors_reqd = rev(brewer.pal(length(unique(bins_meanrad))-1, "RdBu"))
colors_reqd = brewer.pal(length(unique(bins_meanrad))-1, "RdBu")

plot_cluster_mean_radius = ggplot()+
  geom_polygon(data=shp_df, aes(x=long,y=lat,group = group,fill=Cluster_Radius_Mean))+
  geom_polygon(data=shp_df, aes(x=long,y=lat,group=group),fill=NA,colour="black",size=0.1/2)+
  #geom_text(data=shp_centroid,aes(label=hybas_id, x=long, y=lat),size=2.5)+
  #coord_cartesian(xlim=c(-80, -50), ylim = c(-20, 20))+
  #coord_cartesian(xlim=c(5, 40), ylim = c(55, 65))+
  #scale_fill_viridis_c()+
  #scale_fill_brewer(type="qual",palette = "Paired")+
  #scale_fill_distiller(type = "div", palette ="RdBu")
  #scale_fill_gradient2(low="gold",mid="white",high='firebrick2',na.value="lightgrey")+
  scale_fill_manual(breaks=as.character(levels(shp_df$Cluster_Radius_Mean)),values=colors_reqd,na.value="lightgrey",name="Cluster Radius (Mean)")+
  theme(axis.line = element_line(color="black",size=0.1/2))+
  theme(legend.position = "bottom")+
  theme(text = element_text(family="Times", size= 12))+
  theme(legend.text = element_text(size=10))+
  theme(legend.title = element_text(face="bold"))+
  theme(axis.text = element_text(size=12))+
  guides(fill=guide_legend(title="Cluster Radius (Mean):", face="bold",nrow=1,byrow = TRUE))+
  theme(strip.background = element_blank())+
  theme(panel.background = element_rect(fill=alpha("lightgrey",0.3/2)))+
  theme(panel.grid.major = element_line(colour="black",size=0.1/2),
        panel.grid.minor = element_line(colour="black",size=0.1/2),
        panel.border = element_rect(colour="black",fill=NA))+
  labs(x = "Longitude", y="Latitude")
  #guides(fill=guide_legend(title="Cluster Radius (Mean)"))
ggsave(plot_cluster_mean_radius,filename = "Graphs_For_Manuscript_Revised/Figure_3/Cluster_Radius_Mean_Normalized_noSM2RAIN_withBoundary.png", height=115, width=190,dpi=300,units="mm")
#ggsave(plot_cluster_mean_radius,filename = "Output/Cluster/test.png", height=5, width=10,dpi=300)

# Standard Deviation of Radius
levels_reqd_stdrad = c(paste0("0.0","-", round(quantiles_stdrad[1],2)),
                       paste0(round(quantiles_stdrad[1],2),"-", round(quantiles_stdrad[2],2)),
                       paste0(round(quantiles_stdrad[2],2),"-", round(quantiles_stdrad[3],2)),
                       paste0(round(quantiles_stdrad[3],2),"-", round(quantiles_stdrad[4],2)),
                       paste0(round(quantiles_stdrad[4],2),"-", round(quantiles_stdrad[5],2)),
                       paste0("> ",round(quantiles_stdrad[5],2)))
shp_df$Cluster_Radius_SD = factor(shp_df$Cluster_Radius_SD, levels=levels_reqd_stdrad)
shp_df$Cluster_Radius_SD = factor(shp_df$Cluster_Radius_SD, levels=rev(levels(shp_df$Cluster_Radius_SD)))

plot_cluster_sd_radius = ggplot()+
  geom_polygon(data=shp_df, aes(x=long,y=lat,group=group,fill=Cluster_Radius_SD),colour = "black",size=0.1/2)+
  scale_fill_manual(breaks=as.character(levels(shp_df$Cluster_Radius_SD)),values=colors_reqd,na.value="lightgrey",name="Cluster Radius (Mean)")+
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
  guides(fill=guide_legend(title="Cluster Radius (SD):", face="bold",nrow=1,byrow = TRUE))+
  labs(x = "Longitude", y="Latitude")
  #guides(fill=guide_legend(title="Cluster Radius (Standard Deviation)"))
ggsave(plot_cluster_sd_radius,filename = "Graphs_For_Manuscript_Revised/Figure_3/Cluster_Radius_SD_Normalized_noSM2RAIN_withBoundary.png", height=115, width=190,dpi=300,units="mm")

flag=0
if (flag>1){
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
  guides(fill=guide_legend(title="Cluster Radius (Maximum)"))
#ggsave(plot_cluster_max_radius,filename = "Output/Cluster/Cluster_Radius_Maximum_Normalized.png", height=6, width=10,dpi=300)
#ggsave(plot_cluster_max_radius,filename = "test.png", height=6, width=10,dpi=300)

# Min Radius
levels_reqd_minrad = c(paste0("0.0"," - ", round(quantiles_minrad[1],1)),
                       paste0(round(quantiles_minrad[1],3)," - ", round(quantiles_minrad[2],3)),
                       paste0(round(quantiles_minrad[2],3)," - ", round(quantiles_minrad[3],3)),
                       paste0(round(quantiles_minrad[3],3)," - ", round(quantiles_minrad[4],3)),
                       paste0(round(quantiles_minrad[4],3)," - ", round(quantiles_minrad[5],3)),
                       paste0("> ",round(quantiles_minrad[5],3)))
shp_df$Cluster_Radius_Min = factor(shp_df$Cluster_Radius_Min, levels=levels_reqd_minrad)

plot_cluster_min_radius = ggplot()+
  geom_polygon(data=shp_df, aes(x=long,y=lat,group=group,fill=Cluster_Radius_Min))+
#  scale_fill_viridis()+
scale_fill_brewer(type="qual",palette = "Paired")+
  theme(axis.line = element_line(color="black"))+
  theme(strip.background = element_blank())+
  theme(panel.background = element_rect(fill=alpha("#C8AD7F",0.5)))+
  theme(panel.grid.major = element_line(colour="black"),
        panel.grid.minor = element_line(colour="grey"),
        panel.border = element_rect(colour="black",fill=NA))+
  guides(fill=guide_legend(title="Cluster Radius (Minimum)"))
ggsave(plot_cluster_min_radius,filename = "Output/Cluster/Cluster_Radius_Minimum_Normalized.png", height=6, width=10,dpi=300)

# Standard Deviation of Radius
levels_reqd_stdrad = c(paste0("0.0"," - ", round(quantiles_stdrad[1],3)),
                       paste0(round(quantiles_stdrad[1],3)," - ", round(quantiles_stdrad[2],3)),
                       paste0(round(quantiles_stdrad[2],3)," - ", round(quantiles_stdrad[3],3)),
                       paste0(round(quantiles_stdrad[3],3)," - ", round(quantiles_stdrad[4],3)),
                       paste0(round(quantiles_stdrad[4],3)," - ", round(quantiles_stdrad[5],3)),
                       paste0("> ",round(quantiles_stdrad[5],3)))
shp_df$Cluster_Radius_SD = factor(shp_df$Cluster_Radius_SD, levels=levels_reqd_stdrad)

plot_cluster_sd_radius = ggplot()+
  geom_polygon(data=shp_df, aes(x=long,y=lat,group=group,fill=Cluster_Radius_SD))+
#  scale_fill_viridis()+
scale_fill_brewer(type="qual",palette = "Paired")+
  theme(axis.line = element_line(color="black"))+
  theme(strip.background = element_blank())+
  theme(panel.background = element_rect(fill=alpha("#C8AD7F",0.5)))+
  theme(panel.grid.major = element_line(colour="black"),
        panel.grid.minor = element_line(colour="grey"),
      panel.border = element_rect(colour="black",fill=NA))+
  guides(fill=guide_legend(title="Cluster Radius (Standard Deviation)"))
ggsave(plot_cluster_sd_radius,filename = "Output/Cluster/Cluster_Radius_SD_Normalized.png", height=6, width=10,dpi=300)

# Density (max radius)
levels_reqd_dmaxrad = c(paste0("0.0"," - ", round(quantiles_dmaxrad[1],1)),
                       paste0(round(quantiles_dmaxrad[1],1)," - ", round(quantiles_dmaxrad[2],1)),
                       paste0(round(quantiles_dmaxrad[2],1)," - ", round(quantiles_dmaxrad[3],1)),
                       paste0(round(quantiles_dmaxrad[3],1)," - ", round(quantiles_dmaxrad[4],1)),
                       paste0(round(quantiles_dmaxrad[4],1)," - ", round(quantiles_dmaxrad[5],1)),
                       paste0("> ",round(quantiles_dmaxrad[5],1)))
shp_df$Density_Max_Radius = factor(shp_df$Density_Max_Radius, levels=levels_reqd_dmaxrad)

plot_cluster_density_max_radius = ggplot()+
  geom_polygon(data=shp_df, aes(x=long,y=lat,group=group,fill=Density_Max_Radius))+
#  scale_fill_viridis()+
scale_fill_brewer(type="qual",palette = "Paired")+
  theme(axis.line = element_line(color="black"))+
  theme(strip.background = element_blank())+
  theme(panel.background = element_rect(fill=alpha("#C8AD7F",0.5)))+
  theme(panel.grid.major = element_line(colour="black"),
        panel.grid.minor = element_line(colour="grey"),
        panel.border = element_rect(colour="black",fill=NA))+
  guides(fill=guide_legend(title="Density (Max Radius)"))
ggsave(plot_cluster_density_max_radius,filename = "Output/Cluster/Cluster_Density_MaxRadius_Normalized.png", height=6, width=10,dpi=300)

# Density (mean radius)
levels_reqd_dmeanrad = c(paste0("0.0"," - ", round(quantiles_dmeanrad[1],1)),
                        paste0(round(quantiles_dmeanrad[1],1)," - ", round(quantiles_dmeanrad[2],1)),
                        paste0(round(quantiles_dmeanrad[2],1)," - ", round(quantiles_dmeanrad[3],1)),
                        paste0(round(quantiles_dmeanrad[3],1)," - ", round(quantiles_dmeanrad[4],1)),
                        paste0(round(quantiles_dmeanrad[4],1)," - ", round(quantiles_dmeanrad[5],1)),
                        paste0("> ",round(quantiles_dmeanrad[5],1)))
shp_df$Density_Mean_Radius = factor(shp_df$Density_Mean_Radius, levels=levels_reqd_dmeanrad)

plot_cluster_density_mean_radius = ggplot()+
  geom_polygon(data=shp_df, aes(x=long,y=lat,group=group,fill=Density_Mean_Radius))+
  #scale_fill_viridis()+
  scale_fill_brewer(type="qual",palette = "Paired")+
  theme(axis.line = element_line(color="black"))+
  theme(strip.background = element_blank())+
  theme(panel.background = element_rect(fill=alpha("#C8AD7F",0.5)))+
  theme(panel.grid.major = element_line(colour="black"),
        panel.grid.minor = element_line(colour="grey"),
        panel.border = element_rect(colour="black",fill=NA))+
  guides(fill=guide_legend(title="Density (Mean Radius)"))
ggsave(plot_cluster_density_mean_radius,filename = "Output/Cluster/Cluster_Density_MeanRadius_Normalized.png", height=6, width=10,dpi=300)
}
















