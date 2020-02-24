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

# read in the cluster metric results
cluster_metric = read.table("Output/Budyko_Cluster_Metrics_RS_noSM2RAIN.txt",header=TRUE)

# separate the area and AI data
hydrosheds_id = as.character(dist_metric$HydroSHEDS_ID)
area = dist_metric$Area
ai = dist_metric$AI
dist_metric = dist_metric[,-c(1:3)]

# normalize the cluster metrics by the aridity index
dist_metric = dist_metric/ai
cluster_metric = cluster_metric/ai
cluster_metric = cluster_metric[,-c(1:3)]

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
dist_temp_metric = matrix(data = NA, nrow=nrow(dist_metric),ncol = 4)
dist_temp_metric = data.frame(dist_temp_metric)
colnames(dist_temp_metric) = c("Max_Radius","Mean_Radius","Min_Radius","Sd_Radius")
dist_temp_metric$Max_Radius  = apply(dist_metric,1,max,na.rm=TRUE)
dist_temp_metric$Max_Radius[which(dist_temp_metric$Max_Radius==Inf)] = NA
dist_temp_metric$Max_Radius[which(dist_temp_metric$Max_Radius==-Inf)] = NA
dist_temp_metric$Mean_Radius = apply(dist_metric,1,mean,na.rm=TRUE)
dist_temp_metric$Min_Radius  = apply(dist_metric,1,min,na.rm=TRUE)
dist_temp_metric$Min_Radius[which(dist_temp_metric$Min_Radius==Inf)] = NA
dist_temp_metric$Min_Radius[which(dist_temp_metric$Min_Radius==-Inf)] = NA
dist_temp_metric$Sd_Radius   = apply(dist_metric,1,sd,na.rm=TRUE)

# create different classes for each of the metric (maybe five or six classes each) and then assign colors to them.
# I think that will create a better map than what I have now. 
# The bins are created based on quantiles
bins_reqd = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)
# Mean Radius
ratio_CR_CE = (cluster_metric$Mean_Radius/(cluster_metric$Mean_Radius + dist_temp_metric$Mean_Radius))
#quantiles_cr2ce = quantile(ratio_CR_CE, probs=c(0.1,0.25,0.5,0.75,0.9),na.rm=TRUE)
bins_cr2ce = rep(NA,nrow(cluster_metric))
#bins_cr2ce[which(ratio_CR_CE <= quantiles_cr2ce[1])]                                                     = paste0("0.0","-", round(quantiles_cr2ce[1],2))
#bins_cr2ce[which(ratio_CR_CE <= quantiles_cr2ce[2] & ratio_CR_CE > quantiles_cr2ce[1])] = paste0(round(quantiles_cr2ce[1],2),"-", round(quantiles_cr2ce[2],2))
#bins_cr2ce[which(ratio_CR_CE <= quantiles_cr2ce[3] & ratio_CR_CE > quantiles_cr2ce[2])] = paste0(round(quantiles_cr2ce[2],2),"-", round(quantiles_cr2ce[3],2))
#bins_cr2ce[which(ratio_CR_CE <= quantiles_cr2ce[4] & ratio_CR_CE > quantiles_cr2ce[3])] = paste0(round(quantiles_cr2ce[3],2),"-", round(quantiles_cr2ce[4],2))
#bins_cr2ce[which(ratio_CR_CE <= quantiles_cr2ce[5] & ratio_CR_CE > quantiles_cr2ce[4])] = paste0(round(quantiles_cr2ce[4],2),"-", round(quantiles_cr2ce[5],2))
#bins_cr2ce[which(ratio_CR_CE >  quantiles_cr2ce[5])]                                                     = paste0("> ",round(quantiles_cr2ce[5],2))

bins_cr2ce[which(ratio_CR_CE <= bins_reqd[2] & ratio_CR_CE > bins_reqd[1])] = "0.0-0.2"
bins_cr2ce[which(ratio_CR_CE <= bins_reqd[3] & ratio_CR_CE > bins_reqd[2])] = "0.2-0.4"
bins_cr2ce[which(ratio_CR_CE <= bins_reqd[4] & ratio_CR_CE > bins_reqd[3])] = "0.4-0.6"
bins_cr2ce[which(ratio_CR_CE <= bins_reqd[5] & ratio_CR_CE > bins_reqd[4])] = "0.6-0.8"
bins_cr2ce[which(ratio_CR_CE <= bins_reqd[6] & ratio_CR_CE > bins_reqd[5])] = "0.8-1.0"


# ratio of mean P to mean ET
ratio_p2e = data_cluster_sens$P/(data_cluster_sens$ET + data_cluster_sens$P)
#quantiles_sens_p2et = quantile(ratio_p2e, probs=c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm=TRUE)
bins_sens_p2et = rep(NA,nrow(data_cluster_sens))
#bins_sens_p2et[which(ratio_p2e <= 1)]                                                     = paste0("0.0","-", "1.00")
#bins_sens_p2et[which(ratio_p2e <= quantiles_sens_p2et[2] & ratio_p2e > 1)] = paste0("1.00","-", round(quantiles_sens_p2et[2],2))
#bins_sens_p2et[which(ratio_p2e <= quantiles_sens_p2et[3] & ratio_p2e > quantiles_sens_p2et[2])] = paste0(round(quantiles_sens_p2et[2],2),"-", round(quantiles_sens_p2et[3],2))
#bins_sens_p2et[which(ratio_p2e <= quantiles_sens_p2et[4] & ratio_p2e > quantiles_sens_p2et[3])] = paste0(round(quantiles_sens_p2et[3],2),"-", round(quantiles_sens_p2et[4],2))
#bins_sens_p2et[which(ratio_p2e <= quantiles_sens_p2et[5] & ratio_p2e > quantiles_sens_p2et[4])] = paste0(round(quantiles_sens_p2et[4],2),"-", round(quantiles_sens_p2et[5],2))
#bins_sens_p2et[which(ratio_p2e >  quantiles_sens_p2et[5])]                                                     = paste0("> ",round(quantiles_sens_p2et[5],2))

bins_sens_p2et[which(ratio_p2e <= bins_reqd[2] & ratio_p2e > bins_reqd[1])] = "0.0-0.2"
bins_sens_p2et[which(ratio_p2e <= bins_reqd[3] & ratio_p2e > bins_reqd[2])] = "0.2-0.4"
bins_sens_p2et[which(ratio_p2e <= bins_reqd[4] & ratio_p2e > bins_reqd[3])] = "0.4-0.6"
bins_sens_p2et[which(ratio_p2e <= bins_reqd[5] & ratio_p2e > bins_reqd[4])] = "0.6-0.8"
bins_sens_p2et[which(ratio_p2e <= bins_reqd[6] & ratio_p2e > bins_reqd[5])] = "0.8-1.0"

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
shp$CR2CE = bins_cr2ce
shp$Sensitivity_P2ET    = bins_sens_p2et
#shp$CR2CE1 = ratio_CR_CE

shp@data$id = rownames(shp@data)
shp_df = fortify(shp)
shp_df = join(shp_df,shp@data,by="id")

flag = 1
if (flag>0){

  
# color for the six bins (may be we need more bins)
colors_reqd = brewer.pal(5, "RdBu")
colors_reqd = viridis(5)

# Ration of CR to CE
levels_reqd_cr2ce = c("0.0-0.2","0.2-0.4","0.4-0.6","0.6-0.8")
shp_df$CR2CE = factor(shp_df$CR2CE, levels=rev(levels_reqd_cr2ce))

plot_cr2ce = ggplot()+
  geom_polygon(data=shp_df, aes(x=long,y=lat,group=group,fill=CR2CE),colour="black",size=0.1/2)+
  scale_fill_manual(breaks=as.character(levels(shp_df$CR2CE)),values=colors_reqd[c(2:5)],na.value="lightgrey",name="Closure Radius (Mean ET)")+
  #scale_fill_viridis()+
  theme(axis.line = element_line(color="black",size=0.1/2))+
  theme(legend.position = "bottom")+
  theme(text = element_text(family="Times", size= 12))+
  theme(legend.text = element_text(size=10))+
  theme(legend.title = element_text(face="bold"))+
  theme(axis.text = element_text(size=12))+
  guides(fill=guide_legend(title="R(CR, CE):", face="bold",nrow=1,byrow = TRUE))+
  theme(strip.background = element_blank())+
  theme(panel.background = element_rect(fill=alpha("lightgrey",0.3/2)))+
  theme(panel.grid.major = element_line(colour="black",size=0.1/2),
        panel.grid.minor = element_line(colour="black",size=0.1/2),
        panel.border = element_rect(colour="black",fill=NA))+
  labs(x = "Longitude", y="Latitude")
#guides(fill=guide_legend(title="Cluster Radius (Mean ET)"))
ggsave(plot_cr2ce,filename = "Graphs_For_Manuscript_Revised/Figure_8/Cluster_CR2CE_noSM2RAIN_withBoundary.png", height=115, width=205,dpi=300, units="mm")

# Ratio of P to ET
levels_reqd_sens_p2et = c("0.0-0.2","0.2-0.4","0.4-0.6","0.6-0.8", "0.8-1.0")
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
  guides(fill=guide_legend(title=expression(bold("R"*"(CR"["P"]*","~"CR"["ET"]*"):")), face="bold",nrow=1,byrow = TRUE))+
  theme(strip.background = element_blank())+
  theme(panel.background = element_rect(fill=alpha("lightgrey",0.3/2)))+
  theme(panel.grid.major = element_line(colour="black",size=0.1/2),
        panel.grid.minor = element_line(colour="black",size=0.1/2),
        panel.border = element_rect(colour="black",fill=NA))+
  labs(x = "Longitude", y="Latitude")
#guides(fill=guide_legend(title="Cluster Radius (Mean ET)"))
ggsave(plot_cluster_sens_p2et,filename = "Graphs_For_Manuscript_Revised/Figure_8/Cluster_Sensitivity_Mean_P2ET_noSM2RAIN_withBoundary.png", height=115, width=205,dpi=300, units="mm")


}

