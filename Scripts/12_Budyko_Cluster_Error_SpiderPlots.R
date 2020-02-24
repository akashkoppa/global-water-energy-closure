# Script to plot spider plots of CR and CE for different classifications
# Author - Akash Koppa
# Date - 2019-01-08
# Description - Focus on three catchment characteristics - 1) Aridity, 2) CTI, and 3) NDVI. Others come later. 

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
library(ggradar)
library(dplyr)
library(scales)

# set working directory 
setwd("/home/akash/Documents/FifthPaper/Budyko_Analysis_New/")

# read in the shapefile
#shp = shapefile("Input/HydroSHEDS_All_Level05/hybas_all_lev05_v1c.shp")

# read in the clustering metric results
cluster_metric = read.table("Output/Budyko_Cluster_Metrics_RS_noSM2RAIN.txt",header=TRUE)

# separate the area and AI data
hydrosheds_id = as.character(cluster_metric$HydroSHEDS_ID)
area = cluster_metric$Area
ai = cluster_metric$AI
cluster_metric = cluster_metric[,-c(1:3)]

# normalize the cluster metrics by the aridity index
cluster_metric = cluster_metric/ai
cluster_metric[cluster_metric > 10] = NA
cluster_metric$Mean_Radius[which(is.nan(cluster_metric$Mean_Radius))] = NA

# read in the distance metric results
dist_metric = read.table("Output/Budyko_Distance_HydroSHEDS_Lvl05_Rn-CERES_noSM2RAIN.txt",header=TRUE)

# separate the area and AI data
hydrosheds_id = as.character(dist_metric$HydroSHEDS_ID)
area = dist_metric$Area
ai = dist_metric$AI
dist_metric = dist_metric[,-c(1:3)]


# normalize the cluster metrics by the aridity index
dist_metric = dist_metric/ai

# calculate mean of all data combinations for each catchment
dist_metric = data.frame(apply(dist_metric, 1, mean, na.rm=TRUE))
colnames(dist_metric) = "Mean_CE"
dist_metric[dist_metric > 10] = NA
dist_metric$Mean_CE[which(is.nan(dist_metric$Mean_CE))] = NA

# load the required classification criteria file
class_criteria = read.table("Input/Predictors_For_HydrobasinsLvl05_Final.txt",header=TRUE)
colnames(class_criteria) = c("HydroSHEDS_ID","Latitude","Longitude","Area","Elevation","Slope","CTI","NDVI")

# derive data of catchment distribution according to different criteria
# Aridity
# in this I follow the UNEP classification (UNEP (United Nations Environment Programme), 1997. World atlas of desertification 2ED. UNEP, London):
# 1) Hyper Arid:    > 33.3 (<0.03)
# 2) Arid:          33.3 to 5 (0.03 - 0.2)
# 3) Semi-Arid:     5 to 2 (0.2 - 0.5)
# 4) Dry Sub-humid: 2 to 1.5 (0.5 - 0.65)
# 5) Humid:         < 1.5 (>0.65)
index_hyperarid = which(ai>33.33)
index_arid = which(ai<=33.33 & ai>5)
index_semiarid = which(ai<=5 & ai>2)
index_drysubhumid = which(ai<=2 & ai>1.5)
index_humid = which(ai<1.5)
catchments_ai_missing = hydrosheds_id[which(is.na(ai))]
catchments_ai_hyperarid = hydrosheds_id[which(ai>33.33)]
catchmets_ai_arid = hydrosheds_id[which(ai<=33.33 & ai>5)]
catchments_ai_semiarid = hydrosheds_id[which(ai<=5 & ai>2)]
catchments_ai_drysubhumid = hydrosheds_id[which(ai<=2 & ai>1.5)]
catchments_ai_humid = hydrosheds_id[which(ai<1.5)]

# Calculate mean CR and CE for each of the catchment characteristics 
# decompose the statistics according to the aridity indices
# Hyper Arid
mean_cr_hyperarid = mean(cluster_metric$Mean_Radius[index_hyperarid], na.rm=TRUE)
mean_ce_hyperarid = mean(dist_metric$Mean_CE[index_hyperarid], na.rm=TRUE)

# Arid
mean_cr_arid = mean(cluster_metric$Mean_Radius[index_arid], na.rm=TRUE)
mean_ce_arid = mean(dist_metric$Mean_CE[index_arid], na.rm=TRUE)

# Semi Arid
mean_cr_semiarid = mean(cluster_metric$Mean_Radius[index_semiarid], na.rm=TRUE)
mean_ce_semiarid = mean(dist_metric$Mean_CE[index_semiarid], na.rm=TRUE)

# Dry Subhumid
mean_cr_drysubhumid = mean(cluster_metric$Mean_Radius[index_drysubhumid], na.rm=TRUE)
mean_ce_drysubhumid = mean(dist_metric$Mean_CE[index_drysubhumid], na.rm=TRUE)

# Humid
mean_cr_humid = mean(cluster_metric$Mean_Radius[index_humid], na.rm=TRUE)
mean_ce_humid = mean(dist_metric$Mean_CE[index_humid], na.rm=TRUE)

# data frame for spider plot (aridity)
data_aridity_spider = data.frame(Metric = c("Cluster Radius (CR)", "Closure Error (CE)"),
                                 Hyper_Arid = c(mean_cr_hyperarid, mean_ce_hyperarid),
                                 Arid = c(mean_cr_arid, mean_ce_arid),
                                 Semi_Arid = c(mean_cr_arid, mean_ce_arid),
                                 Dry_Subhumid = c(mean_cr_drysubhumid, mean_ce_drysubhumid),
                                 Humid = c(mean_cr_humid, mean_ce_humid))
colnames(data_aridity_spider) = c("Metric", "Hyper Arid", "Arid", "Semi Arid", "Dry Sub-Humid", "Humid")

# rescale the dataset for ggradar
data_aridity_spider[1,-1] = rescale(c(unlist(data_aridity_spider[1,-1])))
data_aridity_spider[2,-1] = rescale(c(unlist(data_aridity_spider[2,-1])))

# change levels 
data_aridity_spider$Metric = factor(data_aridity_spider$Metric, levels=c("Cluster Radius (CR)", "Closure Error (CE)"))

# CTI
cti = class_criteria$CTI
cti[is.nan(cti)] = NA
cti_class = quantile(cti, probs=c(0,0.1,0.25,0.5,0.75,0.9,1.0),na.rm=TRUE)
index_cti1 = which(cti>=cti_class[1] & cti<cti_class[2])
index_cti2 = which(cti>=cti_class[2] & cti<cti_class[3])
index_cti3 = which(cti>=cti_class[3] & cti<cti_class[4])
index_cti4 = which(cti>=cti_class[4] & cti<cti_class[5])
index_cti5 = which(cti>=cti_class[5] & cti<cti_class[6])
index_cti6 = which(cti>=cti_class[6] & cti<cti_class[7])

# decompose the statistics according to different cti classes
# CTI class 1
mean_cr_cti1 = mean(cluster_metric$Mean_Radius[index_cti1], na.rm=TRUE)
mean_ce_cti1 = mean(dist_metric$Mean_CE[index_cti1], na.rm=TRUE)

# CTI class 2
mean_cr_cti2 = mean(cluster_metric$Mean_Radius[index_cti2], na.rm=TRUE)
mean_ce_cti2 = mean(dist_metric$Mean_CE[index_cti2], na.rm=TRUE)

# CTI class 3
mean_cr_cti3 = mean(cluster_metric$Mean_Radius[index_cti3], na.rm=TRUE)
mean_ce_cti3 = mean(dist_metric$Mean_CE[index_cti3], na.rm=TRUE)

# CTI class 4
mean_cr_cti4 = mean(cluster_metric$Mean_Radius[index_cti4], na.rm=TRUE)
mean_ce_cti4 = mean(dist_metric$Mean_CE[index_cti4], na.rm=TRUE)

# CTI class 5
mean_cr_cti5 = mean(cluster_metric$Mean_Radius[index_cti5], na.rm=TRUE)
mean_ce_cti5 = mean(dist_metric$Mean_CE[index_cti5], na.rm=TRUE)

# CTI class 6
mean_cr_cti6 = mean(cluster_metric$Mean_Radius[index_cti6], na.rm=TRUE)
mean_ce_cti6 = mean(dist_metric$Mean_CE[index_cti6], na.rm=TRUE)

# data frame for spider plot (CTI)
data_cti_spider = data.frame(Metric = c("Cluster Radius (CR)", "Closure Error (CE)"),
                             CTI1 = c(mean_cr_cti1, mean_ce_cti1),
                             CTI2 = c(mean_cr_cti2, mean_ce_cti2),
                             CTI3 = c(mean_cr_cti3, mean_ce_cti3),
                             CTI4 = c(mean_cr_cti4, mean_ce_cti4),
                             CTI5 = c(mean_cr_cti5, mean_ce_cti5),
                             CTI6 = c(mean_cr_cti6, mean_ce_cti6))
colnames(data_cti_spider) = c("Metric", "Less Than 9.3", "9.3 - 10.2","10.2 - 11.2",
                              "11.2 - 12.1", "12.1 - 13.1", "13.1 - 47.2")

# rescale data (required before using ggradar)
data_cti_spider[1,-1] = rescale(c(unlist(data_cti_spider[1,-1])))
data_cti_spider[2,-1] = rescale(c(unlist(data_cti_spider[2,-1])))

# change levels 
data_cti_spider$Metric = factor(data_cti_spider$Metric, levels=c("Cluster Radius (CR)", "Closure Error (CE)"))

# NDVI
ndvi = class_criteria$NDVI
ndvi[is.nan(ndvi)] = NA
ndvi_class = quantile(ndvi, probs=c(0,0.1,0.25,0.5,0.75,0.9,1.0),na.rm=TRUE)
index_ndvi1 = which(ndvi>=ndvi_class[1] & ndvi<ndvi_class[2])
index_ndvi2 = which(ndvi>=ndvi_class[2] & ndvi<ndvi_class[3])
index_ndvi3 = which(ndvi>=ndvi_class[3] & ndvi<ndvi_class[4])
index_ndvi4 = which(ndvi>=ndvi_class[4] & ndvi<ndvi_class[5])
index_ndvi5 = which(ndvi>=ndvi_class[5] & ndvi<ndvi_class[6])
index_ndvi6 = which(ndvi>=ndvi_class[6] & ndvi<ndvi_class[7])

# decompose the statistics according to different cti classes
# NDVI class 1
mean_cr_ndvi1 = mean(cluster_metric$Mean_Radius[index_ndvi1], na.rm=TRUE)
mean_ce_ndvi1 = mean(dist_metric$Mean_CE[index_ndvi1], na.rm=TRUE)

# NDVI class 2
mean_cr_ndvi2 = mean(cluster_metric$Mean_Radius[index_ndvi2], na.rm=TRUE)
mean_ce_ndvi2 = mean(dist_metric$Mean_CE[index_ndvi2], na.rm=TRUE)

# NDVI class 3
mean_cr_ndvi3 = mean(cluster_metric$Mean_Radius[index_ndvi3], na.rm=TRUE)
mean_ce_ndvi3 = mean(dist_metric$Mean_CE[index_ndvi3], na.rm=TRUE)

# NDVI class 4
mean_cr_ndvi4 = mean(cluster_metric$Mean_Radius[index_ndvi4], na.rm=TRUE)
mean_ce_ndvi4 = mean(dist_metric$Mean_CE[index_ndvi4], na.rm=TRUE)

# NDVI class 5
mean_cr_ndvi5 = mean(cluster_metric$Mean_Radius[index_ndvi5], na.rm=TRUE)
mean_ce_ndvi5 = mean(dist_metric$Mean_CE[index_ndvi5], na.rm=TRUE)

# NDVI class 6
mean_cr_ndvi6 = mean(cluster_metric$Mean_Radius[index_ndvi6], na.rm=TRUE)
mean_ce_ndvi6 = mean(dist_metric$Mean_CE[index_ndvi6], na.rm=TRUE)

data_ndvi_spider = data.frame(Metric = c("Cluster Radius (CR)", "Closure Error (CE)"),
                             NDVI1 = c(mean_cr_ndvi1, mean_ce_ndvi1),
                             NDVI2 = c(mean_cr_ndvi2, mean_ce_ndvi2),
                             NDVI3 = c(mean_cr_ndvi3, mean_ce_ndvi3),
                             NDVI4 = c(mean_cr_ndvi4, mean_ce_ndvi4),
                             NDVI5 = c(mean_cr_ndvi5, mean_ce_ndvi5),
                             NDVI6 = c(mean_cr_ndvi6, mean_ce_ndvi6))

colnames(data_ndvi_spider) = c("Metric", "Less than 0.09", "0.09 - 0.20","0.20 - 0.38",
                              "0.38 - 0.55", "0.55 - 0.67", "0.67 - 0.85")

# rescale data (required before using ggradar)
data_ndvi_spider[1,-1] = rescale(c(unlist(data_ndvi_spider[1,-1])))
data_ndvi_spider[2,-1] = rescale(c(unlist(data_ndvi_spider[2,-1])))

# change levels 
data_ndvi_spider$Metric = factor(data_ndvi_spider$Metric, levels=c("Cluster Radius (CR)", "Closure Error (CE)"))

# Elevation
elev = class_criteria$Elevation
elev[is.nan(elev)] = NA
elev_class = quantile(elev, probs=c(0,0.1,0.25,0.5,0.75,0.9,1.0),na.rm=TRUE)
index_elev1 = which(elev>=elev_class[1] & elev<elev_class[2])
index_elev2 = which(elev>=elev_class[2] & elev<elev_class[3])
index_elev3 = which(elev>=elev_class[3] & elev<elev_class[4])
index_elev4 = which(elev>=elev_class[4] & elev<elev_class[5])
index_elev5 = which(elev>=elev_class[5] & elev<elev_class[6])
index_elev6 = which(elev>=elev_class[6] & elev<elev_class[7])

# decompose the statistics according to different elev classes
# elev class 1
mean_cr_elev1 = mean(cluster_metric$Mean_Radius[index_elev1], na.rm=TRUE)
mean_ce_elev1 = mean(dist_metric$Mean_CE[index_elev1], na.rm=TRUE)

# elev class 2
mean_cr_elev2 = mean(cluster_metric$Mean_Radius[index_elev2], na.rm=TRUE)
mean_ce_elev2 = mean(dist_metric$Mean_CE[index_elev2], na.rm=TRUE)

# elev class 3
mean_cr_elev3 = mean(cluster_metric$Mean_Radius[index_elev3], na.rm=TRUE)
mean_ce_elev3 = mean(dist_metric$Mean_CE[index_elev3], na.rm=TRUE)

# elev class 4
mean_cr_elev4 = mean(cluster_metric$Mean_Radius[index_elev4], na.rm=TRUE)
mean_ce_elev4 = mean(dist_metric$Mean_CE[index_elev4], na.rm=TRUE)

# elev class 5
mean_cr_elev5 = mean(cluster_metric$Mean_Radius[index_elev5], na.rm=TRUE)
mean_ce_elev5 = mean(dist_metric$Mean_CE[index_elev5], na.rm=TRUE)

# elev class 6
mean_cr_elev6 = mean(cluster_metric$Mean_Radius[index_elev6], na.rm=TRUE)
mean_ce_elev6 = mean(dist_metric$Mean_CE[index_elev6], na.rm=TRUE)

# data frame for spider plot (elev)
data_elev_spider = data.frame(Metric = c("Cluster Radius (CR)", "Closure Error (CE)"),
                             elev1 = c(mean_cr_elev1, mean_ce_elev1),
                             elev2 = c(mean_cr_elev2, mean_ce_elev2),
                             elev3 = c(mean_cr_elev3, mean_ce_elev3),
                             elev4 = c(mean_cr_elev4, mean_ce_elev4),
                             elev5 = c(mean_cr_elev5, mean_ce_elev5),
                             elev6 = c(mean_cr_elev6, mean_ce_elev6))
colnames(data_elev_spider) = c("Metric", "Less Than 59", "59 - 155","155 - 360",
                              "360 - 748", "748 - 1286", "1286 - 5300")

# rescale data (required before using ggradar)
data_elev_spider[1,-1] = rescale(c(unlist(data_elev_spider[1,-1])))
data_elev_spider[2,-1] = rescale(c(unlist(data_elev_spider[2,-1])))

# change levels 
data_elev_spider$Metric = factor(data_elev_spider$Metric, levels=c("Cluster Radius (CR)", "Closure Error (CE)"))



# plot ggradar plots

png(filename = "Output/Spider_Aridity.png",width = 8,height = 6, units="in", res = 300)
ggradar(data_aridity_spider,
        base.size = 10, 
        font.radar="sans",
        #plot.title="c) Elevation",
        group.point.size = rep(c(3,4.5),each=6),
        group.line.width = rep(c(1,1.5),each=6))
dev.off()

png(filename = "Output/Spider_Elev.png",width = 8,height = 6, units="in", res = 300)
ggradar(data_elev_spider,
        base.size = 10, 
        font.radar="sans",
        #plot.title="c) Elevation",
        group.point.size = rep(c(3,4.5),each=7),
        group.line.width = rep(c(1,1.5),each=7))
dev.off()

png(filename = "Output/Spider_CTI.png",width = 8,height = 6, units="in", res = 300)
ggradar(data_cti_spider,
        base.size = 10, 
        font.radar="sans",
        #plot.title="c) Elevation",
        group.point.size = rep(c(3,4.5),each=7),
        group.line.width = rep(c(1,1.5),each=7))
dev.off()

png(filename = "Output/Spider_NDVI.png",width = 8,height = 6, units="in", res = 300)
ggradar(data_ndvi_spider,
        base.size = 10, 
        font.radar="sans",
        #plot.title="c) Elevation",
        group.point.size = rep(c(3,4.5),each=7),
        group.line.width = rep(c(1,1.5),each=7))
dev.off()




