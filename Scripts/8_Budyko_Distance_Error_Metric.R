# Script to summarize and plot error metrics based on the Budyko distances
# Author - Akash Koppa
# Date - 2019-06-13
# Summarize the error metrics according to: 
# 1) Continents, 2) Latitude, 3) Longitude, 4) NDVI, 5) Elevation, 6) CTI, 7) Land cover
# Only RMSE is considered (Do not complicate this too much)
# I think the best plot type is a heatmap. The second option is a bar plot.

# clear workspace
rm(list=ls())

# set working directory
setwd("/home/akash/Documents/FifthPaper/Budyko_Analysis_New/")

# load required libraries
library(ggplot2)
library(raster)
library(stringr)
library(viridis)
library(scales)
library(data.table)
library(reshape2)
library(ggthemes)
library(gridExtra)
library(ggridges)
library(maptools)
library(plyr)
library(RColorBrewer)

# load the required classification criteria file
class_criteria = read.table("Input/Predictors_For_HydrobasinsLvl05_Final.txt",header=TRUE)
colnames(class_criteria) = c("HydroSHEDS_ID","Latitude","Longitude","Area","Elevation","Slope","CTI","NDVI")

# read in the shapefile
shp = shapefile("Input/HydroSHEDS_All_Level05/hybas_all_lev05_v1c.shp")

# read in the distance metric results
dist_metric = read.table("Output/Budyko_Distance_HydroSHEDS_Lvl05_Rn-CERES_noSM2RAIN.txt",header=TRUE)

# read in the parameter 

# separate the area and AI data
hydrosheds_id = as.character(dist_metric$HydroSHEDS_ID)
area = dist_metric$Area
ai = dist_metric$AI
dist_metric = dist_metric[,-c(1:3)]

# get the name of the combinations of P and ET datasets
p_et_dataset = colnames(dist_metric)

# read in the long-term average p, et, and pet datasets
p = read.table("Input/Precipitation_Long_Term.txt", header=TRUE)
p = p[,-1]
p_data = colnames(p)
p = p[,-which(p_data %in% c("CPC.Unifiedv1.0","CRU.TSv4.03","ERA5.Land","ERA5","GPCCv7.0","PREC.Land","UDELv5.0","SM2RAIN.CCI"))]
p_data = p_data[-which(p_data %in% c("CPC.Unifiedv1.0","CRU.TSv4.03","ERA5.Land","ERA5","GPCCv7.0","PREC.Land","UDELv5.0","SM2RAIN.CCI"))]

et = read.table("Input/Evapotranspiration_Long_Term.txt", header=TRUE)
et = et[,-1]
et_data = colnames(et)
et = et[,-which(et_data %in% c("ERA5.Land","FLDAS","FluxCom.RSM","GLDASv2.1"))]
et_data = et_data[-which(et_data %in% c("ERA5.Land","FLDAS","FluxCom.RSM","GLDASv2.1"))]

p_dataset = p_data
et_dataset = et_data

# some preliminary analysis of the catchments
# distribution of catchments according to ai
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

# distribution of catchments according to continents
temp = str_sub(hydrosheds_id,start=1,end=1)
catchments_northamerica = hydrosheds_id[which(temp=="7")]
catchments_southamerica = hydrosheds_id[which(temp=="6")]
catchments_eumideast = hydrosheds_id[which(temp=="2")]
catchments_africa = hydrosheds_id[which(temp=="1")]
catchments_asia = hydrosheds_id[which(temp=="4")]
catchments_australia = hydrosheds_id[which(temp=="5")]
catchments_arctic = hydrosheds_id[which(temp=="8")]
catchments_greenland = hydrosheds_id[which(temp=="9")]
catchments_siberia = hydrosheds_id[which(temp=="3")]
index_na = which(temp=="7")
index_sa = which(temp=="6")
index_eu = which(temp=="2")
index_af = which(temp=="1")
index_as = which(temp=="4")
index_au = which(temp=="5")
index_ar = which(temp=="8")
index_gr = which(temp=="9")
index_si = which(temp=="3")
remove(temp)

# distribution of catchments according to elevation classes
elev = class_criteria$Elevation
elev_class = quantile(elev, probs=c(0,0.1,0.25,0.5,0.75,0.9,1.0))
index_elev1 = which(elev>=elev_class[1] & elev<elev_class[2])
index_elev2 = which(elev>=elev_class[2] & elev<elev_class[3])
index_elev3 = which(elev>=elev_class[3] & elev<elev_class[4])
index_elev4 = which(elev>=elev_class[4] & elev<elev_class[5])
index_elev5 = which(elev>=elev_class[5] & elev<elev_class[6])
index_elev6 = which(elev>=elev_class[6] & elev<elev_class[7])

# distribution of catchments according to NDVI
ndvi = class_criteria$NDVI
ndvi[is.nan(ndvi)] = NA
ndvi_class = quantile(ndvi, probs=c(0,0.1,0.25,0.5,0.75,0.9,1.0),na.rm=TRUE)
index_ndvi1 = which(ndvi>=ndvi_class[1] & ndvi<ndvi_class[2])
index_ndvi2 = which(ndvi>=ndvi_class[2] & ndvi<ndvi_class[3])
index_ndvi3 = which(ndvi>=ndvi_class[3] & ndvi<ndvi_class[4])
index_ndvi4 = which(ndvi>=ndvi_class[4] & ndvi<ndvi_class[5])
index_ndvi5 = which(ndvi>=ndvi_class[5] & ndvi<ndvi_class[6])
index_ndvi6 = which(ndvi>=ndvi_class[6] & ndvi<ndvi_class[7])

# distribution of catchments according to Latitude
# distribution of catchments according to longitude
# distribution of catchments according to CTI
cti = class_criteria$CTI
cti[is.nan(cti)] = NA
cti_class = quantile(cti, probs=c(0,0.1,0.25,0.5,0.75,0.9,1.0),na.rm=TRUE)
index_cti1 = which(cti>=cti_class[1] & cti<cti_class[2])
index_cti2 = which(cti>=cti_class[2] & cti<cti_class[3])
index_cti3 = which(cti>=cti_class[3] & cti<cti_class[4])
index_cti4 = which(cti>=cti_class[4] & cti<cti_class[5])
index_cti5 = which(cti>=cti_class[5] & cti<cti_class[6])
index_cti6 = which(cti>=cti_class[6] & cti<cti_class[7])

# distribution of catchments according to Land cover classification
landcover = read.table("Input/LandCover/landcover.txt",sep=",",header=TRUE)
unique_landcover = as.character(unlist(c(unique(landcover))))
landcover = as.character(unlist(c(landcover)))
index_lc1 = which(landcover == unique_landcover[1])
index_lc2 = which(landcover == unique_landcover[2])
index_lc3 = which(landcover == unique_landcover[3])
index_lc4 = which(landcover == unique_landcover[4])
index_lc5 = which(landcover == unique_landcover[5])
index_lc6 = which(landcover == unique_landcover[6])
index_lc7 = which(landcover == unique_landcover[7])
index_lc8 = which(landcover == unique_landcover[8])
index_lc9 = which(landcover == unique_landcover[9])
index_lc10 = which(landcover == unique_landcover[10])
index_lc11 = which(landcover == unique_landcover[11])
index_lc12 = which(landcover == unique_landcover[12])
index_lc13 = which(landcover == unique_landcover[13])
index_lc14 = which(landcover == unique_landcover[14])
index_lc15 = which(landcover == unique_landcover[15])
index_lc16 = which(landcover == unique_landcover[16])
index_lc17 = which(landcover == unique_landcover[17])
index_lc18 = which(landcover == unique_landcover[18])
index_lc19 = which(landcover == unique_landcover[19])
index_lc20 = which(landcover == unique_landcover[20])
index_lc21 = which(landcover == unique_landcover[21])
index_lc22 = which(landcover == unique_landcover[22])
index_lc23 = which(landcover == unique_landcover[23])
index_lc24 = which(landcover == unique_landcover[24])
index_lc25 = which(landcover == unique_landcover[25])
index_lc26 = which(landcover == unique_landcover[26])
index_lc27 = which(landcover == unique_landcover[27])
index_lc28 = which(landcover == unique_landcover[28])
index_lc29 = which(landcover == unique_landcover[29])

# carefully go through each combination and determine the catchments which have anomalously high distance metrics.
# And store their catchment ID, AI and distance in a list
# also set these values = 5000 (or NA) in the main dist_metric matrix (this is to prevent removing the effect of bad datasets)
#dist_anomalous = list()
#for (i in 1:ncol(dist_metric)){
#  temp = dist_metric[,i]
#  index_reqd = which(temp>1000)
#  data_temp = data.frame(HydroSHEDS_ID = hydrosheds_id[index_reqd],
#                         Area = area[index_reqd],
#                         AI = ai[index_reqd],
#                         Distance = temp[index_reqd],
#                         Index = index_reqd)
#  dist_anomalous[[i]] = data_temp
#  dist_metric[index_reqd,i] = NA # may be change to NA later on but not now
#  dist_metric[which(is.infinite(dist_metric[,i]))] = NA
#}
#names(dist_anomalous) = p_et_dataset
#remove(data_temp)

## all the non-spatial analyses
# calculate some statistics for all the catchments together
rmse_all = sqrt(apply(dist_metric^2,2,mean,na.rm=TRUE))

# decompose the statistics according to the aridity indices
# Hyper Arid
rmse_hyperarid = sqrt(apply(dist_metric[index_hyperarid,]^2,2,mean,na.rm=TRUE))

# Arid
rmse_arid = sqrt(apply(dist_metric[index_arid,]^2,2,mean,na.rm=TRUE))

# Semi Arid
rmse_semiarid = sqrt(apply(dist_metric[index_semiarid,]^2,2,mean,na.rm=TRUE))

# Dry Subhumid
rmse_drysubhumid = sqrt(apply(dist_metric[index_drysubhumid,]^2,2,mean,na.rm=TRUE))

# Humid
rmse_humid = sqrt(apply(dist_metric[index_humid,]^2,2,mean,na.rm=TRUE))

# create a data frame for ggplots later
rmse_aridity = data.frame(Aridity = c("Hyper Arid","Arid","Semi Arid","Dry Sub-humid","Humid"))
temp = as.data.frame(cbind(rmse_hyperarid,rmse_arid,rmse_semiarid,rmse_drysubhumid,rmse_humid))
#temp = apply(temp,1,rescale)
#p_temp = unlist(strsplit(rownames(temp),"[.]"))[seq(from=1,to=length(p_et_dataset)*2,by=2)]
#et_temp = unlist(strsplit(rownames(temp),"[.]"))[seq(from=2,to=length(p_et_dataset)*2,by=2)]
p_temp = rep(p_dataset,each=length(et_dataset))
et_temp = rep(et_dataset,times = length(p_dataset))

p_rmse_aridity = data.frame(P_Dataset = as.character(p_temp),
                          Hyper_Arid = temp[,1],
                          Arid = temp[,2],
                          Semi_Arid = temp[,3],
                          Dry_SubHumid = temp[,4],
                          Humid = temp[,5])
p_rmse_aridity = data.table(p_rmse_aridity)

et_rmse_aridity = data.frame(ET_Dataset = et_temp,
                            Hyper_Arid = temp[,1],
                            Arid = temp[,2],
                            Semi_Arid = temp[,3],
                            Dry_SubHumid = temp[,4],
                            Humid = temp[,5])
et_rmse_aridity = data.table(et_rmse_aridity)

# average according to P and ET datasets
p_rmse_aridity = p_rmse_aridity[,lapply(.SD, mean),by=P_Dataset]
p_rmse_aridity = data.frame(p_rmse_aridity)
et_rmse_aridity = et_rmse_aridity[,lapply(.SD, mean),by=ET_Dataset]
et_rmse_aridity = data.frame(et_rmse_aridity)

# scale the data between 0 and 1
p_rmse_aridity[-1] = lapply(p_rmse_aridity[-1],rescale)
et_rmse_aridity[-1] = lapply(et_rmse_aridity[-1],rescale)

# decompose the statistics according to continents
# North America
rmse_na = sqrt(apply(dist_metric[index_na,]^2,2,mean,na.rm=TRUE))
rmse_na[is.nan(rmse_na)] = NA

# South America
rmse_sa = sqrt(apply(dist_metric[index_sa,]^2,2,mean,na.rm=TRUE))
rmse_na[is.nan(rmse_na)] = NA

# Europe and Middle East
rmse_eu = sqrt(apply(dist_metric[index_eu,]^2,2,mean,na.rm=TRUE))
rmse_sa[is.nan(rmse_sa)] = NA

# Africa
rmse_af = sqrt(apply(dist_metric[index_af,]^2,2,mean,na.rm=TRUE))
rmse_af[is.nan(rmse_af)] = NA

# Asia
rmse_as = sqrt(apply(dist_metric[index_as,]^2,2,mean,na.rm=TRUE))
rmse_as[is.nan(rmse_as)] = NA

# Australia
rmse_au = sqrt(apply(dist_metric[index_au,]^2,2,mean,na.rm=TRUE))
rmse_au[is.nan(rmse_au)] = NA

# Arctic
rmse_ar = sqrt(apply(dist_metric[index_ar,]^2,2,mean,na.rm=TRUE))
rmse_ar[is.nan(rmse_ar)] = NA

# Greenland
rmse_gr = sqrt(apply(dist_metric[index_gr,]^2,2,mean,na.rm=TRUE))
rmse_gr[is.nan(rmse_gr)] = NA

# Siberia
rmse_si = sqrt(apply(dist_metric[index_si,]^2,2,mean,na.rm=TRUE))
rmse_si[is.nan(rmse_si)] = NA

# create a data frame for ggplots later
rmse_continent = data.frame(Aridity = c("North America","South America","Europe and Middle East","Asia","Australia","Arctic","Greenland","Siberia"))
temp = as.data.frame(cbind(rmse_na,rmse_sa,rmse_eu,rmse_eu,rmse_af,rmse_as,rmse_au,rmse_ar,rmse_gr,rmse_si))
#temp[is.nan(temp)] = NA
#temp = apply(temp,1,rescale)
#p_temp = unlist(strsplit(rownames(temp),"[.]"))[seq(from=1,to=length(p_et_dataset)*2,by=2)]
#et_temp = unlist(strsplit(rownames(temp),"[.]"))[seq(from=2,to=length(p_et_dataset)*2,by=2)]
p_rmse_continent = data.frame(P_Dataset = as.character(p_temp),
                            North_America = temp[,1],
                            South_America = temp[,2],
                            Europe_ME = temp[,3],
                            Asia = temp[,4],
                            Australia = temp[,5],
                            Arctic = temp[,6],
                            Greenland = temp[,7],
                            Siberia = temp[,8])
p_rmse_continent = data.table(p_rmse_continent)

et_rmse_continent = data.frame(ET_Dataset = as.character(et_temp),
                              North_America = temp[,1],
                              South_America = temp[,2],
                              Europe_ME = temp[,3],
                              Asia = temp[,4],
                              Australia = temp[,5],
                              Arctic = temp[,6],
                              Greenland = temp[,7],
                              Siberia = temp[,8])
et_rmse_continent = data.table(et_rmse_continent)

# average according to P and ET datasets
p_rmse_continent = p_rmse_continent[,lapply(.SD, mean,na.rm=TRUE),by=P_Dataset]
p_rmse_continent = data.frame(p_rmse_continent)
et_rmse_continent = et_rmse_continent[,lapply(.SD, mean,na.rm=TRUE),by=ET_Dataset]
et_rmse_continent = data.frame(et_rmse_continent)

# scale the data between 0 and 1
p_rmse_continent[-1] = lapply(p_rmse_continent[-1],rescale)
et_rmse_continent[-1] = lapply(et_rmse_continent[-1],rescale)

## Decompose the RMSE metric according to elevation
# Elev Class 1
rmse_elev1 = sqrt(apply(dist_metric[index_elev1,]^2,2,mean,na.rm=TRUE))
rmse_elev1[is.nan(rmse_elev1)] = NA

# Elev Class 2
rmse_elev2 = sqrt(apply(dist_metric[index_elev2,]^2,2,mean,na.rm=TRUE))
rmse_elev2[is.nan(rmse_elev2)] = NA

# Elev Class 3
rmse_elev3 = sqrt(apply(dist_metric[index_elev3,]^2,2,mean,na.rm=TRUE))
rmse_elev3[is.nan(rmse_elev3)] = NA

# Elev Class 4
rmse_elev4 = sqrt(apply(dist_metric[index_elev4,]^2,2,mean,na.rm=TRUE))
rmse_elev4[is.nan(rmse_elev4)] = NA

# Elev Class 5
rmse_elev5 = sqrt(apply(dist_metric[index_elev5,]^2,2,mean,na.rm=TRUE))
rmse_elev5[is.nan(rmse_elev5)] = NA

# Elev Class 6
rmse_elev6 = sqrt(apply(dist_metric[index_elev6,]^2,2,mean,na.rm=TRUE))
rmse_elev6[is.nan(rmse_elev6)] = NA

# create a data frame for ggplots later
rmse_elevation = data.frame(Elevation = c("Elevation1","Elevation2","Elevation3","Elevation4","Elevation5","Elevation6"))
temp = as.data.frame(cbind(rmse_elev1,rmse_elev2,rmse_elev3,rmse_elev4,rmse_elev5,rmse_elev6))
#temp[is.nan(temp)] = NA
#temp = apply(temp,1,rescale)
#p_temp = unlist(strsplit(rownames(temp),"[.]"))[seq(from=1,to=length(p_et_dataset)*2,by=2)]
#et_temp = unlist(strsplit(rownames(temp),"[.]"))[seq(from=2,to=length(p_et_dataset)*2,by=2)]
p_rmse_elevation = data.frame(P_Dataset = as.character(p_temp),
                              Elev_0_To_60 = temp[,1],
                              Elev_60_To_156 = temp[,2],
                              Elev_156_To_360 = temp[,3],
                              Elev_360_To_748 = temp[,4],
                              Elev_748_To_1287 = temp[,5],
                              Elev_1287_To_5300 = temp[,6])
p_rmse_elevation = data.table(p_rmse_elevation)

et_rmse_elevation = data.frame(ET_Dataset = as.character(et_temp),
                               Elev_0_To_60 = temp[,1],
                               Elev_60_To_156 = temp[,2],
                               Elev_156_To_360 = temp[,3],
                               Elev_360_To_748 = temp[,4],
                               Elev_748_To_1287 = temp[,5],
                               Elev_1287_To_5300 = temp[,6])
et_rmse_elevation = data.table(et_rmse_elevation)

# average according to P and ET datasets
p_rmse_elevation = p_rmse_elevation[,lapply(.SD, mean,na.rm=TRUE),by=P_Dataset]
p_rmse_elevation = data.frame(p_rmse_elevation)
et_rmse_elevation = et_rmse_elevation[,lapply(.SD, mean,na.rm=TRUE),by=ET_Dataset]
et_rmse_elevation = data.frame(et_rmse_elevation)

# scale the data between 0 and 1
p_rmse_elevation[-1] = lapply(p_rmse_elevation[-1],rescale)
et_rmse_elevation[-1] = lapply(et_rmse_elevation[-1],rescale)

## Decompose the RMSE metric according to NDVI
# NDVI Class 1
rmse_ndvi1 = sqrt(apply(dist_metric[index_ndvi1,]^2,2,mean,na.rm=TRUE))
rmse_ndvi1[is.nan(rmse_ndvi1)] = NA

# ndvi Class 2
rmse_ndvi2 = sqrt(apply(dist_metric[index_ndvi2,]^2,2,mean,na.rm=TRUE))
rmse_ndvi2[is.nan(rmse_ndvi2)] = NA

# ndvi Class 3
rmse_ndvi3 = sqrt(apply(dist_metric[index_ndvi3,]^2,2,mean,na.rm=TRUE))
rmse_ndvi3[is.nan(rmse_ndvi3)] = NA

# ndvi Class 4
rmse_ndvi4 = sqrt(apply(dist_metric[index_ndvi4,]^2,2,mean,na.rm=TRUE))
rmse_ndvi4[is.nan(rmse_ndvi4)] = NA

# ndvi Class 5
rmse_ndvi5 = sqrt(apply(dist_metric[index_ndvi5,]^2,2,mean,na.rm=TRUE))
rmse_ndvi5[is.nan(rmse_ndvi5)] = NA

# ndvi Class 6
rmse_ndvi6 = sqrt(apply(dist_metric[index_ndvi6,]^2,2,mean,na.rm=TRUE))
rmse_ndvi6[is.nan(rmse_ndvi6)] = NA

# create a data frame for ggplots later
rmse_ndvi = data.frame(NDVI = c("NDVI1","NDVI2","NDVI3","NDVI4","NDVI5","NDVI6"))
temp = as.data.frame(cbind(rmse_ndvi1,rmse_ndvi2,rmse_ndvi3,rmse_ndvi4,rmse_ndvi5,rmse_ndvi6))
#temp[is.nan(temp)] = NA
#temp = apply(temp,1,rescale)
#p_temp = unlist(strsplit(rownames(temp),"[.]"))[seq(from=1,to=length(p_et_dataset)*2,by=2)]
#et_temp = unlist(strsplit(rownames(temp),"[.]"))[seq(from=2,to=length(p_et_dataset)*2,by=2)]
p_rmse_ndvi = data.frame(P_Dataset = as.character(p_temp),
                              NDVI_Less_Than_0.09 = temp[,1],
                              NDVI_0.09_To_0.20 = temp[,2],
                              NDVI_0.20_To_0.38 = temp[,3],
                              NDVI_0.38_To_0.55 = temp[,4],
                              NDVI_0.55_To_0.67 = temp[,5],
                              NDVI_0.67_To_0.855 = temp[,6])
p_rmse_ndvi = data.table(p_rmse_ndvi)

et_rmse_ndvi = data.frame(ET_Dataset = as.character(et_temp),
                          NDVI_Less_Than_0.09 = temp[,1],
                          NDVI_0.09_To_0.20 = temp[,2],
                          NDVI_0.20_To_0.38 = temp[,3],
                          NDVI_0.38_To_0.55 = temp[,4],
                          NDVI_0.55_To_0.67 = temp[,5],
                          NDVI_0.67_To_0.855 = temp[,6])
et_rmse_ndvi = data.table(et_rmse_ndvi)

# average according to P and ET datasets
p_rmse_ndvi = p_rmse_ndvi[,lapply(.SD, mean,na.rm=TRUE),by=P_Dataset]
p_rmse_ndvi = data.frame(p_rmse_ndvi)
et_rmse_ndvi = et_rmse_ndvi[,lapply(.SD, mean,na.rm=TRUE),by=ET_Dataset]
et_rmse_ndvi = data.frame(et_rmse_ndvi)

# scale the data between 0 and 1
p_rmse_ndvi[-1] = lapply(p_rmse_ndvi[-1],rescale)
et_rmse_ndvi[-1] = lapply(et_rmse_ndvi[-1],rescale)

# CTI 
## Decompose the RMSE metric according to cti
# cti Class 1
rmse_cti1 = sqrt(apply(dist_metric[index_cti1,]^2,2,mean,na.rm=TRUE))
rmse_cti1[is.nan(rmse_cti1)] = NA

# cti Class 2
rmse_cti2 = sqrt(apply(dist_metric[index_cti2,]^2,2,mean,na.rm=TRUE))
rmse_cti2[is.nan(rmse_cti2)] = NA

# cti Class 3
rmse_cti3 = sqrt(apply(dist_metric[index_cti3,]^2,2,mean,na.rm=TRUE))
rmse_cti3[is.nan(rmse_cti3)] = NA

# cti Class 4
rmse_cti4 = sqrt(apply(dist_metric[index_cti4,]^2,2,mean,na.rm=TRUE))
rmse_cti4[is.nan(rmse_cti4)] = NA

# cti Class 5
rmse_cti5 = sqrt(apply(dist_metric[index_cti5,]^2,2,mean,na.rm=TRUE))
rmse_cti5[is.nan(rmse_cti5)] = NA

# cti Class 6
rmse_cti6 = sqrt(apply(dist_metric[index_cti6,]^2,2,mean,na.rm=TRUE))
rmse_cti6[is.nan(rmse_cti6)] = NA

# create a data frame for ggplots later
rmse_cti = data.frame(cti = c("cti1","cti2","cti3","cti4","cti5","cti6"))
temp = as.data.frame(cbind(rmse_cti1,rmse_cti2,rmse_cti3,rmse_cti4,rmse_cti5,rmse_cti6))
#temp[is.nan(temp)] = NA
#temp = apply(temp,1,rescale)
#p_temp = unlist(strsplit(rownames(temp),"[.]"))[seq(from=1,to=length(p_et_dataset)*2,by=2)]
#et_temp = unlist(strsplit(rownames(temp),"[.]"))[seq(from=2,to=length(p_et_dataset)*2,by=2)]
p_rmse_cti = data.frame(P_Dataset = as.character(p_temp),
                         CTI_Less_Than_9.3 = temp[,1],
                         CTI_9.3_To_10.2 = temp[,2],
                         CTI_10.2_To_11.2 = temp[,3],
                         CTI_11.2_To_12.1 = temp[,4],
                         CTI_12.1_To_13.1 = temp[,5],
                         CTI_13.1_To_47.1 = temp[,6])
p_rmse_cti = data.table(p_rmse_cti)

et_rmse_cti = data.frame(ET_Dataset = as.character(et_temp),
                         CTI_Less_Than_9.3 = temp[,1],
                         CTI_9.3_To_10.2 = temp[,2],
                         CTI_10.2_To_11.2 = temp[,3],
                         CTI_11.2_To_12.1 = temp[,4],
                         CTI_12.1_To_13.1 = temp[,5],
                         CTI_13.1_To_47.1 = temp[,6])
et_rmse_cti = data.table(et_rmse_cti)

# average according to P and ET datasets
p_rmse_cti = p_rmse_cti[,lapply(.SD, mean,na.rm=TRUE),by=P_Dataset]
p_rmse_cti = data.frame(p_rmse_cti)
et_rmse_cti = et_rmse_cti[,lapply(.SD, mean,na.rm=TRUE),by=ET_Dataset]
et_rmse_cti = data.frame(et_rmse_cti)

# scale the data between 0 and 1
p_rmse_cti[-1] = lapply(p_rmse_cti[-1],rescale)
et_rmse_cti[-1] = lapply(et_rmse_cti[-1],rescale)

# Land use 
## RMSE according to landuse class
rmse_lc1 = sqrt(apply(dist_metric[index_lc1,]^2,2,mean,na.rm=TRUE)); rmse_lc1[is.nan(rmse_lc1)] = NA
rmse_lc2 = sqrt(apply(dist_metric[index_lc2,]^2,2,mean,na.rm=TRUE)); rmse_lc1[is.nan(rmse_lc2)] = NA
rmse_lc3 = sqrt(apply(dist_metric[index_lc3,]^2,2,mean,na.rm=TRUE)); rmse_lc1[is.nan(rmse_lc3)] = NA
rmse_lc4 = sqrt(apply(dist_metric[index_lc4,]^2,2,mean,na.rm=TRUE)); rmse_lc1[is.nan(rmse_lc4)] = NA
rmse_lc5 = sqrt(apply(dist_metric[index_lc5,]^2,2,mean,na.rm=TRUE)); rmse_lc1[is.nan(rmse_lc5)] = NA
rmse_lc6 = sqrt(apply(dist_metric[index_lc6,]^2,2,mean,na.rm=TRUE)); rmse_lc1[is.nan(rmse_lc6)] = NA
rmse_lc7 = sqrt(apply(dist_metric[index_lc7,]^2,2,mean,na.rm=TRUE)); rmse_lc1[is.nan(rmse_lc7)] = NA
rmse_lc8 = sqrt(apply(dist_metric[index_lc8,]^2,2,mean,na.rm=TRUE)); rmse_lc1[is.nan(rmse_lc8)] = NA
rmse_lc9 = sqrt(apply(dist_metric[index_lc9,]^2,2,mean,na.rm=TRUE)); rmse_lc1[is.nan(rmse_lc9)] = NA
rmse_lc10 = sqrt(apply(dist_metric[index_lc10,]^2,2,mean,na.rm=TRUE)); rmse_lc1[is.nan(rmse_lc10)] = NA
rmse_lc11 = sqrt(apply(dist_metric[index_lc11,]^2,2,mean,na.rm=TRUE)); rmse_lc1[is.nan(rmse_lc11)] = NA
rmse_lc12 = sqrt(apply(dist_metric[index_lc12,]^2,2,mean,na.rm=TRUE)); rmse_lc1[is.nan(rmse_lc12)] = NA
rmse_lc13 = sqrt(apply(dist_metric[index_lc13,]^2,2,mean,na.rm=TRUE)); rmse_lc1[is.nan(rmse_lc13)] = NA
rmse_lc14 = sqrt(apply(dist_metric[index_lc14,]^2,2,mean,na.rm=TRUE)); rmse_lc1[is.nan(rmse_lc14)] = NA
rmse_lc15 = sqrt(apply(dist_metric[index_lc15,]^2,2,mean,na.rm=TRUE)); rmse_lc1[is.nan(rmse_lc15)] = NA
rmse_lc16 = sqrt(apply(dist_metric[index_lc16,]^2,2,mean,na.rm=TRUE)); rmse_lc1[is.nan(rmse_lc16)] = NA
rmse_lc17 = sqrt(apply(dist_metric[index_lc17,]^2,2,mean,na.rm=TRUE)); rmse_lc1[is.nan(rmse_lc17)] = NA
rmse_lc18 = sqrt(apply(dist_metric[index_lc18,]^2,2,mean,na.rm=TRUE)); rmse_lc1[is.nan(rmse_lc18)] = NA
rmse_lc19 = sqrt(apply(dist_metric[index_lc19,]^2,2,mean,na.rm=TRUE)); rmse_lc1[is.nan(rmse_lc19)] = NA
rmse_lc20 = sqrt(apply(dist_metric[index_lc20,]^2,2,mean,na.rm=TRUE)); rmse_lc1[is.nan(rmse_lc20)] = NA
rmse_lc21 = sqrt(apply(dist_metric[index_lc21,]^2,2,mean,na.rm=TRUE)); rmse_lc1[is.nan(rmse_lc21)] = NA
rmse_lc22 = sqrt(apply(dist_metric[index_lc22,]^2,2,mean,na.rm=TRUE)); rmse_lc1[is.nan(rmse_lc22)] = NA
rmse_lc23 = sqrt(apply(dist_metric[index_lc23,]^2,2,mean,na.rm=TRUE)); rmse_lc1[is.nan(rmse_lc23)] = NA
rmse_lc24 = sqrt(apply(dist_metric[index_lc24,]^2,2,mean,na.rm=TRUE)); rmse_lc1[is.nan(rmse_lc24)] = NA
rmse_lc25 = sqrt(apply(dist_metric[index_lc25,]^2,2,mean,na.rm=TRUE)); rmse_lc1[is.nan(rmse_lc25)] = NA
rmse_lc26 = sqrt(apply(dist_metric[index_lc26,]^2,2,mean,na.rm=TRUE)); rmse_lc1[is.nan(rmse_lc26)] = NA
rmse_lc27 = sqrt(apply(dist_metric[index_lc27,]^2,2,mean,na.rm=TRUE)); rmse_lc1[is.nan(rmse_lc27)] = NA
rmse_lc28 = sqrt(apply(dist_metric[index_lc28,]^2,2,mean,na.rm=TRUE)); rmse_lc1[is.nan(rmse_lc28)] = NA
rmse_lc29 = sqrt(apply(dist_metric[index_lc29,]^2,2,mean,na.rm=TRUE)); rmse_lc1[is.nan(rmse_lc29)] = NA

# create a data frame for ggplots later
rmse_lc = data.frame(cti = unique_landcover)
temp = as.data.frame(cbind(rmse_lc1, rmse_lc2, rmse_lc3, rmse_lc4, rmse_lc5, rmse_lc6, rmse_lc7, rmse_lc8, rmse_lc9, rmse_lc10,
                           rmse_lc11, rmse_lc12, rmse_lc13, rmse_lc14, rmse_lc15, rmse_lc16, rmse_lc17, rmse_lc18, rmse_lc19, rmse_lc20,
                           rmse_lc21, rmse_lc22, rmse_lc23, rmse_lc24, rmse_lc25, rmse_lc26, rmse_lc27, rmse_lc28, rmse_lc29))
#temp[is.nan(temp)] = NA
#temp = apply(temp,1,rescale)
#p_temp = unlist(strsplit(rownames(temp),"[.]"))[seq(from=1,to=length(p_et_dataset)*2,by=2)]
#et_temp = unlist(strsplit(rownames(temp),"[.]"))[seq(from=2,to=length(p_et_dataset)*2,by=2)]
p_rmse_lc = data.frame(P_Dataset = as.character(p_temp))
p_rmse_lc = cbind(p_rmse_lc,temp)
colnames(p_rmse_lc) = c("P_Dataset",unique_landcover)
p_rmse_lc = data.table(p_rmse_lc)

et_rmse_lc = data.frame(ET_Dataset = as.character(et_temp))
et_rmse_lc = cbind(et_rmse_lc,temp)
colnames(et_rmse_lc) = c("ET_Dataset",unique_landcover)
et_rmse_lc = data.table(et_rmse_lc)

# average according to P and ET datasets
p_rmse_lc = p_rmse_lc[,lapply(.SD, mean,na.rm=TRUE),by=P_Dataset]
p_rmse_lc = data.frame(p_rmse_lc)
et_rmse_lc = et_rmse_lc[,lapply(.SD, mean,na.rm=TRUE),by=ET_Dataset]
et_rmse_lc = data.frame(et_rmse_lc)

# scale the data between 0 and 1
p_rmse_lc[-1] = lapply(p_rmse_lc[-1],rescale)
et_rmse_lc[-1] = lapply(et_rmse_lc[-1],rescale)

# Data for the ridge or box plot of distances of all the precipitation and evapotranpiration
p_distance_all = NULL
for (i in 1:length(p_dataset)){
  test = dist_metric[,which(p_temp==p_dataset[i])]
  test = apply(test,1,mean,na.rm=TRUE)
  p_distance_all = cbind(p_distance_all, test)
}
colnames(p_distance_all) = p_dataset

et_distance_all = NULL
for (i in 1:length(et_dataset)){
  test = dist_metric[,which(et_temp==et_dataset[i])]
  test = apply(test,1,mean,na.rm=TRUE)
  et_distance_all = cbind(et_distance_all, test)
}
colnames(et_distance_all) = et_dataset
remove(test)


flag = 2
if(flag>0){
## All heat map (or bar) plots according to different criteria
# Aridity
p_temp_rmse_aridity = melt(data = p_rmse_aridity,id.vars = "P_Dataset")
p_temp_rmse_aridity = cbind(p_temp_rmse_aridity,as.data.frame("Precipitation"))
colnames(p_temp_rmse_aridity) = c("Dataset","Classification","RMSE","Variable")
et_temp_rmse_aridity = melt(data = et_rmse_aridity,id.vars = "ET_Dataset")
et_temp_rmse_aridity = cbind(et_temp_rmse_aridity,as.data.frame("Evapotranspiration"))
colnames(et_temp_rmse_aridity) = c("Dataset","Classification","RMSE","Variable")
data_rmse_aridity = rbind(p_temp_rmse_aridity,et_temp_rmse_aridity)
#remove(p_temp_rmse_aridity,et_temp_rmse_aridity)

# plot the aridity heat map
plot_p_rmse_aridity = ggplot(data = p_temp_rmse_aridity, aes(x=Classification, y=Dataset,fill=RMSE))+
                    geom_tile(color="white") +
                    scale_fill_viridis()+
                    geom_point(aes(size=ifelse(RMSE==0.0, "dot","no_dot")),color="red",shape=8) +
                    scale_size_manual(values=c(dot=2, no_dot=NA), guide="none")+
                    theme_tufte(base_family="Helvetica")+
                    theme(axis.ticks=element_blank())+
                    theme(axis.text=element_text(size=6))+
  theme(axis.text.x=element_text(size=6))+
                    theme(panel.background = element_blank())+
                    theme(plot.background=element_blank())+
                    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
                    xlab("Aridity")+ggtitle("a) Precipitation")
  #coord_flip()

plot_et_rmse_aridity = ggplot(data = et_temp_rmse_aridity, aes(x=Classification, y=Dataset,fill=RMSE))+
                      geom_tile(color="white") +
                      scale_fill_viridis()+
                      geom_point(aes(size=ifelse(RMSE==0.0, "dot","no_dot")),color="red",shape=8) +
                      scale_size_manual(values=c(dot=2, no_dot=NA), guide="none")+
                      theme_tufte(base_family="Helvetica")+
                      theme(axis.ticks=element_blank())+
                      theme(axis.text=element_text(size=6))+
  theme(axis.text.x=element_text(size=6))+
  #theme(legend.position = "none" )+
                      theme(panel.background = element_blank())+
                      theme(plot.background=element_blank())+
                      theme(axis.text.x = element_text(angle = 90, hjust = 1))+
                      xlab("Aridity")+ggtitle("b) Evapotranspiration")
  #coord_flip()
grid_aridity = grid.arrange(grobs = list(plot_p_rmse_aridity,plot_et_rmse_aridity), layout_matrix=rbind(c(1,2)))
ggsave("Output/RMSE_Aridity_New.png",grid_aridity,width=8, height = 4, dpi=300)

# Continents
p_temp_rmse_continent = melt(data = p_rmse_continent,id.vars = "P_Dataset")
p_temp_rmse_continent = cbind(p_temp_rmse_continent,as.data.frame("Precipitation"))
colnames(p_temp_rmse_continent) = c("Dataset","Classification","RMSE","Variable")
et_temp_rmse_continent = melt(data = et_rmse_continent,id.vars = "ET_Dataset")
et_temp_rmse_continent = cbind(et_temp_rmse_continent,as.data.frame("Evapotranspiration"))
colnames(et_temp_rmse_continent) = c("Dataset","Classification","RMSE","Variable")
data_rmse_continent = rbind(p_temp_rmse_continent,et_temp_rmse_continent)
#remove(p_temp_rmse_aridity,et_temp_rmse_aridity)

# plot the aridity heat map
plot_p_rmse_continent = ggplot(data = p_temp_rmse_continent, aes(x=Classification, y=Dataset,fill=RMSE))+
  geom_tile(color="white") +
  scale_fill_viridis()+
  geom_point(aes(size=ifelse(RMSE==0.0, "dot","no_dot")),color="red",shape=8) +
  scale_size_manual(values=c(dot=2, no_dot=NA), guide="none")+
  theme_tufte(base_family="Helvetica")+
  theme(axis.ticks=element_blank())+
  theme(axis.text=element_text(size=6))+
  theme(axis.text.x=element_text(size=6))+
  theme(panel.background = element_blank())+
  theme(plot.background=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("Regions")+ggtitle("a) Precipitation")
  #coord_flip()

plot_et_rmse_continent = ggplot(data = et_temp_rmse_continent, aes(x=Classification, y=Dataset,fill=RMSE))+
  geom_tile(color="white") +
  scale_fill_viridis()+
  geom_point(aes(size=ifelse(RMSE==0.0, "dot","no_dot")),color="red",shape=8) +
  scale_size_manual(values=c(dot=2, no_dot=NA), guide="none")+
  theme_tufte(base_family="Helvetica")+
  theme(axis.ticks=element_blank())+
  theme(axis.text=element_text(size=6))+
  theme(axis.text.x=element_text(size=6))+
  #theme(legend.position = "none" )+
  theme(panel.background = element_blank())+
  theme(plot.background=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("Regions")+ggtitle("b) Evapotranspiration")
  #coord_flip()
grid_continent = grid.arrange(grobs = list(plot_p_rmse_continent,plot_et_rmse_continent), layout_matrix=rbind(c(1,2)))
ggsave("Output/RMSE_Continent_New.png",grid_continent,width=8, height = 4, dpi=300)

# Elevation
p_temp_rmse_elevation = melt(data = p_rmse_elevation,id.vars = "P_Dataset")
p_temp_rmse_elevation = cbind(p_temp_rmse_elevation,as.data.frame("Precipitation"))
colnames(p_temp_rmse_elevation) = c("Dataset","Classification","RMSE","Variable")
et_temp_rmse_elevation = melt(data = et_rmse_elevation,id.vars = "ET_Dataset")
et_temp_rmse_elevation= cbind(et_temp_rmse_elevation,as.data.frame("Evapotranspiration"))
colnames(et_temp_rmse_elevation) = c("Dataset","Classification","RMSE","Variable")
data_rmse_elevation = rbind(p_temp_rmse_elevation,et_temp_rmse_elevation)
#remove(p_temp_rmse_aridity,et_temp_rmse_aridity)

# plot the aridity heat map
plot_p_rmse_elevation = ggplot(data = p_temp_rmse_elevation, aes(x=Classification, y=Dataset,fill=RMSE))+
  geom_tile(color="white") +
  scale_fill_viridis()+
  geom_point(aes(size=ifelse(RMSE==0.0, "dot","no_dot")),color="red",shape=8) +
  scale_size_manual(values=c(dot=2, no_dot=NA), guide="none")+
  theme_tufte(base_family="Helvetica")+
  theme(axis.ticks=element_blank())+
  theme(axis.text=element_text(size=6))+
  theme(axis.text.x=element_text(size=6))+
  theme(panel.background = element_blank())+
  theme(plot.background=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("Elevation")+ggtitle("a) Precipitation")
  #coord_flip()

plot_et_rmse_elevation = ggplot(data = et_temp_rmse_elevation, aes(x=Classification, y=Dataset,fill=RMSE))+
  geom_tile(color="white") +
  scale_fill_viridis()+
  geom_point(aes(size=ifelse(RMSE==0.0, "dot","no_dot")),color="red",shape=8) +
  scale_size_manual(values=c(dot=2, no_dot=NA), guide="none")+
  theme_tufte(base_family="Helvetica")+
  theme(axis.ticks=element_blank())+
  theme(axis.text=element_text(size=6))+
  theme(axis.text.x=element_text(size=6))+
  #theme(legend.position = "none" )+
  theme(panel.background = element_blank())+
  theme(plot.background=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("Elevation")+ggtitle("b) Evapotranspiration")
  #coord_flip()
grid_elevation = grid.arrange(grobs = list(plot_p_rmse_elevation,plot_et_rmse_elevation), layout_matrix=rbind(c(1,2)))
ggsave("Output/RMSE_Elevation_New.png",grid_elevation,width=8, height = 4, dpi=300)

# NDVI
p_temp_rmse_ndvi = melt(data = p_rmse_ndvi,id.vars = "P_Dataset")
p_temp_rmse_ndvi = cbind(p_temp_rmse_ndvi,as.data.frame("Precipitation"))
colnames(p_temp_rmse_ndvi) = c("Dataset","Classification","RMSE","Variable")
et_temp_rmse_ndvi = melt(data = et_rmse_ndvi,id.vars = "ET_Dataset")
et_temp_rmse_ndvi= cbind(et_temp_rmse_ndvi,as.data.frame("Evapotranspiration"))
colnames(et_temp_rmse_ndvi) = c("Dataset","Classification","RMSE","Variable")
data_rmse_ndvi = rbind(p_temp_rmse_ndvi,et_temp_rmse_ndvi)
#remove(p_temp_rmse_aridity,et_temp_rmse_aridity)

# plot the aridity heat map
plot_p_rmse_ndvi = ggplot(data = p_temp_rmse_ndvi, aes(x=Classification, y=Dataset,fill=RMSE))+
  geom_tile(color="white") +
  scale_fill_viridis()+
  geom_point(aes(size=ifelse(RMSE==0.0, "dot","no_dot")),color="red",shape=8) +
  scale_size_manual(values=c(dot=2, no_dot=NA), guide="none")+
  theme_tufte(base_family="Helvetica")+
  theme(axis.ticks=element_blank())+
  theme(axis.text=element_text(size=6))+
  theme(axis.text.x=element_text(size=6))+
  theme(panel.background = element_blank())+
  theme(plot.background=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("NDVI")+ggtitle("a) Precipitation")
  #coord_flip()

plot_et_rmse_ndvi = ggplot(data = et_temp_rmse_ndvi, aes(x=Classification, y=Dataset,fill=RMSE))+
  geom_tile(color="white") +
  scale_fill_viridis()+
  geom_point(aes(size=ifelse(RMSE==0.0, "dot","no_dot")),color="red",shape=8) +
  scale_size_manual(values=c(dot=2, no_dot=NA), guide="none")+
  theme_tufte(base_family="Helvetica")+
  theme(axis.ticks=element_blank())+
  theme(axis.text=element_text(size=6))+
  theme(axis.text.x=element_text(size=6))+
  #theme(legend.position = "none" )+
  theme(panel.background = element_blank())+
  theme(plot.background=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("NDVI")+ggtitle("b) Evapotranspiration")
  #coord_flip()
grid_ndvi = grid.arrange(grobs = list(plot_p_rmse_ndvi,plot_et_rmse_ndvi), layout_matrix=rbind(c(1,2)))
ggsave("Output/RMSE_NDVI_New.png",grid_ndvi,width=8, height = 4, dpi=300)

# CTI
p_temp_rmse_cti = melt(data = p_rmse_cti,id.vars = "P_Dataset")
p_temp_rmse_cti = cbind(p_temp_rmse_cti,as.data.frame("Precipitation"))
colnames(p_temp_rmse_cti) = c("Dataset","Classification","RMSE","Variable")
et_temp_rmse_cti = melt(data = et_rmse_cti,id.vars = "ET_Dataset")
et_temp_rmse_cti= cbind(et_temp_rmse_cti,as.data.frame("Evapotranspiration"))
colnames(et_temp_rmse_cti) = c("Dataset","Classification","RMSE","Variable")
data_rmse_cti = rbind(p_temp_rmse_cti,et_temp_rmse_cti)
#remove(p_temp_rmse_aridity,et_temp_rmse_aridity)

# plot the aridity heat map
plot_p_rmse_cti = ggplot(data = p_temp_rmse_cti, aes(x=Classification, y=Dataset,fill=RMSE))+
  geom_tile(color="white") +
  scale_fill_viridis()+
  geom_point(aes(size=ifelse(RMSE==0.0, "dot","no_dot")),color="red",shape=8) +
  scale_size_manual(values=c(dot=2, no_dot=NA), guide="none")+
  theme_tufte(base_family="Helvetica")+
  theme(axis.ticks=element_blank())+
  theme(axis.text=element_text(size=6))+
  theme(axis.text.x=element_text(size=6))+
  theme(panel.background = element_blank())+
  theme(plot.background=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("CTI")+ggtitle("a) Precipitation")
  #coord_flip()

plot_et_rmse_cti = ggplot(data = et_temp_rmse_cti, aes(x=Classification, y=Dataset,fill=RMSE))+
  geom_tile(color="white") +
  scale_fill_viridis()+
  geom_point(aes(size=ifelse(RMSE==0.0, "dot","no_dot")),color="red",shape=8) +
  scale_size_manual(values=c(dot=2, no_dot=NA), guide="none")+
  theme_tufte(base_family="Helvetica")+
  theme(axis.ticks=element_blank())+
  theme(axis.text=element_text(size=6))+
  theme(axis.text.x=element_text(size=6))+
  #theme(legend.position = "none" )+
  theme(panel.background = element_blank())+
  theme(plot.background=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("CTI")+ggtitle("b) Evapotranspiration")
  #coord_flip()
grid_cti = grid.arrange(grobs = list(plot_p_rmse_cti,plot_et_rmse_cti), layout_matrix=rbind(c(1,2)))
ggsave("Output/RMSE_CTI_New.png",grid_cti,width=8, height = 4, dpi=300)

# Land cover (Vegetation)
# change the land cover names accordingly
colnames(p_rmse_lc) = c("P_Dataset",
                        "Bare Areas",
                        "Shrubland",
                        "Sparse Vegetation",
                        "Grassland", 
                        "Cropland (Rainfed)",
                        "Mosaic (Herbaceous (>50) and Tree-Shrub (<50))",
                        "Mosaic (Herbaceous (<50) and Tree-Shrub (>50))",
                        "Tree Cover (Broad and Needle)",
                        "Herbaceous",
                        "Tree Cover (Broad and Deciduous (15-40))",
                        "Mosaic (Crop (>50) and Natural (<50)",
                        "Tree Cover (Broad and Deciduous (>15))",
                        "Shrub or Herbaceous (Flooded)",
                        "Tree Cover (Flooded)",
                        "Water Bodies",
                        "Cropland (Irrigated)",
                        "Tree Cover (Broad and Evergreen (>15))",
                        "Tree Cover (Flooded)",
                        "Herbaceous (<15)",
                        "Mosaic (Natural (>50) and Cropland (<50)",
                        "Tree or Shrub",
                        "Permanent snow and ice",
                        "Tree cover (Needle and Evergreen (>40)",
                        "Tree cover (Needle and Evergreen (>15)",
                        "Lichens and mosses",
                        "Urban",
                        "Tree Cover (Broad and Deciduous (>40)",
                        "Deciduous Shrubland")
colnames(et_rmse_lc) = c("ET_Dataset",
                        "Bare Areas",
                        "Shrubland",
                        "Sparse Vegetation",
                        "Grassland", 
                        "Cropland (Rainfed)",
                        "Mosaic (Herbaceous (>50) and Tree-Shrub (<50))",
                        "Mosaic (Herbaceous (<50) and Tree-Shrub (>50))",
                        "Tree Cover (Broad and Needle)",
                        "Herbaceous",
                        "Tree Cover (Broad and Deciduous (15-40))",
                        "Mosaic (Crop (>50) and Natural (<50)",
                        "Tree Cover (Broad and Deciduous (>15))",
                        "Shrub or Herbaceous (Flooded)",
                        "Tree Cover (Flooded)",
                        "Water Bodies",
                        "Cropland (Irrigated)",
                        "Tree Cover (Broad and Evergreen (>15))",
                        "Tree Cover (Flooded)",
                        "Herbaceous (<15)",
                        "Mosaic (Natural (>50) and Cropland (<50)",
                        "Tree or Shrub",
                        "Permanent snow and ice",
                        "Tree cover (Needle and Evergreen (>40)",
                        "Tree cover (Needle and Evergreen (>15)",
                        "Lichens and mosses",
                        "Urban",
                        "Tree Cover (Broad and Deciduous (>40)",
                        "Deciduous Shrubland")

p_temp_rmse_lc = melt(data = p_rmse_lc,id.vars = "P_Dataset")
p_temp_rmse_lc = cbind(p_temp_rmse_lc,as.data.frame("Precipitation"))
colnames(p_temp_rmse_lc) = c("Dataset","Classification","RMSE","Variable")
et_temp_rmse_lc = melt(data = et_rmse_lc,id.vars = "ET_Dataset")
et_temp_rmse_lc= cbind(et_temp_rmse_lc,as.data.frame("Evapotranspiration"))
colnames(et_temp_rmse_lc) = c("Dataset","Classification","RMSE","Variable")
data_rmse_lc = rbind(p_temp_rmse_lc,et_temp_rmse_lc)
#remove(p_temp_rmse_aridity,et_temp_rmse_aridity)

# plot the aridity heat map
plot_p_rmse_lc = ggplot(data = p_temp_rmse_lc, aes(x=Classification, y=Dataset,fill=RMSE))+
  geom_tile(color="white") +
  scale_fill_viridis()+
  geom_point(aes(size=ifelse(RMSE==0.0, "dot","no_dot")),color="red",shape=8) +
  scale_size_manual(values=c(dot=2, no_dot=NA), guide="none")+
  theme_tufte(base_family="Helvetica")+
  theme(axis.ticks=element_blank())+
  theme(axis.text=element_text(size=6))+
  theme(axis.text.x=element_text(size=6))+
  theme(text = element_text(family="Times"))+
  theme(panel.background = element_blank())+
  theme(plot.background=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("Land Cover")+ggtitle("a) Precipitation")
 # coord_flip()

plot_et_rmse_lc = ggplot(data = et_temp_rmse_lc, aes(x=Classification, y=Dataset,fill=RMSE))+
  geom_tile(color="white") +
  scale_fill_viridis()+
  geom_point(aes(size=ifelse(RMSE==0.0, "dot","no_dot")),color="red",shape=8) +
  scale_size_manual(values=c(dot=2, no_dot=NA), guide="none")+
  theme_tufte(base_family="Helvetica")+
  theme(axis.ticks=element_blank())+
  theme(axis.text=element_text(size=6))+
  theme(axis.text.x=element_text(size=6))+
  theme(text = element_text(family="Times"))+
  #theme(legend.position = "none" )+
  theme(panel.background = element_blank())+
  theme(plot.background=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("Land Cover")+ggtitle("b) Evapotranspiration")
 # coord_flip()
grid_lc = grid.arrange(grobs = list(plot_p_rmse_lc,plot_et_rmse_lc), layout_matrix=rbind(c(1,2)))
ggsave("Output/RMSE_LandCover_New.png",grid_lc,width=10, height = 4.5, dpi=300)


# Ridge plots
# prepare the distance values for ridge plots in ggplot
# rescale the p_distance_all metric for each catchment
p_distance_all_rescaled = apply(p_distance_all,1,rescale,na.rm=TRUE)
p_distance_all_rescaled = t(p_distance_all_rescaled)
p_data_ridge = melt(p_distance_all_rescaled)
p_data_ridge = p_data_ridge[,-1]
p_data_ridge[which(p_data_ridge[,2]>25),2] = 25
colnames(p_data_ridge) = c("Precipitation","Distance")
#plot_p_ridge = ggplot(p_data_ridge, aes(x=Distance,y=Precipitation,fill=0.5-abs(0.5-..ecdf..)))+
#               stat_density_ridges(geom="density_ridges_gradient",calc_ecdf=TRUE)+
#               scale_fill_viridis(name="Probability",direction=-1)
 
#plot_p_ridge = ggplot(p_data_ridge, aes(x=Distance,y=Precipitation))+
#               geom_density_ridges(scale=1)
#ggsave("P_Ridge.png",plot_p_ridge,width=5, height = 12, dpi=300)   


p_data_ridge$Precipitation = factor(p_data_ridge$Precipitation, 
                                    levels = c("CHIRPSv2.0", "CMORPHv0.x.RAW","","PERSIANN", "PERSIANN.CCS", "PERSIANN.CDR",
                                               "TRMM.3B42RT", "TRMM.3B43"))
plot_p_bar = ggplot()+
             geom_boxplot(data = p_data_ridge, aes(x=Precipitation,y=Distance,fill=Precipitation))+
             scale_fill_brewer(type="qual",palette = "Paired")+
             theme(axis.line = element_line(color="black"))+
             theme(strip.background = element_blank())+
             theme(panel.background = element_rect(fill="white"))+
             theme(panel.grid.major = element_line(colour="lightgrey"),
             panel.grid.minor = element_line(colour="lightgrey"),
             panel.border = element_rect(colour="black",fill=NA))+
             theme(axis.text.x = element_blank())
             #theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("Output/P_BoxPlot_New.png",plot_p_bar,width=8, height = 4, dpi=300)

# rescale the p_distance_all metric for each catchment
et_distance_all_rescaled = apply(et_distance_all,1,rescale,na.rm=TRUE)
et_distance_all_rescaled = t(et_distance_all_rescaled)
et_data_ridge = melt(et_distance_all_rescaled)
et_data_ridge = et_data_ridge[,-1]
#p_data_ridge[which(p_data_ridge[,2]>25),2] = 25
colnames(et_data_ridge) = c("Evapotranspiration","Distance")
#plot_p_ridge = ggplot(p_data_ridge, aes(x=Distance,y=Precipitation,fill=0.5-abs(0.5-..ecdf..)))+
#               stat_density_ridges(geom="density_ridges_gradient",calc_ecdf=TRUE)+
#               scale_fill_viridis(name="Probability",direction=-1)

#plot_p_ridge = ggplot(p_data_ridge, aes(x=Distance,y=Precipitation))+
#               geom_density_ridges(scale=1)
#ggsave("P_Ridge.png",plot_p_ridge,width=5, height = 12, dpi=300)   

# bar plot
plot_et_bar = ggplot()+
  geom_boxplot(data = et_data_ridge, aes(x=Evapotranspiration,y=Distance,fill=Evapotranspiration))+
  scale_fill_brewer(type="qual",palette = "Paired")+
  theme(axis.line = element_line(color="black"))+
  theme(strip.background = element_blank())+
  theme(panel.background = element_rect(fill="white"))+
  theme(panel.grid.major = element_line(colour="lightgrey"),
        panel.grid.minor = element_line(colour="lightgrey"),
        panel.border = element_rect(colour="black",fill=NA))+
  theme(axis.text.x = element_blank())
  #theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("Output/ET_BoxPlot_New.png",plot_et_bar,width=8, height = 4, dpi=300)
}

flag=0
if (flag>0){
# plot the maps of different variables 
# CTI
#colnames(p_rmse_cti) = c("Metric", "< 9.3", "9.3-10.2","10.2-11.2",
#                         "11.2-12.1", "12.1-13.1", "13.1-47.2")
#cti_map = rep(NA, length(shp))
#cti_map[index_cti1] = colnames(p_rmse_cti)[2]
#cti_map[index_cti2] = colnames(p_rmse_cti)[3]
#cti_map[index_cti3] = colnames(p_rmse_cti)[4]
#cti_map[index_cti4] = colnames(p_rmse_cti)[5]
#cti_map[index_cti5] = colnames(p_rmse_cti)[6]
#cti_map[index_cti6] = colnames(p_rmse_cti)[7]
#shp$CTI = cti_map

#shp@data$id = rownames(shp@data)
#shp_df = fortify(shp)
#shp_df = join(shp_df,shp@data,by="id")

# reorder factors
#shp_df$CTI = factor(shp_df$CTI, levels = rev(colnames(p_rmse_cti)[2:7]))

#plot_cti_map = ggplot()+
#  geom_polygon(data=shp_df, aes(x=long,y=lat,group=group,fill=factor(CTI)),colour="black",size=0.1/2)+
#  #scale_fill_manual(values = c("#d11d1d","#b2642c","#37842a","#3faec1","#1f5b93"))+
#  #scale_fill_manual(values = c("#d11d1d","#b2642c","#37842a","#3faec1","#1f5b93"))+
#  scale_fill_brewer(breaks=as.character(levels(shp_df$CTI)),type="qual",palette = "Paired",direction = -1)+
#  theme(axis.line = element_line(color="black"))+
#  theme(axis.text = element_text(size=12))+
#  theme(strip.background = element_blank())+
#  theme(legend.position = "bottom")+
#  theme(legend.text = element_text(size=12))+
#  theme(legend.title = element_text(face="bold"))+
#  theme(text = element_text(family="Times", size= 12))+
#  theme(panel.background = element_rect(fill=alpha("lightgrey",0.3/2)))+
#  theme(panel.grid.major = element_line(colour="black",size=0.1/2),
#        panel.grid.minor = element_line(colour="black",size=0.1/2),
#        panel.border = element_rect(colour="black",fill=NA))+
#  guides(fill=guide_legend(title="CTI:", face="bold", nrow=1, byrow=TRUE))+
#labs(x="Longitude", y="Latitude")
#ggsave(plot_cti_map,filename = "/home/akash/Documents/FifthPaper/Budyko_Analysis_New/Graphs_For_Manuscript_Revised/Figure_1/HydroBasins_Lvl05_CTI_Study_Area_Final.png", height=115, width=190,dpi=300,units="mm")
#remove(shp_df)

# NDVI
#colnames(p_rmse_ndvi) = c("Metric", "< 0.09", "0.09-0.20","0.20-0.38",
#                          "0.38-0.55", "0.55-0.67", "0.67-0.85")
#ndvi_map = rep(NA, length(shp))
#ndvi_map[index_ndvi1] = colnames(p_rmse_ndvi)[2]
#ndvi_map[index_ndvi2] = colnames(p_rmse_ndvi)[3]
#ndvi_map[index_ndvi3] = colnames(p_rmse_ndvi)[4]
#ndvi_map[index_ndvi4] = colnames(p_rmse_ndvi)[5]
#ndvi_map[index_ndvi5] = colnames(p_rmse_ndvi)[6]
#ndvi_map[index_ndvi6] = colnames(p_rmse_ndvi)[7]
#shp$NDVI = ndvi_map

#shp@data$id = rownames(shp@data)
#shp_df = fortify(shp)
#shp_df = join(shp_df,shp@data,by="id")

# reorder factors
#shp_df$NDVI = factor(shp_df$NDVI, levels = rev(colnames(p_rmse_ndvi)[2:7]))

#plot_ndvi_map = ggplot()+
#  geom_polygon(data=shp_df, aes(x=long,y=lat,group=group,fill=factor(NDVI)),colour="black",size=0.1/2)+
#  #scale_fill_manual(values = c("#d11d1d","#b2642c","#37842a","#3faec1","#1f5b93"))+
#  scale_fill_brewer(breaks=as.character(levels(shp_df$NDVI)),type="qual",palette = "Paired", direction = -1)+
#    theme(axis.line = element_line(color="black"))+
#    theme(strip.background = element_blank())+
#    theme(legend.position = "bottom")+
#    theme(legend.text = element_text(size=12))+
#    theme(legend.title = element_text(face="bold"))+
#    theme(axis.text = element_text(size=12))+
#    theme(text = element_text(family="Times", size= 12))+
#    theme(panel.background = element_rect(fill=alpha("lightgrey",0.3/2)))+
#    theme(panel.grid.major = element_line(colour="black",size=0.1/2),
#          panel.grid.minor = element_line(colour="black",size=0.1/2),
#          panel.border = element_rect(colour="black",fill=NA))+
#  guides(fill=guide_legend(title="NDVI:", nrow=1, byrow=TRUE))+
#  labs(x="Longitude", y="Latitude")
#ggsave(plot_ndvi_map,filename = "/home/akash/Documents/FifthPaper/Budyko_Analysis_New/Graphs_For_Manuscript_Revised/Figure_1/HydroBasins_Lvl05_NDVI_Study_Area_Final.png", height=115, width=190,dpi=300,units="mm")
#remove(shp_df)


# Land Cover
shp$Land_Cover = landcover

shp@data$id = rownames(shp@data)
shp_df = fortify(shp)
shp_df = join(shp_df,shp@data,by="id")

getPalette = colorRampPalette(brewer.pal(9, "RdYlGn"))
plot_landcover_map = ggplot()+
  geom_polygon(data=shp_df, aes(x=long,y=lat,group=group,fill=factor(Land_Cover)),colour="black",size=0.1/2)+
#scale_fill_manual(values = c("#d11d1d","#b2642c","#37842a","#3faec1","#1f5b93"))+
  #scale_fill_brewer(type="qual",palette = "Paired")+
  scale_fill_manual(values=getPalette(29))+
  theme(axis.line = element_line(color="black"))+
  theme(strip.background = element_blank())+
  theme(panel.background = element_rect(fill=alpha("lightgrey",0.3/2)))+
  theme(panel.grid.major = element_line(colour="black",size=0.1/2),
        panel.grid.minor = element_line(colour="black",size=0.1/2),
        panel.border = element_rect(colour="black",fill=NA),
        legend.position="right")+
  theme(legend.text = element_text(size=8))+
  theme(legend.title = element_text(face="bold"))+
  theme(axis.text = element_text(size=12))+
  theme(text = element_text(family="Times", size= 12))+
  guides(fill=guide_legend(title="Land Cover"))+
  labs(x="Longitude", y="Latitude")
ggsave(plot_landcover_map,filename = "Graphs_For_Manuscript_Revised/Supplementary/LandCover_Study_Area_Legend.png", height=115, width=210,dpi=300,units="mm")
#remove(shp_df)



# Elevation
#colnames(p_rmse_elevation) = c("Metric", "< 59", "59-155","155-360",
#                               "360-748", "748-1286", "1286-5300")
#elev_map = rep(NA, length(shp))
#elev_map[index_elev1] = colnames(p_rmse_elevation)[2]
#elev_map[index_elev2] = colnames(p_rmse_elevation)[3]
#elev_map[index_elev3] = colnames(p_rmse_elevation)[4]
#elev_map[index_elev4] = colnames(p_rmse_elevation)[5]
#elev_map[index_elev5] = colnames(p_rmse_elevation)[6]
#elev_map[index_elev6] = colnames(p_rmse_elevation)[7]
#shp$ELEVATION = elev_map

#shp@data$id = rownames(shp@data)
#shp_df = fortify(shp)
#shp_df = join(shp_df,shp@data,by="id")

# reorder factors
#shp_df$ELEVATION = factor(shp_df$ELEVATION, levels = rev(colnames(p_rmse_elevation)[2:7]))

#plot_elevation_map = ggplot()+
#  geom_polygon(data=shp_df, aes(x=long,y=lat,group=group,fill=factor(ELEVATION)),colour="black",size=0.1/2)+
  #scale_fill_manual(values = c("#d11d1d","#b2642c","#37842a","#3faec1","#1f5b93"))+
#  scale_fill_brewer(breaks=as.character(levels(shp_df$ELEVATION)),type="qual",palette = "Paired", direction = -1)+
#  theme(axis.line = element_line(color="black"))+
#  theme(strip.background = element_blank())+
#  theme(legend.position = "bottom")+
#  theme(legend.text = element_text(size=12))+
# theme(legend.title = element_text(face="bold"))+
#  theme(axis.text = element_text(size=12))+
#  theme(text = element_text(family="Times", size= 12))+
#  theme(panel.background = element_rect(fill=alpha("lightgrey",0.3/2)))+
#  theme(panel.grid.major = element_line(colour="black",size=0.1/2),
#        panel.grid.minor = element_line(colour="black",size=0.1/2),
#        panel.border = element_rect(colour="black",fill=NA))+
#  guides(fill=guide_legend(title="Elevation (in m):", nrow=1, byrow=TRUE))+
#  labs(x="Longitude", y="Latitude")
#ggsave(plot_elevation_map,filename = "/home/akash/Documents/FifthPaper/Budyko_Analysis_New/Graphs_For_Manuscript_Revised/Figure_1/HydroBasins_Lvl05_Elevation_Study_Area_Final.png", height=115, width=190,dpi=300,units="mm")
#remove(shp_df)

}








