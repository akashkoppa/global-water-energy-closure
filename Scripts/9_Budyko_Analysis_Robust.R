# Script to analyze the distance metrics based on the Budyko hypothesis
# Author - Akash Koppa 
# Date - 2019/04/016
# Description - This will be a preliminary script for analyzing the distance metrics

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

# set working directory
setwd("/home/akash/Documents/FifthPaper/Budyko_Analysis_New/")

# read in the shapefile
shp = shapefile("Input/HydroSHEDS_All_Level05/hybas_all_lev05_v1c.shp")

# read in the distance metric results
dist_metric = read.table("Output/Budyko_Distance_HydroSHEDS_Lvl05_Rn-CERES_noSM2RAIN.txt",header=TRUE)

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
mean_all = apply(dist_metric,2,mean,na.rm=TRUE)
std_all = sqrt(apply(dist_metric,2,var,na.rm=TRUE))
rmse_all = sqrt(apply(dist_metric^2,2,mean,na.rm=TRUE))

# decompose the statistics according to the aridity indices
# Hyper Arid
mean_hyperarid = apply(dist_metric[index_hyperarid,],2,mean,na.rm=TRUE)
std_hyperarid = sqrt(apply(dist_metric[index_hyperarid,],2,var,na.rm=TRUE))
rmse_hyperarid = sqrt(apply(dist_metric[index_hyperarid,]^2,2,mean,na.rm=TRUE))
# Arid
mean_arid = apply(dist_metric[index_arid,],2,mean,na.rm=TRUE)
std_arid = sqrt(apply(dist_metric[index_arid,],2,var,na.rm=TRUE))
rmse_arid = sqrt(apply(dist_metric[index_arid,]^2,2,mean,na.rm=TRUE))
# Semi Arid
mean_semiarid = apply(dist_metric[index_semiarid,],2,mean,na.rm=TRUE)
std_semiarid = sqrt(apply(dist_metric[index_semiarid,],2,var,na.rm=TRUE))
rmse_semiarid = sqrt(apply(dist_metric[index_semiarid,]^2,2,mean,na.rm=TRUE))
# Dry Subhumid
mean_drysubhumid = apply(dist_metric[index_drysubhumid,],2,mean,na.rm=TRUE)
std_drysubhumid = sqrt(apply(dist_metric[index_drysubhumid,],2,var,na.rm=TRUE))
rmse_drysubhumid = sqrt(apply(dist_metric[index_drysubhumid,]^2,2,mean,na.rm=TRUE))
# Humid
mean_humid = apply(dist_metric[index_humid,],2,mean,na.rm=TRUE)
std_humid = sqrt(apply(dist_metric[index_humid,],2,var,na.rm=TRUE))
rmse_humid = sqrt(apply(dist_metric[index_humid,]^2,2,mean,na.rm=TRUE))

# decompose the statistics according to continents
# North America
mean_na = apply(dist_metric[index_na,],2,mean,na.rm=TRUE)
std_na = sqrt(apply(dist_metric[index_na,],2,var,na.rm=TRUE))
rmse_na = sqrt(apply(dist_metric[index_na,]^2,2,mean,na.rm=TRUE))

# South America
mean_sa = apply(dist_metric[index_sa,],2,mean,na.rm=TRUE)
std_sa = sqrt(apply(dist_metric[index_sa,],2,var,na.rm=TRUE))
rmse_sa = sqrt(apply(dist_metric[index_sa,]^2,2,mean,na.rm=TRUE))

# Europe and Middle East
mean_eu = apply(dist_metric[index_eu,],2,mean,na.rm=TRUE)
std_eu = sqrt(apply(dist_metric[index_eu,],2,var,na.rm=TRUE))
rmse_eu = sqrt(apply(dist_metric[index_eu,]^2,2,mean,na.rm=TRUE))

# Africa
mean_af = apply(dist_metric[index_af,],2,mean,na.rm=TRUE)
std_af = sqrt(apply(dist_metric[index_af,],2,var,na.rm=TRUE))
rmse_af = sqrt(apply(dist_metric[index_af,]^2,2,mean,na.rm=TRUE))

# Asia
mean_as = apply(dist_metric[index_as,],2,mean,na.rm=TRUE)
std_as = sqrt(apply(dist_metric[index_as,],2,var,na.rm=TRUE))
rmse_as = sqrt(apply(dist_metric[index_as,]^2,2,mean,na.rm=TRUE))

# Australia
mean_au = apply(dist_metric[index_au,],2,mean,na.rm=TRUE)
std_au = sqrt(apply(dist_metric[index_au,],2,var,na.rm=TRUE))
rmse_au = sqrt(apply(dist_metric[index_au,]^2,2,mean,na.rm=TRUE))

# Arctic
mean_ar = apply(dist_metric[index_ar,],2,mean,na.rm=TRUE)
std_ar = sqrt(apply(dist_metric[index_ar,],2,var,na.rm=TRUE))
rmse_ar = sqrt(apply(dist_metric[index_ar,]^2,2,mean,na.rm=TRUE))

# Greenland
mean_gr = apply(dist_metric[index_gr,],2,mean,na.rm=TRUE)
std_gr = sqrt(apply(dist_metric[index_gr,],2,var,na.rm=TRUE))
rmse_gr = sqrt(apply(dist_metric[index_gr,]^2,2,mean,na.rm=TRUE))

# Siberia
mean_si = apply(dist_metric[index_si,],2,mean,na.rm=TRUE)
std_si = sqrt(apply(dist_metric[index_si,],2,var,na.rm=TRUE))
rmse_si = sqrt(apply(dist_metric[index_si,]^2,2,mean,na.rm=TRUE))

# create a final data frame of summary metrics 
summary_metric_dist = data.frame(        Mean_all = mean_all,
                                         SD_all = std_all,
                                         RMSE_all = rmse_all,
                                         Mean_hyperarid = mean_hyperarid,
                                         SD_hyperarid = std_hyperarid,
                                         RMSE_hyperarid = rmse_hyperarid,
                                         Mean_arid = mean_arid,
                                         SD_arid = std_arid,
                                         RMSE_arid = rmse_arid,
                                         Mean_semiarid = mean_semiarid,
                                         SD_semiarid = std_semiarid,
                                         RMSE_semiarid = rmse_semiarid,
                                         Mean_drysubhumid = mean_drysubhumid,
                                         SD_drysubhumid = std_drysubhumid,
                                         RMSE_drysubhumid = rmse_drysubhumid,
                                         Mean_humid = mean_humid,
                                         SD_humid = std_humid,
                                         RMSE_humid = rmse_humid,
                                         Mean_na = mean_na,
                                         SD_na = std_na,
                                         RMSE_na = rmse_na,
                                         Mean_sa = mean_sa,
                                         SD_sa = std_sa,
                                         RMSE_sa = rmse_sa,
                                         Mean_eu = mean_eu,
                                         SD_eu = std_eu,
                                         RMSE_eu = rmse_eu,
                                         Mean_af = mean_af,
                                         SD_af = std_af,
                                         RMSE_af = rmse_af,
                                         Mean_as = mean_as,
                                         SD_as = std_as,
                                         RMSE_as = rmse_as,
                                         Mean_au = mean_au,
                                         SD_au = std_au,
                                         RMSE_au = rmse_au,
                                         Mean_ar = mean_ar,
                                         SD_ar = std_ar,
                                         RMSE_ar = rmse_ar,
                                         Mean_gr = mean_gr,
                                         SD_gr = std_gr,
                                         RMSE_gr = rmse_gr,
                                         Mean_si = mean_si,
                                         SD_si = std_si,
                                         RMSE_si = rmse_si)
summary_metric_dist = round(summary_metric_dist,2)

# More detailed analyses - All combinations 
# scatter or line plots 
flag = 0
if (flag>1){
dist_metric_final = read.table("Output/Budyko_Distance_HydroSHEDS_Lvl05_Rn-CERES.txt",header=TRUE)

# order the matrix according to aridity index 
dist_metric_ai = dist_metric_final[order(dist_metric_final$AI),]
dist_metric_ai$HydroSHEDS_ID = c(1:nrow(dist_metric_ai))

plot_dist_ai = melt(dist_metric_ai,id.vars = c("HydroSHEDS_ID","AI","Area"))
p_ai = unlist(strsplit(as.character(plot_dist_ai$variable),"[.]"))[seq(from=1,to=length(plot_dist_ai$variable)*2,by=2)]
e_ai = unlist(strsplit(as.character(plot_dist_ai$variable),"[.]"))[seq(from=2,to=length(plot_dist_ai$variable)*2,by=2)]
plot_dist_ai = cbind(plot_dist_ai,p_ai,e_ai)
colnames(plot_dist_ai) = c("Hydrosheds","AI","Area","P.ET_Combination","Distance","P_Dataset","ET_Dataset")
plot_dist_ai$Distance[plot_dist_ai$Distance>100] = 100
#plot_dist_ai$AI[plot_dist_ai$AI > 75] = 75
#basin_ordered = as.character(dist_metric_ai$HydroSHEDS_ID)
#plot_dist_ai$Hydrosheds = factor(plot_dist_ai$Hydrosheds,levels=basin_ordered)

# plot line (or point) for each catchment
plot_dist_point = ggplot()+
                  geom_point(data = plot_dist_ai,aes(x=Hydrosheds, y=Distance, color=P_Dataset, shape=ET_Dataset, alpha=P_Dataset),size=1.5)+
                  scale_color_brewer(type="qual",palette = "Paired")+
                  scale_alpha_manual(values = c(1.0,0.9,0.85,0.8,0.75,0.7,0.65,0.6,0.55,0.4))+
                  scale_shape_manual(values = c(0,1,2,12,13,5,6,7,11,9))+
  theme(axis.line = element_line(color="black"))+
  theme(strip.background = element_blank())+
  theme(panel.background = element_rect(fill="white"))+
  theme(panel.grid.major = element_line(colour="lightgrey"),
        panel.grid.minor = element_line(colour="lightgrey"),
        panel.border = element_rect(colour="black",fill=NA))
ggsave(filename = "test_point_ai_20.png",plot = plot_dist_point, width = 15, height = 7, dpi = 200)

#plot_dist_line  = ggplot()+
#                  geom_path(data = plot_dist_ai,aes(x=Hydrosheds, y=Distance, color=P_Dataset, linetype=ET_Dataset))+
#                  scale_color_brewer(type="qual",palette = "Paired")
#ggsave(filename = "test_line.png",plot = plot_dist_line, width = 15, height = 8, dpi= 200)

# loop through the P datasets and plot the distances for one P and other ET dataset combinations
#plot_dist_ai_temp = plot_dist_ai
#for (i in 1:length(p_dataset)){
#  index_reqd = which(plot_dist_ai_temp$P_Dataset==p_dataset[i])
#  plot_dist_ai = plot_dist_ai_temp[index_reqd,]
  
#plot_dist_point = ggplot()+
#                  geom_point(data = plot_dist_ai,aes(x=Hydrosheds, y=Distance, color=ET_Dataset,shape=ET_Dataset),size=1.7)+
#                  scale_color_brewer(type="qual",palette = "Paired")+
#                  #scale_alpha_manual(values = c(1.0,0.9,0.85,0.8,0.75,0.7,0.65,0.6,0.55,0.4))+
#                  scale_shape_manual(values = c(0,1,2,12,13,5,6,7,11,9))+
#                  ggtitle(label = p_dataset[i])
#ggsave(filename = paste0("Output/Distance/Distance_All_Combinations_P-",p_dataset[i],".png"),plot = plot_dist_point, width = 10, height = 6, dpi = 200)
#}

# loop through the ET datasets and plot the distances for one ET and other P dataset combinations
#for (i in 1:length(et_dataset)){
#  index_reqd = which(plot_dist_ai_temp$ET_Dataset==et_dataset[i])
#  plot_dist_ai = plot_dist_ai_temp[index_reqd,]
  
#  plot_dist_point = ggplot()+
#    geom_point(data = plot_dist_ai,aes(x=Hydrosheds, y=Distance, color=P_Dataset,shape=P_Dataset),size=1.7)+
#    scale_color_brewer(type="qual",palette = "Paired")+
    #scale_alpha_manual(values = c(1.0,0.9,0.85,0.8,0.75,0.7,0.65,0.6,0.55,0.4))+
#    scale_shape_manual(values = c(0,1,2,12,13,5,6,7,11,9))+
#    ggtitle(label = et_dataset[i])
#  ggsave(filename = paste0("Output/Distance/Distance_All_Combinations_ET-",et_dataset[i],".png"),plot = plot_dist_point, width = 10, height = 6, dpi = 200)
#}

}

flag = 2
if (flag>1){
# Ranking analysis
# In contrast to the other ranking analysis, this is a more robust algorithm
# here the best precipitation and evapotranspiration datasets are selected based on the following steps:
# For each catchment and for each P or ET dataset calculate the distance metric for all other ET or P dataset respectively. 
# Then calculate the minimum of these distances to select the best P and ET dataset (for each catchment)
best_p  = NULL
best_et = NULL
missing_datasets = NULL

for (i in 1:nrow(dist_metric)){
  
  temp = dist_metric[i,]
  ai_temp = ai[i]
  # find the length of NA values in each row
  length_NA = length(which(is.na(temp)))
  missing_datasets = c(missing_datasets,length_NA)
  
  # precipitation
  mean_p  = NULL
  for (j in 1:length(p_dataset)){
    index_p = grep(p_dataset[j],p_et_dataset)
    temp_data = temp[,index_p]
    mean_p = c(mean_p,mean(unlist(c(temp_data)),na.rm=TRUE))
  }
  mean_p[which(is.nan(mean_p))] = NA
  index_min_p = which(mean_p == min(mean_p,na.rm =TRUE))[1]
  if(length(index_min_p)<1){
    best_p = c(best_p,NA)
  }else{
  best_p = c(best_p,p_dataset[index_min_p])
  }
  
  # evapotranspiration
  mean_et = NULL
  for (k in 1:length(et_dataset)){
    index_et = grep(et_dataset[k],p_et_dataset)
    temp_data = temp[,index_et]
    mean_et = c(mean_et,mean(unlist(c(temp_data)),na.rm=TRUE))
  }
  mean_et[which(is.nan(mean_et))] = NA
  index_min_et = which(mean_et == min(mean_et,na.rm=TRUE))[1]
  if(length(index_min_et)<1){
    best_et = c(best_et,NA)
  }else{
  best_et = c(best_et,et_dataset[index_min_et])
  }
}

# Plot the best P and ET dataset map
shp$BEST_P = best_p
shp$BEST_ET = best_et

shp@data$id = rownames(shp@data)
shp_df = fortify(shp)
shp_df = join(shp_df,shp@data,by="id")
shp_df$BEST_P = factor(shp_df$BEST_P)
shp_df$BEST_ET = factor(shp_df$BEST_ET)

plot_best_p = ggplot()+
              geom_polygon(data=shp_df, aes(x=long,y=lat,group=group,fill=factor(BEST_P)),colour="black",size=0.1/2)+
              scale_fill_brewer(breaks= levels(shp_df$BEST_P),type="qual",palette = "Paired")+
              theme(axis.line = element_line(color="black"))+
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
              guides(fill=guide_legend(title="P Datasets:", face="bold",nrow=2,byrow = TRUE))+
              #guides(fill=guide_legend(title="P Datasets"))+
              labs(x="Longitude",y="Latitude")
ggsave(plot_best_p,filename = "Graphs_For_Manuscript_Revised/Figure_6/Best_P_Global_noSM2RAIN_withBoundary.png", height=120, width=190,dpi=300, units="mm")

plot_best_et = ggplot()+
  geom_polygon(data=shp_df, aes(x=long,y=lat,group=group,fill=factor(BEST_ET)),colour="black",size=0.1/2)+
  scale_fill_brewer(breaks=levels(shp_df$BEST_ET),type="qual",palette = "Paired")+
  theme(axis.line = element_line(color="black"))+
  theme(strip.background = element_blank())+
 # theme(panel.background = element_rect(fill="white"))+
#  theme(panel.grid.major = element_line(colour="lightgrey"),
#        panel.grid.minor = element_line(colour="lightgrey"),
#        panel.border = element_rect(colour="black",fill=NA))+
  theme(panel.background = element_rect(fill=alpha("lightgrey",0.3/2)))+
  theme(panel.grid.major = element_line(colour="black",size=0.1/2),
        panel.grid.minor = element_line(colour="black",size=0.1/2),
        panel.border = element_rect(colour="black",fill=NA))+
  theme(legend.position = "bottom")+
  theme(text = element_text(family="Times", size= 12))+
  theme(legend.text = element_text(size=10))+
  theme(legend.title = element_text(face="bold"))+
  theme(axis.text = element_text(size=12))+
  guides(fill=guide_legend(title="ET Datasets:", face="bold",nrow=2,byrow = TRUE))+
  #guides(fill=guide_legend(title="ET Datasets"))+
  labs(x="Longitude",y="Latitude")
ggsave(plot_best_et,filename = "Graphs_For_Manuscript_Revised/Figure_6/Best_ET_Global_noSM2RAIN_withBoundary.png", height=120, width=190,dpi=300, units="mm")

}





