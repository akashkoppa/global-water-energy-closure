# Script to perform sensitivity analysis for the Budyko results
# Author - Akash Koppa
# Date - 2019-08-27
# Description - The main objective will be to use sensitivty analysis to find which variable contributes most to the uncertainties and errors
# I should focus on two main things
# 1) Sources of Uncertainty - Where does the maximum uncertainty come from (P, ET, and ET0)
# 2) Quantifying the Absolute Uncertainty - What is the absolute magnitude of uncertainty in P, ET, and ET0 for each catchment?


# clear workspace
rm(list=ls())

# load required libraries
library(ggplot2)

# set working directory
setwd("/home/akash/Documents/FifthPaper/Budyko_Analysis_New/")

# load the EI and AI data
ei = read.table("Output/Evaporative_Index_Rn-CERES_RS_noSM2RAIN.txt",header=TRUE)
ai = read.table("Output/Aridity_Index_Rn-CERES_RS_noSM2RAIN.txt", header=TRUE)

# actual aridity 
ai_actual = ei$AI_Actual
hydrosheds_id = ei$HydroSHEDS_ID
area = ei$Area

ei = ei[,-c(1:3)]
ai = ai[,-c(1:3)]

# get the name of the combinations of P and ET datasets
p_et_dataset = colnames(ei)

# read in the long-term average p, et, and pet datasets
p = read.table("Input/Precipitation_Long_Term.txt", header=TRUE)
p = p[,-1]
p_dataset = colnames(p)
p = p[,-which(p_dataset %in% c("CPC.Unifiedv1.0","CRU.TSv4.03","ERA5.Land","ERA5","GPCCv7.0","PREC.Land","UDELv5.0","SM2RAIN.CCI"))]
p_dataset = p_dataset[-which(p_dataset %in% c("CPC.Unifiedv1.0","CRU.TSv4.03","ERA5.Land","ERA5","GPCCv7.0","PREC.Land","UDELv5.0","SM2RAIN.CCI"))]


et = read.table("Input/Evapotranspiration_Long_Term.txt", header=TRUE)
et = et[,-1]
et_dataset = colnames(et)
et = et[,-which(et_dataset %in% c("ERA5.Land","FLDAS","FluxCom.RSM","GLDASv2.1"))]
et_dataset = et_dataset[-which(et_dataset %in% c("ERA5.Land","FLDAS","FluxCom.RSM","GLDASv2.1"))]


# Sources of Uncertainty (in the form of cluster density)
# P (one ET and ensemble P)
cluster_density_p = NULL
for (k in 1:length(et_dataset)){
  ei_temp1 = ei[,grep(pattern = et_dataset[k],p_et_dataset)]
  ai_temp1 = ai[,grep(pattern = et_dataset[k],p_et_dataset)]
  
  # calculate the centroid of ei and ai
  ei_centroid = apply(ei_temp1,1,mean,na.rm=TRUE)
  ai_centroid = apply(ai_temp1,1,mean,na.rm=TRUE)
  
  cluster_distance = matrix(data = NA, nrow=nrow(ei_temp1),ncol = ncol(ei_temp1))
  for(i in 1:nrow(cluster_distance)){
    ei_centre_temp = ei_centroid[i]
    ai_centre_temp = ai_centroid[i]
    for(j in 1:ncol(cluster_distance)){
      ei_temp = ei_temp1[i,j]
      ai_temp = ai_temp1[i,j]
      dist_temp = sqrt(((ei_temp - ei_centre_temp)^2) + ((ai_temp - ai_centre_temp)^2))
      cluster_distance[i,j] = dist_temp
    }
  }
  mean_radius_cluster = apply(cluster_distance,1,mean,na.rm=TRUE)
  mean_radius_cluster[is.infinite(mean_radius_cluster)] = NA
  #density_mean_radius = ncol(cluster_distance)/(pi*(mean_radius_cluster^2))
  density_mean_radius = mean_radius_cluster # hack by Akash for changing from density to radius
  density_mean_radius[which(is.nan(density_mean_radius))] = NA
  cluster_density_p = cbind(cluster_density_p,density_mean_radius)
}
colnames(cluster_density_p) = et_dataset

# calculate mean and sd of cluster_density sensitivities
cluster_density_p_mean = apply(cluster_density_p,1,mean,na.rm=TRUE)
cluster_density_p_sd   = apply(cluster_density_p,1,sd,na.rm=TRUE)

# ET (one P and ensemble ET)
cluster_density_et = NULL
for (k in 1:length(p_dataset)){
  ei_temp1 = ei[,grep(pattern = p_dataset[k],p_et_dataset)]
  ai_temp1 = ai[,grep(pattern = p_dataset[k],p_et_dataset)]
  
  # calculate the centroid of ei and ai
  ei_centroid = apply(ei_temp1,1,mean,na.rm=TRUE)
  ai_centroid = apply(ai_temp1,1,mean,na.rm=TRUE)
  
  cluster_distance = matrix(data = NA, nrow=nrow(ei_temp1),ncol = ncol(ei_temp1))
  for(i in 1:nrow(cluster_distance)){
    ei_centre_temp = ei_centroid[i]
    ai_centre_temp = ai_centroid[i]
    for(j in 1:ncol(cluster_distance)){
      ei_temp = ei_temp1[i,j]
      ai_temp = ai_temp1[i,j]
      dist_temp = sqrt(((ei_temp - ei_centre_temp)^2) + ((ai_temp - ai_centre_temp)^2))
      cluster_distance[i,j] = dist_temp
    }
  }
  mean_radius_cluster = apply(cluster_distance,1,mean,na.rm=TRUE)
  mean_radius_cluster[is.infinite(mean_radius_cluster)] = NA
  #density_mean_radius = ncol(cluster_distance)/(pi*(mean_radius_cluster^2))
  density_mean_radius = mean_radius_cluster
  density_mean_radius[which(is.nan(density_mean_radius))] = NA
  cluster_density_et = cbind(cluster_density_et,density_mean_radius)
}
colnames(cluster_density_et) = p_dataset

# calculate mean and sd of cluster_density sensitivities
cluster_density_et_mean = apply(cluster_density_et,1,mean,na.rm=TRUE)
cluster_density_et_sd   = apply(cluster_density_et,1,sd,na.rm=TRUE)

# Write the results to a data table
final_sensitivity = data.frame(CD_P_Mean  = cluster_density_p_mean,
                               CD_P_SD    = cluster_density_p_sd,
                               CD_ET_Mean = cluster_density_et_mean,
                               CD_ET_SD   = cluster_density_et_sd)
write.table(final_sensitivity, file="Output/Budyko_Cluster_Density_Sensitivity_noSM2RAIN.txt",row.names=FALSE, col.names = TRUE)

# Absolute Uncertainty in P, ET, and ET0
# The uncertainty here is fairly straight forward - All the years considered will be averaged
# P
#p_all = NULL
#for (i in 1:length(p_dataset)){
#  p = read.table(paste0("Input/P-",p_dataset[i],"/P-",p_dataset[i],".txt"),header = TRUE)
#  p = p[,-1]
#  p_all = cbind(p_all, apply(p,1,mean,na.rm=TRUE))
#}

# calculate mean and standard deviation
p_all_mean = apply(p,1,mean,na.rm=TRUE)
p_all_sd   = apply(p,1,sd,na.rm=TRUE)

# ET
#et_all = NULL
#for (i in 1:length(et_dataset)){
#  et = read.table(paste0("Input/ET-",et_dataset[i],"/ET-",et_dataset[i],".txt"),header = TRUE)
#  et = et[,-1]
#  et_all = cbind(et_all, apply(et,1,mean,na.rm=TRUE))
#}

# calculate mean and standard deviation
et_all_mean = apply(et,1,mean,na.rm=TRUE)
et_all_sd   = apply(et,1,sd,na.rm=TRUE)

# write the results to a table
final_abs_sensitivity = data.frame(P_Mean  = p_all_mean,
                                   P_SD    = p_all_sd,
                                   ET_Mean = et_all_mean,
                                   ET_SD   = et_all_sd)
write.table(final_abs_sensitivity, file="Output/Absolute_P_ET_Uncertainty_noSM2RAIN.txt",row.names=FALSE, col.names = TRUE)



