# Script to analyze  clusters of P and ET datasets in the Budyko space
# Author - Akash Koppa
# Date - 2019-06-27

# clear workspace
rm(list=ls())

# load required libraries
library(ggplot2)

# set working directory
setwd("/home/akash/Documents/FifthPaper/Budyko_Analysis_New/")

# load the EI and AI data
ei = read.table("Output/Evaporative_Index_Rn-CERES_RS_noSM2RAIN.txt",header=TRUE)
ai = read.table("Output/Aridity_Index_Rn-CERES_RS_noSM2RAIN.txt", header=TRUE)

# subset the problem catchments in Amazonia
#catch_reqd = c(6050207670, 6050156880, 6050207660, 6050287830, 6050266740, 6050266880, 6050294360, 6050294270, 6050247260, 6050247270)
#index_reqd = which(ei$HydroSHEDS_ID %in% catch_reqd)

# subset ei and ai data for the Amazon region
#ei_reqd = ei[index_reqd,]
#ai_reqd = ei[index_reqd,]

# subset a few catchments with higher uncertainty in Amazonia
#catch_reqd_high = c(6050315910, 6050315920, 6050029730)
#index_reqd_high = which(ei$HydroSHEDS_ID %in% catch_reqd_high)
#ei_reqd_high = ei[index_reqd_high,]
#ai_reqd_high = ei[index_reqd_high,]

# actual aridity 
ai_actual = ei$AI_Actual
hydrosheds_id = ei$HydroSHEDS_ID
area = ei$Area

ei = ei[,-c(1:3)]
ai = ai[,-c(1:3)]

# remove inf values
ei[ei==Inf] = NA
ei[ei==-Inf] = NA
ai[ai==Inf] = NA
ai[ai==-Inf] = NA

# plot the number of EI and AI for each catchment (basically a world map)
number_ei = NULL
number_ai = NULL
for (i in 1:nrow(ai)){
  number_ei = c(number_ei, ncol(ei)-length(which(is.na(ei[i,]))))
  number_ai = c(number_ai, ncol(ai)-length(which(is.na(ai[i,]))))
}


# calculate the centroid of ei and ai
ei_centroid = apply(ei,1,mean,na.rm=TRUE)
ai_centroid = apply(ai,1,mean,na.rm=TRUE)


# perform the clustering analysis
# For each catchment calculate the radius and density of the clusters
# calculate the euclidean distance from the centroid to each point
cluster_distance = matrix(data = NA, nrow=nrow(ei),ncol = ncol(ei))
for(i in 1:nrow(cluster_distance)){
  ei_centre_temp = ei_centroid[i]
  ai_centre_temp = ai_centroid[i]
  for(j in 1:ncol(cluster_distance)){
    ei_temp = ei[i,j]
    ai_temp = ai[i,j]
    dist_temp = sqrt(((ei_temp - ei_centre_temp)^2) + ((ai_temp - ai_centre_temp)^2))
    cluster_distance[i,j] = dist_temp
  }
}

# plot the number of EI and AI for each catchment (basically a world map)
number_cd= NULL
for (i in 1:nrow(ai)){
  number_cd = c(number_cd, ncol(cluster_distance)-length(which(is.na(cluster_distance[i,]))))
}

#cluster_distance[is.infinite(cluster_distance)] = NA
# calculate maximum distance
max_radius_cluster = apply(cluster_distance,1,max,na.rm=TRUE)
max_radius_cluster[is.infinite(max_radius_cluster)] = NA
mean_radius_cluster = apply(cluster_distance,1,mean,na.rm=TRUE)
mean_radius_cluster[is.infinite(mean_radius_cluster)] = NA
min_radius_cluster = apply(cluster_distance,1,min,na.rm=TRUE)
min_radius_cluster[is.infinite(min_radius_cluster)] = NA
std_radius_cluster = apply(cluster_distance,1,sd,na.rm=TRUE)
std_radius_cluster[is.infinite(std_radius_cluster)] = NA

# calculate density (points per unit area)
density_max_radius = ncol(cluster_distance)/(pi*(max_radius_cluster^2))
density_mean_radius = ncol(cluster_distance)/(pi*(mean_radius_cluster^2))

# write all the cluster metrics to a data frame
data_cluster_final = data.frame(HydroSHEDS_ID = hydrosheds_id,
                                AI = ai_actual,
                                Area = area,
                                Max_Radius = max_radius_cluster,
                                Mean_Radius = mean_radius_cluster,
                                Min_Radius = min_radius_cluster,
                                Sd_Radius = std_radius_cluster,
                                Density_Max_Radius = density_max_radius,
                                Density_Mean_Radius = density_mean_radius)
write.table(data_cluster_final,file = "Output/Budyko_Cluster_Metrics_RS_noSM2RAIN.txt",row.names=FALSE,quote = FALSE)
