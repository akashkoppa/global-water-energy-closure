# A script to determine the distance metric for all the HydroSHEDS level 05 basins
# Author - Akash Koppa
# Date - 2019-03-31
# Description - 1) the script is as modular as possible. 
#               2) only the user defined settings are required to be changed
#               3) this involves adding more names to the p_data and et_data variables
# Note: 1) As the PET dataset from GLEAM stretches from 1998 to 2012, the current implementation only considers data from 1998 to 2015. 
#       2) This needs to be discussed as it is a long-term water-energy balance, keeping the years consistent across all combinations may not be needed.
# Updated (2019-05-13): 1) replace PET with CERES net radiation, 2) replace Global Aridity index with SRB and WorldClim data

# clear workspace
rm(list=ls())

# set working directory
setwd("/home/akash/Documents/FifthPaper/Budyko_Analysis_New/")

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

rn = read.table("Input/NetRadiation_Long_Term.txt", header=TRUE)
rn = rn[,-1]
rn_data = colnames(rn)
rn = data.frame(rn[,"CERES"])
rn_data = "CERES"

# load the Budyko parameter (w) data
w_temp      = read.table("Input/Budyko_Parameter_HydroBasinsLvl05_All.txt",header = TRUE)
w_param     = w_temp$W_Combined # currently considering the parameter derived from the combined equation in Xu et al 2013
area        = w_temp$Area
hyd_basins  = as.character(w_temp$HydroSHEDS_ID)
remove(w_temp)

# load the global aridity index data
ai = read.table("Input/Global_Aridity_Index_SRB.txt", header=TRUE)
ai = ai$AI

for (k in 1:length(rn_data)){
rn_temp = rn[,k]

# loop through different combinations of P and E datasets 
distance_final = NULL
combination = NULL # for column name of final distance dataset
  for (i in 1:length(p_data)){
    p_temp = p[,i]
  
    for (j in 1:length(et_data)){
      et_temp = et[,j]  
    
# determine ai_budyko, ei_budyko, ai_data, ei_data (ai_data and ei_data are from the datasets under consideration)
  ai_budyko = ai
  ei_budyko = 1 + (ai_budyko) - (1 + (ai_budyko)^w_param)^(1/w_param)
  # remove inf values
  ai_data = rn_temp/p_temp
  ei_data = et_temp/p_temp
  ei_data[ei_data==Inf] = NA
  ei_data[ei_data==-Inf] = NA
  ai_data[ai_data==Inf] = NA
  ai_data[ai_data==-Inf] = NA
  
  # determine the distance metric
  distance_temp = sqrt(((ei_budyko - ei_data))^2 + ((ai_budyko - ai_data)^2))
  
  # append the distance metric data to the final dataset
  distance_final = cbind(distance_final, distance_temp)
  # append the combination names to the final vector of combinations
  combination = c(combination, paste0(p_data[i],".",et_data[j]))
  
  }
}

# prepare the final dataset for output
distance_final = round( distance_final,4)
distance_final = cbind(hyd_basins,area,round(ai,4),distance_final)
colnames(distance_final) = c("HydroSHEDS_ID","Area","AI",combination)
distance_final = gsub("NaN", "NA", distance_final)

# write the distance table to a text file
write.table(distance_final,paste0("Output/Budyko_Distance_HydroSHEDS_Lvl05_Rn-",rn_data[k],"_noSM2RAIN.txt"),row.names=FALSE, quote = FALSE)

}



