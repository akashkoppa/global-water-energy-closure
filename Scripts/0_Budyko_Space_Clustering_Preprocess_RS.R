# Script to preprocess the clusters of P and ET datasets in the Budyko space
# Author - Akash Koppa
# Date - 2019-09-22

# clear workspace
rm(list=ls())

# load required libraries

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

# subset ei and ai data for the Amazon region
#p_reqd = p[index_reqd,]

# subset a few catchments with higher uncertainty in Amazonia
#catch_reqd_high = c(6050315910, 6050315920, 6050029730)
#index_reqd_high = which(ei$HydroSHEDS_ID %in% catch_reqd_high)
#p_reqd_high = ei[index_reqd_high,]


# load the global aridity index data
ai = read.table("Input/Global_Aridity_Index_SRB.txt", header=TRUE)
ai = ai$AI

# load the Budyko parameter (w) data
w_temp      = read.table("Input/Budyko_Parameter_HydroBasinsLvl05_All.txt",header = TRUE)
w_param     = w_temp$W_Combined # currently considering the parameter derived from the combined equation in Xu et al 2013
area        = w_temp$Area
hyd_basins  = as.character(w_temp$HydroSHEDS_ID)
remove(w_temp)

# plot the number of precipitation and evapotranspiration for each catchment (basically a world map)
number_p = NULL
number_et = NULL
for (i in 1:nrow(p)){
  number_p  = c(number_p, ncol(p)-length(which(is.na(p[i,]))))
  number_et = c(number_et, ncol(et)-length(which(is.na(et[i,]))))
}

# loop through different combinations of P and E datasets 
EI_final = NULL
AI_final = NULL
combination = NULL
  for (i in 1:length(p_data)){
    p_temp = p[,i]
    
    for (j in 1:length(et_data)){
      et_temp = et[,j]
      rn_temp = rn[,1]
      
      # determine ai_data and ei_data
      ai_data = rn_temp/p_temp
      ei_data = et_temp/p_temp
      
      EI_final = cbind(EI_final, ei_data)
      AI_final = cbind(AI_final, ai_data)
      
      combination = c(combination, paste0(p_data[i],".",et_data[j]))
    }
  }

# plot the number of EI and AI for each catchment (basically a world map)
number_ei = NULL
number_ai = NULL
for (i in 1:nrow(AI_final)){
  number_ei = c(number_ei, ncol(EI_final)-length(which(is.na(EI_final[i,]))))
  number_ai = c(number_ai, ncol(AI_final)-length(which(is.na(AI_final[i,]))))
}

# prepare the final dataset for output
EI_final = round(EI_final,4)
EI_final = cbind(hyd_basins,area,round(ai,4),EI_final)
colnames(EI_final) = c("HydroSHEDS_ID","Area","AI_Actual",combination)
EI_final = gsub("NaN", "NA", EI_final)

AI_final = round(AI_final,4)
AI_final = cbind(hyd_basins,area,round(ai,4),AI_final)
colnames(AI_final) = c("HydroSHEDS_ID","Area","AI_Actual",combination)
AI_final = gsub("NaN", "NA", AI_final)

# write the data into text files
write.table(EI_final,paste0("Output/Evaporative_Index_Rn-",rn_data[1],"_RS_noSM2RAIN.txt"),row.names=FALSE, quote = FALSE)
write.table(AI_final,paste0("Output/Aridity_Index_Rn-",rn_data[1],"_RS_noSM2RAIN.txt"),row.names=FALSE, quote = FALSE)

  