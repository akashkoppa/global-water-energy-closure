# Script to perform some fine detailed analysis and lock into some final results
# Author - Akash Koppa
# Date - 2019 - 08 - 21
# Description - Separate the catchments according to aridity and carry out the scatter and regression with different variables.

# clear workspace
rm(list=ls())

# load required libraries
library(quantreg)
library(ggplot2)

# set working directory
setwd("/home/akash/Documents/FifthPaper/Budyko_Analysis_New/")

# load the error data 
cluster_metric_all = read.table("Output/Budyko_Cluster_Metrics_RS.txt",header=TRUE)
cluster_metric_all = cluster_metric_all[,-c(1:3)]

# load the uncertainty data
dist_metric_all = read.table("Output/Budyko_Distance_HydroSHEDS_Lvl05_Rn-CERES.txt",header=TRUE)

# load the predictors data
class_criteria_all = read.table("Input/Predictors_For_HydrobasinsLvl05_Final.txt",header=TRUE)

# separate the area and AI data
hydrosheds_id = as.character(dist_metric_all$HydroSHEDS_ID)
area = dist_metric_all$Area
ai = dist_metric_all$AI
dist_metric_all = dist_metric_all[,-c(1:3)]
ai_all = ai

# normalize the cluster metrics by the aridity index
dist_metric_all = dist_metric_all/ai
cluster_metric_all = cluster_metric_all/ai

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
# create a list of indices
list_indices = list(index_hyperarid,index_arid,index_semiarid,index_drysubhumid,index_humid)

catchments_ai_missing = hydrosheds_id[which(is.na(ai))]
catchments_ai_hyperarid = hydrosheds_id[which(ai>33.33)]
catchmets_ai_arid = hydrosheds_id[which(ai<=33.33 & ai>5)]
catchments_ai_semiarid = hydrosheds_id[which(ai<=5 & ai>2)]
catchments_ai_drysubhumid = hydrosheds_id[which(ai<=2 & ai>1.5)]
catchments_ai_humid = hydrosheds_id[which(ai<1.5)]

# create a vector of catchment aridity
aridity = rep(NA,length(ai))
aridity[index_hyperarid] = "Hyper Arid"
aridity[index_arid] = "Arid"
aridity[index_semiarid] = "Semi-Arid"
aridity[index_drysubhumid] = "Dry Sub-Humid"
aridity[index_humid] = "Humid"
aridity_all = aridity

# get the name of the combinations of P and ET datasets
p_et_dataset = colnames(dist_metric_all)
p_dataset = unique(unlist(strsplit(p_et_dataset,"[.]"))[seq(from=1,to=length(p_et_dataset)*2,by=2)])
et_dataset = unique(unlist(strsplit(p_et_dataset,"[.]"))[seq(from=2,to=length(p_et_dataset)*2,by=2)]) 

arid_unique = c("Hyper_Arid", "Arid", "Semi-Arid", "Dry_Sub-Humid", "Humid")
colors_reqd = c("#d11d1d","#b2642c","#37842a","#3faec1","#1f5b93")
              
slope_25_all = NULL
slope_50_all = NULL
slope_75_all = NULL
# loop through the different aridity classes and perform the analyses
for (x in 1:length(arid_unique)){
  # subset the cluster and distance metric data frames, and all the other required datasets
  cluster_metric = cluster_metric_all[list_indices[[x]],]
  dist_metric = dist_metric_all[list_indices[[x]],]
  class_criteria = class_criteria_all[list_indices[[x]],]
  ai = ai_all[list_indices[[x]]]
  aridity = aridity_all[list_indices[[x]]]
  
  # 1) Aridity 
  print("Processing Regression 1")
  Y = cbind(cluster_metric$Density_Mean_Radius)
  X = cbind(ai)
  
  # remove NAs
  index_na = which(is.na(Y[,1]))
  aridity_temp = aridity
  if(length(index_na>0)){
  Y = cbind(Y[-index_na,1])
  X = cbind(X[-index_na,1])
  aridity_temp = aridity_temp[-index_na]
  }
  
  # remove Inf
  index_inf = which(is.infinite(X[,1]))
  if(length(index_inf>0)){
  X = cbind(X[-index_inf,1])
  Y = cbind(Y[-index_inf,1])
  aridity_temp = aridity_temp[-index_inf]
  }
  
  colnames(Y) = "Y"
  colnames(X) = "X"
  # normalize the data 
  #Y = (Y - min(Y))/(max(Y) - min(Y))
  #X = (X - min(X))/(max(X) - min(X))
  QR_arid = rq(Y~X, tau=seq(from=0.05, to=0.95, by=0.1))
  sumQR_arid = summary(QR_arid)
  
  # Determine if the 25% quantile and 75% quantile coefficients are different
  qreg_arid_25 = rq(Y~X, tau=0.25)
  qreg_arid_50 = rq(Y~X, tau=0.50)
  qreg_arid_75 = rq(Y~X, tau=0.75)
  sign_arid    = anova(qreg_arid_25, qreg_arid_75)
  
  # 2) Scatter Plots with Regression Lines
  data_scatter = data.frame(Y = Y,
                            X = X,
                            Aridity=aridity_temp)
  #data_scatter$Aridity = factor(data_scatter$Aridity, levels=c("Hyper Arid","Arid","Semi-Arid","Dry Sub-Humid","Humid"))
  data_coeff_25 = data.frame(x = X,
                             y = qreg_arid_25$coefficients[2]*X + qreg_arid_25$coefficients[1],
                             Quantile = "0.25 Quantile") 
  data_coeff_50 = data.frame(x = X,
                             y = qreg_arid_50$coefficients[2]*X + qreg_arid_50$coefficients[1],
                             Quantile = "0.50 Quantile") 
  data_coeff_75 = data.frame(x = X,
                             y = qreg_arid_75$coefficients[2]*X + qreg_arid_75$coefficients[1],
                             Quantile = "0.75 Quantile")
  data_regline = rbind(data_coeff_25,data_coeff_50,data_coeff_75)
  colnames(data_regline) = c("X","Y","Quantile")
  # plot the scatter diagram
  plot_arid_scatter = ggplot()+
    geom_point(data=data_scatter, aes(x=X, y=Y),shape=1,color=colors_reqd[x])+
    #scale_color_manual(values = c("#d11d1d","#b2642c","#37842a","#3faec1","#1f5b93"))+
    geom_line(data = data_regline, aes(x=X,y=Y,linetype = Quantile))+
    coord_cartesian(xlim=c(quantile(X,probs = 0.0),quantile(X,probs = 0.9)), ylim=c(quantile(Y,probs = 0.0),quantile(Y,probs=0.9)))+
    theme(axis.line = element_line(color="black"))+
    theme(strip.background = element_blank())+
    theme(panel.background = element_rect(fill=alpha("#C8AD7F",0.5)))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour="black",fill=NA))+
    theme(legend.key = element_blank())+
    labs(x="Aridity",y= "Cluster Density")+
    ggtitle("Scatter of Cluster Density vs Aridity")
  ggsave(plot = plot_arid_scatter, file=paste0("Output/Cluster/Aridity_Class/Scatter_Arid_Uncertainty_",arid_unique[x],".png"),width=6,height=4,dpi=300)
  ##----------------------------------------------------------------------------------------------------------##
  
  ##----------------------------------------------------------------------------------------------------------##
  # 2) Elevation
  print("Processing Regression 2")
  Y = cbind(cluster_metric$Density_Mean_Radius)
  X = cbind(class_criteria$Elevation)
  # remove NAs
  index_na = which(is.na(Y[,1]))
  aridity_temp = aridity
  aridity_temp = aridity_temp[-index_na]
  Y = cbind(Y[-index_na,1])
  X = cbind(X[-index_na,1])
  colnames(Y) = "Y"
  colnames(X) = "X"
  # normalize the data 
  #Y = (Y - min(Y))/(max(Y) - min(Y))
  #X = (X - min(X))/(max(X) - min(X))
  QR_elev = rq(Y~X, tau=seq(from=0.05, to=0.95, by=0.1))
  sumQR_elev = summary(QR_elev, se="rank")
  
  # Determine if the 25% quantile and 75% quantile coefficients are different
  qreg_elev_25 = rq(Y~X, tau=0.25)
  qreg_elev_50 = rq(Y~X, tau=0.50)
  qreg_elev_75 = rq(Y~X, tau=0.75)
  sign_elev    = anova(qreg_elev_25, qreg_elev_75)
  
  # Plot the required quantile regression plots
  # 1) Regression Slopes
  data_plot_slope = data.frame(Quantile = seq(from=0.05, to=0.95, by=0.1),
                               Slope = QR_elev$coefficients[2,])
  
  # 2) Scatter Plots with Regression Lines
  data_scatter = data.frame(Y = Y,
                            X = X,
                            Aridity = aridity_temp)
  #data_scatter$Aridity = factor(data_scatter$Aridity, levels=c("Hyper Arid","Arid","Semi-Arid","Dry Sub-Humid","Humid"))
  data_coeff_25 = data.frame(x = X,
                             y = qreg_elev_25$coefficients[2]*X + qreg_elev_25$coefficients[1],
                             Quantile = "0.25 Quantile") 
  data_coeff_50 = data.frame(x = X,
                             y = qreg_elev_50$coefficients[2]*X + qreg_elev_50$coefficients[1],
                             Quantile = "0.50 Quantile") 
  data_coeff_75 = data.frame(x = X,
                             y = qreg_elev_75$coefficients[2]*X + qreg_elev_75$coefficients[1],
                             Quantile = "0.75 Quantile")
  data_regline = rbind(data_coeff_25,data_coeff_50,data_coeff_75)
  colnames(data_regline) = c("X","Y","Quantile")
  # plot the scatter diagram
  plot_elev_scatter = ggplot()+
    geom_point(data=data_scatter, aes(x=X, y=Y), shape=1,color=colors_reqd[x])+
    #scale_color_manual(values = c("#d11d1d","#b2642c","#37842a","#3faec1","#1f5b93"))+
    geom_line(data = data_regline, aes(x=X,y=Y,linetype = Quantile))+
    coord_cartesian(xlim=c(quantile(X,probs = 0.0),quantile(X,probs = 0.9)), ylim=c(quantile(Y,probs = 0.0),quantile(Y,probs=0.9)))+
    theme(axis.line = element_line(color="black"))+
    theme(strip.background = element_blank())+
    theme(panel.background = element_rect(fill=alpha("#C8AD7F",0.5)))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour="black",fill=NA))+
    theme(legend.key = element_blank())+
    labs(x="Elevation",y= "Cluster Density")+
    ggtitle("Scatter of Cluster Density vs Elevation")
  ggsave(plot = plot_elev_scatter, file=paste0("Output/Cluster/Aridity_Class/Scatter_Elevation_Uncertainty_",arid_unique[x],".png"),width=6,height=4,dpi=300)
  ##----------------------------------------------------------------------------------------------------------##
  
  ##----------------------------------------------------------------------------------------------------------##
  # 3) CTI
  print("Processing Regression 3")
  Y = cbind(cluster_metric$Density_Mean_Radius)
  X = cbind(class_criteria$CTI)
  # remove NAs
  index_na = which(is.na(Y[,1]))
  Y = cbind(Y[-index_na,1])
  X = cbind(X[-index_na,1])
  aridity_temp = aridity
  aridity_temp = aridity_temp[-index_na]
  colnames(Y) = "Y"
  colnames(X) = "X"
  # normalize the data 
  #Y = (Y - min(Y))/(max(Y) - min(Y))
  #X = (X - min(X))/(max(X) - min(X))
  QR_cti = rq(Y~X, tau=seq(from=0.05, to=0.95, by=0.1))
  sumQR_cti = summary(QR_cti)
  
  # Determine if the 25% quantile and 75% quantile coefficients are different
  qreg_cti_25 = rq(Y~X, tau=0.25)
  qreg_cti_50 = rq(Y~X, tau=0.50)
  qreg_cti_75 = rq(Y~X, tau=0.75)
  sign_cti    = anova(qreg_cti_25, qreg_cti_75)
  
  # Plot the required quantile regression plots
  # 1) Regression Slopes
  data_plot_slope = data.frame(Quantile = seq(from=0.05, to=0.95, by=0.1),
                               Slope = QR_cti$coefficients[2,])
  
  # 2) Scatter Plots with Regression Lines
  data_scatter = data.frame(Y = Y,
                            X = X,
                            Aridity = aridity_temp)
  #data_scatter$Aridity = factor(data_scatter$Aridity, levels=c("Hyper Arid","Arid","Semi-Arid","Dry Sub-Humid","Humid"))
  data_coeff_25 = data.frame(x = X,
                             y = qreg_cti_25$coefficients[2]*X + qreg_cti_25$coefficients[1],
                             Quantile = "0.25 Quantile") 
  data_coeff_50 = data.frame(x = X,
                             y = qreg_cti_50$coefficients[2]*X + qreg_cti_50$coefficients[1],
                             Quantile = "0.50 Quantile") 
  data_coeff_75 = data.frame(x = X,
                             y = qreg_cti_75$coefficients[2]*X + qreg_cti_75$coefficients[1],
                             Quantile = "0.75 Quantile")
  data_regline = rbind(data_coeff_25,data_coeff_50,data_coeff_75)
  colnames(data_regline) = c("X","Y","Quantile")
  # plot the scatter diagram
  plot_cti_scatter = ggplot()+
    geom_point(data=data_scatter, aes(x=X, y=Y), shape=1, color=colors_reqd[x])+
    #scale_color_manual(values = c("#d11d1d","#b2642c","#37842a","#3faec1","#1f5b93"))+
    geom_line(data = data_regline, aes(x=X,y=Y,linetype = Quantile))+
    coord_cartesian(xlim=c(quantile(X,probs = 0.0),quantile(X,probs = 0.9)), ylim=c(quantile(Y,probs = 0.0),quantile(Y,probs=0.9)))+
    theme(axis.line = element_line(color="black"))+
    theme(strip.background = element_blank())+
    theme(panel.background = element_rect(fill=alpha("#C8AD7F",0.5)))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour="black",fill=NA))+
    theme(legend.key = element_blank())+
    labs(x="CTI",y= "Cluster Density")+
    ggtitle("Scatter of Cluster Density vs CTI")
  ggsave(plot = plot_cti_scatter, file=paste0("Output/Cluster/Aridity_Class/Scatter_CTI_Uncertainty_",arid_unique[x],".png"),width=6,height=4,dpi=300)
  ##----------------------------------------------------------------------------------------------------------##
  
  ##----------------------------------------------------------------------------------------------------------##
  # 4) NDVI
  print("Processing Regression 4")
  Y = cbind(cluster_metric$Density_Mean_Radius)
  X = cbind(class_criteria$NDVI)
  # remove NAs
  index_na = which(is.na(Y[,1]))
  Y = cbind(Y[-index_na,1])
  X = cbind(X[-index_na,1])
  aridity_temp = aridity
  aridity_temp = aridity_temp[-index_na]
  colnames(Y) = "Y"
  colnames(X) = "X"
  index_na = which(is.na(X[,1]))
  Y = cbind(Y[-index_na,1])
  X = cbind(X[-index_na,1])
  aridity_temp = aridity_temp[-index_na]
  # normalize the data 
  #Y = (Y - min(Y))/(max(Y) - min(Y))
  #X = (X - min(X))/(max(X) - min(X))
  QR_ndvi = rq(Y~X, tau=seq(from=0.05, to=0.95, by=0.1))
  sumQR_ndvi = summary(QR_ndvi)
  
  # Determine if the 25% quantile and 75% quantile coefficients are different
  qreg_ndvi_25 = rq(Y~X, tau=0.25)
  qreg_ndvi_50 = rq(Y~X, tau=0.50)
  qreg_ndvi_75 = rq(Y~X, tau=0.75)
  sign_ndvi    = anova(qreg_ndvi_25, qreg_ndvi_75)
  
  # 2) Scatter Plots with Regression Lines
  data_scatter = data.frame(Y = Y,
                            X = X,
                            Aridity = aridity_temp)
  #data_scatter$Aridity = factor(data_scatter$Aridity, levels=c("Hyper Arid","Arid","Semi-Arid","Dry Sub-Humid","Humid"))
  data_coeff_25 = data.frame(x = X,
                             y = qreg_ndvi_25$coefficients[2]*X + qreg_ndvi_25$coefficients[1],
                             Quantile = "0.25 Quantile") 
  data_coeff_50 = data.frame(x = X,
                             y = qreg_ndvi_50$coefficients[2]*X + qreg_ndvi_50$coefficients[1],
                             Quantile = "0.50 Quantile") 
  data_coeff_75 = data.frame(x = X,
                             y = qreg_ndvi_75$coefficients[2]*X + qreg_ndvi_75$coefficients[1],
                             Quantile = "0.75 Quantile")
  data_regline = rbind(data_coeff_25,data_coeff_50,data_coeff_75)
  colnames(data_regline) = c("X","Y","Quantile")
  # plot the scatter diagram
  plot_ndvi_scatter = ggplot()+
    geom_point(data=data_scatter, aes(x=X, y=Y), shape=1, color=colors_reqd[x])+
    #scale_color_manual(values = c("#d11d1d","#b2642c","#37842a","#3faec1","#1f5b93"))+
    geom_line(data = data_regline, aes(x=X,y=Y,linetype = Quantile))+
    coord_cartesian(xlim=c(quantile(X,probs = 0.0),quantile(X,probs = 0.9)), ylim=c(quantile(Y,probs = 0.0),quantile(Y,probs=0.9)))+
    theme(axis.line = element_line(color="black"))+
    theme(strip.background = element_blank())+
    theme(panel.background = element_rect(fill=alpha("#C8AD7F",0.5)))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour="black",fill=NA))+
    theme(legend.key = element_blank())+
    labs(x="NDVI",y= "Cluster Density")+
    ggtitle("Scatter of Cluster Density vs NDVI")
  ggsave(plot = plot_ndvi_scatter, file=paste0("Output/Cluster/Aridity_Class/Scatter_NDVI_Uncertainty_",arid_unique[x],".png"),width=6,height=4,dpi=300)
  
  
  # For Budyko errors
  dist_mean = apply(dist_metric,1,mean,na.rm=TRUE)
  dist_mean[which(is.nan(dist_mean))] = NA
  dist_mean[which(is.infinite(dist_mean))] = NA
  
  ##----------------------------------------------------------------------------------------------------------##
  # 1) Aridity
  print("Processing Regression 5")
  Y = cbind(dist_mean)
  Y[Y>1000] = NA
  X = cbind(ai)
  # remove NAs
  index_na = which(is.na(Y[,1]))
  Y = cbind(Y[-index_na,1])
  X = cbind(X[-index_na,1])
  aridity_temp = aridity
  aridity_temp = aridity_temp[-index_na]
  colnames(Y) = "Y"
  colnames(X) = "X"
  # normalize the data 
  #Y = (Y - min(Y))/(max(Y) - min(Y))
  #X = (X - min(X))/(max(X) - min(X))
  QR_arid_err = rq(Y~X, tau=seq(from=0.05, to=0.95, by=0.1))
  sumQR_arid_err = summary(QR_arid_err, se="boot")
  
  # Determine if the 25% quantile and 75% quantile coefficients are different
  qreg_arid_err_25 = rq(Y~X, tau=0.25)
  qreg_arid_err_50 = rq(Y~X, tau=0.50)
  qreg_arid_err_75 = rq(Y~X, tau=0.75)
  sign_arid_err    = anova(qreg_arid_err_25, qreg_arid_err_75)
  
  # 2) Scatter Plots with Regression Lines
  data_scatter = data.frame(Y = Y,
                            X = X,
                            Aridity = aridity_temp)
 # data_scatter$Aridity = factor(data_scatter$Aridity, levels=c("Hyper Arid","Arid","Semi-Arid","Dry Sub-Humid","Humid"))
  data_coeff_25 = data.frame(x = X,
                             y = qreg_arid_err_25$coefficients[2]*X + qreg_arid_err_25$coefficients[1],
                             Quantile = "0.25 Quantile") 
  data_coeff_50 = data.frame(x = X,
                             y = qreg_arid_err_50$coefficients[2]*X + qreg_arid_err_50$coefficients[1],
                             Quantile = "0.50 Quantile") 
  data_coeff_75 = data.frame(x = X,
                             y = qreg_arid_err_75$coefficients[2]*X + qreg_arid_err_75$coefficients[1],
                             Quantile = "0.75 Quantile")
  data_regline = rbind(data_coeff_25,data_coeff_50,data_coeff_75)
  colnames(data_regline) = c("X","Y","Quantile")
  # plot the scatter diagram
  plot_arid_err_scatter = ggplot()+
    geom_point(data=data_scatter, aes(x=X, y=Y), shape=1, color=colors_reqd[x])+
    scale_color_manual(values = c("#d11d1d","#b2642c","#37842a","#3faec1","#1f5b93"))+
    geom_line(data = data_regline, aes(x=X,y=Y,linetype = Quantile))+
    coord_cartesian(xlim=c(quantile(X,probs = 0.0),quantile(X,probs = 0.9)), ylim=c(quantile(Y,probs = 0.0),quantile(Y,probs=0.9)))+
    theme(axis.line = element_line(color="black"))+
    theme(strip.background = element_blank())+
    theme(panel.background = element_rect(fill=alpha("#C8AD7F",0.5)))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour="black",fill=NA))+
    theme(legend.key = element_blank())+
    labs(x="Aridity",y= "Budyko Error")+
    ggtitle("Scatter of Budyko Error vs Aridity")
  ggsave(plot = plot_arid_err_scatter, file=paste0("Output/Cluster/Aridity_Class/Scatter_Arid_Error_",arid_unique[x],".png"),width=6,height=4,dpi=300)
  ##----------------------------------------------------------------------------------------------------------##
  
  ##----------------------------------------------------------------------------------------------------------##
  # 2) Elevation
  print("Processing Regression 6")
  Y = cbind(dist_mean)
  Y[Y>1000] = NA
  X = cbind(class_criteria$Elevation)
  # remove NAs
  index_na = which(is.na(Y[,1]))
  Y = cbind(Y[-index_na,1])
  X = cbind(X[-index_na,1])
  aridity_temp = aridity
  aridity_temp = aridity_temp[-index_na]
  colnames(Y) = "Y"
  colnames(X) = "X"
  # normalize the data 
  #Y = (Y - min(Y))/(max(Y) - min(Y))
  #X = (X - min(X))/(max(X) - min(X))
  QR_elev_err = rq(Y~X, tau=seq(from=0.05, to=0.95, by=0.1))
  sumQR_elev_err = summary(QR_elev_err, se="boot")
  
  # Determine if the 25% quantile and 75% quantile coefficients are different
  qreg_elev_err_25 = rq(Y~X, tau=0.25)
  qreg_elev_err_50 = rq(Y~X, tau=0.50)
  qreg_elev_err_75 = rq(Y~X, tau=0.75)
  sign_elev_err    = anova(qreg_elev_err_25, qreg_elev_err_75)
  
  # Plot the required quantile regression plots
  # 2) Scatter Plots with Regression Lines
  data_scatter = data.frame(Y = Y,
                            X = X,
                            Aridity = aridity_temp)
  #data_scatter$Aridity = factor(data_scatter$Aridity, levels=c("Hyper Arid","Arid","Semi-Arid","Dry Sub-Humid","Humid"))
  data_coeff_25 = data.frame(x = X,
                             y = qreg_elev_err_25$coefficients[2]*X + qreg_elev_err_25$coefficients[1],
                             Quantile = "0.25 Quantile") 
  data_coeff_50 = data.frame(x = X,
                             y = qreg_elev_err_50$coefficients[2]*X + qreg_elev_err_50$coefficients[1],
                             Quantile = "0.50 Quantile") 
  data_coeff_75 = data.frame(x = X,
                             y = qreg_elev_err_75$coefficients[2]*X + qreg_elev_err_75$coefficients[1],
                             Quantile = "0.75 Quantile")
  data_regline = rbind(data_coeff_25,data_coeff_50,data_coeff_75)
  colnames(data_regline) = c("X","Y","Quantile")
  # plot the scatter diagram
  plot_elev_err_scatter = ggplot()+
    geom_point(data=data_scatter, aes(x=X, y=Y), shape=1, color=colors_reqd[x])+
    scale_color_manual(values = c("#d11d1d","#b2642c","#37842a","#3faec1","#1f5b93"))+
    geom_line(data = data_regline, aes(x=X,y=Y,linetype = Quantile))+
    coord_cartesian(xlim=c(quantile(X,probs = 0.0),quantile(X,probs = 0.9)), ylim=c(quantile(Y,probs = 0.0),quantile(Y,probs=0.9)))+
    theme(axis.line = element_line(color="black"))+
    theme(strip.background = element_blank())+
    theme(panel.background = element_rect(fill=alpha("#C8AD7F",0.5)))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour="black",fill=NA))+
    theme(legend.key = element_blank())+
    labs(x="Elevation",y= "Budyko Error")+
    ggtitle("Scatter of Budyko Error vs Elevation")
  ggsave(plot = plot_elev_err_scatter, file=paste0("Output/Cluster/Aridity_Class/Scatter_Elevation_Error_",arid_unique[x],".png"),width=6,height=4,dpi=300)
  
  ##----------------------------------------------------------------------------------------------------------##
  
  ##----------------------------------------------------------------------------------------------------------##
  
  # 3) CTI
  print("Processing Regression 7")
  Y = cbind(dist_mean)
  Y[Y>1000] = NA
  X = cbind(class_criteria$CTI)
  # remove NAs
  index_na = which(is.na(Y[,1]))
  Y = cbind(Y[-index_na,1])
  X = cbind(X[-index_na,1])
  aridity_temp = aridity
  aridity_temp = aridity_temp[-index_na]
  colnames(Y) = "Y"
  colnames(X) = "X"
  # normalize the data 
  #Y = (Y - min(Y))/(max(Y) - min(Y))
  #X = (X - min(X))/(max(X) - min(X))
  QR_cti_err = rq(Y~X, tau=seq(from=0.05, to=0.95, by=0.1))
  sumQR_cti_err = summary(QR_cti_err, se="boot")
  
  # Determine if the 25% quantile and 75% quantile coefficients are different
  qreg_cti_err_25 = rq(Y~X, tau=0.25)
  qreg_cti_err_50 = rq(Y~X, tau=0.50)
  qreg_cti_err_75 = rq(Y~X, tau=0.75)
  sign_cti_err    = anova(qreg_cti_err_25, qreg_cti_err_75)
  
  # Plot the required quantile regression plots
   # 2) Scatter Plots with Regression Lines
  data_scatter = data.frame(Y = Y,
                            X = X,
                            Aridity = aridity_temp)
#  data_scatter$Aridity = factor(data_scatter$Aridity, levels=c("Hyper Arid","Arid","Semi-Arid","Dry Sub-Humid","Humid"))
  data_coeff_25 = data.frame(x = X,
                             y = qreg_cti_err_25$coefficients[2]*X + qreg_cti_err_25$coefficients[1],
                             Quantile = "0.25 Quantile") 
  data_coeff_50 = data.frame(x = X,
                             y = qreg_cti_err_50$coefficients[2]*X + qreg_cti_err_50$coefficients[1],
                             Quantile = "0.50 Quantile") 
  data_coeff_75 = data.frame(x = X,
                             y = qreg_cti_err_75$coefficients[2]*X + qreg_cti_err_75$coefficients[1],
                             Quantile = "0.75 Quantile")
  data_regline = rbind(data_coeff_25,data_coeff_50,data_coeff_75)
  colnames(data_regline) = c("X","Y","Quantile")
  # plot the scatter diagram
  plot_cti_err_scatter = ggplot()+
    geom_point(data=data_scatter, aes(x=X, y=Y), shape=1, color=colors_reqd[x])+
    scale_color_manual(values = c("#d11d1d","#b2642c","#37842a","#3faec1","#1f5b93"))+
    geom_line(data = data_regline, aes(x=X,y=Y,linetype = Quantile))+
    coord_cartesian(xlim=c(quantile(X,probs = 0.0),quantile(X,probs = 0.9)), ylim=c(quantile(Y,probs = 0.0),quantile(Y,probs=0.9)))+
    theme(axis.line = element_line(color="black"))+
    theme(strip.background = element_blank())+
    theme(panel.background = element_rect(fill=alpha("#C8AD7F",0.5)))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour="black",fill=NA))+
    theme(legend.key = element_blank())+
    labs(x="CTI",y= "Budyko Error")+
    ggtitle("Scatter of Budyko Error vs CTI")
  ggsave(plot = plot_cti_err_scatter, file=paste0("Output/Cluster/Aridity_Class/Scatter_CTI_Error_",arid_unique[x],".png"),width=6,height=4,dpi=300)
  
  
  ##----------------------------------------------------------------------------------------------------------##
  
  ##----------------------------------------------------------------------------------------------------------##
  # 4) NDVI
  print("Processing Regression 8")
  Y = cbind(dist_mean)
  Y[Y>1000] = NA
  X = cbind(class_criteria$NDVI)
  # remove NAs
  index_na = which(is.na(Y[,1]))
  Y = cbind(Y[-index_na,1])
  X = cbind(X[-index_na,1])
  aridity_temp = aridity
  aridity_temp = aridity_temp[-index_na]
  colnames(Y) = "Y"
  colnames(X) = "X"
  #index_na = which(is.na(X[,1]))
  #Y = cbind(Y[-index_na,1])
  #X = cbind(X[-index_na,1])
  # normalize the data 
  #Y = (Y - min(Y))/(max(Y) - min(Y))
  #X = (X - min(X))/(max(X) - min(X))
  QR_ndvi_err = rq(Y~X, tau=seq(from=0.05, to=0.95, by=0.1))
  sumQR_ndvi_err = summary(QR_ndvi_err, se="boot")
  
  # Determine if the 25% quantile and 75% quantile coefficients are different
  qreg_ndvi_err_25 = rq(Y~X, tau=0.25)
  qreg_ndvi_err_50 = rq(Y~X, tau=0.50)
  qreg_ndvi_err_75 = rq(Y~X, tau=0.75)
  sign_ndvi_err    = anova(qreg_ndvi_err_25, qreg_ndvi_err_75)
  
  # Plot the required quantile regression plots
  # 2) Scatter Plots with Regression Lines
  data_scatter = data.frame(Y = Y,
                            X = X,
                            Aridity = aridity_temp)
  #data_scatter$Aridity = factor(data_scatter$Aridity, levels=c("Hyper Arid","Arid","Semi-Arid","Dry Sub-Humid","Humid"))
  data_coeff_25 = data.frame(x = X,
                             y = qreg_ndvi_err_25$coefficients[2]*X + qreg_ndvi_err_25$coefficients[1],
                             Quantile = "0.25 Quantile") 
  data_coeff_50 = data.frame(x = X,
                             y = qreg_ndvi_err_50$coefficients[2]*X + qreg_ndvi_err_50$coefficients[1],
                             Quantile = "0.50 Quantile") 
  data_coeff_75 = data.frame(x = X,
                             y = qreg_ndvi_err_75$coefficients[2]*X + qreg_ndvi_err_75$coefficients[1],
                             Quantile = "0.75 Quantile")
  data_regline = rbind(data_coeff_25,data_coeff_50,data_coeff_75)
  colnames(data_regline) = c("X","Y","Quantile")
  # plot the scatter diagram
  plot_ndvi_err_scatter = ggplot()+
    geom_point(data=data_scatter, aes(x=X, y=Y), shape=1, color=colors_reqd[x])+
    scale_color_manual(values = c("#d11d1d","#b2642c","#37842a","#3faec1","#1f5b93"))+
    geom_line(data = data_regline, aes(x=X,y=Y,linetype = Quantile))+
    coord_cartesian(xlim=c(quantile(X,probs = 0.0),quantile(X,probs = 0.9)), ylim=c(quantile(Y,probs = 0.0),quantile(Y,probs=0.9)))+
    theme(axis.line = element_line(color="black"))+
    theme(strip.background = element_blank())+
    theme(panel.background = element_rect(fill=alpha("#C8AD7F",0.5)))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour="black",fill=NA))+
    theme(legend.key = element_blank())+
    labs(x="NDVI",y= "Budyko Error")+
    ggtitle("Scatter of Budyko Error vs NDVI")
  ggsave(plot = plot_ndvi_err_scatter, file=paste0("Output/Cluster/Aridity_Class/Scatter_NDVI_Error_",arid_unique[x],".png"),width=6,height=4,dpi=300)
  
  ##----------------------------------------------------------------------------------------------------------#
  
# store the slope and intercept values 
slope_25_all = rbind(slope_25_all, c(qreg_arid_25$coefficients[2],qreg_elev_25$coefficients[2],qreg_cti_25$coefficients[2],qreg_ndvi_25$coefficients[2]))
slope_50_all = rbind(slope_50_all, c(qreg_arid_50$coefficients[2],qreg_elev_50$coefficients[2],qreg_cti_50$coefficients[2],qreg_ndvi_50$coefficients[2]))
slope_75_all = rbind(slope_75_all, c(qreg_arid_75$coefficients[2],qreg_elev_75$coefficients[2],qreg_cti_75$coefficients[2],qreg_ndvi_75$coefficients[2]))
}

slope_25_all = data.frame(slope_25_all)
slope_50_all = data.frame(slope_50_all) 
slope_75_all = data.frame(slope_75_all)
colnames(slope_25_all) = c("Aridity","Elevation","CTI","NDVI")
colnames(slope_50_all) = c("Aridity","Elevation","CTI","NDVI")
colnames(slope_75_all) = c("Aridity","Elevation","CTI","NDVI")

rownames(slope_25_all) = arid_unique
rownames(slope_50_all) = arid_unique
rownames(slope_75_all) = arid_unique

write.table(slope_25_all,"Output/Cluster/Aridity_Class/Slope_25_All.txt",row.names=TRUE, col.names=TRUE)
write.table(slope_50_all,"Output/Cluster/Aridity_Class/Slope_50_All.txt",row.names=TRUE, col.names=TRUE)
write.table(slope_75_all,"Output/Cluster/Aridity_Class/Slope_75_All.txt",row.names=TRUE, col.names=TRUE)

