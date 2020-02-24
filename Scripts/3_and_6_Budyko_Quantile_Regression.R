# Script to perform some fine detailed analysis and lock into some final results
# Author - Akash Koppa
# Date - 2019 - 08 - 21

# clear workspace
rm(list=ls())

# load required libraries
library(quantreg)
library(ggplot2)

# set working directory 
setwd("/home/akash/Documents/FifthPaper/Budyko_Analysis_New/")

# load the error data 
cluster_metric = read.table("Output/Budyko_Cluster_Metrics_RS.txt",header=TRUE)
cluster_metric = cluster_metric[,-c(1:3)]

# load the uncertainty data
dist_metric = read.table("Output/Budyko_Distance_HydroSHEDS_Lvl05_Rn-CERES.txt",header=TRUE)

# load the predictors data
class_criteria = read.table("Input/Predictors_For_HydrobasinsLvl05_Final.txt",header=TRUE)
#colnames(class_criteria) = c("HydroSHEDS_ID","Latitude","Longitude","Area","Elevation","Slope","CTI","NDVI")

# separate the area and AI data
hydrosheds_id = as.character(dist_metric$HydroSHEDS_ID)
area = dist_metric$Area
ai = dist_metric$AI
dist_metric = dist_metric[,-c(1:3)]

# normalize the cluster metrics by the aridity index
dist_metric = dist_metric/ai
cluster_metric = cluster_metric/ai

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

# create a vector of catchment aridity
aridity = rep(NA,length(ai))
aridity[index_hyperarid] = "Hyper Arid"
aridity[index_arid] = "Arid"
aridity[index_semiarid] = "Semi-Arid"
aridity[index_drysubhumid] = "Dry Sub-Humid"
aridity[index_humid] = "Humid"

# get the name of the combinations of P and ET datasets
p_et_dataset = colnames(dist_metric)
p_dataset = unique(unlist(strsplit(p_et_dataset,"[.]"))[seq(from=1,to=length(p_et_dataset)*2,by=2)])
et_dataset = unique(unlist(strsplit(p_et_dataset,"[.]"))[seq(from=2,to=length(p_et_dataset)*2,by=2)])     

# Implement quantile regression for the uncertainties/errors and the predictor variables
# For uncertainties (cluster density)

# 1) Aridity 
print("Processing Regression 1")
Y = cbind(cluster_metric$Mean_Radius)
X = cbind(ai)

# remove NAs
index_na = which(is.na(Y[,1]))
Y = cbind(Y[-index_na,1])
X = cbind(X[-index_na,1])
aridity_temp = aridity
aridity_temp = aridity_temp[-index_na]
# remove Inf
index_inf = which(is.infinite(X[,1]))
X = cbind(X[-index_inf,1])
Y = cbind(Y[-index_inf,1])
aridity_temp = aridity_temp[-index_inf]
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

# Plot the required quantile regression plots
# 1) Regression Slopes
data_plot_slope = data.frame(Quantile = seq(from=0.05, to=0.95, by=0.1),
                             Slope = QR_arid$coefficients[2,])
# collect all the standard errors
std_err_arid = NULL
for (i in 1:length(sumQR_arid)){
  std_err_arid = c(std_err_arid, sumQR_arid[[i]]$coefficients[2,2])
}
data_plot_ribbon = data.frame(Quantile = seq(from=0.05, to=0.95, by=0.1),
                              Ymin = QR_arid$coefficients[2,] - std_err_arid,
                              Ymax = QR_arid$coefficients[2,] + std_err_arid)
plot_arid = ggplot()+
            geom_ribbon(data = data_plot_ribbon, aes(x=Quantile, ymin=Ymin, ymax=Ymax),fill="#77916e",alpha=0.6)+
            geom_line(data = data_plot_slope, aes(x=Quantile, y=Slope))+
            geom_point(data = data_plot_slope, aes(x=Quantile, y=Slope))+
            theme(axis.line = element_line(color="black"))+
            theme(strip.background = element_blank())+
            theme(panel.background = element_rect(fill=alpha("#C8AD7F",0.7)))+
            theme(panel.grid.major = element_line(colour="#C8AD7F"),
                  panel.grid.minor = element_line(colour="#C8AD7F"),
                  panel.border = element_rect(colour="black",fill=NA))+
            theme(legend.key = element_blank())+
            labs(x="Quantiles",y= "Slope")+
            ggtitle("Slope of QR (Cluster Radius with Aridity) ")
#ggsave(plot = plot_arid, file="Output/Cluster/Slope_Arid_Uncertainty.png",width=7,height=4,dpi=300)

# 2) Scatter Plots with Regression Lines
data_scatter = data.frame(Y = Y,
                          X = X,
                          Aridity=aridity_temp)
data_scatter$Aridity = factor(data_scatter$Aridity, levels=c("Hyper Arid","Arid","Semi-Arid","Dry Sub-Humid","Humid"))
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
                    geom_point(data=data_scatter, aes(x=X, y=Y,color=Aridity),shape=1)+
                    scale_color_manual(values = c("#d11d1d","#b2642c","#37842a","#3faec1","#1f5b93"))+
                    geom_line(data = data_regline, aes(x=X,y=Y,linetype = Quantile))+
                    coord_cartesian(xlim=c(0,quantile(X,probs = 0.9)), ylim=c(0,quantile(Y,probs=0.9)))+
                    theme(axis.line = element_line(color="black"))+
                    theme(strip.background = element_blank())+
                    theme(panel.background = element_rect(fill=alpha("#C8AD7F",0.5)))+
                    theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_rect(colour="black",fill=NA))+
                    theme(legend.key = element_blank())+
                    labs(x="Aridity",y= "Cluster Radius")+
                    ggtitle("Scatter of Cluster Radius vs Aridity")
#ggsave(plot = plot_arid_scatter, file="Output/Cluster/Scatter_Arid_Uncertainty.png",width=6,height=4,dpi=300)
##----------------------------------------------------------------------------------------------------------##

##----------------------------------------------------------------------------------------------------------##
# 2) Elevation
print("Processing Regression 2")
Y = cbind(cluster_metric$Mean_Radius)
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
# collect all the standard errors
std_err_elev = NULL
for (i in 1:length(sumQR_elev)){
  std_err_elev = c(std_err_elev, sumQR_elev[[i]]$coefficients[2,2])
}
data_plot_ribbon = data.frame(Quantile = seq(from=0.05, to=0.95, by=0.1),
                              Ymin = QR_elev$coefficients[2,] - std_err_elev,
                              Ymax = QR_elev$coefficients[2,] + std_err_elev)
plot_elev = ggplot()+
  geom_ribbon(data = data_plot_ribbon, aes(x=Quantile, ymin=Ymin, ymax=Ymax),fill="#77916e",alpha=0.6)+
  geom_line(data = data_plot_slope, aes(x=Quantile, y=Slope))+
  geom_point(data = data_plot_slope, aes(x=Quantile, y=Slope))+
  theme(axis.line = element_line(color="black"))+
  theme(strip.background = element_blank())+
  theme(panel.background = element_rect(fill=alpha("#C8AD7F",0.7)))+
  theme(panel.grid.major = element_line(colour="#C8AD7F"),
        panel.grid.minor = element_line(colour="#C8AD7F"),
        panel.border = element_rect(colour="black",fill=NA))+
  theme(legend.key = element_blank())+
  labs(x="Quantiles",y= "Slope")+
  ggtitle("Slope of QR (Cluster Radius with Elevation) ")
#ggsave(plot = plot_elev, file="Output/Cluster/Slope_Elevation_Uncertainty.png",width=7,height=4,dpi=300)

# 2) Scatter Plots with Regression Lines
data_scatter = data.frame(Y = Y,
                          X = X,
                          Aridity = aridity_temp)
data_scatter$Aridity = factor(data_scatter$Aridity, levels=c("Hyper Arid","Arid","Semi-Arid","Dry Sub-Humid","Humid"))
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
  geom_point(data=data_scatter, aes(x=X, y=Y, color=Aridity), shape=1)+
  scale_color_manual(values = c("#d11d1d","#b2642c","#37842a","#3faec1","#1f5b93"))+
  geom_line(data = data_regline, aes(x=X,y=Y,linetype = Quantile))+
  coord_cartesian(xlim=c(0,quantile(X,probs = 0.9)), ylim=c(0,quantile(Y,probs=0.9)))+
  theme(axis.line = element_line(color="black"))+
  theme(strip.background = element_blank())+
  theme(panel.background = element_rect(fill=alpha("#C8AD7F",0.5)))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour="black",fill=NA))+
  theme(legend.key = element_blank())+
  labs(x="Elevation",y= "Cluster Radius")+
  ggtitle("Scatter of Cluster Radius vs Elevation")
#ggsave(plot = plot_elev_scatter, file="Output/Cluster/Scatter_Elevation_Uncertainty.png",width=6,height=4,dpi=300)
##----------------------------------------------------------------------------------------------------------##

##----------------------------------------------------------------------------------------------------------##
# 3) CTI
print("Processing Regression 3")
Y = cbind(cluster_metric$Mean_Radius)
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
# collect all the standard errors
std_err_cti = NULL
for (i in 1:length(sumQR_elev)){
  std_err_cti = c(std_err_cti, sumQR_cti[[i]]$coefficients[2,2])
}
data_plot_ribbon = data.frame(Quantile = seq(from=0.05, to=0.95, by=0.1),
                              Ymin = QR_cti$coefficients[2,] - std_err_cti,
                              Ymax = QR_cti$coefficients[2,] + std_err_cti)
plot_cti = ggplot()+
  geom_ribbon(data = data_plot_ribbon, aes(x=Quantile, ymin=Ymin, ymax=Ymax),fill="#77916e",alpha=0.6)+
  geom_line(data = data_plot_slope, aes(x=Quantile, y=Slope))+
  geom_point(data = data_plot_slope, aes(x=Quantile, y=Slope))+
  theme(axis.line = element_line(color="black"))+
  theme(strip.background = element_blank())+
  theme(panel.background = element_rect(fill=alpha("#C8AD7F",0.7)))+
  theme(panel.grid.major = element_line(colour="#C8AD7F"),
        panel.grid.minor = element_line(colour="#C8AD7F"),
        panel.border = element_rect(colour="black",fill=NA))+
  theme(legend.key = element_blank())+
  labs(x="Quantiles",y= "Slope")+
  ggtitle("Slope of QR (Cluster Radius with CTI) ")
#ggsave(plot = plot_cti, file="Output/Cluster/Slope_CTI_Uncertainty.png",width=7,height=4,dpi=300)

# 2) Scatter Plots with Regression Lines
data_scatter = data.frame(Y = Y,
                          X = X,
                          Aridity = aridity_temp)
data_scatter$Aridity = factor(data_scatter$Aridity, levels=c("Hyper Arid","Arid","Semi-Arid","Dry Sub-Humid","Humid"))
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
X = X[-which(is.na(X))]
# plot the scatter diagram
plot_cti_scatter = ggplot()+
  geom_point(data=data_scatter, aes(x=X, y=Y, color=Aridity), shape=1)+
  scale_color_manual(values = c("#d11d1d","#b2642c","#37842a","#3faec1","#1f5b93"))+
  geom_line(data = data_regline, aes(x=X,y=Y,linetype = Quantile))+
  coord_cartesian(xlim=c(quantile(X,probs = 0.0),quantile(X,probs = 0.9)), ylim=c(0,quantile(Y,probs=0.9)))+
  theme(axis.line = element_line(color="black"))+
  theme(strip.background = element_blank())+
  theme(panel.background = element_rect(fill=alpha("#C8AD7F",0.5)))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour="black",fill=NA))+
  theme(legend.key = element_blank())+
  labs(x="CTI",y= "Cluster Radius")+
  ggtitle("Scatter of Cluster Radius vs CTI")
#ggsave(plot = plot_cti_scatter, file="Output/Cluster/Scatter_CTI_Uncertainty.png",width=6,height=4,dpi=300)
##----------------------------------------------------------------------------------------------------------##

##----------------------------------------------------------------------------------------------------------##
# 4) NDVI
print("Processing Regression 4")
Y = cbind(cluster_metric$Mean_Radius)
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

# Plot the required quantile regression plots
# 1) Regression Slopes
data_plot_slope = data.frame(Quantile = seq(from=0.05, to=0.95, by=0.1),
                             Slope = QR_ndvi$coefficients[2,])
# collect all the standard errors
std_err_ndvi = NULL
for (i in 1:length(sumQR_elev)){
  std_err_ndvi = c(std_err_ndvi, sumQR_ndvi[[i]]$coefficients[2,2])
}
data_plot_ribbon = data.frame(Quantile = seq(from=0.05, to=0.95, by=0.1),
                              Ymin = QR_ndvi$coefficients[2,] - std_err_ndvi,
                              Ymax = QR_ndvi$coefficients[2,] + std_err_ndvi)
plot_ndvi = ggplot()+
  geom_ribbon(data = data_plot_ribbon, aes(x=Quantile, ymin=Ymin, ymax=Ymax),fill="#77916e",alpha=0.6)+
  geom_line(data = data_plot_slope, aes(x=Quantile, y=Slope))+
  geom_point(data = data_plot_slope, aes(x=Quantile, y=Slope))+
  theme(axis.line = element_line(color="black"))+
  theme(strip.background = element_blank())+
  theme(panel.background = element_rect(fill=alpha("#C8AD7F",0.7)))+
  theme(panel.grid.major = element_line(colour="#C8AD7F"),
        panel.grid.minor = element_line(colour="#C8AD7F"),
        panel.border = element_rect(colour="black",fill=NA))+
  theme(legend.key = element_blank())+
  labs(x="Quantiles",y= "Slope")+
  ggtitle("Slope of QR (Cluster Radius with NDVI) ")
#ggsave(plot = plot_ndvi, file="Output/Cluster/Slope_NDVI_Uncertainty.png",width=7,height=4,dpi=300)

# 2) Scatter Plots with Regression Lines
data_scatter = data.frame(Y = Y,
                          X = X,
                          Aridity = aridity_temp)
data_scatter$Aridity = factor(data_scatter$Aridity, levels=c("Hyper Arid","Arid","Semi-Arid","Dry Sub-Humid","Humid"))
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
  geom_point(data=data_scatter, aes(x=X, y=Y, color=Aridity), shape=1)+
  scale_color_manual(values = c("#d11d1d","#b2642c","#37842a","#3faec1","#1f5b93"))+
  geom_line(data = data_regline, aes(x=X,y=Y,linetype = Quantile))+
  coord_cartesian(xlim=c(0,quantile(X,probs = 0.9)), ylim=c(0,quantile(Y,probs=0.9)))+
  theme(axis.line = element_line(color="black"))+
  theme(strip.background = element_blank())+
  theme(panel.background = element_rect(fill=alpha("#C8AD7F",0.5)))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour="black",fill=NA))+
  theme(legend.key = element_blank())+
  labs(x="NDVI",y= "Cluster Radius")+
  ggtitle("Scatter of Cluster Radius vs NDVI")
#ggsave(plot = plot_ndvi_scatter, file="Output/Cluster/Scatter_NDVI_Uncertainty.png",width=6,height=4,dpi=300)


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

# Plot the required quantile regression plots
# 1) Regression Slopes
data_plot_slope = data.frame(Quantile = seq(from=0.05, to=0.95, by=0.1),
                             Slope = QR_arid_err$coefficients[2,])
# collect all the standard errors
std_err_arid_err = NULL
for (i in 1:length(sumQR_arid_err)){
  std_err_arid_err = c(std_err_arid_err, sumQR_arid_err[[i]]$coefficients[2,2])
}
data_plot_ribbon = data.frame(Quantile = seq(from=0.05, to=0.95, by=0.1),
                              Ymin = QR_arid_err$coefficients[2,] - std_err_arid_err,
                              Ymax = QR_arid_err$coefficients[2,] + std_err_arid_err)
plot_arid_err = ggplot()+
  geom_ribbon(data = data_plot_ribbon, aes(x=Quantile, ymin=Ymin, ymax=Ymax),fill="#77916e",alpha=0.6)+
  geom_line(data = data_plot_slope, aes(x=Quantile, y=Slope))+
  geom_point(data = data_plot_slope, aes(x=Quantile, y=Slope))+
  theme(axis.line = element_line(color="black"))+
  theme(strip.background = element_blank())+
  theme(panel.background = element_rect(fill=alpha("#C8AD7F",0.7)))+
  theme(panel.grid.major = element_line(colour="#C8AD7F"),
        panel.grid.minor = element_line(colour="#C8AD7F"),
        panel.border = element_rect(colour="black",fill=NA))+
  theme(legend.key = element_blank())+
  labs(x="Quantiles",y= "Slope")+
  ggtitle("Slope of QR (Budyko Error with Aridity) ")
#ggsave(plot = plot_arid_err, file="Output/Cluster/Slope_Arid_Error.png",width=7,height=4,dpi=300)

# 2) Scatter Plots with Regression Lines
data_scatter = data.frame(Y = Y,
                          X = X,
                          Aridity = aridity_temp)
data_scatter$Aridity = factor(data_scatter$Aridity, levels=c("Hyper Arid","Arid","Semi-Arid","Dry Sub-Humid","Humid"))
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
  geom_point(data=data_scatter, aes(x=X, y=Y, color=Aridity), shape=1)+
  scale_color_manual(values = c("#d11d1d","#b2642c","#37842a","#3faec1","#1f5b93"))+
  geom_line(data = data_regline, aes(x=X,y=Y,linetype = Quantile))+
  coord_cartesian(xlim=c(0,quantile(X,probs = 0.9)), ylim=c(0,quantile(Y,probs=0.9)))+
  theme(axis.line = element_line(color="black"))+
  theme(strip.background = element_blank())+
  theme(panel.background = element_rect(fill=alpha("#C8AD7F",0.5)))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour="black",fill=NA))+
  theme(legend.key = element_blank())+
  labs(x="Aridity",y= "Budyko Error")+
  ggtitle("Scatter of Budyko Error vs Aridity")
#ggsave(plot = plot_arid_err_scatter, file="Output/Cluster/Scatter_Arid_Error.png",width=6,height=4,dpi=300)
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
# 1) Regression Slopes
data_plot_slope = data.frame(Quantile = seq(from=0.05, to=0.95, by=0.1),
                             Slope = QR_elev_err$coefficients[2,])
# collect all the standard errors
std_err_elev_err = NULL
for (i in 1:length(sumQR_elev_err)){
  std_err_elev_err = c(std_err_elev_err, sumQR_elev_err[[i]]$coefficients[2,2])
}
data_plot_ribbon = data.frame(Quantile = seq(from=0.05, to=0.95, by=0.1),
                              Ymin = QR_elev_err$coefficients[2,] - std_err_elev_err,
                              Ymax = QR_elev_err$coefficients[2,] + std_err_elev_err)
plot_elev_err = ggplot()+
  geom_ribbon(data = data_plot_ribbon, aes(x=Quantile, ymin=Ymin, ymax=Ymax),fill="#77916e",alpha=0.6)+
  geom_line(data = data_plot_slope, aes(x=Quantile, y=Slope))+
  geom_point(data = data_plot_slope, aes(x=Quantile, y=Slope))+
  theme(axis.line = element_line(color="black"))+
  theme(strip.background = element_blank())+
  theme(panel.background = element_rect(fill=alpha("#C8AD7F",0.7)))+
  theme(panel.grid.major = element_line(colour="#C8AD7F"),
        panel.grid.minor = element_line(colour="#C8AD7F"),
        panel.border = element_rect(colour="black",fill=NA))+
  theme(legend.key = element_blank())+
  labs(x="Quantiles",y= "Slope")+
  ggtitle("Slope of QR (Budyko Error with Elevation) ")
#ggsave(plot = plot_elev_err, file="Output/Cluster/Slope_Elevation_Error.png",width=7,height=4,dpi=300)

# 2) Scatter Plots with Regression Lines
data_scatter = data.frame(Y = Y,
                          X = X,
                          Aridity = aridity_temp)
data_scatter$Aridity = factor(data_scatter$Aridity, levels=c("Hyper Arid","Arid","Semi-Arid","Dry Sub-Humid","Humid"))
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
  geom_point(data=data_scatter, aes(x=X, y=Y, color=Aridity), shape=1)+
  scale_color_manual(values = c("#d11d1d","#b2642c","#37842a","#3faec1","#1f5b93"))+
  geom_line(data = data_regline, aes(x=X,y=Y,linetype = Quantile))+
  coord_cartesian(xlim=c(0,quantile(X,probs = 0.9)), ylim=c(0,quantile(Y,probs=0.9)))+
  theme(axis.line = element_line(color="black"))+
  theme(strip.background = element_blank())+
  theme(panel.background = element_rect(fill=alpha("#C8AD7F",0.5)))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour="black",fill=NA))+
  theme(legend.key = element_blank())+
  labs(x="Elevation",y= "Budyko Error")+
  ggtitle("Scatter of Budyko Error vs Elevation")
#ggsave(plot = plot_elev_err_scatter, file="Output/Cluster/Scatter_Elevation_Error.png",width=6,height=4,dpi=300)

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
# 1) Regression Slopes
data_plot_slope = data.frame(Quantile = seq(from=0.05, to=0.95, by=0.1),
                             Slope = QR_cti_err$coefficients[2,])
# collect all the standard errors
std_err_cti_err = NULL
for (i in 1:length(sumQR_cti_err)){
  std_err_cti_err = c(std_err_cti_err, sumQR_cti_err[[i]]$coefficients[2,2])
}
data_plot_ribbon = data.frame(Quantile = seq(from=0.05, to=0.95, by=0.1),
                              Ymin = QR_cti_err$coefficients[2,] - std_err_cti_err,
                              Ymax = QR_cti_err$coefficients[2,] + std_err_cti_err)
plot_cti_err = ggplot()+
  geom_ribbon(data = data_plot_ribbon, aes(x=Quantile, ymin=Ymin, ymax=Ymax),fill="#77916e",alpha=0.6)+
  geom_line(data = data_plot_slope, aes(x=Quantile, y=Slope))+
  geom_point(data = data_plot_slope, aes(x=Quantile, y=Slope))+
  theme(axis.line = element_line(color="black"))+
  theme(strip.background = element_blank())+
  theme(panel.background = element_rect(fill=alpha("#C8AD7F",0.7)))+
  theme(panel.grid.major = element_line(colour="#C8AD7F"),
        panel.grid.minor = element_line(colour="#C8AD7F"),
        panel.border = element_rect(colour="black",fill=NA))+
  theme(legend.key = element_blank())+
  labs(x="Quantiles",y= "Slope")+
  ggtitle("Slope of QR (Budyko Error with CTI) ")
#ggsave(plot = plot_cti_err, file="Output/Cluster/Slope_CTI_Error.png",width=7,height=4,dpi=300)

# 2) Scatter Plots with Regression Lines
data_scatter = data.frame(Y = Y,
                          X = X,
                          Aridity = aridity_temp)
data_scatter$Aridity = factor(data_scatter$Aridity, levels=c("Hyper Arid","Arid","Semi-Arid","Dry Sub-Humid","Humid"))
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
  geom_point(data=data_scatter, aes(x=X, y=Y, color=Aridity), shape=1)+
  scale_color_manual(values = c("#d11d1d","#b2642c","#37842a","#3faec1","#1f5b93"))+
  geom_line(data = data_regline, aes(x=X,y=Y,linetype = Quantile))+
  coord_cartesian(xlim=c(quantile(X,probs = 0.0),quantile(X,probs = 0.9)), ylim=c(0,quantile(Y,probs=0.9)))+
  theme(axis.line = element_line(color="black"))+
  theme(strip.background = element_blank())+
  theme(panel.background = element_rect(fill=alpha("#C8AD7F",0.5)))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour="black",fill=NA))+
  theme(legend.key = element_blank())+
  labs(x="CTI",y= "Budyko Error")+
  ggtitle("Scatter of Budyko Error vs CTI")
#ggsave(plot = plot_cti_err_scatter, file="Output/Cluster/Scatter_CTI_Error.png",width=6,height=4,dpi=300)


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
# 1) Regression Slopes
data_plot_slope = data.frame(Quantile = seq(from=0.05, to=0.95, by=0.1),
                             Slope = QR_ndvi_err$coefficients[2,])
# collect all the standard errors
std_err_ndvi_err = NULL
for (i in 1:length(sumQR_ndvi_err)){
  std_err_ndvi_err = c(std_err_ndvi_err, sumQR_ndvi_err[[i]]$coefficients[2,2])
}
data_plot_ribbon = data.frame(Quantile = seq(from=0.05, to=0.95, by=0.1),
                              Ymin = QR_ndvi_err$coefficients[2,] - std_err_ndvi_err,
                              Ymax = QR_ndvi_err$coefficients[2,] + std_err_ndvi_err)
plot_ndvi_err = ggplot()+
  geom_ribbon(data = data_plot_ribbon, aes(x=Quantile, ymin=Ymin, ymax=Ymax),fill="#77916e",alpha=0.6)+
  geom_line(data = data_plot_slope, aes(x=Quantile, y=Slope))+
  geom_point(data = data_plot_slope, aes(x=Quantile, y=Slope))+
  theme(axis.line = element_line(color="black"))+
  theme(strip.background = element_blank())+
  theme(panel.background = element_rect(fill=alpha("#C8AD7F",0.7)))+
  theme(panel.grid.major = element_line(colour="#C8AD7F"),
        panel.grid.minor = element_line(colour="#C8AD7F"),
        panel.border = element_rect(colour="black",fill=NA))+
  theme(legend.key = element_blank())+
  labs(x="Quantiles",y= "Slope")+
  ggtitle("Slope of QR (Budyko Error with NDVI) ")
#ggsave(plot = plot_ndvi_err, file="Output/Cluster/Slope_NDVI_Error.png",width=7,height=4,dpi=300)

# 2) Scatter Plots with Regression Lines
data_scatter = data.frame(Y = Y,
                          X = X,
                          Aridity = aridity_temp)
data_scatter$Aridity = factor(data_scatter$Aridity, levels=c("Hyper Arid","Arid","Semi-Arid","Dry Sub-Humid","Humid"))
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
  geom_point(data=data_scatter, aes(x=X, y=Y, color=Aridity), shape=1)+
  scale_color_manual(values = c("#d11d1d","#b2642c","#37842a","#3faec1","#1f5b93"))+
  geom_line(data = data_regline, aes(x=X,y=Y,linetype = Quantile))+
  coord_cartesian(xlim=c(0,quantile(X,probs = 0.9)), ylim=c(0,quantile(Y,probs=0.9)))+
  theme(axis.line = element_line(color="black"))+
  theme(strip.background = element_blank())+
  theme(panel.background = element_rect(fill=alpha("#C8AD7F",0.5)))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour="black",fill=NA))+
  theme(legend.key = element_blank())+
  labs(x="NDVI",y= "Budyko Error")+
  ggtitle("Scatter of Budyko Error vs NDVI")
#ggsave(plot = plot_ndvi_err_scatter, file="Output/Cluster/Scatter_NDVI_Error.png",width=6,height=4,dpi=300)

##----------------------------------------------------------------------------------------------------------##

# Plot of ET cluster density vs P cluster density with aridity
cluster_sens = read.table("Output/Budyko_Cluster_Density_Sensitivity_noSM2RAIN.txt",header=TRUE)
cluster_sens = cluster_sens/ai
cluster_sens$CD_P_Mean[which(cluster_sens$CD_P_Mean==0)]= NA
cluster_sens$CD_ET_Mean[which(cluster_sens$CD_ET_Mean==0)]= NA

data_cluster_sens = data.frame(P = cluster_sens$CD_P_Mean,
                               ET = cluster_sens$CD_ET_Mean,
                               Aridity = aridity)
data_cluster_sens$Aridity = factor(data_cluster_sens$Aridity, levels=c("Hyper Arid","Arid","Semi-Arid","Dry Sub-Humid","Humid"))

plot_cluster_sens_scatter = ggplot()+
  geom_point(data=data_cluster_sens, aes(x=ET, y=P, color=Aridity), shape=1)+
  geom_abline()+
  scale_color_manual(values = c("#d11d1d","#b2642c","#37842a","#3faec1","#1f5b93"))+
  coord_cartesian(xlim=c(0,quantile(data_cluster_sens$ET,probs = 0.90,na.rm=TRUE)), ylim=c(0,quantile(data_cluster_sens$P,probs=0.90,na.rm=TRUE)))+
  theme(axis.line = element_line(color="black"))+
  theme(strip.background = element_blank())+
  theme(panel.background = element_rect(fill=alpha("lightgrey",0.3)))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour="black",fill=NA))+
  theme(legend.key = element_blank())+
  labs(x="Cluster Radius (Mean ET)",y= "Cluster Radius (Mean P)")
  #ggtitle("Scatter of Budyko Error vs NDVI")
ggsave(plot = plot_cluster_sens_scatter, file="Output/Cluster/Scatter_Cluster_Sensitivity.png",width=6,height=4,dpi=300)


# Plot the coefficient of variation of P and ET calculated from absolute values of uncertainty
abs_uncertainty = read.table("Output/Absolute_P_ET_Uncertainty_noSM2RAIN.txt", header=TRUE)

data_abs_uncertainty = data.frame(P  = abs_uncertainty$P_SD/abs_uncertainty$P_Mean,
                                  ET = abs_uncertainty$ET_SD/abs_uncertainty$ET_Mean,
                                  Aridity = aridity)
data_abs_uncertainty$Aridity = factor(data_abs_uncertainty$Aridity, levels=c("Hyper Arid","Arid","Semi-Arid","Dry Sub-Humid","Humid"))
data_abs_uncertainty[which(data_abs_uncertainty$ET>10),"ET"] = NA

plot_abs_uncer = ggplot()+
  geom_point(data=data_abs_uncertainty, aes(x=ET, y=P, color=Aridity), shape=1)+
  geom_abline()+
  scale_color_manual(breaks=as.character(c("Hyper Arid","Arid","Semi-Arid","Dry Sub-Humid","Humid")),values = c("#d11d1d","#b2642c","#37842a","#3faec1","#1f5b93"))+
  coord_cartesian(xlim = c(0,2.5), ylim=c(0,2.5))+
  #coord_cartesian(xlim=c(0,quantile(data_cluster_sens$P,probs = 0.90,na.rm=TRUE)), ylim=c(0,quantile(data_cluster_sens$ET,probs=0.90,na.rm=TRUE)))+
  theme(axis.line = element_line(color="black"))+
  theme(text = element_text(family="Times", size=12))+
  theme(strip.background = element_blank())+
  theme(panel.background = element_rect(fill=alpha("lightgrey",0.3)))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour="black",fill=NA))+
  theme(legend.key = element_blank())+
  theme(legend.position = "bottom")+
  theme(legend.title = element_text(face="bold"))+
  guides(fill=guide_legend(title="Aridity:", face="bold",nrow=1,byrow = TRUE))+
  labs(x="ET",y= "P")
  #ggtitle("Coefficient of Variation")
ggsave(plot = plot_abs_uncer, file="Output/Cluster/Scatter_Absolute_Uncertainty.png",width=150,height=160,dpi=300, units="mm")

# Plot the absolute standard deviation
data_std_uncertainty = data.frame(P  = abs_uncertainty$P_SD,
                                  ET = abs_uncertainty$ET_SD,
                                  Aridity = aridity)
data_std_uncertainty$Aridity = factor(data_std_uncertainty$Aridity, levels=c("Hyper Arid","Arid","Semi-Arid","Dry Sub-Humid","Humid"))

plot_std_uncer = ggplot()+
  geom_point(data=data_std_uncertainty, aes(x=ET, y=P, color=Aridity), shape=1)+
  geom_abline()+
  scale_color_manual(values = c("#d11d1d","#b2642c","#37842a","#3faec1","#1f5b93"))+
  #coord_cartesian(xlim=c(0,quantile(data_cluster_sens$P,probs = 0.90,na.rm=TRUE)), ylim=c(0,quantile(data_cluster_sens$ET,probs=0.90,na.rm=TRUE)))+
  theme(axis.line = element_line(color="black"))+
  theme(strip.background = element_blank())+
  theme(panel.background = element_rect(fill=alpha("#C8AD7F",0.5)))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour="black",fill=NA))+
  theme(legend.key = element_blank())+
  labs(x="ET",y= "P")+
  ggtitle("Standard Deviation")
#ggsave(plot = plot_std_uncer, file="Output/Cluster/Scatter_SD_Uncertainty.png",width=6,height=4,dpi=300)

# Plot of ET cluster density vs P cluster density with aridity
cluster_sens = read.table("Output/Budyko_Cluster_Density_Sensitivity.txt",header=TRUE)
cluster_sens = cluster_sens/ai

data_cluster_sens = data.frame(P = cluster_sens$CD_P_Mean,
                               ET = cluster_sens$CD_ET_Mean,
                               Aridity = aridity)
data_cluster_sens$ET[which(data_cluster_sens$ET>0.3)] = NA
data_cluster_sens$ET[which(data_cluster_sens$P>0.3)] = NA
data_cluster_sens$Aridity = factor(data_cluster_sens$Aridity, levels=c("Hyper Arid","Arid","Semi-Arid","Dry Sub-Humid","Humid"))

plot_cluster_sens_scatter = ggplot()+
  geom_point(data=data_cluster_sens, aes(x=ET, y=P, color=Aridity), shape=1)+
  geom_abline()+
  scale_color_manual(breaks=as.character(c("Hyper Arid","Arid","Semi-Arid","Dry Sub-Humid","Humid")),values = c("#d11d1d","#b2642c","#37842a","#3faec1","#1f5b93"))+
  #coord_cartesian(xlim=c(0,quantile(data_cluster_sens$P,probs = 0.90,na.rm=TRUE)), ylim=c(0,quantile(data_cluster_sens$ET,probs=0.90,na.rm=TRUE)))+
  coord_cartesian(xlim=c(0,0.4),ylim=c(0,0.4))+
  theme(axis.line = element_line(color="black"))+
  theme(text = element_text(family="Times", size=12))+
  theme(strip.background = element_blank())+
  theme(panel.background = element_rect(fill=alpha("lightgrey",0.3)))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour="black",fill=NA))+
  theme(legend.key = element_blank())+
  theme(legend.position = "bottom")+
  theme(legend.title = element_text(face="bold"))+
  guides(fill=guide_legend(title="Aridity:", face="bold",nrow=1,byrow = TRUE))+
  labs(x="Cluster Radius ET",y= "Cluster Radius P")
#ggtitle("Scatter of Budyko Error vs NDVI")
ggsave(plot = plot_cluster_sens_scatter, file="Output/Cluster/Scatter_Cluster_Sensitivity.png",width=150,height=140,dpi=300,units="mm")


# Plot the coefficient of variation of P and ET calculated from absolute values of uncertainty
abs_uncertainty = read.table("Output/Absolute_P_ET_Uncertainty.txt", header=TRUE)

data_abs_uncertainty = data.frame(P  = abs_uncertainty$P_SD/abs_uncertainty$P_Mean,
                                  ET = abs_uncertainty$ET_SD/abs_uncertainty$ET_Mean,
                                  Aridity = aridity)
data_abs_uncertainty$Aridity = factor(data_abs_uncertainty$Aridity, levels=c("Hyper Arid","Arid","Semi-Arid","Dry Sub-Humid","Humid"))

plot_abs_uncer = ggplot()+
  geom_point(data=data_abs_uncertainty, aes(x=ET, y=P, color=Aridity), shape=1)+
  geom_abline()+
  scale_color_manual(values = c("#d11d1d","#b2642c","#37842a","#3faec1","#1f5b93"))+
  #coord_cartesian(xlim=c(0,quantile(data_cluster_sens$P,probs = 0.90,na.rm=TRUE)), ylim=c(0,quantile(data_cluster_sens$ET,probs=0.90,na.rm=TRUE)))+
  theme(axis.line = element_line(color="black"))+
  theme(strip.background = element_blank())+
  theme(panel.background = element_rect(fill=alpha("#C8AD7F",0.5)))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour="black",fill=NA))+
  theme(legend.key = element_blank())+
  labs(x="ET",y= "P")+
  ggtitle("Coefficient of Variation")
ggsave(plot = plot_abs_uncer, file="Output/Cluster/Scatter_Absolute_Uncertainty.png",width=6,height=4,dpi=300)

# Plot the absolute standard deviation
data_std_uncertainty = data.frame(P  = abs_uncertainty$P_SD,
                                  ET = abs_uncertainty$ET_SD,
                                  Aridity = aridity)
data_std_uncertainty$Aridity = factor(data_std_uncertainty$Aridity, levels=c("Hyper Arid","Arid","Semi-Arid","Dry Sub-Humid","Humid"))

plot_std_uncer = ggplot()+
  geom_point(data=data_std_uncertainty, aes(x=ET, y=P, color=Aridity), shape=1)+
  geom_abline()+
  scale_color_manual(values = c("#d11d1d","#b2642c","#37842a","#3faec1","#1f5b93"))+
  #coord_cartesian(xlim=c(0,quantile(data_cluster_sens$P,probs = 0.90,na.rm=TRUE)), ylim=c(0,quantile(data_cluster_sens$ET,probs=0.90,na.rm=TRUE)))+
  theme(axis.line = element_line(color="black"))+
  theme(strip.background = element_blank())+
  theme(panel.background = element_rect(fill=alpha("#C8AD7F",0.5)))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour="black",fill=NA))+
  theme(legend.key = element_blank())+
  labs(x="ET",y= "P")+
  ggtitle("Standard Deviation")
ggsave(plot = plot_std_uncer, file="Output/Cluster/Scatter_SD_Uncertainty.png",width=6,height=4,dpi=300)









