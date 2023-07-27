################################################################################
#################### eurasia scale analyses revision ###########################
################################################################################
require(quantreg)
require(Qtools)
require(vegan)
require(ggplot2)
require(qgam)

setwd("C:/Users/yaqchang.D/OneDrive - ETH Zurich/PhD-chapter 1")
load("01-data/cleaned data/RData/revision_eurasia_0313.RData")
sr_eurasia_revision <- read.csv("01-data/cleaned data/table/sr_eurasia_revision.csv")
HDS_template <- raster("01-data/cleaned data/raster/HDS_template_1k.tif")
sr_eurasia_revision <- merge(sr_eurasia_revision,kew_eurasia_hete_shannon2,by="kew")
kew_family <- colnames(sr_eurasia_revision[4:100])
kew_region_tot <- read.csv("01-data/cleaned data/table/kew_region_tot.csv")
sr_eurasia_hetero <- sr_eurasia_revision[,c("area","LEVEL3_COD","gdd","rock","swb","compound","hetero_shannon")] 
sr_eurasia_revision <- merge(kew_region_tot,sr_eurasia_hetero,by="LEVEL3_COD")
kew_family <- colnames(sr_eurasia_revision[7:114])

#--------------------------------------------------------
# making exploration points to compare different methods
#--------------------------------------------------------

# method 1: sr~residual(heterogeneity~area)

plot(sr_eurasia_revision$area,sr_eurasia_revision$compound)
kew_family <- colnames(sr_eurasia_revision[3:99])

for (i in 1:length(kew_family)) {
  print(i)
  
dat_sub <- sr_eurasia_revision[,c(kew_family[i],"area","gdd","swb","rock","compound")]
chc <- which(dat_sub$compound>4.3284058&dat_sub$compound<4.328406)
dat_sub <- na.omit(dat_sub)
head(dat_sub)

tiff(paste("03-results/figure_2/revision/",kew_family[i],"_area_resi",".tif",sep = ""),width = 40,height = 20,units = "cm",res=300,compression = "lzw")

par(mfrow=c(2,2))

nqr_50 <- nlrq(compound~SSarrhenius(area, k,z),data = dat_sub,tau=0.50,trace=TRUE)
plot(residuals(nqr_50),dat_sub[,1],main=paste(kew_family[i],"compound residual"))
points(x=residuals(nqr_50)[chc],y=dat_sub[chc,1],col="red")

nqr_50 <- nlrq(gdd~SSarrhenius(area, k,z),data = dat_sub,tau=0.50,trace=TRUE)
plot(residuals(nqr_50),dat_sub[,1],main=paste(kew_family[i],"gdd residual"))
points(x=residuals(nqr_50)[chc],y=dat_sub[chc,1],col="red")

nqr_50 <- nlrq(swb~SSarrhenius(area, k,z),data = dat_sub,tau=0.50,trace=TRUE)
plot(residuals(nqr_50),dat_sub[,1],main=paste(kew_family[i],"swb residual"))
points(x=residuals(nqr_50)[chc],y=dat_sub[chc,1],col="red")

nqr_50 <- nlrq(rock~SSarrhenius(area, k,z),data = dat_sub,tau=0.50,trace=TRUE)
plot(residuals(nqr_50),dat_sub[,1],main=paste(kew_family[i],"rock residual"))
points(x=residuals(nqr_50)[chc],y=dat_sub[chc,1],col="red")
dev.off()
}

# method 2: residual (sr~area) ~ heterogeneity

for (i in 1:length(kew_family)) {
  print(i)
  
  dat_sub <- sr_eurasia_revision[,c(kew_family[i],"area","gdd","swb","rock","compound")]
  dat_sub <- na.omit(dat_sub)
  chc <- which(dat_sub$compound>4.3284058&dat_sub$compound<4.328406)
  names(dat_sub)[1] <- "resp"
  
  area_50 <- nlrq(resp~SSarrhenius(area, k,z),data = dat_sub,tau=0.50,trace=TRUE)
  hist(residuals(area_50))
  
  tiff(paste("03-results/figure_2/revision/",kew_family[i],"_richness_resi",".tif",sep = ""),width = 40,height = 20,units = "cm",res=300,compression = "lzw")
  par(mfrow=c(2,2))
  
  resi=log(residuals(area_50)+abs(min(residuals(area_50)))+1)
  plot(dat_sub$compound,resi,main=paste(kew_family[i],"compound to area residual"))
  points(x=dat_sub[chc,"compound"],y=resi[chc],col="red")
  plot(dat_sub$gdd,resi,main=paste(kew_family[i],"gdd to area residual"))
  points(x=dat_sub[chc,"gdd"],y=resi[chc],col="red")
  plot(dat_sub$swb,resi,main=paste(kew_family[i],"swb to area residual"))
  points(x=dat_sub[chc,"swb"],y=resi[chc],col="red")
  plot(dat_sub$rock,resi,main=paste(kew_family[i],"rock to area residual"))
  points(x=dat_sub[chc,"rock"],y=resi[chc],col="red")
dev.off()

}

#---------------------------------
# non-linear quantile regression
#---------------------------------
compound_qnt_sum <- data.frame()
gdd_qnt_sum <- data.frame()
swb_qnt_sum <- data.frame()
rock_qnt_sum <- data.frame()
full_qnt_sum <- data.frame()
for (i in 1:length(kew_family)) {
  print(kew_family[i])
  compound_kew_qnt <- c(family=kew_family[i],kew_qnt_mod(sr_eurasia_revision,loc = "LEVEL3_COD",pred = "compound",resp = kew_family[i],mod="idiv"))
  compound_qnt_sum <- rbind(compound_qnt_sum,compound_kew_qnt)
  names(compound_qnt_sum) <- c("family","qnt025","qnt050","qnt975","SR","diff","std_diff","pseudo_r_05","pseudo_r_50","pseudo_r_95","RSS_05","RSS_50","RSS_95","null_resi_05","null_resi_50","null_resi_95")
  gdd_kew_qnt <- c(family=kew_family[i],kew_qnt_mod(sr_eurasia_revision,loc = "LEVEL3_COD",pred = "gdd",resp = kew_family[i],mod="idiv"))
  gdd_qnt_sum <- rbind(gdd_qnt_sum,gdd_kew_qnt)
  names(gdd_qnt_sum) <- c("family","qnt025","qnt050","qnt975","SR","diff","std_diff","pseudo_r_05","pseudo_r_50","pseudo_r_95","RSS_05","RSS_50","RSS_95","null_resi_05","null_resi_50","null_resi_95")
  swb_kew_qnt <- c(family=kew_family[i],kew_qnt_mod(sr_eurasia_revision,loc = "LEVEL3_COD",pred = "swb",resp = kew_family[i],mod="idiv"))
  swb_qnt_sum <- rbind(swb_qnt_sum,swb_kew_qnt)
  names(swb_qnt_sum) <- c("family","qnt025","qnt050","qnt975","SR","diff","std_diff","pseudo_r_05","pseudo_r_50","pseudo_r_95","RSS_05","RSS_50","RSS_95","null_resi_05","null_resi_50","null_resi_95")
  rock_kew_qnt <- c(family=kew_family[i],kew_qnt_mod(sr_eurasia_revision,loc = "LEVEL3_COD",pred = "rock",resp = kew_family[i],mod="idiv"))
  rock_qnt_sum <- rbind(rock_qnt_sum,rock_kew_qnt)
  names(rock_qnt_sum) <- c("family","qnt025","qnt050","qnt975","SR","diff","std_diff","pseudo_r_05","pseudo_r_50","pseudo_r_95","RSS_05","RSS_50","RSS_95","null_resi_05","null_resi_50","null_resi_95")
  full_kew_qnt <- c(family=kew_family[i],kew_qnt_mod(sr_eurasia_revision,loc = "LEVEL3_COD",pred = c("gdd","swb","rock"),resp = kew_family[i],mod="full"))
  full_qnt_sum <- rbind(full_qnt_sum,full_kew_qnt)
  names(full_qnt_sum) <- c("family","qnt025","qnt050","qnt975","SR","diff","std_diff","pseudo_r_05","pseudo_r_50","pseudo_r_95","RSS_05","RSS_50","RSS_95","null_resi_05","null_resi_50","null_resi_95")
  
}
write.csv(compound_qnt_sum,"01-data/cleaned data/table/compound_qnt_sum.csv")
write.csv(gdd_qnt_sum,"01-data/cleaned data/table/gdd_qnt_sum.csv")
write.csv(swb_qnt_sum,"01-data/cleaned data/table/swb_qnt_sum.csv")
write.csv(rock_qnt_sum,"01-data/cleaned data/table/rock_qnt_sum.csv")
write.csv(full_qnt_sum,"01-data/cleaned data/table/full_qnt_sum.csv")

full_qnt_outlier <- as.data.frame(cbind(family=full_qnt_sum[as.numeric(full_qnt_sum$diff)>0,"family"],full=1))
compound_qnt_outlier <- as.data.frame(cbind(family=compound_qnt_sum[as.numeric(compound_qnt_sum$diff)>0,"family"],compound=1))
gdd_qnt_outlier <- as.data.frame(cbind(family=gdd_qnt_sum[as.numeric(gdd_qnt_sum$diff)>0,"family"],gdd=1))
swb_qnt_outlier <- as.data.frame(cbind(family=swb_qnt_sum[as.numeric(swb_qnt_sum$diff)>0,"family"],swb=1))
rock_qnt_outlier <- as.data.frame(cbind(family=rock_qnt_sum[as.numeric(rock_qnt_sum$diff)>0,"family"],rock=1))

total_qnt_outlier <- full_join(full_qnt_outlier,compound_qnt_outlier,by="family")%>%
  full_join(.,gdd_qnt_outlier,by="family")%>%
  full_join(.,swb_qnt_outlier,by="family")%>%
  full_join(.,rock_qnt_outlier,by="family")
total_qnt_outlier$full <- as.numeric(total_qnt_outlier$full)
total_qnt_outlier$compound <- as.numeric(total_qnt_outlier$compound)
total_qnt_outlier$gdd <- as.numeric(total_qnt_outlier$gdd)
total_qnt_outlier$swb <- as.numeric(total_qnt_outlier$swb)
total_qnt_outlier$rock <- as.numeric(total_qnt_outlier$rock)

total_qnt_outlier$sum <- rowSums(total_qnt_outlier[,2:6],na.rm = T)
family_sr <- list.files("01-data/cleaned data/jpeg/fam_richness",pattern = ".tif")
file.copy(paste("01-data/cleaned data/jpeg/fam_richness/",total_qnt_outlier[which(total_qnt_outlier$sum==4),"family"],".tif",sep = ""),
          paste("01-data/cleaned data/jpeg/fam_richness/outlier_4"))
file.copy(paste("01-data/cleaned data/jpeg/fam_richness/",total_qnt_outlier[which(total_qnt_outlier$sum==3),"family"],".tif",sep = ""),
          paste("01-data/cleaned data/jpeg/fam_richness/outlier_3"))
file.copy(paste("01-data/cleaned data/jpeg/fam_richness/",total_qnt_outlier[which(total_qnt_outlier$sum==2),"family"],".tif",sep = ""),
          paste("01-data/cleaned data/jpeg/fam_richness/outlier_2"))
file.copy(paste("01-data/cleaned data/jpeg/fam_richness/",total_qnt_outlier[which(total_qnt_outlier$sum==1),"family"],".tif",sep = ""),
          paste("01-data/cleaned data/jpeg/fam_richness/outlier_1"))

file.copy(paste("01-data/cleaned data/jpeg/fam_richness/",total_qnt_outlier[which(total_qnt_outlier$full==1),"family"],".tif",sep = ""),
          paste("01-data/cleaned data/jpeg/fam_richness/outlier_full"))
file.copy(paste("01-data/cleaned data/jpeg/fam_richness/",total_qnt_outlier[which(total_qnt_outlier$compound==1),"family"],".tif",sep = ""),
          paste("01-data/cleaned data/jpeg/fam_richness/outlier_compound"))
file.copy(paste("01-data/cleaned data/jpeg/fam_richness/",total_qnt_outlier[which(total_qnt_outlier$gdd==1),"family"],".tif",sep = ""),
          paste("01-data/cleaned data/jpeg/fam_richness/outlier_gdd"))
file.copy(paste("01-data/cleaned data/jpeg/fam_richness/",total_qnt_outlier[which(total_qnt_outlier$swb==1),"family"],".tif",sep = ""),
          paste("01-data/cleaned data/jpeg/fam_richness/outlier_swb"))
file.copy(paste("01-data/cleaned data/jpeg/fam_richness/",total_qnt_outlier[which(total_qnt_outlier$rock==1),"family"],".tif",sep = ""),
          paste("01-data/cleaned data/jpeg/fam_richness/outlier_rock"))


# bootstrapping
compound_bootstrapping_sum <- data.frame()
gdd_bootstrapping_sum <- data.frame()
swb_bootstrapping_sum <- data.frame()
rock_bootstrapping_sum <- data.frame()
full_bootstrapping_sum <- data.frame()

for (i in 1:length(kew_family)) {
  print(kew_family[i])
  compound_kew_bootstrapping <- c(family=kew_family[i],kew_bootstrapping(sr_eurasia_revision,loc = "LEVEL3_COD",pred = "compound",resp = kew_family[i],mod="idiv"))
  compound_bootstrapping_sum <- rbind(compound_bootstrapping_sum,compound_kew_bootstrapping)
  names(compound_bootstrapping_sum) <- c("family","qnt025","qnt050","qnt975","SR","diff_upper","diff_lower","coefficient","p","explained_deviance","adj_R_square","AIC")
  gdd_kew_bootstrapping <- c(family=kew_family[i],kew_bootstrapping(sr_eurasia_revision,loc = "LEVEL3_COD",pred = "gdd",resp = kew_family[i],mod="idiv"))
  gdd_bootstrapping_sum <- rbind(gdd_bootstrapping_sum,gdd_kew_bootstrapping)
  names(gdd_bootstrapping_sum) <- c("family","qnt025","qnt050","qnt975","SR","diff_upper","diff_lower","coefficient","p","explained_deviance","adj_R_square","AIC")
  swb_kew_bootstrapping <- c(family=kew_family[i],kew_bootstrapping(sr_eurasia_revision,loc = "LEVEL3_COD",pred = "swb",resp = kew_family[i],mod="idiv"))
  swb_bootstrapping_sum <- rbind(swb_bootstrapping_sum,swb_kew_bootstrapping)
  names(swb_bootstrapping_sum) <- c("family","qnt025","qnt050","qnt975","SR","diff_upper","diff_lower","coefficient","p","explained_deviance","adj_R_square","AIC")
  rock_kew_bootstrapping <- c(family=kew_family[i],kew_bootstrapping(sr_eurasia_revision,loc = "LEVEL3_COD",pred = "rock",resp = kew_family[i],mod="idiv"))
  rock_bootstrapping_sum <- rbind(rock_bootstrapping_sum,rock_kew_bootstrapping)
  names(rock_bootstrapping_sum) <- c("family","qnt025","qnt050","qnt975","SR","diff_upper","diff_lower","coefficient","p","explained_deviance","adj_R_square","AIC")
  full_kew_bootstrapping <- c(family=kew_family[i],kew_bootstrapping(sr_eurasia_revision,loc = "LEVEL3_COD",pred = c("gdd","swb","rock"),resp = kew_family[i],mod="full"))
  full_bootstrapping_sum <- rbind(full_bootstrapping_sum,full_kew_bootstrapping)
  names(full_bootstrapping_sum) <- c("family","qnt025","qnt050","qnt975","SR","diff_upper","diff_lower","coefficient","p","explained_deviance","adj_R_square","AIC")
  
}
write.csv(compound_bootstrapping_sum,"01-data/cleaned data/table/compound_bootstrapping_sum.csv")
write.csv(gdd_bootstrapping_sum,"01-data/cleaned data/table/gdd_bootstrapping_sum.csv")
write.csv(swb_bootstrapping_sum,"01-data/cleaned data/table/swb_bootstrapping_sum.csv")
write.csv(rock_bootstrapping_sum,"01-data/cleaned data/table/rock_bootstrapping_sum.csv")
write.csv(full_bootstrapping_sum,"01-data/cleaned data/table/full_bootstrapping_sum.csv")

compound_bootstrapping_sum_ord <- compound_bootstrapping_sum[order(compound_bootstrapping_sum$diff,decreasing = T),]
file.copy(paste("01-data/cleaned data/jpeg/fam_richness/",compound_bootstrapping_sum_ord[1:24,"family"],".tif",sep = ""),
          paste("01-data/cleaned data/jpeg/fam_richness/outlier_compound"))

gdd_bootstrapping_sum_ord <- gdd_bootstrapping_sum[order(gdd_bootstrapping_sum$diff,decreasing = T),]
file.copy(paste("01-data/cleaned data/jpeg/fam_richness/",gdd_bootstrapping_sum_ord[1:24,"family"],".tif",sep = ""),
          paste("01-data/cleaned data/jpeg/fam_richness/outlier_gdd"))

swb_bootstrapping_sum_ord <- swb_bootstrapping_sum[order(swb_bootstrapping_sum$diff,decreasing = T),]
file.copy(paste("01-data/cleaned data/jpeg/fam_richness/",swb_bootstrapping_sum_ord[1:24,"family"],".tif",sep = ""),
          paste("01-data/cleaned data/jpeg/fam_richness/outlier_swb"))

rock_bootstrapping_sum_ord <- rock_bootstrapping_sum[order(rock_bootstrapping_sum$diff,decreasing = T),]
file.copy(paste("01-data/cleaned data/jpeg/fam_richness/",rock_bootstrapping_sum_ord[1:24,"family"],".tif",sep = ""),
          paste("01-data/cleaned data/jpeg/fam_richness/outlier_rock"))

full_bootstrapping_sum_ord <- full_bootstrapping_sum[order(full_bootstrapping_sum$diff,decreasing = T),]
file.copy(paste("01-data/cleaned data/jpeg/fam_richness/",full_bootstrapping_sum_ord[1:24,"family"],".tif",sep = ""),
          paste("01-data/cleaned data/jpeg/fam_richness/outlier_full"))

options(rasterMaxMemory = 2e10)

SR_mapping <- function(hetero,order){
SR= HDS_template
if(order==1){
  family_list <- hetero[1:24,"family"]
}else if (order==2){
  family_list <- hetero[25:48,"family"]
}else if (order==3){
  family_list <- hetero[49:72,"family"]
}else{
  family_list <- hetero[73:97,"family"]
}
# see how fast it goes
for(m in 1:length(family_list)){
  rst_m=readAll(raster(paste("01-data/cleaned data/raster/tif_family_update/",family_list[m],".tif",sep = "")))
  SR= SR + rst_m
  print(m)
} 
return(SR)
}

for (i in 1:4) {
  family_compound_sr <- SR_mapping(compound_bootstrapping_sum_ord,i)
  writeRaster(family_compound_sr,paste("01-data/cleaned data/raster/quantile_residual/compound_",i,".tif",sep = ""))
  tiff(paste("01-data/cleaned data/jpeg/fam_richness/quantile_residual/compound_",i,".tif",sep = ""),
       width = 30,height = 20,units = "cm",res=300,compression = "lzw")
  plot(family_compound_sr)
  dev.off()
  
  family_full_sr <- SR_mapping(full_bootstrapping_sum_ord,i)
  writeRaster(family_full_sr,paste("01-data/cleaned data/raster/quantile_residual/full_",i,".tif",sep = ""))
  tiff(paste("01-data/cleaned data/jpeg/fam_richness/quantile_residual/full_",i,".tif",sep = ""),
       width = 30,height = 20,units = "cm",res=300,compression = "lzw")
  plot(family_full_sr)
  dev.off()
  
  family_gdd_sr <- SR_mapping(gdd_bootstrapping_sum_ord,i)
  writeRaster(family_gdd_sr,paste("01-data/cleaned data/raster/quantile_residual/gdd_",i,".tif",sep = ""))
  tiff(paste("01-data/cleaned data/jpeg/fam_richness/quantile_residual/gdd_",i,".tif",sep = ""),
       width = 30,height = 20,units = "cm",res=300,compression = "lzw")
  plot(family_gdd_sr)
  dev.off()
  
  family_swb_sr <- SR_mapping(swb_bootstrapping_sum_ord,i)
  writeRaster(family_swb_sr,paste("01-data/cleaned data/raster/quantile_residual/swb_",i,".tif",sep = ""))
  tiff(paste("01-data/cleaned data/jpeg/fam_richness/quantile_residual/swb_",i,".tif",sep = ""),
       width = 30,height = 20,units = "cm",res=300,compression = "lzw")
  plot(family_swb_sr)
  dev.off()
  
  family_rock_sr <- SR_mapping(rock_bootstrapping_sum_ord,i)
  writeRaster(family_rock_sr,paste("01-data/cleaned data/raster/quantile_residual/rock_",i,".tif",sep = ""))
  tiff(paste("01-data/cleaned data/jpeg/fam_richness/quantile_residual/rock_",i,".tif",sep = ""),
       width = 30,height = 20,units = "cm",res=300,compression = "lzw")
  plot(family_rock_sr)
  dev.off()
}

writeRaster(SR,paste("output/combined/richness/paleorichness_",m,"_Ma.tif",sep = ""))


###########################################################################
########### HOPE this is the last version of eurasia models #############
#########################################################################
# qgam
compound_qgam_sum <- data.frame()
hetero2_qgam_sum <- data.frame()
gdd_qgam_sum <- data.frame()
swb_qgam_sum <- data.frame()
rock_qgam_sum <- data.frame()
full_qgam_sum <- data.frame()

compound_resi_05 <- data.frame(loc=sr_eurasia_revision$LEVEL3_COD)
hetero2_resi_05 <- data.frame(loc=sr_eurasia_revision$LEVEL3_COD)
gdd_resi_05 <- data.frame(loc=sr_eurasia_revision$LEVEL3_COD)
swb_resi_05 <- data.frame(loc=sr_eurasia_revision$LEVEL3_COD)
rock_resi_05 <- data.frame(loc=sr_eurasia_revision$LEVEL3_COD)
full_resi_05 <- data.frame(loc=sr_eurasia_revision$LEVEL3_COD)

compound_resi_50 <- data.frame(loc=sr_eurasia_revision$LEVEL3_COD)
hetero2_resi_50 <- data.frame(loc=sr_eurasia_revision$LEVEL3_COD)
gdd_resi_50 <- data.frame(loc=sr_eurasia_revision$LEVEL3_COD)
swb_resi_50 <- data.frame(loc=sr_eurasia_revision$LEVEL3_COD)
rock_resi_50 <- data.frame(loc=sr_eurasia_revision$LEVEL3_COD)
full_resi_50 <- data.frame(loc=sr_eurasia_revision$LEVEL3_COD)

compound_resi_95 <- data.frame(loc=sr_eurasia_revision$LEVEL3_COD)
hetero2_resi_95 <- data.frame(loc=sr_eurasia_revision$LEVEL3_COD)
gdd_resi_95 <- data.frame(loc=sr_eurasia_revision$LEVEL3_COD)
swb_resi_95 <- data.frame(loc=sr_eurasia_revision$LEVEL3_COD)
rock_resi_95 <- data.frame(loc=sr_eurasia_revision$LEVEL3_COD)
full_resi_95 <- data.frame(loc=sr_eurasia_revision$LEVEL3_COD)

for (i in 1:length(kew_family)) {
  print(kew_family[i])
  compound_kew_qgam <- kew_qgam(sr_eurasia_revision,loc = "LEVEL3_COD",pred = "compound",resp = kew_family[i],mod="idiv")
  compound_qgam_sum <- rbind(compound_qgam_sum,c(family=kew_family[i],compound_kew_qgam[["eva_matrics"]]))
  names(compound_qgam_sum) <- c("family","qnt05","qnt050","qnt95","SR","resi_upper","resi_lower","resi_fit","coefficient","p","explained_deviance","adj_R_square","AIC","BIC")
  compound_resi_05 <- merge(compound_resi_05,compound_kew_qgam[["resi_05"]],all.x = T)
  compound_resi_50 <- merge(compound_resi_50,compound_kew_qgam[["resi_50"]],all.x = T)
  compound_resi_95 <- merge(compound_resi_95,compound_kew_qgam[["resi_95"]],all.x = T)
  
  hetero2_kew_qgam <- kew_qgam(sr_eurasia_revision,loc = "LEVEL3_COD",pred = "hetero_shannon",resp = kew_family[i],mod="idiv")
  hetero2_qgam_sum <- rbind(hetero2_qgam_sum,c(family=kew_family[i],hetero2_kew_qgam[["eva_matrics"]]))
  names(hetero2_qgam_sum) <- c("family","qnt05","qnt050","qnt95","SR","resi_upper","resi_lower","resi_fit","coefficient","p","explained_deviance","adj_R_square","AIC","BIC")
  hetero2_resi_05 <- merge(hetero2_resi_05,hetero2_kew_qgam[["resi_05"]],all.x = T)
  hetero2_resi_50 <- merge(hetero2_resi_50,hetero2_kew_qgam[["resi_50"]],all.x = T)
  hetero2_resi_95 <- merge(hetero2_resi_95,hetero2_kew_qgam[["resi_95"]],all.x = T)
  
  gdd_kew_qgam <- kew_qgam(sr_eurasia_revision,loc = "LEVEL3_COD",pred = "gdd",resp = kew_family[i],mod="idiv")
  gdd_qgam_sum <- rbind(gdd_qgam_sum,c(family=kew_family[i],gdd_kew_qgam[["eva_matrics"]]))
  names(gdd_qgam_sum) <- c("family","qnt05","qnt050","qnt95","SR","resi_upper","resi_lower","resi_fit","coefficient","p","explained_deviance","adj_R_square","AIC","BIC")
  gdd_resi_05 <- merge(gdd_resi_05,gdd_kew_qgam[["resi_05"]],all.x = T)
  gdd_resi_50 <- merge(gdd_resi_50,gdd_kew_qgam[["resi_50"]],all.x = T)
  gdd_resi_95 <- merge(gdd_resi_95,gdd_kew_qgam[["resi_95"]],all.x = T)
  
  swb_kew_qgam <- kew_qgam(sr_eurasia_revision,loc = "LEVEL3_COD",pred = "swb",resp = kew_family[i],mod="idiv")
  swb_qgam_sum <- rbind(swb_qgam_sum,c(family=kew_family[i],swb_kew_qgam[["eva_matrics"]]))
  names(swb_qgam_sum) <- c("family","qnt05","qnt050","qnt95","SR","resi_upper","resi_lower","resi_fit","coefficient","p","explained_deviance","adj_R_square","AIC","BIC")
  swb_resi_05 <- merge(swb_resi_05,swb_kew_qgam[["resi_05"]],all.x = T)
  swb_resi_50 <- merge(swb_resi_50,swb_kew_qgam[["resi_50"]],all.x = T)
  swb_resi_95 <- merge(swb_resi_95,swb_kew_qgam[["resi_95"]],all.x = T)
  
  rock_kew_qgam <- kew_qgam(sr_eurasia_revision,loc = "LEVEL3_COD",pred = "rock",resp = kew_family[i],mod="idiv")
  rock_qgam_sum <- rbind(rock_qgam_sum,c(family=kew_family[i],rock_kew_qgam[["eva_matrics"]]))
  names(rock_qgam_sum) <- c("family","qnt05","qnt050","qnt95","SR","resi_upper","resi_lower","resi_fit","coefficient","p","explained_deviance","adj_R_square","AIC","BIC")
  rock_resi_05 <- merge(rock_resi_05,rock_kew_qgam[["resi_05"]],all.x = T)
  rock_resi_50 <- merge(rock_resi_50,rock_kew_qgam[["resi_50"]],all.x = T)
  rock_resi_95 <- merge(rock_resi_95,rock_kew_qgam[["resi_95"]],all.x = T)
  
  full_kew_qgam <- kew_qgam(sr_eurasia_revision,loc = "LEVEL3_COD",pred = c("gdd","swb","rock"),resp = kew_family[i],mod="full")
  full_qgam_sum <- rbind(full_qgam_sum,c(family=kew_family[i],full_kew_qgam[["eva_matrics"]]))
  names(full_qgam_sum) <- c("family","qnt05","qnt050","qnt95","SR","resi_upper","resi_lower","resi_fit","coefficient","p","explained_deviance","adj_R_square","AIC","BIC")
  full_resi_05 <- merge(full_resi_05,full_kew_qgam[["resi_05"]],all.x = T)
  full_resi_50 <- merge(full_resi_50,full_kew_qgam[["resi_50"]],all.x = T)
  full_resi_95 <- merge(full_resi_95,full_kew_qgam[["resi_95"]],all.x = T)
}

compound_qgam_sum_num <- as.data.frame(sapply(compound_qgam_sum, as.numeric))
compound_qgam_sum_num$family <- compound_qgam_sum$family
hetero2_qgam_sum_num <- as.data.frame(sapply(hetero2_qgam_sum, as.numeric))
hetero2_qgam_sum_num$family <- hetero2_qgam_sum$family
full_qgam_sum_num <- as.data.frame(sapply(full_qgam_sum, as.numeric))
full_qgam_sum_num$family <- full_qgam_sum$family
gdd_qgam_sum_num <- as.data.frame(sapply(gdd_qgam_sum, as.numeric))
gdd_qgam_sum_num$family <- gdd_qgam_sum$family
swb_qgam_sum_num <- as.data.frame(sapply(swb_qgam_sum, as.numeric))
swb_qgam_sum_num$family <- swb_qgam_sum$family
rock_qgam_sum_num <- as.data.frame(sapply(rock_qgam_sum, as.numeric))
rock_qgam_sum_num$family <- rock_qgam_sum$family

reorganizing <- function(resi_50){
resi_50 <- as.data.frame(t(resi_50))
colnames(resi_50) <- resi_50[1,]
resi_50 <- resi_50[-1,]
resi_50_num <- as.data.frame(sapply(resi_50, as.numeric))
rownames(resi_50_num) <- rownames(resi_50)
sort(colMeans(resi_50_num,na.rm = T))

resi_50_sub <- dplyr::select(resi_50_num,
                                 names(sort(colMeans(resi_50_num,na.rm = T))[78:87]) )
resi_50_sub$family=rownames(resi_50_sub)
resi_50_sub <- pivot_longer(resi_50_sub,cols = 1:10,values_to = "residuals",names_to = "region")
return(resi_50_sub)
}

compound_resi_50_sub=reorganizing(compound_resi_50)
compound_resi_50_sub$hetero <- "compound"
hetero2_resi_50_sub=reorganizing(hetero2_resi_50)
hetero2_resi_50_sub$hetero <- "hetero2"
gdd_resi_50_sub <- reorganizing(gdd_resi_50)
gdd_resi_50_sub$hetero <- "gdd"
swb_resi_50_sub <- reorganizing(swb_resi_50)
swb_resi_50_sub$hetero <- "swb"
rock_resi_50_sub <- reorganizing(rock_resi_50)
rock_resi_50_sub$hetero <- "rock"
full_resi_50_sub <- reorganizing(full_resi_50)
full_resi_50_sub$hetero <- "full"

# try 95th 
compound_resi_50_sub=reorganizing(compound_resi_95)
compound_resi_50_sub$hetero <- "compound"
hetero2_resi_50_sub=reorganizing(hetero2_resi_95)
hetero2_resi_50_sub$hetero <- "hetero2"
gdd_resi_50_sub <- reorganizing(gdd_resi_95)
gdd_resi_50_sub$hetero <- "gdd"
swb_resi_50_sub <- reorganizing(swb_resi_95)
swb_resi_50_sub$hetero <- "swb"
rock_resi_50_sub <- reorganizing(rock_resi_95)
rock_resi_50_sub$hetero <- "rock"
full_resi_50_sub <- reorganizing(full_resi_95)
full_resi_50_sub$hetero <- "full"
########################
# Figure 2 plotting
########################
resi_50_tot <- rbind(compound_resi_50_sub,gdd_resi_50_sub,swb_resi_50_sub,rock_resi_50_sub,full_resi_50_sub)
resi_50_tot[resi_50_tot$hetero=="full","hetero"]="multivariate"
resi_50_tot[resi_50_tot$hetero=="rock","hetero"]="bedrock"
resi_50_tot$hetero <- factor(resi_50_tot$hetero, levels=c('compound','gdd','bedrock','swb','multivariate'))
resi_50_tot$region <- as.factor(resi_50_tot$region)
resi_50_tot_num <- na.omit(resi_50_tot)

resi_50_tot_num$region <- factor(resi_50_tot_num$region,levels = c("CHC","CHS","MYA","EHM","ASS","CHT","JAP","NEP","TUR","CHN","IRN","FRA","KAZ","YUG","SPA"))
tiff("03-results/figure_2/regions_bp_full.tif",width = 40,height = 20,units = "cm",res=300,compression = "lzw")
p <- ggplot(data = resi_50_tot_num,aes(x=region,y=residuals,fill=hetero))+
  geom_boxplot()+
  scale_fill_manual(values = c("#B33039","#CC9966","#CDCD66","#62AC5D","#65A8F5"))+
  labs(x="Region",y="Residuals")+
  geom_hline(yintercept =0,linetype="dashed",color="grey")+
  theme(legend.position="bottom")+
  theme_classic(base_size = 24)
# p <- p+theme(legend.position="bottom",
#         legend.justification="center",
#         legend.margin=margin(5,0,20,-70))
print(p)
dev.off()

resi_50_tot_num_sub <- dplyr::filter(resi_50_tot_num,region %in% c("CHC","CHS","FRA","KAZ","TUR"))
resi_50_tot_num_sub$region <- factor(resi_50_tot_num_sub$region,levels = c("CHC","CHS","KAZ","TUR","FRA"))
tiff("03-results/figure_2/regions_bp_sub.tif",width = 30,height = 20,units = "cm",res=600,compression = "lzw")
p <- ggplot(data = resi_50_tot_num_sub,aes(x=region,y=residuals,fill=hetero))+
  geom_boxplot(show.legend = FALSE)+
  scale_fill_manual(values = c("#B33039","#FAB168","#F2F265","#9FCBFD","#65A8F5"))+
  labs(x="Region",y="Residuals")+
  geom_hline(yintercept =0,linetype="dashed",color="grey")+
#  theme(legend.position="bottom")+
  theme_classic(base_size = 24)
# p <- p+theme(legend.position="bottom",
#              legend.justification="center",
#              legend.margin=margin(5,0,20,-70))
print(p)
dev.off()

adj_R_square <- as.data.frame(cbind(family=compound_qgam_sum_num$family,
                                    compound=compound_qgam_sum_num$adj_R_square,
                                    gdd=gdd_qgam_sum_num$adj_R_square,
                                    swb=swb_qgam_sum_num$adj_R_square,
                                    multivariate=full_qgam_sum_num$adj_R_square,
                                    bedrock=rock_qgam_sum_num$adj_R_square))
adj_R_square_long <- pivot_longer(adj_R_square,cols = 2:6,values_to = "adj_R2",names_to = "heterogeneity")
str(adj_R_square_long)
adj_R_square_long$adj_R2 <- as.numeric(adj_R_square_long$adj_R2)
adj_R_square_long$heterogeneity <- factor(adj_R_square_long$heterogeneity,levels=c('compound','gdd','bedrock','swb','multivariate'))

tiff("03-results/figure_2/hetero_adjR_bp.tif",width = 30,height = 20,units = "cm",res=300,compression = "lzw")
p <- ggplot(data = adj_R_square_long,aes(x=heterogeneity,y=adj_R2,fill=heterogeneity))+
  geom_boxplot()+
  scale_fill_manual(values = c("#B33039","#FAB168","#F2F265","#9FCBFD","#65A8F5"))+
  labs(x="Heterogeneity",y="adjusted R square")+
  theme(legend.position="bottom")+
  theme_classic(base_size = 24)
p <- p+theme(legend.position="bottom",
             legend.justification="center",
             legend.margin=margin(5,0,20,-70))
print(p)
dev.off()

# significance test for the heterogeneity model
mod <- lm(adj_R2~heterogeneity,data = adj_R_square_long)
print(summary(anova(mod)))
print(TukeyHSD(aov(mod)))

# try all the outlier combination
compound_outlier <- compound_qgam_sum_num[compound_qgam_sum_num$resi_upper>0,"family"]
gdd_outlier <- gdd_qgam_sum_num[gdd_qgam_sum_num$resi_upper>0,"family"]
swb_outlier <- swb_qgam_sum_num[swb_qgam_sum_num$resi_upper>0,"family"]
rock_outlier <- rock_qgam_sum_num[rock_qgam_sum_num$resi_upper>0,"family"]
full_outlier <- full_qgam_sum_num[full_qgam_sum_num$resi_upper>0,"family"]

outlier_all <- intersect(intersect(intersect(intersect(compound_outlier,gdd_outlier),swb_outlier),rock_outlier),full_outlier)
write.csv(outlier_all,"01-data/cleaned data/table/outlier_all.csv")
kew_final_50 <- kew_final[kew_final$spp_number>50,"family"]
outlier_intersect <- intersect(outlier_all,kew_final_50)
SR_outlier_all <- sr_mapping(outlier_all)
plot(SR_outlier_all)
writeRaster(SR_outlier_all,"01-data/cleaned data/raster/revision_outlier.tif")

file.copy(paste("01-data/cleaned data/jpeg/fam_richness/",outlier_all,".tif",sep = ""),
          paste("01-data/cleaned data/jpeg/fam_richness/outlier_all"))

non_outlier <- kew_family[!kew_family %in% outlier_all]
SR_non_outlier <- sr_mapping(non_outlier)
plot(SR_non_outlier)
writeRaster(SR_non_outlier,"01-data/cleaned data/raster/revision_non_outlier.tif")

# annova analyses for the different region
for (het in c('compound','gdd','bedrock','swb','multivariate')) {
  print(het)
  region_sub <- resi_50_tot_num_sub[resi_50_tot_num_sub$hetero==het,]
  mod <- lm(residuals~region,data = region_sub)
  print(anova(mod))
  print(TukeyHSD(aov(mod)))
}
resi_CHC_compound <- dplyr::filter(resi_50_tot_num_sub,region=="CHC"&hetero=="compound")
mean(resi_CHC_compound$residuals)
sd(resi_CHC_compound$residuals)

resi_CHC_gdd <- dplyr::filter(resi_50_tot_num_sub,region=="CHC"&hetero=="gdd")
mean(resi_CHC_gdd$residuals)
sd(resi_CHC_gdd$residuals)

resi_CHC_multivariate <- dplyr::filter(resi_50_tot_num_sub,region=="CHC"&hetero=="multivariate")
mean(resi_CHC_multivariate$residuals)
sd(resi_CHC_multivariate$residuals)

# compute mean of the outlier non outlier families adjusted r square
outlier_label <- as.data.frame(cbind(family=outlier_all,cat="outlier"))
compound_qgam_sum_num<- merge(compound_qgam_sum_num,outlier_label,by="family",all.x=T)
compound_qgam_sum_num[is.na(compound_qgam_sum_num$cat)==T,"cat"]="Non-outlier"
mean(dplyr::filter(compound_qgam_sum_num,cat=="outlier")$adj_R_square)
sd(dplyr::filter(compound_qgam_sum_num,cat=="outlier")$adj_R_square)
mean(dplyr::filter(compound_qgam_sum_num,cat=="Non-outlier")$adj_R_square)
sd(dplyr::filter(compound_qgam_sum_num,cat=="Non-outlier")$adj_R_square)
wilcox.test(dplyr::filter(compound_qgam_sum_num,cat=="Non-outlier")$adj_R_square,dplyr::filter(compound_qgam_sum_num,cat=="outlier")$adj_R_square)

outlier_label <- as.data.frame(cbind(family=outlier_all,cat="outlier"))
gdd_qgam_sum_num<- merge(gdd_qgam_sum_num,outlier_label,by="family",all.x=T)
gdd_qgam_sum_num[is.na(gdd_qgam_sum_num$cat)==T,"cat"]="Non-outlier"
mean(dplyr::filter(gdd_qgam_sum_num,cat=="outlier")$adj_R_square)
sd(dplyr::filter(gdd_qgam_sum_num,cat=="outlier")$adj_R_square)
mean(dplyr::filter(gdd_qgam_sum_num,cat=="Non-outlier")$adj_R_square)
sd(dplyr::filter(gdd_qgam_sum_num,cat=="Non-outlier")$adj_R_square)
wilcox.test(dplyr::filter(gdd_qgam_sum_num,cat=="Non-outlier")$adj_R_square,dplyr::filter(gdd_qgam_sum_num,cat=="outlier")$adj_R_square)

outlier_label <- as.data.frame(cbind(family=outlier_all,cat="outlier"))
swb_qgam_sum_num<- merge(swb_qgam_sum_num,outlier_label,by="family",all.x=T)
swb_qgam_sum_num[is.na(swb_qgam_sum_num$cat)==T,"cat"]="Non-outlier"
mean(dplyr::filter(swb_qgam_sum_num,cat=="outlier")$adj_R_square)
sd(dplyr::filter(swb_qgam_sum_num,cat=="outlier")$adj_R_square)
mean(dplyr::filter(swb_qgam_sum_num,cat=="Non-outlier")$adj_R_square)
sd(dplyr::filter(swb_qgam_sum_num,cat=="Non-outlier")$adj_R_square)
wilcox.test(dplyr::filter(swb_qgam_sum_num,cat=="Non-outlier")$adj_R_square,dplyr::filter(swb_qgam_sum_num,cat=="outlier")$adj_R_square)

outlier_label <- as.data.frame(cbind(family=outlier_all,cat="outlier"))
rock_qgam_sum_num<- merge(rock_qgam_sum_num,outlier_label,by="family",all.x=T)
rock_qgam_sum_num[is.na(rock_qgam_sum_num$cat)==T,"cat"]="Non-outlier"
mean(dplyr::filter(rock_qgam_sum_num,cat=="outlier")$adj_R_square)
sd(dplyr::filter(rock_qgam_sum_num,cat=="outlier")$adj_R_square)
mean(dplyr::filter(rock_qgam_sum_num,cat=="Non-outlier")$adj_R_square)
sd(dplyr::filter(rock_qgam_sum_num,cat=="Non-outlier")$adj_R_square)
t.test(dplyr::filter(rock_qgam_sum_num,cat=="Non-outlier")$adj_R_square,dplyr::filter(rock_qgam_sum_num,cat=="outlier")$adj_R_square)


outlier_label <- as.data.frame(cbind(family=outlier_all,cat="outlier"))
full_qgam_sum_num<- merge(full_qgam_sum_num,outlier_label,by="family",all.x=T)
full_qgam_sum_num[is.na(full_qgam_sum_num$cat)==T,"cat"]="Non-outlier"
mean(dplyr::filter(full_qgam_sum_num,cat=="outlier")$adj_R_square)
sd(dplyr::filter(full_qgam_sum_num,cat=="outlier")$adj_R_square)
mean(dplyr::filter(full_qgam_sum_num,cat=="Non-outlier")$adj_R_square)
sd(dplyr::filter(full_qgam_sum_num,cat=="Non-outlier")$adj_R_square)
wilcox.test(dplyr::filter(full_qgam_sum_num,cat=="Non-outlier")$adj_R_square,dplyr::filter(full_qgam_sum_num,cat=="outlier")$adj_R_square)

# last last try: different residuals based on the response variables
compound_qgam_sum_num$resi_response_fit <- compound_qgam_sum_num$SR-compound_qgam_sum_num$qnt050
compound_qgam_sum_num$resi_response_upper <- compound_qgam_sum_num$SR-compound_qgam_sum_num$qnt95
compound_qgam_sum_num_fit_order <- compound_qgam_sum_num[order(compound_qgam_sum_num$resi_response_fit,decreasing = T),]
compound_qgam_sum_num_upper_order <- compound_qgam_sum_num[order(compound_qgam_sum_num$resi_response_upper,decreasing = T),]
compound_q1_response_upper <- compound_qgam_sum_num_upper_order[1:24,"family"]
compound_q1_response <- compound_qgam_sum_num_fit_order[1:25,"family"]
  
# decided to keep the outlier information but match whether endemic families are fully covered
endemic_spp <- read.csv("01-data/cleaned data/table/Endemic_final_common_removed.csv")
endemic_genus <- distinct(endemic_spp,Genus,.keep_all = T)
names(endemic_genus)[2] <- "genus"
endemic_genus$genus <- str_to_title(endemic_genus$genus)
endemic_family<- merge(endemic_genus,kew_accept_spp,by="genus",all.x=T)
endemic_family_list <- unique(endemic_family$family)

SR_outlier_all <- sr_mapping(outlier_intersect)
plot(SR_outlier_all)

gdd_qgam_sum_num <- as.data.frame(sapply(gdd_qgam_sum, as.numeric))
gdd_qgam_sum_num$family <- gdd_qgam_sum$family
family_gdd_gam_outlier <- gdd_qgam_sum_num[gdd_qgam_sum_num$diff_upper>0,"family"]
mean(gdd_qgam_sum_num[gdd_qgam_sum_num$diff_upper>0,"explained_deviance"])
family_gdd_gam_nonoutlier <- gdd_qgam_sum_num[gdd_qgam_sum_num$diff_upper<0,"family"]


HDS_template <- raster("01-data/cleaned data/raster/HDS_template_1k.tif")

sr_mapping <- function(fam_list){
  SR=HDS_template
for (m in 1:length(fam_list)) {
  rst_m=readAll(raster(paste("01-data/cleaned data/raster/tif_family_update/",fam_list[m],".tif",sep = "")))
  SR= SR + rst_m
  print(m)
}
  return(SR)
}

