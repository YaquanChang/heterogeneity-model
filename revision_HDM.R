################################################################################
# this is the revise Hengduan scale analyses, the model would be the linear 
# quantile regression model, for every single family
################################################################################
#install.packages("tidyverse")
#install.packages("MuMIn")
require(raster)
require(Qtools)
require(tidyr)
require(tidyverse)
require(MuMIn)
require(rgdal)

if(interactive()){
  setwd("C:/Users/yaqchang.D/OneDrive - ETH Zurich/PhD-chapter 1")
  dir <- "01-data/cleaned data/"

} else {
  #########       CLUSTER   #########
  dir <- "input/"
  ###################################
}
source("02-codes/function/Conservation et services2_Hotspots_Congruence.R")
source("02-codes/function/Map_the_world.R")
load("01-data/cleaned data/RData/revision_HDS.RData")

HDS_template <- raster(paste(dir,"raster/HDS_template_1k.tif",sep = ""))
kew_family <- read.csv(paste(dir,"table/kew_family.csv",sep=""))
kew_family <- kew_family$x
HDS_hillshade <- raster("01-data/cleaned data/raster/HDS_hillshade.tif")
HDS_river <- readOGR("01-data/cleaned data/shape/HDS_river.shp")
HDS_river <- spTransform(HDS_river,CRSobj = crs(HDS_hillshade))
HDS_template_df <- na.omit(data.frame(coordinates(HDS_template),HDS_template=getValues(HDS_template)))
scale_vec <- seq(from = 5, to = 200, by = 20)
mod_sum <- data.frame()

# mapping non-outlier patterns
kew_outlier <- read.csv("01-data/cleaned data/table/outlier_all.csv")
kew_outlier <- as.data.frame(cbind(family=kew_outlier$x,cat="outlier"))
kew_fam_tot <- as.data.frame(cbind(family=kew_family))
kew_fam_tot <- merge(kew_outlier,kew_fam_tot,by="family",all.y=T)
kew_fam_tot[is.na(kew_fam_tot$cat)==T,2]="non_outlier"
kew_non_outlier <- kew_fam_tot[kew_fam_tot$cat=="non_outlier","family"]
sr_non_outlier <- HDS_template
HDM_fam_tot <- list.files("01-data/cleaned data/raster/tif_family_update/")
HDM_fam_tot <- gsub(".tif","",HDM_fam_tot)
HDM_fam_tot <- HDM_fam_tot[HDM_fam_tot!="Thumbs.db"]
HDM_fam_tot <- as.data.frame(cbind(family=HDM_fam_tot))
HDM_fam_tot <- merge(HDM_fam_tot,kew_fam_tot,by="family",all.x=T)
HDM_fam_tot[is.na(HDM_fam_tot$cat)==T,"cat"]="rest"
HDM_rest <- HDM_fam_tot[HDM_fam_tot$cat=="rest","family"]
for (i in 1:length(kew_non_outlier)) {
  print(i)
  sr_fam <- raster(paste("01-data/cleaned data/raster/tif_family_update/",kew_non_outlier[i],".tif",sep=""))
  sr_non_outlier <- sr_non_outlier+sr_fam
}
writeRaster(sr_non_outlier,"01-data/cleaned data/raster/revision_non_outlier.tif")

sr_HDM_rest <- HDS_template
for (i in 1:length(HDM_rest)) {
  print(i)
  sr_fam <- raster(paste("01-data/cleaned data/raster/tif_family_update/",HDM_rest[i],".tif",sep=""))
  sr_HDM_rest <- sr_HDM_rest+sr_fam
}
writeRaster(sr_HDM_rest,"01-data/cleaned data/raster/revision_HDM_rest.tif")
#
# HDM_mod <- function(family){
#   scale_sum <- data.frame()
#   for (scale in scale_vec) {
#     print(paste("scale is",scale))
#     # prepare the data
#     compound_rst <- raster(paste(dir,"raster/hetero/HDS_hetero_shan_",scale,".tif",sep = ""))
#     print("work!")
#     gdd_rst <- raster(paste(dir,"raster/hetero/HDS_gdd0_shan_",scale,".tif",sep = ""))
#     print("work!")
#     swb_rst <- raster(paste(dir,"raster/hetero/HDS_swb_shan_",scale,".tif",sep = ""))
#     print("work!")
#     rock_rst <- raster(paste(dir,"raster/hetero/HDS_rock_shan_",scale,".tif",sep = ""))
#     print("work!")
#     resp <- raster(paste(dir,"raster/tif_family_update/",family,".tif",sep = ""))
#     print("work!")
#     resp_rst <- focal(resp,matrix(1,scale,scale),fun=mean,pad=T)
#     print("work!")
#     writeRaster(resp_rst,paste(dir,"raster/family_scale/",family,"_",scale,".tif",sep = ""),overwrite=T)
#     mod_stk <- stack(compound_rst,gdd_rst,swb_rst,rock_rst,resp_rst)
#     pseudo_r_matric <- data.frame(matrix(NA,100,5))
#     colnames(pseudo_r_matric) <- c("compound","gdd","swb","rock","full")
#     scale_m <- scale*1000
#     for (j in 1:100) {
#       print(j)
#       center <- HDS_template_df[sample(nrow(HDS_template_df),size = 1),1:2]
#       xmin <- center[1]-scale_m*ceiling(1372000/scale_m)
#       xmax <- center[1]+scale_m*ceiling(1372000/scale_m)
#       ymin <- center[2]-scale_m*ceiling(1240000/scale_m)
#       ymax <- center[2]+scale_m*ceiling(1240000/scale_m)
#       ex <- c(xmin$x,xmax$x,ymin$y,ymax$y)
#       samples <- raster(xmn=xmin$x,xmx=xmax$x,ymn=ymin$y,ymx=ymax$y,res=scale_m)
#       values(samples) <- NA
#       HDS_samples <- crop(samples,HDS_template)
#       HDS_samples_df <- data.frame(coordinates(HDS_samples))
#       HDS_samples_dat <- na.omit(data.frame(HDS_samples_df,extract(mod_stk,HDS_samples_df)))
#       names(HDS_samples_dat) <- c("x","y","compound","gdd","swb","rock","SR")
#       
#       # build quantile regression model at 50th quantile and extract the explained deviance
#       
#       # pseudo
#       
#       compound_50 <- rq.counts(SR~compound,data = HDS_samples_dat,tau=0.50, M = 50, zeta = 1e-05)
#       gdd_50 <- rq.counts(SR~gdd,data = HDS_samples_dat,tau=0.50, M = 50, zeta = 1e-05)
#       swb_50 <- rq.counts(SR~swb,data = HDS_samples_dat,tau=0.50, M = 50, zeta = 1e-05)
#       rock_50 <- rq.counts(SR~rock,data = HDS_samples_dat,tau=0.50, M = 50, zeta = 1e-05)
#       full_50 <- rq.counts(SR~gdd+swb+rock,data = HDS_samples_dat,tau=0.50, M = 50, zeta = 1e-05)
#       
#       null_50 <- midrq(SR~1,data = HDS_samples_dat,tau=0.50)
#       null_resi_50 <- sum(residuals.rq.counts(null_50)^2)/nrow(HDS_samples_dat)
#       RSS_compound_50 <- sum(residuals.rq.counts(compound_50)^2)/nrow(HDS_samples_dat)
#       pseudo_r_compound_50 <- 1-(RSS_compound_50/null_resi_50)
#       
#       RSS_gdd_50 <- sum(residuals.rq.counts(gdd_50)^2)/nrow(HDS_samples_dat)
#       pseudo_r_gdd_50 <- 1-(RSS_gdd_50/null_resi_50)
#       
#       RSS_swb_50 <- sum(residuals.rq.counts(swb_50)^2)/nrow(HDS_samples_dat)
#       pseudo_r_swb_50 <- 1-(RSS_swb_50/null_resi_50)
#       
#       RSS_rock_50 <- sum(residuals.rq.counts(rock_50)^2)/nrow(HDS_samples_dat)
#       pseudo_r_rock_50 <- 1-(RSS_rock_50/null_resi_50)
#       
#       RSS_full_50 <- sum(residuals.rq.counts(full_50)^2)/nrow(HDS_samples_dat)
#       pseudo_r_full_50 <- 1-(RSS_full_50/null_resi_50)
#       
#       pseudo_r <- c(pseudo_r_compound_50,pseudo_r_gdd_50,pseudo_r_swb_50,pseudo_r_rock_50,pseudo_r_full_50)
#       pseudo_r_matric[j,] <- pseudo_r
#     }
#     pseudo_r_matric_mean <- colMeans(pseudo_r_matric)
#     mod_sum <- as.data.frame(cbind(scale,compound=pseudo_r_matric_mean["compound"],
#                                    gdd=pseudo_r_matric_mean["gdd"],
#                                    swb=pseudo_r_matric_mean["swb"],
#                                    rock=pseudo_r_matric_mean["rock"],
#                                    full=pseudo_r_matric_mean["full"]))
#     scale_sum <- rbind(mod_sum,scale_sum)
#   }
#   
#   return(mod_sum)
# }
sr_mapping <- function(raster_layer,path){
  colourCount = length(unique(getValues(raster_layer)))
  tiff(path,width=12,height=12,units="cm",compression="lzw",res=600)
  plot(HDS_hillshade,
       col = grey(0:100/100),
       legend = FALSE,
       bty="n",xaxt = "n", yaxt = "n",box=F,main=scale_vec[i],adj=0,cex.main = 2)
  plot(raster_layer,col=getPalette(colourCount),alpha=0.5,add=T,axis.args=list(cex.axis=1.2,line=0))
  plot(HDS_river,col="#9191AC",add=T)
  dev.off()
}

HDM_mod_gam <- function(family){
  scale_sum <- data.frame()
  rsquare_raw <- list()
  for (i in 1:length(scale_vec)) {
    print(paste("scale is",scale_vec[i]))
    
    # prepare the data
    compound_rst <- raster(paste(dir,"raster/hetero/HDS_hetero_shan_",scale_vec[i],".tif",sep = ""))
    compound2_rst <- raster(paste(dir,"raster/hetero/output/HDS_hetero2_",scale_vec[i],".tif",sep = ""))
    print("work!")
    gdd_rst <- raster(paste(dir,"raster/hetero/HDS_gdd0_shan_",scale_vec[i],".tif",sep = ""))
    print("work!")
    swb_rst <- raster(paste(dir,"raster/hetero/HDS_swb_shan_",scale_vec[i],".tif",sep = ""))
    print("work!")
    rock_rst <- raster(paste(dir,"raster/hetero/HDS_rock_shan_",scale_vec[i],".tif",sep = ""))
    print("work!")
    
    # # this is the heterogeneity for each family
    #resp <- raster(paste(dir,"raster/tif_family_update/",family,".tif",sep = ""))
    #  resp <- raster(paste(dir,"raster/",family,".tif",sep = ""))
    #  
    #  print("work!")
    #  resp_rst <- focal(resp,matrix(1,scale_vec[i],scale_vec[i]),fun=mean,na.rm=T,pad=T)
    #  resp_rst <- mask(resp_rst,HDS_template)
    #  print("work!")
    #  writeRaster(resp_rst,paste(dir,"raster/revision_scale/",family,"_",scale_vec[i],".tif",sep = ""),overwrite=T)
    # 
    # this is the heterogeneity for outlier richness
    resp_rst <- raster(paste(dir,"raster/revision_scale/",family,"_",scale_vec[i],".tif",sep = ""))

    mod_stk <- stack(compound_rst,compound2_rst,gdd_rst,swb_rst,rock_rst,resp_rst)
    r_square_matric <- data.frame(matrix(NA,100,6))
    colnames(r_square_matric) <- c("compound","compound2","gdd","swb","rock","full")
    exp_dev_matric <- data.frame(matrix(NA,100,6))
    colnames(exp_dev_matric) <- c("compound","compound2","gdd","swb","rock","full")
    coef_matric <- data.frame(matrix(NA,100,6))
    colnames(coef_matric) <- c("compound","compound2","gdd","swb","rock","full")
    p_matric <- data.frame(matrix(NA,100,6))
    colnames(p_matric) <- c("compound","compound2","gdd","swb","rock","full")
    aic_matric <- data.frame(matrix(NA,100,6))
    colnames(aic_matric) <- c("compound","compound2","gdd","swb","rock","full")
    sample_matric <- matrix(NA,1,100)
    scale_m <- 185*1000 # 185 is the scale number can get around 25 samples
    compound_md <- list()
    for (j in 1:200) {
      print(j)
      center <- HDS_template_df[sample(nrow(HDS_template_df),size = 1),1:2]
      xmin <- center[1]-scale_m*ceiling(1372000/scale_m)
      xmax <- center[1]+scale_m*ceiling(1372000/scale_m)
      ymin <- center[2]-scale_m*ceiling(1240000/scale_m)
      ymax <- center[2]+scale_m*ceiling(1240000/scale_m)
      ex <- c(xmin$x,xmax$x,ymin$y,ymax$y)
      samples <- raster(xmn=xmin$x,xmx=xmax$x,ymn=ymin$y,ymx=ymax$y,res=scale_m)
      values(samples) <- NA
      HDS_samples <- crop(samples,HDS_template)
      HDS_samples_df <- data.frame(coordinates(HDS_samples))
      HDS_samples_dat <- na.omit(data.frame(HDS_samples_df,raster::extract(mod_stk,HDS_samples_df)))
      names(HDS_samples_dat) <- c("x","y","compound","compound2","gdd","swb","rock","SR")
      #HDS_samples_dat <- HDS_samples_dat[HDS_samples_dat$SR>0,]
      
      # build quantile regression model at 50th quantile and extract the explained deviance
      
      # pseudo
      
      compound_50 <- mgcv::gam(SR~compound,data = HDS_samples_dat)
      compound2_50 <- mgcv::gam(SR~compound2,data = HDS_samples_dat)
      gdd_50 <- mgcv::gam(SR~gdd,data = HDS_samples_dat)
      swb_50 <- mgcv::gam(SR~swb,data = HDS_samples_dat)
      rock_50 <- mgcv::gam(SR~rock,data = HDS_samples_dat)
      full_50 <- mgcv::gam(SR~gdd+swb+rock,data = HDS_samples_dat)
      
      compound_md[[j]] <- compound_50
      
      r_square <- c(summary(compound_50)$r.sq,summary(compound2_50)$r.sq,summary(gdd_50)$r.sq,summary(swb_50)$r.sq,
                    summary(rock_50)$r.sq,summary(full_50)$r.sq)
      r_square_matric[j,] <- r_square
      
      exp_dev <- c(summary(compound_50)$dev.expl,summary(compound2_50)$dev.expl,summary(gdd_50)$dev.expl,summary(swb_50)$dev.expl,
                    summary(rock_50)$dev.expl,summary(full_50)$dev.expl)
      exp_dev_matric[j,] <- exp_dev
      
      coef <- c(compound_50$coefficients[2],compound2_50$coefficients[2],gdd_50$coefficients[2],swb_50$coefficients[2],
                rock_50$coefficients[2],full_50$coefficients[2])
      coef_matric[j,] <- coef
      
      p <- c(summary(compound_50)$p.pv[2],summary(compound2_50)$p.pv[2],summary(gdd_50)$p.pv[2],summary(swb_50)$p.pv[2],
                summary(rock_50)$p.pv[2],summary(full_50)$p.pv[2])
      p_matric[j,] <- p
      
      AIC <- c(AIC(compound_50),AIC(compound2_50),AIC(gdd_50),AIC(swb_50),
               AIC(rock_50),AIC(full_50))
      aic_matric[j,] <- AIC
      
      sample_matric[j] <-nrow(HDS_samples_dat)
      print(paste("sample size is",nrow(HDS_samples_dat)))
    }
    r_square_matric_mean <- colMeans(r_square_matric)
    exp_dev_matric_mean <- colMeans(exp_dev_matric)
    coef_matric_mean <- colMeans(coef_matric)
    p_matric_mean <- colMeans(p_matric)
    aic_matric_mean <- colMeans(aic_matric)
    sample_matric_mean <- mean(sample_matric,na.rm = T)
    mod_sum <- as.data.frame(rbind(r_square_matric_mean,
                                   exp_dev_matric_mean,
                                   coef_matric_mean,
                                   p_matric_mean,
                                   aic_matric_mean))
    mod_sum$scale=scale_vec[i]
    mod_sum$sample=sample_matric_mean
    scale_sum <- rbind(mod_sum,scale_sum)
    rsquare_raw[[i]] <- r_square_matric
  }
  
  return(list(scale_sum=scale_sum,rsquare_raw=rsquare_raw))
}
se <- function(a){sqrt(sum((a-mean(a))^2/(length(a)-1)))/sqrt(length(a))}

scale_fun <- function(non_outlier_fam,label){
  non_outlier_fam_matrix <- non_outlier_fam$scale_sum
  non_outlier_fam_matrix$X <- rownames(non_outlier_fam_matrix)
  non_outlier_fam_matrix$X <-  gsub("[[:digit:]]","",non_outlier_fam_matrix$X)
  matric_sum_non_outlier <- non_outlier_fam_matrix
  
  non_outlier_rsquare_raw <- non_outlier_fam$rsquare_raw
  non_outlier_sd_tot <- data.frame()
  for (i in 1:15) {
    sd_sub <- sapply(non_outlier_rsquare_raw[[i]], se)
    non_outlier_sd_tot <- rbind(non_outlier_sd_tot,sd_sub)
  }
  names(non_outlier_sd_tot) <- c("compound","compound2","gdd","swb","rock","full")
  non_outlier_sd_tot$scale <- seq(from = 285, to = 5, by = -20) # change to 285
  non_outlier_sd_tot <- pivot_longer(non_outlier_sd_tot,cols = c("compound","compound2","gdd","swb","rock","full"),
                                     names_to = "heterogeneity",
                                     values_to = "sd")
  
  # compute some matric
  # file.list <- list.files("01-data/cleaned data/table/HDS_try/",pattern = ".csv",full.names = T)
  # family_name <- gsub("\\.csv$","",basename(file.list))
  # matric_sum <- data.frame()
  # for (i in 1:length(file.list)) {
  #   print(i)
  #   file <- read.csv(file.list[i])
  #   file$family=family_name[i]
  #   matric_sum <- rbind(matric_sum,file)
  # }
  # non outlier
  non_outlier_matric_long <- pivot_longer(matric_sum_non_outlier,cols = c("compound","gdd","swb","rock","full"),
                                          names_to = "heterogeneity",
                                          values_to = "matric")
  
  non_outlier_matric_long$X <- gsub("[[:digit:]]","",non_outlier_matric_long$X)
  non_outlier_matric_vec <- unique(non_outlier_matric_long$X)
  
  matric <- non_outlier_matric_long
  sd <- non_outlier_sd_tot
  
  
  #for (i in 1:length(matric_vec,matric,sd)) {
  dat_sub <- matric[matric$X=="r_square_matric_mean",]
  dat_sub <- merge(dat_sub,sd,by=c("heterogeneity","scale"))
  dat_sub[dat_sub$heterogeneity=="full","heterogeneity"]="multivariate"
  dat_sub[dat_sub$heterogeneity=="rock","heterogeneity"]="bedrock"
  dat_sub$heterogeneity <- factor(dat_sub$heterogeneity,levels=c('compound','gdd','bedrock','swb','multivariate'))
  #tiff(paste("03-results/supp/family_matrics/outlier_",matric_vec[i],".tif",sep = ""),width = 30,height = 20,units = "cm",res=300,compression = "lzw")
  tiff(paste("03-results/figure_3/revision/",label,"_r_square_matric_mean",".tif",sep = ""),width = 30,height = 20,units = "cm",res=300,compression = "lzw")
  p <- ggplot(data = dat_sub,aes(x=scale,y=matric,color=heterogeneity))+
    geom_line(size=1)+
    scale_color_manual(values = c("#B33039","#CC9966","#CDCD66","#9FCBFD","#0066CC"))+
    geom_ribbon(aes(scale, ymin=matric-sd,ymax=matric+sd,fill=heterogeneity),alpha=0.3,outline.type = "full")+
    scale_fill_manual(values = c("#B33039","#CC9966","#CDCD66","#9FCBFD","#0066CC"))+
    labs(y="Adjusted R square",x="Window sizes")+
    theme_classic(base_size = 24)+
    theme(legend.position="bottom",
          legend.justification="center",
          legend.margin=margin(0,0,0,-70))
  print(p)
  dev.off()
}





for (i in 1:97) {
  print(kew_family[i])
  test <- HDM_mod_gam(kew_family[i])
  write.csv(test,paste(dir,"table/HDS_try/",kew_family[i],".csv",sep = ""))
}

outlier_fam  <- HDM_mod_gam("revision_outlier")
outlier_fam_matrix <- outlier_fam$scale_sum
outlier_fam_matrix$X <- rownames(outlier_fam_matrix)
outlier_fam_matrix$X <-  gsub("[[:digit:]]","",outlier_fam_matrix$X)
matric_sum_outlier <- outlier_fam_matrix
write.csv(filter(outlier_fam_matrix,X=="r_square_matric_mean"),
          "01-data/cleaned data/table/outlier_r_square.csv")
write.csv(filter(outlier_fam_matrix,X=="exp_dev_matric_mean"),
          "01-data/cleaned data/table/exp_dev_matric_mean.csv")
write.csv(filter(outlier_fam_matrix,X=="coef_matric_mean"),
          "01-data/cleaned data/table/coef_matric_mean.csv")
write.csv(filter(outlier_fam_matrix,X=="p_matric_mean"),
          "01-data/cleaned data/table/p_matric_mean.csv")
write.csv(filter(outlier_fam_matrix,X=="r_square_matric_mean"),
          "01-data/cleaned data/table/outlier_r_square.csv")
write.csv(outlier_fam_matrix,
          "01-data/cleaned data/table/outlier_fam_matrix.csv")

outlier_rsquare_raw <- outlier_fam$rsquare_raw
se <- function(a){sqrt(sum((a-mean(a))^2/(length(a)-1)))/sqrt(length(a))}
sd_tot <- data.frame()
for (i in 1:15) {
  sd_sub <- sapply(outlier_rsquare_raw[[i]], se)
  sd_tot <- rbind(sd_tot,sd_sub)
}
names(sd_tot) <- c("compound","compound2","gdd","swb","rock","full")
sd_tot$scale <- seq(from = 285, to = 5, by = -20)
sd_tot <- pivot_longer(sd_tot,cols = c("compound","compound2","gdd","swb","rock","full"),
                       names_to = "heterogeneity",
                       values_to = "sd")

non_outlier_fam <- HDM_mod_gam("revision_non_outlier")
non_outlier_fam$X <- rownames(non_outlier_fam)
non_outlier_fam$X <-  gsub("[[:digit:]]","",non_outlier_fam$X)
matric_sum_nonoutlier <- non_outlier_fam

# compute some matric
file.list <- list.files("01-data/cleaned data/table/HDS_try/",pattern = ".csv",full.names = T)
family_name <- gsub("\\.csv$","",basename(file.list))
matric_sum <- data.frame()
for (i in 1:length(file.list)) {
  print(i)
  file <- read.csv(file.list[i])
  file$family=family_name[i]
  matric_sum <- rbind(matric_sum,file)
}

matric_long <- pivot_longer(matric_sum_outlier,cols = c("compound","gdd","swb","rock","full"),
                            names_to = "heterogeneity",
                            values_to = "matric")

matric_long$X <- gsub("[[:digit:]]","",matric_long$X)
matric_vec <- unique(matric_long$X)



for (i in 1:length(matric_vec)) {
  dat_sub <- matric_long[matric_long$X==matric_vec[i],]
  dat_sub <- merge(dat_sub,sd_tot,by=c("heterogeneity","scale"))
  dat_sub[dat_sub$heterogeneity=="full","heterogeneity"]="multivariate"
  dat_sub[dat_sub$heterogeneity=="rock","heterogeneity"]="bedrock"
  dat_sub$heterogeneity <- factor(dat_sub$heterogeneity,levels=c('compound','gdd','bedrock','swb','multivariate'))
  #tiff(paste("03-results/supp/family_matrics/outlier_",matric_vec[i],".tif",sep = ""),width = 30,height = 20,units = "cm",res=300,compression = "lzw")
  tiff(paste("03-results/figure_3/revision/outlier_",matric_vec[i],".tif",sep = ""),width = 30,height = 20,units = "cm",res=300,compression = "lzw")
  p <- ggplot(data = dat_sub,aes(x=scale,y=matric,color=heterogeneity))+
    geom_line(size=1)+
    scale_color_manual(values = c("#B33039","#CC9966","#CDCD66","#62AC5D","#65A8F5"))+
    geom_ribbon(aes(scale, ymin=matric-sd,ymax=matric+sd,fill=heterogeneity),alpha=0.3,outline.type = "full")+
    scale_fill_manual(values = c("#B33039","#CC9966","#CDCD66","#62AC5D","#65A8F5"))+
    labs(y="Adjusted R square",x="Window sizes")+
    theme_classic(base_size = 24)+
    theme(legend.position="bottom",
            legend.justification="center",
            legend.margin=margin(0,0,0,-70))
  print(p)
  dev.off()
}

tiff(paste("03-results/supp/family_matrics/samples.tif",sep = ""),width = 30,height = 20,units = "cm",res=300,compression = "lzw")
p <- ggplot(data = dat_sub,aes(x=as.factor(scale),y=sample))+
  geom_boxplot()+
  labs(title="samples")
print(p)
dev.off()

# # congruence analyses 
# HDS_con <- function (hetero,family){
#   out_all <- data_frame()
#   for (i in 1:length(scale_vec)) {
#   print(c(hetero,family))
#   spp_raster <- raster(paste(dir,"raster/revision_scale/",family,"_",scale_vec[i],".tif",sep = ""))
#   pred_rst <- raster(paste(dir,"raster/hetero/",hetero,scale_vec[i],".tif",sep = ""))
#   Cong_data <- data.frame(coordinates(HDS_template),pred=getValues(pred_rst),SR=getValues(spp_raster))
#   Cong_data <-  na.omit(Cong_data)
#   plot(rasterFromXYZ(Cong_data[,c(1,2,3:4)]))
#   Data_end <- Cong_data[,-c(1:2)]
#   res_hot_spot <- Cons_Serv_Permut(Data_end,999,seuil=10/100)
#   print(res_hot_spot$Resultats)
#   Hotspot <- cbind(res_hot_spot$Mat.Hotspot,X=Cong_data$x,Y=Cong_data$y)
#   
#   
#   tiff(filename=paste("03-results/supp/congruence/",hetero,"_",family,"_",scale_vec[i],".tif",sep = ""),width=12,height=12,units="cm",compression="lzw",res=600)
#   
#   carto_function(y=Hotspot,x="pred",index2="SR",names_fig=paste(family,"_",scale[i],sep = ""),cex_point=0.1,cex_text_fig=0.3,cex_legend=0.3,
#                  text=c("Non-hotspot values",paste("Values in top 10% for",family,scale[i]),"Values in top 10% for shannon heterogeneity index","Congruence zone"),
#                  xlim=c(-248500,1123500),ylim=c(2563500,3803500),pos_legY=88,pos_legX=0.5,cex.axis=0.4,ncol=1,pos_leg="bottomleft",pt.cexleg=1)
#   dev.off()
#   out <- t(as.data.frame(unlist(c(family,hetero,scale,res_hot_spot$Resultats))))
#   names(out) <- c("pattern","heterogeneity","scale","pred","SR","Nb.hotspots.communs.Obs",
#                   "Espere=NiNj/Nt","Theo>=Obs","p_value","sign")
#   out_all <- rbind(out_all,out)
#   }
#   out_all <- as.data.frame(out_all)
#   return(out_all)
# }
# hetero <- c("HDS_hetero_shan_","HDS_gdd0_shan_","HDS_swb_shan_","HDS_rock_shan_")

statistic_all <- list()
for (i in hetero) {
  statistic <- HDS_con(i,"revision_outlier")
  statistic_all <- c(statistic_all,statistic)
}
statistic_sum[,1] <- scale_vec
statistic_sum$compound <- statistic_all[[6]]
statistic_sum$gdd <- statistic_all[[16]]
statistic_sum$swb <- statistic_all[[26]]
statistic_sum$rock <- statistic_all[[36]]
View(statistic_sum)
statistic_sum <- as.data.frame(sapply(statistic_sum, as.numeric))
statistic_sum <- pivot_longer(statistic_sum,cols = 2:5,values_to = "congruence_pixels",names_to = "heterogeneity")
names(statistic_sum)[1:2] <- c("scale","heterogeneity")
p <- ggplot(statistic_sum,aes(x=scale,y=congruence_pixels,color=heterogeneity))+
  geom_line()
p

# predict spatially and get the residuals
compound <- na.omit(data.frame(cbind(coordinates(compound_rst),compound=getValues(compound_rst),sr=getValues(resp_rst))))
compound$pred_sr <- predict(compound_50,newdata=compound, type = "response")
compound$residuals <- compound$sr-compound$pred_sr 
resi_165 <- rasterFromXYZ(compound[,c(1,2,6)])

par(mfrow=c(1,2))
plot(resp_rst,main="richness at 165km")
plot(resi_165,main="residuals at 165 km")

getPalette = colorRampPalette(rev(brewer.pal(10, "Spectral")))
# a model for the mapping
HDM_gam_mapping <- function(family){
  for (i in 1:length(scale_vec)) {
    print(paste("scale is",scale_vec[i]))
    
    # prepare the data
    compound_rst <- raster(paste(dir,"raster/hetero/HDS_hetero_shan_",scale_vec[i],".tif",sep = ""))
    compound2_rst <- raster(paste(dir,"raster/hetero/output/HDS_hetero2_",scale_vec[i],".tif",sep = ""))
    print("work!")
    gdd_rst <- raster(paste(dir,"raster/hetero/HDS_gdd0_shan_",scale_vec[i],".tif",sep = ""))
    print("work!")
    swb_rst <- raster(paste(dir,"raster/hetero/HDS_swb_shan_",scale_vec[i],".tif",sep = ""))
    print("work!")
    rock_rst <- raster(paste(dir,"raster/hetero/HDS_rock_shan_",scale_vec[i],".tif",sep = ""))
    print("work!")
    
    # # # this is the heterogeneity for each family
    #  resp <- raster(paste(dir,"raster/tif_family_update/",family,".tif",sep = ""))
    #  print("work!")
    #  resp_rst <- focal(resp,matrix(1,scale_vec[i],scale_vec[i]),fun=mean,na.rm=T,pad=T)
    #  resp_rst <- mask(resp_rst,HDS_template)
    #  print("work!")
    #  writeRaster(resp_rst,paste(dir,"raster/family_scale/",family,"_",scale_vec[i],".tif",sep = ""),overwrite=T)
    
    # this is the heterogeneity for outlier richness
    resp_rst <- raster(paste(dir,"raster/revision_scale/",family,"_",scale_vec[i],".tif",sep = ""))
    
    # plot pattern
    colourCount = length(unique(getValues(resp_rst)))
    tiff(paste("03-results/supp/scale_dependency/HDS_outlier_revision_",scale_vec[i],".tif",sep = ""),width=12,height=12,units="cm",compression="lzw",res=600)
    plot(HDS_hillshade,
         col = grey(0:100/100),
         legend = FALSE,
         bty="n",xaxt = "n", yaxt = "n",box=F,main=scale_vec[i],adj=0,cex.main = 2)
    plot(resp_rst,col=getPalette(colourCount),alpha=0.5,add=T,axis.args=list(cex.axis=1.2,line=0))
    plot(HDS_river,col="#9191AC",add=T)
    dev.off()
    
    mod_stk <- stack(compound_rst,compound2_rst,gdd_rst,swb_rst,rock_rst,resp_rst)
    mod_stk_df <- na.omit(as.data.frame(cbind(coordinates(mod_stk),getValues(mod_stk))))
    names(mod_stk_df) <- c("x","y","compound","compound2","gdd","swb","rock","SR")
    scale_m <- 185*1000 # 185 is the scale number can get around 25 samples
    
      center <- HDS_template_df[sample(nrow(HDS_template_df),size = 1),1:2]
      xmin <- center[1]-scale_m*ceiling(1372000/scale_m)
      xmax <- center[1]+scale_m*ceiling(1372000/scale_m)
      ymin <- center[2]-scale_m*ceiling(1240000/scale_m)
      ymax <- center[2]+scale_m*ceiling(1240000/scale_m)
      ex <- c(xmin$x,xmax$x,ymin$y,ymax$y)
      samples <- raster(xmn=xmin$x,xmx=xmax$x,ymn=ymin$y,ymx=ymax$y,res=scale_m)
      values(samples) <- NA
      HDS_samples <- crop(samples,HDS_template)
      HDS_samples_df <- data.frame(coordinates(HDS_samples))
      HDS_samples_dat <- na.omit(data.frame(HDS_samples_df,raster::extract(mod_stk,HDS_samples_df)))
      names(HDS_samples_dat) <- c("x","y","compound","compound2","gdd","swb","rock","SR")
      
      compound_50 <- mgcv::gam(SR~compound,data = HDS_samples_dat)
      compound2_50 <- mgcv::gam(SR~compound2,data = HDS_samples_dat)
      gdd_50 <- mgcv::gam(SR~gdd,data = HDS_samples_dat)
      swb_50 <- mgcv::gam(SR~swb,data = HDS_samples_dat)
      rock_50 <- mgcv::gam(SR~rock,data = HDS_samples_dat)
      full_50 <- mgcv::gam(SR~gdd+swb+rock,data = HDS_samples_dat)
      
      mod_stk_df$pred_compound <- predict(compound_50,newdata=mod_stk_df, type = "response")
      mod_stk_df$pred_compound2 <- predict(compound2_50,newdata=mod_stk_df, type = "response")
      mod_stk_df$pred_gdd <- predict(gdd_50,newdata=mod_stk_df, type = "response")
      mod_stk_df$pred_swb <- predict(swb_50,newdata=mod_stk_df, type = "response")
      mod_stk_df$pred_rock <- predict(rock_50,newdata=mod_stk_df, type = "response")
      mod_stk_df$pred_full <- predict(full_50,newdata=mod_stk_df, type = "response")
      
      mod_stk_df$resi_compound <- mod_stk_df$SR-mod_stk_df$pred_compound
      mod_stk_df$resi_compound2 <- mod_stk_df$SR-mod_stk_df$pred_compound2
      mod_stk_df$resi_gdd <- mod_stk_df$SR-mod_stk_df$pred_gdd
      mod_stk_df$resi_swb <- mod_stk_df$SR-mod_stk_df$pred_swb
      mod_stk_df$resi_rock <- mod_stk_df$SR-mod_stk_df$pred_rock
      mod_stk_df$resi_full <- mod_stk_df$SR-mod_stk_df$pred_full
      
      resi_compound <- rasterFromXYZ(mod_stk_df[,c("x","y","resi_compound")])
      resi_compound2 <- rasterFromXYZ(mod_stk_df[,c("x","y","resi_compound2")])
      resi_gdd <- rasterFromXYZ(mod_stk_df[,c("x","y","resi_gdd")])
      resi_swb <- rasterFromXYZ(mod_stk_df[,c("x","y","resi_swb")])
      resi_rock <- rasterFromXYZ(mod_stk_df[,c("x","y","resi_rock")])
      resi_full <- rasterFromXYZ(mod_stk_df[,c("x","y","resi_full")])
      
      # plot resi compound
      sr_mapping(resi_compound,path = paste("03-results/figure_4/revision/resi_compound_",scale_vec[i],".tif",sep = ""))
      sr_mapping(resi_compound2,path = paste("03-results/figure_4/revision/resi_compound2_",scale_vec[i],".tif",sep = ""))
      sr_mapping(resi_gdd,path = paste("03-results/figure_4/revision/resi_gdd_",scale_vec[i],".tif",sep = ""))
      sr_mapping(resi_swb,path = paste("03-results/figure_4/revision/resi_swb_",scale_vec[i],".tif",sep = ""))
      sr_mapping(resi_rock,path = paste("03-results/figure_4/revision/resi_rock_",scale_vec[i],".tif",sep = ""))
      sr_mapping(resi_full,path = paste("03-results/figure_4/revision/resi_full_",scale_vec[i],".tif",sep = ""))
      
      
  }
}

HDM_gam_mapping("revision_outlier")
sr_mapping(raster("01-data/cleaned data/raster/revision_outlier.tif"),path = 
             "03-results/eurasia/HDS_outlier.tif")

rest_fam <- HDM_mod_gam("revision_HDM_rest") 
scale_fun(non_outlier_fam = rest_fam,label = "rest")
sr_mapping(raster_layer = raster("01-data/cleaned data/raster/revision_scale/revision_HDM_rest_205.tif"),
           path = "03-results/figure_3/revision/rest_205.tif")

sr_mapping(raster_layer = raster("01-data/cleaned data/raster/revision_HDM_rest.tif"),
           path = "03-results/figure_3/revision/rest_original.tif")

sr_mapping(raster_layer = raster("01-data/cleaned data/raster/revision_non_outlier.tif"),
           path = "03-results/figure_3/revision/non_original.tif")

sr_mapping(raster_layer = raster("01-data/cleaned data/raster/elevation_outlier.tif"),
           path = "03-results/figure_3/revision/ele_original.tif")


outlier_fam <- HDM_mod_gam("revision_outlier")
scale_fun(non_outlier_fam = outlier_fam,label = "outlier_185")

scale_fun(non_outlier_fam = outlier_fam,label = "outlier_285")





