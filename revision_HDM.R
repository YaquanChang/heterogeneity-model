######################################################################################################
# this is the Hengduan scale analysis for paper:  10.1111/nph.19206, the model would be the nonlinear 
# quantile regression model, for outlier richness pattern. The raster layer of the heterogeneity predictor
# at different scales can be obtained from 10.16904/envidat.424.
######################################################################################################
#install.packages("tidyverse")
#install.packages("MuMIn")
require(raster)
require(Qtools)
require(tidyr)
require(tidyverse)
require(MuMIn)
require(rgdal)

setwd("C:/Users/yaqchang.D/OneDrive - ETH Zurich/PhD-chapter 1")
load("01-data/cleaned data/RData/revision_HDS.RData")

scale_vec <- seq(from = 5, to = 200, by = 20)

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

