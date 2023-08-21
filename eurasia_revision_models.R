################################################################################
#################### Eurasia scale analyses revision ###########################
################################################################################

# This code is the Eurasia scale analyses for Chang et al (10.1111/nph.19206), including a function to compute quantile regression models and actual calculations for the whole Eurasia.
require(quantreg)
require(vegan)
require(ggplot2)
require(qgam)

setwd("YOUR PATH")
sr_eurasia_revision <- read.csv("YOUR PATH/sr_eurasia_revision.csv")
kew_family <- colnames(sr_eurasia_revision[4:100])

#--------------------------------------------------------
# quantile regression models
#--------------------------------------------------------

kew_qgam <- function(dat,loc,pred,resp,mod){
  if (mod=="idiv"){
    dat_sub <- na.omit(dplyr::select(dat,c(loc,pred,resp,"area"))) # first col is location, second is predictor, and third is response variable
    dat_sub <- dat_sub[order(dat_sub[,2]),]
    names(dat_sub) <- c("loc","pred","resp","area")
    CHC <- which(dat_sub$loc=="CHC")
    
    # create a heterogeneity area model
    area_50 <- nlrq(pred~SSarrhenius(area, k,z),data = dat_sub,tau=0.50,trace=TRUE)
    # calculate heterogeneity area residual
    pred_resi <- residuals(area_50)
    pred_seq <- c(seq(min(pred_resi), max(pred_resi), length.out = 99),pred_resi[CHC])
    dat_sub$pred_resi <- pred_resi
    chc_order=100
    
    # full model
    nqr_05 <- qgam(resp~s(pred_resi,k=3),data = dat_sub,qu = 0.05)
    nqr_50 <- qgam(resp~s(pred_resi,k=3),data = dat_sub,qu = 0.50)
    nqr_95 <- qgam(resp~s(pred_resi,k=3),data = dat_sub,qu = 0.95)
    sr_pred <- as.data.frame(cbind(pred_seq,
                                   Lower=predict(nqr_05,data.frame(pred_resi=pred_seq),link="response"),
                                   Fitted=predict(nqr_50,data.frame(pred_resi=pred_seq),link="response"),
                                   Upper=predict(nqr_95,data.frame(pred_resi=pred_seq),link="response")))
    
    tiff(paste("FIGURE PATH/","pred_",kew_family[i],"_",print(pred),".tif",sep = ""),width = 30,height = 20,units = "cm",res=300,compression = "lzw")
    p=ggplot(dat_sub, aes(x = pred_resi, y = resp)) +
      geom_point(size=2) +
      geom_point(data=dat_sub[CHC,], aes(x = pred_resi, y = resp),col="red",size=2) +
      geom_text(data = dat_sub[CHC,],aes(label=loc),hjust=0, vjust=0,col="red",size=8)+
      geom_ribbon(data = sr_pred, aes(ymin = Lower, ymax = Upper, x = pred_seq),
                  fill = "steelblue2", alpha = 0.2, inherit.aes = FALSE) +
      geom_line(data = sr_pred, aes(y = Fitted, x = pred_seq),size=2) +
      #ggtitle("rq.counts function") +
      #ggtitle(paste(kew_family[i],print(pred),"richness relationship",sep = " ")) +
      # labs(y = "Species richness", x = print(pred))+ # this is for both area and heterogeniety relationship
      labs(y = "Species richness", x = pred)+
      theme_classic(base_size = 27)
    print(p)
    dev.off()
  } else if (mod=="full"){
    
    dat_sub <- na.omit(dplyr::select(dat,c(loc,pred,resp,"area"))) # first col is location, second is predictor, and third is response variable
    dat_sub <- dat_sub[order(dat_sub[,2]),]
    names(dat_sub)[c(1,5)] <- c("loc","resp")
    CHC <- which(dat_sub$loc=="CHC")
    chc_order=CHC
    
    # create a heterogeneity area model
    gdd_50 <- nlrq(gdd~SSarrhenius(area, k,z),data = dat_sub,tau=0.50,trace=TRUE)
    swb_50 <- nlrq(swb~SSarrhenius(area, k,z),data = dat_sub,tau=0.50,trace=TRUE)
    rock_50 <- nlrq(rock~SSarrhenius(area, k,z),data = dat_sub,tau=0.50,trace=TRUE)
    
    dat_sub$gdd_50 <- residuals(gdd_50) 
    dat_sub$swb_50 <- residuals(swb_50)
    dat_sub$rock_50 <- residuals(rock_50)
    
    # full model
    
    nqr_05 <- qgam(resp~s(swb_50+gdd_50+rock_50,k=3),data = dat_sub,qu=0.05)    
    nqr_50 <- qgam(resp~s(swb_50+gdd_50+rock_50,k=3),data = dat_sub,qu=0.50)
    nqr_95 <- qgam(resp~s(swb_50+gdd_50+rock_50,k=3),data = dat_sub,qu=0.95)
    sr_pred <- as.data.frame(cbind(dat_sub[,c(2,3,4)],
                                   Lower=predict(nqr_05, link="response"),
                                   Fitted=predict(nqr_50, link="response"),
                                   Upper=predict(nqr_95, link="response")))
  }
  coef <- summary(nqr_50)$edf
  p_value <- summary(nqr_50)$s.pv
  exp_dev <- summary(nqr_50)$dev.expl
  adj_R2 <- summary(nqr_50)$r.sq
  aic <- AIC(nqr_50)
  bic <- BIC(nqr_50)
  resi_50 <- data.frame(loc=dat_sub$loc,resi=residuals.gam(nqr_50))
  names(resi_50)[2] <- paste(resp,"_resi",sep="")
  resi_05 <- data.frame(loc=dat_sub$loc,resi=residuals.gam(nqr_05))
  names(resi_05)[2] <- paste(resp,"_resi",sep="")
  resi_95 <- data.frame(loc=dat_sub$loc,resi=residuals.gam(nqr_95))
  names(resi_95)[2] <- paste(resp,"_resi",sep="")
  
  eva_matrics <- c(lower=sr_pred[chc_order,"Lower"],mean=sr_pred[chc_order,"Fitted"],upper=sr_pred[chc_order,"Upper"],CHC=dat_sub[CHC,"resp"],
                   resi_upper=resi_95[CHC,2],resi_lower=resi_05[CHC,2],resi_fit=resi_50[CHC,2],
                   coef=coef,p_value=p_value,exp_dev=exp_dev,adj_R2=adj_R2,aic=aic,bic=bic)
  
  out <- list(eva_matrics=eva_matrics,resi_05=resi_05,resi_50=resi_50,resi_95=resi_95)
  return(out)
}



###########################################################################
################### non linear quantile regression model ##################
###########################################################################
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
