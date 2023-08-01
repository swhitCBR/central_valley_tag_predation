library(dplyr)
library(grid)
library(ggplot2)
library(ggplotify)
library(FactoMineR)

source("scripts/predator_filters/ROE_method/phase_1/ROE_clust_fxns.R")
source("scripts/predator_filters/ROE_method/phase_1/clust_analysis_function.R")

load("data/tmp_data/P1_workspace_v1.3.1_.Rdata")


smolt_studs=unique(sdsDF_raw$sdsID)[grep(unique(sdsDF_raw$sdsID),pattern="smolts")]

mean_var_vec <- c('cumul_time_hrz',
                  'migration_rate',
                  'BLPS',
                  'time_since_last_at_station',
                  "log_resid_time_near",
                  'log_resid_time_mid'
                  )


max_var_vec <- c('cumul_resid_time',
                 'cumul_visit_count',
                 'consec_visit_count',
                 'cumul_down_up_switch',
                 'cumul_time_hrz',
                 'cumul_against_flow'
                 )

####################################### #
# Transformations and data cleaning
####################################### #

### Log transfomration
sdsDF_raw$log_resid_time_near=log(as.numeric(sdsDF_raw$resid_time_near)+0.0001) # significance digits as min
sdsDF_raw$log_resid_time_mid=log(as.numeric(sdsDF_raw$resid_time_mid)+0.0001) # significance digits as min

### removed tag
sdsDF_raw=subset(sdsDF_raw,tagID!="sochn_1206018") # dropped because of temporal overlap
my_var_ls <- build_summ_ls(mean_vars=mean_var_vec, max_vars=max_var_vec)
table(is.na(sdsDF_raw$migration_rate))
sdsDF_raw$migration_rate=ifelse(sdsDF_raw$gen_loc=="Medford_SouthDelta_2015" & 
                                  sdsDF_raw$first_event,0,sdsDF_raw$migration_rate)
sdsDF_raw$BLPS=ifelse(sdsDF_raw$gen_loc=="Medford_SouthDelta_2015" & 
                        sdsDF_raw$first_event,0,sdsDF_raw$BLPS)

####################################### #

# creating empty list objects
hist_plt_ls=det_plt_ls=tag_MAT_DF_ls=raw_var_mat_ls=TSS_plt_ls=Sil_plt_ls=PC12_plt_ls=cat_var_mat_ls=PCA_var_comp_ls=list()
PCA_obj=SUM_tab=PCA=PCA_alt=sil_MAT_l_ls=agg_vars=TSS_nclust_ls=outlier_tab_ls=OUT_ls=PCA_varcomp_plt_ls=OUT_ls=list()

# ii=3
for(ii in 1:length(smolt_studs)){
  # ii=3
  tmp=sdsDF_raw[sdsDF_raw$sdsID==smolt_studs[ii],]


  my_var_ls <- build_summ_ls(mean_vars=mean_var_vec, max_vars=max_var_vec)
  
  table(is.na(tmp$migration_rate))
  tmp$migration_rate=ifelse(tmp$gen_loc=="Medford_SouthDelta_2015" & #is.na(tmp$migration_rate) &
                              tmp$first_event,0,tmp$migration_rate)
  tmp$BLPS=ifelse(tmp$gen_loc=="Medford_SouthDelta_2015" & #is.na(tmp$BLPS) &
                    tmp$first_event,0,tmp$BLPS)
  
  
  tag_MAT_DF <- build_tag_summ_MAT(
                     var_grps = my_var_ls,
                     input_DF = tmp,
                     tagID_var = "tagID",
                     sub0s="time_since_last_at_station_MEAN"
                     )
  tag_MAT_DF
  
  

  tag_MAT_DF_ls[[ii]]=tag_MAT_DF
  summ_vers=tag_dat_summ(input_DF = tmp,
                         var_grps = my_var_ls,
                         tag_mat_in=tag_MAT_DF,
                         plot = F)
  
  raw_var_mat_ls[[ii]]=summ_vers$raw_var_mat
  cat_var_mat_ls[[ii]]=summ_vers$var_mat_wCAT
  
  hist_plt_ls[[ii]]=as.grob(summ_vers$plts$raw_hist_plt + ggtitle(smolt_studs[ii]))
  det_plt_ls[[ii]]=as.grob(summ_vers$plts$det_plt + ggtitle(smolt_studs[ii]))
  
  
  # scaled hist plots
  scl_raw_vars=scale(summ_vers$raw_var_mat)
  out_rows=which(scl_raw_vars>5,arr.ind = T)[,1]
  # View(scl_raw_vars[which(scl_raw_vars>5,arr.ind = T)[,1],])
  summ_vers$var_mat_wCAT[out_rows,c("tagID")]
  # out_cols=which(scl_raw_vars[out_rows,]>5,arr.ind = T)[,2]
  out_cols=which(scl_raw_vars>5,arr.ind = T)[,2]
  scl_vals=sapply(1:length(out_rows),function(x) scl_raw_vars[out_rows[x],out_cols[x]] )
  raw_vals=sapply(1:length(out_rows),function(x) summ_vers$raw_var_mat[out_rows[x],out_cols[x]] )

  outlier_tab_ls[[ii]]=data.frame(summ_vers$var_mat_wCAT[out_rows,c("tagID","n_dets")],colnames(scl_raw_vars)[out_cols],raw_vals,scl_vals,out_cols)
  
  # cat_var_mat_ls[[ii]][unique(which(is.na(raw_var_mat_ls[[ii]]),arr.ind = T)[,1]),]
  
  
  OUT=clust_anlys_SW(DF_in = raw_var_mat_ls[[ii]],
                     method_nm="hierWARD2",
                     n_clust=2:20,
                     distance_mat="euclidean")
  
  # cat_var_mat_ls[[ii]]
  
  sil_MAT_l_ls[[ii]]=data.frame(sds_ID=smolt_studs[ii],OUT$sil_MAT_l)
  
  PCA_obj[[ii]]=prcomp(raw_var_mat_ls[[ii]],scale. = F)
  # print(summary(PCA_obj))
  
  PCA_var_comp_ls[[ii]]=data.frame(sds_ID=smolt_studs[ii],
                 var_exp=summary(prcomp(raw_var_mat_ls[[ii]],scale. = F))$"importance"["Cumulative Proportion",1:10],
                 x=1:10)
  
  PCA_varcomp_plt_ls[[ii]] <- ggplot(data=PCA_var_comp_ls[[ii]],aes(x=x,y=var_exp*100)) +
    ggtitle(smolt_studs[ii]) +
    scale_y_continuous(limits=c(40,100))+
    geom_point() + geom_line() + mytheme +
    scale_x_continuous(limits=c(1,10),breaks=seq(1,10),expand=c(0.02,0)) +
    ylab("Variance explained (%)") +# ggtitle("a.") +
    xlab("Principal components") +
    theme(panel.spacing.x = unit(2, "lines"),
          axis.title.x = element_blank(),axis.title.y = element_blank(),
          # axis.title.y=element_text(size = 11),
          plot.title = element_text(size=12))
  
  
  
  
  SUM_tab[[ii]]=summary(PCA_obj[[ii]])
  PCA_alt[[ii]]=PCA(raw_var_mat_ls[[ii]],graph = F,scale.unit = F)
  
  PCA_12_DF=data.frame(sdsID=smolt_studs[ii],
                       PC1=PCA_alt[[ii]]$ind$coord[,1],
                       PC2=PCA_alt[[ii]]$ind$coord[,2])
  
  PC12_plt_ls[[ii]] <-ggplot(data=PCA_12_DF,aes(x=PC1,y=PC2)) +
    geom_point(shape=21,fill="white") + 
    facet_wrap(~sdsID,scales="free",ncol=1) +
    geom_hline(yintercept = 0,linetype="dashed") + 
    geom_vline(xintercept = 0,linetype="dashed") + 
    mytheme + theme(#axis.title.x = element_text(size=11),
      # axis.title.y = element_text(size=11),
      axis.title.x = element_blank(),axis.title.y = element_blank())
  
  clust_lab=as.character(OUT$n_clust)
  TSS=c()
  for(jj in 1:length(OUT$n_clust)){
    agg_vars[[jj]]=aggregate(raw_var_mat_ls[[ii]],by=list(
      OUT$clustID[,clust_lab[jj]]),
      function(x) sum(scale(x,scale=FALSE)^2))
    TSS[jj]=sum(rowSums(agg_vars[[jj]],na.rm = T))
  }
  
  # View(tmp[tmp$tagID %in% c("sochn_1205927", "sochn_1205928", "sochn_1205965", "sochn_1205975"),])
  
  
  TSS_nclustDF=data.frame(sds_ID=smolt_studs[ii],cluster=clust_lab,method="hierWARD2",TSS)
  TSS_nclustDF$prop_var=(1-TSS_nclustDF$TSS/max(TSS_nclustDF$TSS))*100
  
  TSS_nclust_ls[[ii]]=TSS_nclustDF
  
  
  TSS_plt_ls[[ii]]=ggplot(subset(TSS_nclust_ls[[ii]]),aes(x=as.numeric(cluster),y=prop_var)) +ggtitle(smolt_studs[ii]) +
    geom_point() + geom_line() + mytheme +
    scale_y_continuous(limits=c(0,100))+
    scale_x_continuous(limits=c(2,18),breaks=seq(2,18,2),expand=c(0.02,0)) +
    ylab("Variance explained (%)") +# ggtitle("a.") +
    xlab("Number of clusters") +
    theme(panel.spacing.x = unit(2, "lines"),
          axis.title.x = element_blank(),axis.title.y = element_blank(),
          # axis.title.y=element_text(size = 11),
          plot.title = element_text(size=12))
  
    
  # for(ii in 1:4){
  Sil_plt_ls[[ii]]=ggplot(sil_MAT_l_ls[[ii]],aes(x=Clusters,y=mean)) + ggtitle(smolt_studs[ii]) +
               geom_errorbar(aes(ymin=lcl,ymax=ucl),color="gray",linewidth=0.5,width=0) +
               xlab("Number of clusters") +  geom_line(linewidth=0.75) + 
                geom_point(shape=16,color="black")+scale_x_continuous(breaks=seq(2,18,2)) +
               mytheme  + theme(legend.position = "none") + ylab("Silhoutte") +
               theme(panel.spacing.x = unit(2, "lines"),
                     axis.title.x = element_blank(),axis.title.y = element_blank(),
                     # axis.title.y=element_text(size = 11),
                     plot.title = element_text(size=12)#,
                     # strip.text.x = element_text(size=11),
                     # strip.text.y = element_text(size=11)
                     ) 
  
  OUT_ls[[ii]]=OUT
  
  
}

# save.image("data/tmp_data/ROE_res_env.Rdata")


CHOICES

clust_sel_v=c(5,8,3,8) # number of clusters choosen to divide data

# seems like there is a missmatch happening
ROE_tag_assg_ls=list()
# par(mfrow=c(2,2))
for(ii in 1:length(smolt_studs)){
  ROE_tag_assg_ls[[ii]]=data.frame(
    sdsID=smolt_studs[ii],
    cat_var_mat_ls[[ii]],
    clusterID=OUT_ls[[ii]]$"clustID"[,(clust_sel_v[ii]-1)], # -1 here because 2 is the minimum clust number
    PC1=PCA_alt[[ii]]$ind$coord[,1],
    PC2=PCA_alt[[ii]]$ind$coord[,2])
    # Cluster ids are ordered by the number so 1 is the largest group (i.e., "smolts")
    ROE_tag_assg_ls[[ii]]$smolt=ROE_tag_assg_ls[[ii]]$clusterID==1
    ROE_tag_assg_ls[[ii]]$predator=!ROE_tag_assg_ls[[ii]]$smolt
}

ROE_tag_assgDF=do.call(rbind,ROE_tag_assg_ls)

saveRDS(ROE_tag_assgDF,"data/tmp_data/tmpROE_tag_assgDF.rds")

# table(ROE_tag_assgDF$predator,ROE_tag_assgDF$sdsID)
# N_smolts_2015 N_smolts_2016 S_smolts_2015 S_smolts_2016
# SMOLT               556           547            42            75
# PREDATOR             11            19           756          1071
# stat_ellipse(data=subset(DF,filter_res=="smolt"),type = "norm") 

########################################### #

ii=1

PCA_ellip_plt=list()
for(ii in 1:length(smolt_studs)){
  PCA_12_smlt_only=ROE_tag_assgDF[ROE_tag_assgDF$sdsID==smolt_studs[[ii]] & ROE_tag_assgDF$smolt==TRUE,]
  PCA_12_pred_only=ROE_tag_assgDF[ROE_tag_assgDF$sdsID==smolt_studs[[ii]] & ROE_tag_assgDF$smolt==FALSE,]
  mu_s=colMeans(PCA_12_smlt_only[,c("PC1","PC2")])
  Sigma <- cov(PCA_12_smlt_only[,c("PC1","PC2")])
  Zsmolt <- SIBER::pointsToEllipsoid(PCA_12_smlt_only[,c("PC1","PC2")],
                                Sigma = Sigma,mu = mu_s)
  Zpred=SIBER::pointsToEllipsoid(PCA_12_pred_only[,c("PC1","PC2")],
                           Sigma = Sigma,mu = mu_s)

  PCA_ellip_plt[[ii]] <- ggplot(data=PCA_12_smlt_only,aes(x=PC1,y=PC2)) +
  stat_ellipse(type = "norm",level=0.95) + 
  geom_point(shape=21,color="gray20",fill="black") + ggtitle(smolt_studs[ii]) +
  geom_point(data=PCA_12_pred_only,aes(x=PC1,y=PC2),shape=21,fill="red") +
  mytheme+ theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
                 plot.title = element_text(size=12))
}
