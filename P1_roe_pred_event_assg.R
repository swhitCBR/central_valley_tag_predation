library(dplyr)
library(grid)
library(ggplot2)
library(ggplotify)
# library(scatterplot3d)
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
                  'log_resid_time_mid')

max_var_vec <- c('cumul_resid_time',
                 'cumul_visit_count',
                 'consec_visit_count',
                 'cumul_down_up_switch',
                 'cumul_time_hrz',
                 'cumul_against_flow')

####################################################### #
# log trans
####################################################### #

sdsDF_raw$log_resid_time_near=log(as.numeric(sdsDF_raw$resid_time_near)+0.0001) # significance digits as min
sdsDF_raw$log_resid_time_mid=log(as.numeric(sdsDF_raw$resid_time_mid)+0.0001) # significance digits as min
sdsDF_raw=subset(sdsDF_raw,tagID!="sochn_1206018") # dropped because of temporal overlap
my_var_ls <- build_summ_ls(mean_vars=mean_var_vec, max_vars=max_var_vec)
table(is.na(sdsDF_raw$migration_rate))
sdsDF_raw$migration_rate=ifelse(sdsDF_raw$gen_loc=="Medford_SouthDelta_2015" & 
                                  sdsDF_raw$first_event,0,sdsDF_raw$migration_rate)
sdsDF_raw$BLPS=ifelse(sdsDF_raw$gen_loc=="Medford_SouthDelta_2015" & 
                        sdsDF_raw$first_event,0,sdsDF_raw$BLPS)
####################################################### #

options(dplyr.summarise.inform = FALSE)

# plot longing for artifacts of study design and release/receiver locations 
# site specific artifact

# creating empty list objects
hist_plt_ls=det_plt_ls=tag_MAT_DF_ls=raw_var_mat_ls=TSS_plt_ls=Sil_plt_ls=PC12_plt_ls=cat_var_mat_ls=PCA_var_comp_ls=list()
PCA_obj=SUM_tab=PCA=PCA_alt=sil_MAT_l_ls=agg_vars=TSS_nclust_ls=outlier_tab_ls=OUT_ls=PCA_varcomp_plt_ls=list()
SUM_tab=PCA_alt=PCA_obj=tag_MAT_DF_ls=raw_var_mat_ls=cat_var_mat_ls=eval_per_ls=PC12_plt_eve_ls=tms=list()




####################################################### #
# Reading in actual tag assignments
####################################################### #

ROE_tag_assgDF <- readRDS("data/tmp_data/tmpROE_tag_assgDF.rds")

####################################################### #


# ii=1
PC12_plt_eve_ls=tmp_ls=list()
mu_ii_ls=sig_ii_ls=mu_ls=Sigma_ls=Z_ls=ellips_smlt_tags_ls=list()

bt=proc.time()
ii=1;jj=1
for(ii in 1:length(smolt_studs)){
  
  ellips_smlt_tags_ls[[ii]]=ROE_tag_assgDF[ROE_tag_assgDF$sdsID==smolt_studs[ii] & ROE_tag_assgDF$smolt,"tagID"]
  
  mu_ls[[ii]]=list()
  Sigma_ls[[ii]]=list()
  Z_ls[[ii]]=list()
  PC12_plt_eve_ls[[ii]]=list()
  
  # defining the number of resummarizations that must occur
  eval_per=2:max(sdsDF_raw[sdsDF_raw$sdsID==smolt_studs[ii],"det_event_ID"])
  
  for(jj in 1:length(eval_per)){
      
  print(paste(ii,jj,"of",length(smolt_studs),length(eval_per)))
  
  tmp=sdsDF_raw[sdsDF_raw$sdsID==smolt_studs[ii] & sdsDF_raw$det_event_ID<=eval_per[jj],]

  tag_MAT_DF <- build_tag_summ_MAT(
    var_grps = my_var_ls,
    input_DF = tmp,
    tagID_var = "tagID",
    sub0s="time_since_last_at_station_MEAN")

  tag_MAT_DF_ls[[ii]]=tag_MAT_DF
  summ_vers=suppressMessages(tag_dat_summ(input_DF = tmp,
                         var_grps = my_var_ls,
                         tag_mat_in=tag_MAT_DF,
                         plot = F,verbose = F))

  raw_var_mat_ls[[ii]]=summ_vers$raw_var_mat
  cat_var_mat_ls[[ii]]=summ_vers$var_mat_wCAT

  scl_raw_vars=scale(summ_vers$raw_var_mat)
  out_rows=which(scl_raw_vars>5,arr.ind = T)[,1]

  out_cols=which(scl_raw_vars>5,arr.ind = T)[,2]
  scl_vals=sapply(1:length(out_rows),function(x) scl_raw_vars[out_rows[x],out_cols[x]] )
  raw_vals=sapply(1:length(out_rows),function(x) summ_vers$raw_var_mat[out_rows[x],out_cols[x]] )


  PCA_obj[[ii]]=prcomp(raw_var_mat_ls[[ii]],scale. = F)


  SUM_tab[[ii]]=summary(PCA_obj[[ii]])
  PCA_alt[[ii]]=PCA(raw_var_mat_ls[[ii]],graph = F,scale.unit = F)

  PCA_12_DF=data.frame(sdsID_ind=ii,
                       det_ind=jj,
                       tagID=summ_vers$var_mat_wCAT$tagID,
                       sdsID=smolt_studs[ii],
                       PC1=PCA_alt[[ii]]$ind$coord[,1],
                       PC2=PCA_alt[[ii]]$ind$coord[,2])

  PC12_plt_eve_ls[[ii]][[jj]]=PCA_12_DF
  
  
  subsub <- PC12_plt_eve_ls[[ii]][[jj]][PC12_plt_eve_ls[[ii]][[jj]]$tagID  %in%  ellips_smlt_tags_ls[[ii]], ]
  
  if(jj==1){
  print(c(nrow(PCA_12_DF),nrow(subsub)))}
  
  # tmp=ROE_event_assngDF[ROE_event_assngDF$det_ind==jj,]
  mu_ls[[ii]][[jj]]=colMeans(subsub[,c("PC1","PC2")])
  Sigma_ls[[ii]][[jj]] <- cov(subsub[,c("PC1","PC2")])
  Z_ls[[ii]][[jj]]=SIBER::pointsToEllipsoid(X = subsub[,c("PC1","PC2")],
                           Sigma_ls[[ii]][[jj]],
                           mu_ls[[ii]][[jj]])
  }

  mu_ii_ls[[ii]]=do.call(rbind,mu_ls[[ii]])
  sig_ii_ls[[ii]]=do.call(rbind,Sigma_ls[[ii]])
  
  tmp_ls[[ii]] <- do.call(rbind,PC12_plt_eve_ls[[ii]])
  tms[[ii]]=proc.time()-bt
  eval_per_ls[[ii]]=eval_per
}

proc.time()-bt

TAB <- cbind(do.call(rbind,mu_ii_ls),
             do.call(rbind,sig_ii_ls))
ROE_event_assngDF <- do.call(rbind,tmp_ls)

saveRDS(ROE_event_assngDF,"data/tmp_data/tmpROE_event_assgDF.rds")

