library(dplyr)
library(grid)
library(ggplot2)
library(ggplotify)
library(FactoMineR)
library(dendextend)
library(gridExtra)
library(SIBER)

load("data/tmp_data/ROE_res_env.Rdata")
t5_prim_dir="results/Phase_1/task5/primary"

# PCA plots
png(file.path(t5_prim_dir,"ROE_full_PCA_plt.png"),width = 8,height=10.5,units="in",res=450)
grid.arrange(grobs=c(PC12_plt_ls,PCA_varcomp_plt_ls),ncol=2,as.table = FALSE)
dev.off()

# PCA plots
png(file.path(t5_prim_dir,"ROE_PCA_w_assgns_plt.png"),width = 8,height=10.5,units="in",res=450)
grid.arrange(grobs=c(PCA_ellip_plt),ncol=2,as.table = FALSE)
dev.off()

# n_clusters
png(file.path(t5_prim_dir,"ROE_cluster_plt.png"),width = 8,height=10.5,units="in",res=450)
grid.arrange(grobs=c(Sil_plt_ls,TSS_plt_ls),ncol=2,as.table = FALSE)
dev.off()


# Dendrogram
# colors for red and black dendrograms
cols_ls=lapply(smlt_clst_ID,function(x){aa=rep(2,8);aa[x]=1;aa})
smlt_clst_ID=c(3,4,1,1) # order of smolts in the colored dendrogram


png(file.path(t5_prim_dir,"ROE_cluster_dend_plt.png"),width = 10.5,height=8,units="in",res=450)
par(mfrow=c(2,2),mar=c(4,1.5,3,1.5))
for(i in 1:4){
  # plot(OUT_ls[[i]]$fit_ls[[1]],main=smolt_studs[i])
  
  dend2=as.dendrogram(OUT_ls[[i]]$fit_ls[[1]])
  dend2 <- dend2 %>% dendextend::color_branches(dend2,k=clust_sel_v[i],
                                                col = cols_ls[[i]][1:clust_sel_v[i]])
  plot(dend2,horiz=T,leaflab = "none")
  title(paste0(smolt_studs[i]," (k = ",clust_sel_v[i],")"))
  
  cat_var_mat_ls[[1]]
  OUT$clustID
  }
dev.off()


