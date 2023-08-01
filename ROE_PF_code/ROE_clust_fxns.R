

build_summ_ls <- function(mean_vars,max_vars){
  require(dplyr)
  
  tmp_labs=c(
    rep("MEAN",length(mean_vars)),
    rep("MAX",length(max_vars))
    # rep("SUM",length(sum_vars))
    )
  
  list("MEAN"=mean_vars,
       "MAX"=max_vars)#,
       # "SUM"=sum_vars)
}


build_tag_summ_MAT <- function(var_grps,input_DF,tagID_var,sub0s){
  summ_types=c("MEAN","MAX")#,"SUM")#,"LEN","MAXABS","SD","PROP_VIS")
  FN_ls=list(function(x){mean(x,na.rm=T)},
             function(x){ifelse(is.infinite(max(x,na.rm=T)),0,max(x,na.rm=T))}
             # function(x){max(x,na.rm=T)}
             # function(x){ifelse(is.infinite(max(x,na.rm=T)),0,max(x,na.rm=T))}
             # function(x){ifelse(is.infinite(max(x,na.rm=T)),0,max(x,na.rm=T))},
             # function(x){sum(x,na.rm=T)}#,
             # function(x){length(x)},
             # function(x){ifelse(is.infinite(max(abs(x),na.rm=T)),0,max(abs(x),na.rm=T))},
             # function(x){ifelse(is.na(sd(x,na.rm=T)),0,sd(x,na.rm=T))},
             # function(x){sum(x)/length(x)
               # }#, # assuming that there's a T/F summary
    )
  
    names(FN_ls)=summ_types
    col_ind=list()
    MAT=NULL
    
  for(j in 1:length(var_grps)){
    # print(j)
    ncols=length(var_grps[[j]])
    # list of column indices calculated by matching lists of strings for each summary stat grouping
    # col_ind[[j]]=which(names(input_DF) %in% var_grps[[j]])
    tmp=match(names(input_DF),var_grps[[j]])
    tmp_key=data.frame(ls_order=tmp[which(!is.na(tmp))],input_col=which(!is.na(tmp)))

    col_ind[[j]]=tmp_key$input_col
    for(i in 1:ncols){
      MAT=cbind(MAT,tapply(input_DF[,col_ind[[j]][i]],input_DF[,tagID_var],FN_ls[[j]]))
      # labeling each new variable column as it is added
      colnames(MAT)[ncol(MAT)]=paste(names(input_DF)[col_ind[[j]][i]],summ_types[j],sep="_")
    }}
    
    MAT=data.frame(tagID=rownames(MAT),MAT)
    
    if(!missing(sub0s)){
    for(ii in 1:length(sub0s)){
      MAT[,sub0s[ii]]=ifelse(is.na(MAT[,sub0s[ii]]),0,MAT[,sub0s[ii]])
    }}

    tag_MAT_DF= invisible(data.frame(
      input_DF %>% group_by(tagID) %>%
      summarize(
        rel_type=unique(rel_type),
        gen_locs=length(unique(gen_loc)),
        n_visits=length(unique(visit_event)),
        n_dets=length(det_event_ID),
        det_window=lubridate::interval(min(first_detection),max(last_detection)),
        det_window_days=lubridate::time_length(lubridate::as.duration(det_window),"days")) %>%
      select(-det_window) %>%
      left_join(MAT)))
    
    tag_MAT_DF
}


tag_dat_summ <- function(input_DF,var_grps,tag_mat_in,plot=FALSE,verbose=TRUE){
  
  require(ggplot2)
  require(gridExtra)
  require(lubridate)
  require(dplyr)
 
  # look for correlations
  # corr_mat=cor(tag_mat_in[,-c(1:6),])
  # sum_col_nms=which(upper.tri(corr_mat),arr.ind = T)
  # corr_DF=data.frame(X=rownames(corr_mat)[sum_col_nms[,2]],
  #                    Y=rownames(corr_mat)[sum_col_nms[,1]],
  #                    corr=corr_mat[upper.tri(corr_mat)])

  tag_mat_in_l=reshape2::melt(tag_mat_in,id=1:6,variable="tag_sum_met")
    
  var_cols=c("tagID",unique(unlist(my_var_ls)))
  var_TAB_l=reshape2::melt(input_DF[,var_cols],id=1,variable="metric")
  var_sum_TAB = invisible( var_TAB_l %>% group_by(metric,tagID) %>%
    summarize(
      n_tot=length(value),
      n_NAs=length(which(is.na(value))),
      n_Zeros=length(which(!is.na(value) & (value==0)))) %>%
    mutate(n_nonnum=n_NAs+n_Zeros) %>%
    mutate(bad=n_nonnum==n_tot))
  
  var_sum_sparse_TAB <- invisible(var_sum_TAB %>%
    filter(bad==TRUE) %>%
    group_by(metric) %>%
    summarize(n_tags=nrow(tag_mat_in),
              tag_w_NA_or_Zero=length(bad),
              percentage=tag_w_NA_or_Zero/nrow(tag_mat_in)))
  # return(var_sum_sparse_TAB)
  
  if(verbose){
    
  var_sum_sparse_TAB
  det_plt=ggplot(tag_mat_in,aes(x=det_window_days,y=n_dets,fill=tagID,color=tagID)) + 
    geom_point(color="gray10",shape=21) +
    theme(legend.position="none") #+ mytheme
  # log_det_plt=ggplot(tag_mat_in,aes(x=det_window_days,y=log(n_dets))) + geom_point()
  # return(var_sum_TAB)
  
  par(mfrow=c(1,1),mar=c(5,4,7,2))
  
  raw_hist_plt=ggplot(tag_mat_in_l,aes(x=value,group=tag_sum_met)) +
    geom_histogram() + facet_wrap(~tag_sum_met,scales = "free",ncol=4)
  
  # raw_hist_plt=ggplot(tag_mat_in_l,aes(x=log(value+0.001),group=tag_sum_met)) +
  # geom_histogram() + facet_wrap(~tag_sum_met,scales = "free",ncol=4)

  # tags with no observations for a metric
  var_NAs_zeros_plt=ggplot(data=var_sum_sparse_TAB,aes(x=percentage,y=metric)) + geom_bar(stat="identity") + scale_x_continuous(limit=c(0,1))+ labs(x="Proportion of tags with all NAs/0s")
  
  plts=list(
    "raw_hist_plt"=raw_hist_plt,
    "var_NAs_zeros_plt"=var_NAs_zeros_plt,
    "det_plt"=det_plt)#,
    # "log_det_plt"=log_det_plt)
  
  if(plot){
  grid.arrange(
    det_plt,
    var_NAs_zeros_plt,
    # log_det_plt,
    raw_hist_plt,
    layout_matrix = rbind(c(1,3,3),
                          c(2,3,3))
    # ncol=2
    )
    
    }

  out=list(
    "var_mat_wCAT"=tag_mat_in,
    "raw_var_mat"=tag_mat_in[,-c(1:6)],
    "var_sum_TAB"=var_sum_TAB,
    # "corr_DF"=corr_DF,
    "tag_mat_in_l"=tag_mat_in_l,
    "var_sum_sparse_TAB"=var_sum_sparse_TAB,
    "plts"=plts)
 
  return(out) 
  }
  
  out=list(
    "var_mat_wCAT"=tag_mat_in,
    "raw_var_mat"=tag_mat_in[,-c(1:6)])
  
  return(out) 
  
}
  


# message("variables with high correlations")
# 
# 
# corrplot(corr_mat,
#          method="number",
#          # tl.cex = 0.5,number.cex = 0.5,
#          type="upper",
#          mar = c(0,0,3,0),
#          title="Full correlation matrix")
# 

# # look for correlations
# sum_col_nms=which(upper.tri(corr_mat),arr.ind = T)
# corr_DF=data.frame(X=rownames(corr_mat)[sum_col_nms[,2]],
#                    Y=rownames(corr_mat)[sum_col_nms[,1]],
#                    corr=corr_mat[upper.tri(corr_mat)])
# print(corr_DF[abs(corr_DF$corr)>0.8,])
# # 
# 


clust_anlys_SW <- function(DF_in,method_nm="hierWARD2",n_clust,distance_mat="euclidean"){
  require(StatMatch) # for gower.distance
  require(cluster)  # for diana()
  require(utils) # for progress bar
  require(mclust)
  require(reshape2)
  
  par(mfrow=c(1,2))
  bt=proc.time()
  dat_nms=names(DF_in)
  stopifnot(!is.null(dat_nms))
  # Applies all methods
  if(any(method_nm=="all")){method_nm=c("diana","hierWARD2","hierAVG","hierWARD1","hierCOMP","k-means","hierSING")} #,"hierWARD","MCclust"
  fit_ls=d_ls=groups_ls=cls_mat_ls=DF=list()
  n_rws=nrow(DF_in)
  
  # OLD CODE FOR CREATING AN ARRAY
  n_rws=nrow(DF_in)
    # Matrix of silhouette values for all configurations of clusters and methods
  # sil_mat=gap_mat=matrix(nrow=length(DF_in),ncol=length(n_clust))
  # Empty Array of cluster assignments to datasets
  clustID=array(dim=c(n_rws,length(n_clust)),
                dimnames=list(1:n_rws,n_clust))
  # print(str(clustID))
    # Array of silhouette values
  sil_arr=clustID
  
  
  DF_tag=DF_vis=list();
  purity_measDF=NULL  
  
  # pb1 <- txtProgressBar(min = 1, max = length(DF_in)*length(method_nm)*length(n_clust), style = 3)
  # for(i in 1:length(DF_in)){
    # Creates a unique PDF for each file
    # pdf(width = 8,height = 10.5,file = paste(names(DF_in)[i],"_clusters.pdf",sep=""))
    # if(distance_mat=="euclidean"){  d_ls[[i]] <- dist(DF_in[[i]]$for_anly, method = "euclidean")} 
    if(distance_mat=="euclidean"){  d_ls <- dist(DF_in, method = "euclidean") }
    
    # return(d_ls)
    # 
    for(k in 1:length(n_clust)){
      for(j in 1:length(method_nm)){
        grp_asg=c() # temporary vector of cluster assignments

        # Initial clustering method that I used that seemed to intuitively perform Ward clustering
        # However, according to Murtaugh and Legendre (2014), this obvious choice  does not use
        # a squared distance matrix and therefore does not actually do the Ward (1963) clustering,
        # and isn't in the same space as a PCA.
        # if(method_nm[j]=="hierWARD"){
        #   fit_ls[[k]] <- hclust(d_ls, method="ward.D")
        #   groups_ls[[k]] <- cutree(fit_ls[[k]], k=n_clust[k])
        #   grp_asg <- cutree(fit_ls[[k]], k=n_clust[k])
        #   clustID[i,k,j,]=grp_asg
        #   sil_arr[i,k,j,]=silhouette(dist = d_ls,x = cutree(fit_ls[[k]], k=n_clust[k]))[,3]
        # }

        # if(method_nm[j]=="hierWARD2"){
          fit_ls[[k]] <- hclust(d_ls, method="ward.D2")
          # EQUIVALENT TO: fit_ls[[k]] <- hclust(d_ls^2, method="ward.D")
          groups_ls[[k]] <- cutree(fit_ls[[k]], k=n_clust[k])
          grp_asg <- cutree(fit_ls[[k]], k=n_clust[k])
          clustID[,k]=grp_asg
          sil_arr[,k]=silhouette(dist = d_ls,x = cutree(fit_ls[[k]], k=n_clust[k]))[,3]
          # if(k==1){plot(fit_ls[[k]],main=dat_nms[i], cex=0.0001,xlab="",ylim=c(0,100))}

      # }
      }}
  # }
  
  # return(sil_arr)

  # mean silhouettte columns are methods and rows are clusters clusters
  sil_MAT=apply(sil_arr,2,mean)
  sil_MAT_q1=apply(sil_arr,2,function(x) quantile(x,c(0.025)))
  sil_MAT_q2=apply(sil_arr,2,function(x) quantile(x,c(0.5)))
  sil_MAT_q3=apply(sil_arr,2,function(x) quantile(x,c(0.975)))
  
  
  sil_MAT_l=data.frame(n_clust=n_clust,
                       mean=sil_MAT,
                       median=sil_MAT_q2,
                       lcl=sil_MAT_q1,
                       ucl=sil_MAT_q3)
  
  
  names(sil_MAT_l)=c("Clusters","mean","median","lcl","ucl")


  out=list(sil_arr,sil_MAT_l,method_nm,dat_nms,n_clust,clustID,fit_ls)
  names(out)=c("sil_arr","sil_MAT_l","method_nm","dat_nms","n_clust","clustID","fit_ls")

  print(proc.time()-bt)
  return(out)
  
  }
    