# Check correlation
# SCALED AND LOG TRANSFORMED DATA ONLY
check_cors=function(DATA,corr_thresh=0.8){
  cmat=cor(DATA)
  ind_gr.8=which(cor(DATA)>corr_thresh,arr.ind = T)
  ind_gr.8=ind_gr.8[apply(ind_gr.8,1,function(x) x[1]!=x[2]),]
  cor_mat=data.frame(ind_gr.8,NA,NA,corr=NA)
  for(i in 1:nrow(ind_gr.8)){
    cor_mat[i,3]=rownames(cor(DATA))[cor_mat[i,1]]
    cor_mat[i,4]=colnames(cor(DATA))[cor_mat[i,2]]
    cor_mat[i,5]=cor(DATA)[ind_gr.8[i,1],ind_gr.8[i,2]]
  }
  
  up_T=upper.tri(cmat)
  cor_mat=data.frame(which(up_T,arr.ind = T),cmat[up_T])
  cor_mat=data.frame(X1=rownames(cmat)[cor_mat$row],X2=colnames(cmat)[cor_mat$col],corr=cor_mat[,3])
  cor_mat=cor_mat[order(abs(cor_mat$corr),decreasing = T),]
  cor_mat=subset(cor_mat,abs(corr)>corr_thresh)
  
  return(cor_mat)}


# Steve's gini purity 
gini=function(TAB){
  n_grps=nrow(TAB)
  wts=rowSums(TAB)/sum(TAB) # weights for how many obs are in each cluster
  ps=TAB[,2]/rowSums(TAB)   # proportions of successes
  sum(2*ps*(1-ps)*wts)      # weighted average of suceess proportions
}

# Function for applying gap statistic to 
hclusCut <- function(x, k, d.meth = "euclidean", ...)
  list(cluster = cutree(hclust(dist(x, method=d.meth), ...), k=k))
# 
# nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
#                 cex = 0.7, col = "blue")


# FUNCTION THAT RUNS A CLUSTER ANALYSIS

# Performs analysis on alternative data sets, numbers of clusters, methods.Returns a list with many outputs including:
# $ data_sets = list of dataframes used in the analysis ("FULL" used in the paper)
# $ method_nm =  : chr [1:7] "diana"     "hierWARD2" "hierAVG"   "hierWARD1" "hierCOMP"  "k-means"   "hierSING" 
# $ dat_nms   =  : chr [1:2] "FULL" "REDUCED"
# $ n_clust   =  : int [1:19] 2 3 4 5 6 7 8 9 10 11 ...20
# "purity_measDF","fit_ls"  =  Holdover's from past approaches
# $ clustID   =  : int [1:2, 1:19, 1:7, 1:618] (BIG ARRAY WITH ALL CLUSTER ID ASSIGNEMENTS ACROSS METHODS AND DATA SETS)
# $ sil_arr   =  : int [1:2, 1:19, 1:7, 1:618] (BIG ARRAY WITH OBSERVATION-LEVEL SILLHOUTE VALUES ASSIGNEMENTS ACROSS METHODS AND DATA SETS)

clust_anlys=function(data_sets,method_nm,n_clust,distance_mat="euclidean"){
  
  require(StatMatch) # for gower.distance
  require(cluster)  # for diana()
  require(utils) # for progress bar
  require(mclust)
  par(mfrow=c(1,2))
  bt=proc.time()
  dat_nms=names(data_sets)
  stopifnot(!is.null(dat_nms))
  # Applies all methods
  if(any(method_nm=="all")){method_nm=c("diana","hierWARD2","hierAVG","hierWARD1","hierCOMP","k-means","hierSING")} #,"hierWARD","MCclust"
  fit_ls=d_ls=groups_ls=cls_mat_ls=DF=list()
  n_rws=nrow(data_sets[[1]]$for_anly)

  # Matrix of silhouette values for all configurations of clusters and methods
  # sil_mat=gap_mat=matrix(nrow=length(data_sets),ncol=length(n_clust))
  
  # Empty Array of cluster assignments to datasets
  clustID=array(dim=c(length(dat_nms),length(n_clust),length(method_nm),n_rws),
                dimnames=list(dat_nms,n_clust,method_nm,1:n_rws))
  # print(str(clustID))
  # Array of silhouette values
  sil_arr=clustID
  
  DF_tag=DF_vis=list();
  purity_measDF=NULL  

  # pb1 <- txtProgressBar(min = 1, max = length(data_sets)*length(method_nm)*length(n_clust), style = 3)
  for(i in 1:length(data_sets)){
    # Creates a unique PDF for each file
    # pdf(width = 8,height = 10.5,file = paste(names(data_sets)[i],"_clusters.pdf",sep=""))
    
    # distance measures
    if(distance_mat=="gower"){  d_ls[[i]] <- daisy(data_sets[[i]]$for_anly,metric="gower")} 
    if(distance_mat=="euclidean"){  d_ls[[i]] <- dist(data_sets[[i]]$for_anly, method = "euclidean")} 
    
    
    for(k in 1:length(n_clust)){
      for(j in 1:length(method_nm)){
        grp_asg=c() # temporary vector of cluster assignments
        
        
        # Initial clustering method that I used that seemed to intuitively perform Ward clustering
        # However, according to Murtaugh and Legendre (2014), this obvious choice  does not use
        # a squared distance matrix and therefore does not actually do the Ward (1963) clustering,
        # and isn't in the same space as a PCA.
        # if(method_nm[j]=="hierWARD"){
        #   fit_ls[[k]] <- hclust(d_ls[[i]], method="ward.D")
        #   groups_ls[[k]] <- cutree(fit_ls[[k]], k=n_clust[k])
        #   grp_asg <- cutree(fit_ls[[k]], k=n_clust[k])
        #   clustID[i,k,j,]=grp_asg
        #   sil_arr[i,k,j,]=silhouette(dist = d_ls[[i]],x = cutree(fit_ls[[k]], k=n_clust[k]))[,3]
        # } 
        
        if(method_nm[j]=="hierWARD2"){
          fit_ls[[k]] <- hclust(d_ls[[i]], method="ward.D2")
          # EQUIVALENT TO: fit_ls[[k]] <- hclust(d_ls[[i]]^2, method="ward.D")
          groups_ls[[k]] <- cutree(fit_ls[[k]], k=n_clust[k])
          grp_asg <- cutree(fit_ls[[k]], k=n_clust[k])
          clustID[i,k,j,]=grp_asg
          sil_arr[i,k,j,]=silhouette(dist = d_ls[[i]],x = cutree(fit_ls[[k]], k=n_clust[k]))[,3]
          # if(k==1){plot(fit_ls[[k]],main=dat_nms[i], cex=0.0001,xlab="",ylim=c(0,100))}
        } 
        
        if(method_nm[j]=="hierWARD1"){
          fit_ls[[k]] <- hclust(d_ls[[i]], method="ward.D")
          groups_ls[[k]] <- cutree(fit_ls[[k]], k=n_clust[k])
          grp_asg <- cutree(fit_ls[[k]], k=n_clust[k])
          clustID[i,k,j,]=grp_asg
          sil_arr[i,k,j,]=silhouette(dist = d_ls[[i]],x = cutree(fit_ls[[k]], k=n_clust[k]))[,3]
        } 
        
        if(method_nm[j]=="hierSING"){
          fit_ls[[k]] <- hclust(d_ls[[i]], method="single")
          groups_ls[[k]] <- cutree(fit_ls[[k]], k=n_clust[k])
          grp_asg <- cutree(fit_ls[[k]], k=n_clust[k])
          clustID[i,k,j,]=grp_asg
          sil_arr[i,k,j,]=silhouette(dist = d_ls[[i]],x = cutree(fit_ls[[k]], k=n_clust[k]))[,3]
        }

        if(method_nm[j]=="hierAVG"){
          fit_ls[[k]] <- hclust(d_ls[[i]], method="average")
          groups_ls[[k]] <- cutree(fit_ls[[k]], k=n_clust[k])
          grp_asg <- cutree(fit_ls[[k]], k=n_clust[k])
          clustID[i,k,j,]=grp_asg
          sil_arr[i,k,j,]=silhouette(dist = d_ls[[i]],x = cutree(fit_ls[[k]], k=n_clust[k]))[,3]
        }
        
        if(method_nm[j]=="hierCOMP"){
          fit_ls[[k]] <- hclust(d_ls[[i]], method="complete")
          groups_ls[[k]] <- cutree(fit_ls[[k]], k=n_clust[k])
          grp_asg <- cutree(fit_ls[[k]], k=n_clust[k])
          clustID[i,k,j,]=grp_asg
          sil_arr[i,k,j,]=silhouette(dist = d_ls[[i]],x = cutree(fit_ls[[k]], k=n_clust[k]))[,3]
        }
        # Divisive clustering method, as opposed to agglomerative
        if(method_nm[j]=="diana"){
          fit_ls[[k]] <- diana(d_ls[[i]])
          groups_ls[[k]] <- cutree(fit_ls[[k]], k=n_clust[k])
          grp_asg <- cutree(fit_ls[[k]], k=n_clust[k])
          clustID[i,k,j,]=grp_asg
          sil_arr[i,k,j,]=silhouette(dist = d_ls[[i]],x = cutree(fit_ls[[k]], k=n_clust[k]))[,3]
        }
        
        # Non-hierarchical non-parametric clustering based on optimizing the position of principle points
        if(method_nm[j]=="k-means"){
          fit_ls[[k]] <- kmeans(d_ls[[i]],centers=n_clust[k])
          groups_ls[[k]] <- fit_ls[[k]]$cluster
          grp_asg <- fit_ls[[k]]$cluster
          clustID[i,k,j,]=grp_asg
          sil_arr[i,k,j,]=silhouette(dist = d_ls[[i]],x = fit_ls[[k]]$cluster, k=n_clust[k])[,3]
          # print(fit_ls[[k]]$tot.withinss)
        }
        
        # Non-hierarchical model based clustering assuming multivariate normals
        if(method_nm[j]=="MCclust"){
          # NOTE-- The distance matrix calculation occurs inside the Mclust() function
          fit_ls[[k]]=Mclust(data_sets[[i]]$for_anly,G = n_clust[k],verbose = F)
          grp_asg=fit_ls[[k]]$classification
          clustID[i,k,j,]=grp_asg
          sil_arr[i,k,j,]=silhouette(dist = d_ls[[i]],x = grp_asg, k=n_clust[k])[,3]
        }

        # Table of predator assignments
        cls_mat=t(table(data_sets[[i]]$with_cats$pred_score>1,grp_asg))

        DF[[i]]=data_sets[[i]]$with_cats
        DF[[i]]=cbind(DF[[i]],grp_asg)
        names(DF[[i]])[ncol(DF[[i]])]=paste("mtd=",method_nm[j],"ncls=",n_clust[k],sep="")
        
        # purity and gini estimates for evaluating 
        purity_measDF=rbind(purity_measDF,data.frame(dataset=dat_nms[i],Method=method_nm[j],num_clusts=n_clust[k],
                                                     Metric="Purity",
                                                     Value=sum(apply(cls_mat,1,max))/sum(cls_mat)))
        purity_measDF=rbind(purity_measDF,data.frame(dataset=dat_nms[i],Method=method_nm[j],num_clusts=n_clust[k],
                                                     Metric="Gini",
                                                     Value=gini(cls_mat)))
        big_rw=cls_mat[which.max(rowSums(cls_mat)),]
        purity_measDF=rbind(purity_measDF,data.frame(dataset=dat_nms[i],Method=method_nm[j],num_clusts=n_clust[k],
                                                     Metric="P_pred",
                                                     Value=big_rw[2]/sum(big_rw)))
        # setTxtProgressBar(pb1,j+length(method_nm)*(k-1)+(length(n_clust)*length(method_nm))*(i-1))
    }
    }
    }

  out=list(data_sets,DF,purity_measDF,sil_arr,method_nm,dat_nms,n_clust,clustID,fit_ls)
  names(out)=c("data_sets","clust_lab_DFs","purity_measDF","sil_arr","method_nm","dat_nms","n_clust","clustID","fit_ls")
  
  print(proc.time()-bt)
  return(out)}


# FUNCTION THAT RUNS A CLUSTER ANALYSIS
clust_anlys_NO_PURITY=function(data_sets,method_nm,n_clust,distance_mat="euclidean"){
  
  require(StatMatch) # for gower.distance
  require(cluster)  # for diana()
  require(utils) # for progress bar
  require(mclust)
  par(mfrow=c(1,2))
  bt=proc.time()
  dat_nms=names(data_sets)
  stopifnot(!is.null(dat_nms))
  # Applies all methods
  if(any(method_nm=="all")){method_nm=c("diana","hierWARD2","hierAVG","hierWARD1","hierCOMP","k-means","hierSING")} #,"hierWARD","MCclust"
  fit_ls=d_ls=groups_ls=cls_mat_ls=DF=list()
  n_rws=nrow(data_sets[[1]]$for_anly)
  
  # Matrix of silhouette values for all configurations of clusters and methods
  # sil_mat=gap_mat=matrix(nrow=length(data_sets),ncol=length(n_clust))
  
  # Empty Array of cluster assignments to datasets
  clustID=array(dim=c(length(dat_nms),length(n_clust),length(method_nm),n_rws),
                dimnames=list(dat_nms,n_clust,method_nm,1:n_rws))
  # print(str(clustID))
  # Array of silhouette values
  sil_arr=clustID
  
  DF_tag=DF_vis=list();
  purity_measDF=NULL  
  
  # pb1 <- txtProgressBar(min = 1, max = length(data_sets)*length(method_nm)*length(n_clust), style = 3)
  for(i in 1:length(data_sets)){
    # Creates a unique PDF for each file
    # pdf(width = 8,height = 10.5,file = paste(names(data_sets)[i],"_clusters.pdf",sep=""))
    
    # distance measures
    if(distance_mat=="gower"){  d_ls[[i]] <- daisy(data_sets[[i]]$for_anly,metric="gower")} 
    if(distance_mat=="euclidean"){  d_ls[[i]] <- dist(data_sets[[i]]$for_anly, method = "euclidean")} 
    
    
    for(k in 1:length(n_clust)){
      for(j in 1:length(method_nm)){
        grp_asg=c() # temporary vector of cluster assignments
        
        
        # Initial clustering method that I used that seemed to intuitively perform Ward clustering
        # However, according to Murtaugh and Legendre (2014), this obvious choice  does not use
        # a squared distance matrix and therefore does not actually do the Ward (1963) clustering,
        # and isn't in the same space as a PCA.
        # if(method_nm[j]=="hierWARD"){
        #   fit_ls[[k]] <- hclust(d_ls[[i]], method="ward.D")
        #   groups_ls[[k]] <- cutree(fit_ls[[k]], k=n_clust[k])
        #   grp_asg <- cutree(fit_ls[[k]], k=n_clust[k])
        #   clustID[i,k,j,]=grp_asg
        #   sil_arr[i,k,j,]=silhouette(dist = d_ls[[i]],x = cutree(fit_ls[[k]], k=n_clust[k]))[,3]
        # } 
        
        if(method_nm[j]=="hierWARD2"){
          fit_ls[[k]] <- hclust(d_ls[[i]], method="ward.D2")
          # EQUIVALENT TO: fit_ls[[k]] <- hclust(d_ls[[i]]^2, method="ward.D")
          groups_ls[[k]] <- cutree(fit_ls[[k]], k=n_clust[k])
          grp_asg <- cutree(fit_ls[[k]], k=n_clust[k])
          clustID[i,k,j,]=grp_asg
          sil_arr[i,k,j,]=silhouette(dist = d_ls[[i]],x = cutree(fit_ls[[k]], k=n_clust[k]))[,3]
          # if(k==1){plot(fit_ls[[k]],main=dat_nms[i], cex=0.0001,xlab="",ylim=c(0,100))}
        } 
        
        if(method_nm[j]=="hierWARD1"){
          fit_ls[[k]] <- hclust(d_ls[[i]], method="ward.D")
          groups_ls[[k]] <- cutree(fit_ls[[k]], k=n_clust[k])
          grp_asg <- cutree(fit_ls[[k]], k=n_clust[k])
          clustID[i,k,j,]=grp_asg
          sil_arr[i,k,j,]=silhouette(dist = d_ls[[i]],x = cutree(fit_ls[[k]], k=n_clust[k]))[,3]
        } 
        
        if(method_nm[j]=="hierSING"){
          fit_ls[[k]] <- hclust(d_ls[[i]], method="single")
          groups_ls[[k]] <- cutree(fit_ls[[k]], k=n_clust[k])
          grp_asg <- cutree(fit_ls[[k]], k=n_clust[k])
          clustID[i,k,j,]=grp_asg
          sil_arr[i,k,j,]=silhouette(dist = d_ls[[i]],x = cutree(fit_ls[[k]], k=n_clust[k]))[,3]
        }
        
        if(method_nm[j]=="hierAVG"){
          fit_ls[[k]] <- hclust(d_ls[[i]], method="average")
          groups_ls[[k]] <- cutree(fit_ls[[k]], k=n_clust[k])
          grp_asg <- cutree(fit_ls[[k]], k=n_clust[k])
          clustID[i,k,j,]=grp_asg
          sil_arr[i,k,j,]=silhouette(dist = d_ls[[i]],x = cutree(fit_ls[[k]], k=n_clust[k]))[,3]
        }
        
        if(method_nm[j]=="hierCOMP"){
          fit_ls[[k]] <- hclust(d_ls[[i]], method="complete")
          groups_ls[[k]] <- cutree(fit_ls[[k]], k=n_clust[k])
          grp_asg <- cutree(fit_ls[[k]], k=n_clust[k])
          clustID[i,k,j,]=grp_asg
          sil_arr[i,k,j,]=silhouette(dist = d_ls[[i]],x = cutree(fit_ls[[k]], k=n_clust[k]))[,3]
        }
        # Divisive clustering method, as opposed to agglomerative
        if(method_nm[j]=="diana"){
          fit_ls[[k]] <- diana(d_ls[[i]])
          groups_ls[[k]] <- cutree(fit_ls[[k]], k=n_clust[k])
          grp_asg <- cutree(fit_ls[[k]], k=n_clust[k])
          clustID[i,k,j,]=grp_asg
          sil_arr[i,k,j,]=silhouette(dist = d_ls[[i]],x = cutree(fit_ls[[k]], k=n_clust[k]))[,3]
        }
        
        # Non-hierarchical non-parametric clustering based on optimizing the position of principle points
        if(method_nm[j]=="k-means"){
          fit_ls[[k]] <- kmeans(d_ls[[i]],centers=n_clust[k])
          groups_ls[[k]] <- fit_ls[[k]]$cluster
          grp_asg <- fit_ls[[k]]$cluster
          clustID[i,k,j,]=grp_asg
          sil_arr[i,k,j,]=silhouette(dist = d_ls[[i]],x = fit_ls[[k]]$cluster, k=n_clust[k])[,3]
          # print(fit_ls[[k]]$tot.withinss)
        }
        
        # Non-hierarchical model based clustering assuming multivariate normals
        if(method_nm[j]=="MCclust"){
          # NOTE-- The distance matrix calculation occurs inside the Mclust() function
          fit_ls[[k]]=Mclust(data_sets[[i]]$for_anly,G = n_clust[k],verbose = F)
          grp_asg=fit_ls[[k]]$classification
          clustID[i,k,j,]=grp_asg
          sil_arr[i,k,j,]=silhouette(dist = d_ls[[i]],x = grp_asg, k=n_clust[k])[,3]
        }
        
        # Table of predator assignments
        # cls_mat=t(table(data_sets[[i]]$with_cats$pred_score>1,grp_asg))
        # 
        # DF[[i]]=data_sets[[i]]$with_cats
        # DF[[i]]=cbind(DF[[i]],grp_asg)
        # names(DF[[i]])[ncol(DF[[i]])]=paste("mtd=",method_nm[j],"ncls=",n_clust[k],sep="")
        # 
        # # purity and gini estimates for evaluating 
        # purity_measDF=rbind(purity_measDF,data.frame(dataset=dat_nms[i],Method=method_nm[j],num_clusts=n_clust[k],
        #                                              Metric="Purity",
        #                                              Value=sum(apply(cls_mat,1,max))/sum(cls_mat)))
        # purity_measDF=rbind(purity_measDF,data.frame(dataset=dat_nms[i],Method=method_nm[j],num_clusts=n_clust[k],
        #                                              Metric="Gini",
        #                                              Value=gini(cls_mat)))
        # big_rw=cls_mat[which.max(rowSums(cls_mat)),]
        # purity_measDF=rbind(purity_measDF,data.frame(dataset=dat_nms[i],Method=method_nm[j],num_clusts=n_clust[k],
        #                                              Metric="P_pred",
        #                                              Value=big_rw[2]/sum(big_rw)))
        # setTxtProgressBar(pb1,j+length(method_nm)*(k-1)+(length(n_clust)*length(method_nm))*(i-1))
      }
    }
  }
  
  # out=list(data_sets,DF,#purity_measDF,
  #          sil_arr,method_nm,dat_nms,n_clust,clustID,fit_ls)
  # names(out)=c("data_sets","clust_lab_DFs",#"purity_measDF",
  #              "sil_arr","method_nm","dat_nms","n_clust","clustID","fit_ls")
  # 
  # print(proc.time()-bt)
  # return(out)
  
  }

# SOME PLOTTING FUNCTIONS FOR COMBINED DENDROGRAMS AND BARPLOTS

# dendextend:::zero_range()
rescale=function (x, to = c(0, 1), from = range(x, na.rm = TRUE)){
  if (dendextend:::zero_range(from) || dendextend:::zero_range(to)) {
    return(rep(mean(to), length(x)))}
  (x - from[1])/diff(from) * diff(to) + to[1]}

colored_barsALT=function (colors, dend, rowLabels = NULL, cex.rowLabels = 0.9, 
                          add = TRUE, y_scale, y_shift, text_shift = 1, sort_by_labels_order = TRUE, 
                          horiz = FALSE,bar_wid=2, ...){
  n_colors <- if (is.null(dim(colors))) 
    length(colors)
  else nrow(colors)
  n_groups <- if (is.null(dim(colors))) 
    1
  else ncol(colors)
  if (!missing(dend)) {
    if (is.hclust(dend)) 
      dend <- as.dendrogram(dend)
    if (!is.dendrogram(dend)) 
      stop("'dend' should be a dendrogram.")
    dend_labels <- labels(dend)
    dend_order <- order.dendrogram(dend)
  }
  else {
    dend_labels <- rep("W", n_colors)
    dend_order <- seq_len(n_colors)
  }
  if (!sort_by_labels_order) 
    dend_order <- seq_len(n_colors)
  # if (!horiz) {
  #     if (missing(y_shift)) 
  #         y_shift <- -max_labels_height(dend_labels) + par("usr")[3L] - 
  #             strheight("X")
  #     if (missing(y_scale)) 
  #         y_scale <- strheight("X") * n_groups
  # }
  # else {
  #     if (missing(y_shift)) 
  #         y_shift <- -(min(strwidth(dend_labels)) + par("usr")[2L] + 
  #             strwidth("X"))
  #     if (missing(y_scale)) 
  y_scale <- strwidth("X") * n_groups*bar_wid # where I widen the bars
  # }
  y_shift <- y_shift - y_scale
  colors <- as.matrix(colors)
  dimC <- dim(colors)
  if (is.null(rowLabels) & (length(dimnames(colors)[[2]]) == 
                            dimC[2])) 
    rowLabels <- names(as.data.frame(colors))
  op <- options()
  pr <- par(no.readonly = TRUE)
  options(stringsAsFactors = FALSE)
  par(xpd = TRUE)
  if (length(dend_order) != dimC[1]) {
    stop("ERROR: length of colors vector not compatible with number of objects in the hierarchical tree.")
  }
  C <- colors[dend_order, ]
  C <- as.matrix(C)
  step <- 1/(n_colors - 1)
  ystep <- 1/n_groups
  if (!add) {
    barplot(height = bar_wid, col = "white", border = FALSE, 
            space = 0, axes = FALSE, ...)
  }
  charWidth <- strwidth("W")/2
  charHeight <- strheight("W")/2
  for (j in 1:n_groups) {
    ind <- (1:n_colors)
    xl <- (ind - 1.5) * step
    xr <- (ind - 0.5) * step
    yb <- rep(ystep * (j - 1), n_colors)
    yt <- rep(ystep * j, n_colors)
    if (add) {
      xl <- rescale(xl, to = c(1 - 0.5, n_colors - 0.5))
      xr <- rescale(xl, to = c(1 + 0.5, n_colors + 0.5))
      yb <- yb * y_scale + y_shift
      yt <- yt * y_scale + y_shift
    }
    if (horiz) {
      rect(-yb, xl, -yt, xr, col = as.character(C[, j]), 
           border = as.character(C[, j]))
      par(srt = 90)
      if (is.null(rowLabels)) {
        s <- as.character(j)
        text(s, pos = 1, offset = 0.5, y = charHeight * 
               text_shift - dendextend:::rotated_str_dim(s)[2]/2, x = -(ystep * 
                                                                          (j) * y_scale + y_shift), cex = cex.rowLabels)
      }
      else {
        s <- rowLabels[j]
        text(s, pos = 1, offset = 0.5, y = charHeight * 
               text_shift - dendextend:::rotated_str_dim(s)[2]/2, x = -(ystep * 
                                                                          (j) * y_scale + y_shift), cex = cex.rowLabels)
      }
    }
    # else {
    #     rect(xl, yb, xr, yt, col = as.character(C[, j]), 
    #         border = as.character(C[, j]))
    #     if (is.null(rowLabels)) {
    #         text(as.character(j), pos = 2, x = charWidth * 
    #           text_shift, y = ystep * (j - 0.5) * y_scale + 
    #           y_shift, cex = cex.rowLabels)
    #     }
    #     else {
    #         text(rowLabels[j], pos = 2, x = charWidth * text_shift, 
    #           y = ystep * (j - 0.5) * y_scale + y_shift, 
    #           cex = cex.rowLabels)
    #     }
    # }
  }
  for (j in 0:n_groups) {
    the_x <- rescale(c(0, 1), to = c(1 - 0.5, n_colors + 
                                       0.5))
    if (horiz) {
      lines(y = the_x, x = -(c(ystep * j, ystep * j) * 
                               y_scale + y_shift))
    }
    else {
      lines(x = the_x, y = c(ystep * j, ystep * j) * y_scale + 
              y_shift)
    }
  }
  options(op)
  par(pr)
  return(invisible(C))
}


mytheme=theme(axis.line = element_line(colour = "black"),
                      axis.title=element_text(colour = "black",size = 18),
                      axis.text = element_text(colour = "black",size=12),
                      axis.ticks = element_line(colour = 'black'),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border =  element_rect(fill=NA,colour = "black"), #element_blank(),
                      panel.background = element_blank(),
                      plot.margin = unit(c(1,1,1,1), "cm"),
                      strip.text.x = element_text(size = 12),
                      strip.background = element_blank(),
                      axis.title.y=element_text(margin=margin(0,20,0,0)),
                      axis.title.x=element_text(margin=margin(20,0,0,0)),
                      legend.key = element_blank(),
                      legend.title = element_text(size=16),
                      legend.text =element_text(size=12),
                      plot.title = element_text(hjust = 0.5,size=28))