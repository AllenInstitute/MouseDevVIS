
###############################################
### Function: doing clustering using Seurat ###
###############################################

SeuratClustering <- function(obj, nfeatures = 2000, resolution = 1, k.filter = 200, n_dim = 50, min.dist = 0.75){
    
    if(length(table(obj$group))!=1 & all(table(obj$group)>100)){
        
        obj.list <- SplitObject(object = obj, split.by = "group")
        
        for (i in 1:length(x = obj.list)) {
            obj.list[[i]] <- NormalizeData(object = obj.list[[i]], verbose = FALSE)
            obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]], 
                                                  selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
        }
        
        reference.list <- obj.list[names(table(obj$group))]
        obj.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:n_dim, k.filter = k.filter,reduction ="cca")
        
        obj.integrated <- IntegrateData(anchorset = obj.anchors, dims = 1:n_dim)
        
        DefaultAssay(object = obj.integrated) <- "integrated"
        
        obj <- obj.integrated 
        
    } else {
        
        obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
        obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeatures)
        
    }
    
    
    obj <- ScaleData(object = obj, verbose = FALSE)
    obj <- RunPCA(object = obj, npcs = n_dim, verbose = FALSE)
    obj <- FindNeighbors(object = obj, dims = 1:n_dim)
    obj <- FindClusters(object = obj, resolution = resolution, save.SNN = T, do.sparse = T)
    obj <- RunUMAP(object = obj, reduction = "pca", reduction.name = "umap3d",dims = 1:n_dim, n.components = 3, min.dist = min.dist)
    obj <- RunUMAP(object = obj, reduction = "pca", reduction.name = "umap2d",dims = 1:n_dim, n.components = 2, min.dist = min.dist)
    
    obj
}


###############################################
### Function: compute entropy of each cell #####
###############################################

library(sp)
library(rgeos)
#library(scrattch.hicat,lib="/allen/programs/celltypes/workgroups/rnaseqanalysis/yuangao/conda-env/R/lib/R/library")
library(DescTools)
library(BiocNeighbors)
get_cell_entropy <- function(rd.dat,cl,k = 15, replication_times=500,entropy.threshold=0,method = "Annoy.Cosine"){
  
  knn.res = list()
  if(method == "Annoy.Cosine"){
    knn.result = get_knn_batch(dat=rd.dat[names(cl), ], ref.dat = rd.dat[names(cl),], k=k, method = "Annoy.Cosine", batch.size = 500000, mc.cores=20, transposed=FALSE,  return.distance=TRUE)
    knn.res$index = knn.result$index
    knn.res$distance = knn.result$distance
    }
  if(method == "kd_tree"){
    knn.result <- get.knnx(rd.dat, rd.dat, k = k)
    knn.res$index = knn.result$nn.index
    knn.res$distance = knn.result$nn.dist
  }
  
  knn.res.index = knn.res$index
  knn.res.distance = knn.res$distance
  
  tmp1 <- matrix(NA,nrow(knn.res.index),ncol(knn.res.index))
  for(i in 1:k){
    tmp1[,i] <- cl.clean[knn.res.index[,i]]
  }
  
  cell.entropy <- apply(tmp1, 1, function(x){Entropy(table(x), base=exp(1))}) 
  return(list(knn.res = knn.res, cell.entropy = cell.entropy))
  }


#####################################################
### Function: finding ancestor node for each node ###
#####################################################

BuildTrajectory <- function(emb, pd, reduction="umap", replication_times=500, removing_cells_ratio=0.2, method="knn", k_neigh = 5,replace = FALSE,filter.entropy=FALSE,filter.method="IQR",med.outlier.th=2,iqr.outlier.th=1.5,filter.cluster=FALSE , cl.sample.size = 200){
    
    print(paste("Before cell entropy filtering, total cell number is ",dim(emb)[1]," ",dim(emb)[2]))
    if(!"Anno" %in% names(pd) | !"day" %in% names(pd)) {print("Error: no Anno or day in pd before entropy filtering")}
    if(sum(rownames(pd)!=rownames(emb))!=0) {print("Error: rownames are not matched before entropy filtering")}

    if(filter.entropy==TRUE){
      cl.stats = pd %>% group_by(stage,Anno) %>% dplyr::summarize(med = median(entropy),
                                                           mad = mad(entropy),
                                                           max = max(entropy),
                                                           q1 = quantile(entropy,0.25),
                                                           iqr = IQR(entropy),
                                                           q3 = quantile(entropy,0.75))
      
      cl.stats = cl.stats %>% mutate(medth = med + med.outlier.th * mad)
      cl.stats = cl.stats %>% mutate(iqrth = q3 + iqr.outlier.th * iqr)
      
      
      if(filter.method=="IQR"){
        pd.expand = pd %>% left_join(cl.stats) %>% mutate(outlier = (entropy > iqrth))
        include.cells = pd.expand %>% filter(outlier == FALSE) %>% pull(sample_id)
        
      }
      if(filter.method=="MAD"){
        pd.expand = pd %>% left_join(cl.stats) %>% mutate(outlier = (entropy > medth))
        include.cells = pd.expand %>% filter(outlier == FALSE) %>% pull(sample_id)
      }
      pd = pd[include.cells,]
      emb = emb[include.cells,]
    }
  
  if(filter.cluster==TRUE){
    
    new_pd <- pd %>% group_by(stage,Anno) %>% slice_sample(n=cl.sample.size)
    
    include.cells = new_pd %>% pull(sample_id)
    pd = pd[include.cells,]
    emb = emb[include.cells,]
  }
  
    print(paste("After cell entropy filtering, total cell number is ",dim(emb)[1]," ",dim(emb)[2]))
    if(!"Anno" %in% names(pd) | !"day" %in% names(pd)) {print("Error: no Anno or day in pd after entropy filtering")}
    if(sum(rownames(pd)!=rownames(emb))!=0) {print("Error: rownames are not matched after entropy filtering")}
    
    pd$state = pd$Anno
    
    res = list()
    
    rep_i = 1
    
    if(method=="knn"){
      while(rep_i < (replication_times+1)){
        
        print(rep_i)
        #sample to subset the obs
        sampling_index = sample(1:nrow(pd),round(nrow(pd)*(1-removing_cells_ratio)),replace = replace)
        
        emb_sub = emb[sampling_index,]
        pd_sub = pd[sampling_index,]
        # pd_sub1 and pd_sub2 are anno table of pre and nex
        irlba_pca_res_1 <- emb_sub[as.vector(pd_sub$day)=="pre",]
        irlba_pca_res_2 <- emb_sub[as.vector(pd_sub$day)=="nex",]
        pd_sub1 <- pd_sub[pd_sub$day == "pre",]
        pd_sub2 <- pd_sub[pd_sub$day == "nex",]
        
        pre_state_min = min(table(as.vector(pd_sub1$state)))
        
        if (pre_state_min < k_neigh & pre_state_min >= 3){
            k_neigh = pre_state_min
            print(k_neigh)
        }
        
        #if (pre_state_min < 3){
        #    next
        #}
        #search from child to parent
        neighbors <- get.knnx(irlba_pca_res_1, irlba_pca_res_2, k = k_neigh)$nn.index
        #tmp1 is a matrix storing top k_neigh nearest clusters for each cell in nex stage, each row=each cell
        tmp1 <- matrix(NA,nrow(neighbors),ncol(neighbors))
        for(i in 1:k_neigh){
          tmp1[,i] <- as.vector(pd_sub1$state)[neighbors[,i]]
        }
        state1 <- names(table(as.vector(pd_sub1$state)))
        state2 <- names(table(as.vector(pd_sub2$state)))
        
        tmp2 <- matrix(NA,length(state2),length(state1))
        #loop through each cluster of nex, x is the cells connect to state2[2]
        for(i in 1:length(state2)){
          x <- c(tmp1[as.vector(pd_sub2$state)==state2[i],])
          for(j in 1:length(state1)){
            tmp2[i,j] <- sum(x==state1[j])
          }
        }
        tmp2 <- tmp2/apply(tmp2,1,sum)
        tmp2 <- data.frame(tmp2)
        row.names(tmp2) = state2
        names(tmp2) = state1
        
        res[[rep_i]] = tmp2
        
        rep_i = rep_i + 1
        
      }
    }
    
    if(method=="mnn"){
      while(rep_i < (replication_times+1)){
        
        print(rep_i)
        #sample to subset the obs
        sampling_index = sample(1:nrow(pd),round(nrow(pd)*(1-removing_cells_ratio)),replace = replace)
        
        emb_sub = emb[sampling_index,]
        pd_sub = pd[sampling_index,]
        # pd_sub1 and pd_sub2 are anno table of pre and nex
        irlba_pca_res_1 <- emb_sub[as.vector(pd_sub$day)=="pre",]
        irlba_pca_res_2 <- emb_sub[as.vector(pd_sub$day)=="nex",]
        pd_sub1 <- pd_sub[pd_sub$day == "pre",]
        pd_sub2 <- pd_sub[pd_sub$day == "nex",]
        
        pre_state_min = min(table(as.vector(pd_sub1$state)))
        
        #if (pre_state_min < k_neigh & pre_state_min >= 3){
        #    k_neigh = pre_state_min
        #    print(k_neigh)
        #}
        
        #if (pre_state_min < 3){
        #    next
        #}
        #search from both direction:irlba_pca_res_1 is pre and irlba_pca_res_2 is nex
        neighbors <- findMutualNN(irlba_pca_res_1, irlba_pca_res_2, k1 = k_neigh)
        neighbors <- as.data.frame(cbind(neighbors$first,neighbors$second))
        colnames(neighbors) <- c("pre","nex")
        #neighbors$pre_cell <- rownames(irlba_pca_res_1)[neighbors$pre]
        #neighbors$nex_cell <- rownames(irlba_pca_res_2)[neighbors$nex]
        neighbors.index <- do.call(
          rbind.data.frame,
          by(neighbors, neighbors$nex, function(df) { df$index <- seq_len(nrow(df)); df; })
        )
        
        tmp1 <- dcast(setDT(neighbors.index), nex ~ index, value.var = 'pre')
        tmp1 <- as.data.frame(tmp1)
        rownames(tmp1) = rownames(pd_sub2)[tmp1$nex]
        tmp1 = tmp1[,2:ncol(tmp1)]
        #tmp1 is a matrix storing top k_neigh nearest clusters for each cell in nex stage, each row=each cell
        #tmp1 <- matrix(NA,nrow(neighbors),ncol(neighbors))
        for(i in 1:k_neigh){
          tmp1[,i] <- as.vector(pd_sub1$state)[tmp1[,i]]
        }
        state1 <- names(table(as.vector(pd_sub1$state)))
        state2 <- names(table(as.vector(pd_sub2$state)))
        
        tmp2 <- matrix(NA,length(state2),length(state1))
        #loop through each cluster of nex, x is the pre cells connect to nex state2[i]
        tmp1_new <- tmp1
        tmp1_new[is.na(tmp1_new)] <- 0
        pd_sub2_select = pd_sub2[rownames(tmp1_new),]
        for(i in 1:length(state2)){
          x <- tmp1_new[which(pd_sub2_select$state==state2[i]),]
          
          for(j in 1:length(state1)){
            tmp2[i,j] <- sum(na.omit(x==state1[j]))
          }
          
        }
        
        rownames(tmp2) = state2
        colnames(tmp2) = state1
        #tmp2 <- tmp2/apply(tmp2,1,sum)
        tmp2 <- data.frame(tmp2)
        
        res[[rep_i]] = tmp2
        
        rep_i = rep_i + 1
        
      }
    }
    
    return(res)
}


BuildTrajectoryCell <- function(emb, pd, reduction="umap", replication_times=500, removing_cells_ratio=0.2, k_neigh = 5,replace = FALSE,filter.entropy=FALSEs,filter.method="IQR"){
  
  print(dim(emb))
  if(!"Anno" %in% names(pd) | !"day" %in% names(pd)) {print("Error: no Anno or day in pd")}
  if(sum(rownames(pd)!=rownames(emb))!=0) {print("Error: rownames are not matched")}
  pd$state = pd$Anno
  
  res = list()
  
  rep_i = 1
  
  while(rep_i < (replication_times+1)){
    
    print(rep_i)
    #sample to subset the obs
    sampling_index = sample(1:nrow(pd),round(nrow(pd)*(1-removing_cells_ratio)),replace = replace)
    
    emb_sub = emb[sampling_index,]
    pd_sub = pd[sampling_index,]
    # pd_sub1 and pd_sub2 are anno table of pre and nex
    irlba_pca_res_1 <- emb_sub[as.vector(pd_sub$day)=="pre",]
    irlba_pca_res_2 <- emb_sub[as.vector(pd_sub$day)=="nex",]
    pd_sub1 <- pd_sub[pd_sub$day == "pre",]
    pd_sub2 <- pd_sub[pd_sub$day == "nex",]
    
    pre_state_min = min(table(as.vector(pd_sub1$state)))
    
    #search from child to parent
    neighbors <- get.knnx(irlba_pca_res_1, irlba_pca_res_2, k = k_neigh)$nn.index
    #tmp1 is a matrix storing top k_neigh nearest clusters for each cell in nex stage, each row=each cell
    tmp1 <- matrix(NA,nrow(neighbors),ncol(neighbors))
    for(i in 1:k_neigh){
      tmp1[,i] <- as.vector(pd_sub1$state)[neighbors[,i]]
    }
    state1 <- names(table(as.vector(pd_sub1$state)))
    state2 <- names(table(as.vector(pd_sub2$state)))
    
    tmp2 <- matrix(NA,length(state2),length(state1))
    #loop through each cluster of nex, x is the cells connect to state2[2]
    for(i in 1:length(state2)){
      x <- c(tmp1[as.vector(pd_sub2$state)==state2[i],])
      for(j in 1:length(state1)){
        tmp2[i,j] <- sum(x==state1[j])
      }
    }
    tmp2 <- tmp2/apply(tmp2,1,sum)
    tmp2 <- data.frame(tmp2)
    row.names(tmp2) = state2
    names(tmp2) = state1
    
    res[[rep_i]] = tmp2
    
    rep_i = rep_i + 1
    
  }
  
  return(res)
}

######################################################
### Function: create tree structure
######################################################
ConvertToTree <- function(dat,parent.node,child.node,first.time){
  
  if(!first.time %in% dat$age.from){
    first.time.node = dat %>% filter(age.to == first.time) %>% pull(cl.to)
    first.time.node.data = as.data.frame(cbind(cl.from=unique(first.time.node),cl.to = "LineageRoot"))
    dat = as.data.frame(dat)
    dat = rbind(dat[,c(parent.node,child.node)],first.time.node.data)
  }
  library(data.tree)
  dat.sub <- dat[,c(parent.node,child.node)]
  devTree <- FromDataFrameNetwork(dat[,c(parent.node,child.node)])
  devTreeDf <- ToDataFrameTypeCol(devTree)
  return(list(data=dat, devTree = devTree,devTreeDf = devTreeDf))
}


