#This set of functions are designed to work in a 



#Descriptive functions

##Get a vector of sample names for a flowset 
getsample<-function(flowset) {
  sample=c()
  for (name in names(flowset@frames)){
    sample=append(sample,rep(name, nrow(flowset[[name]]@exprs)))
  }
  return(sample)}

##Get all informations about markers
getmarkersinfo<-function(exp){  
  print("Columns")
  print(fsApply(exp,function(x){colnames(x@exprs)}))
  print("Column and associated ID")
  print(fsApply(exp,function(x){x@parameters@data$name}))
  print("Marker")
  print(fsApply(exp,function(x){x@parameters@data$desc}))
}

##Rename channels of a expression dataframe with the associated markers
rename_markers<-function(flowset){
  copy_flowset=flowset[seq(along=flowset)]
  for (i in 1:length(copy_flowset)){
    marker.names=copy_flowset[[i]]@parameters@data$desc
    marker.names=lapply(marker.names,function(x){str_replace_all(x,"-","_")})
    colnames(copy_flowset[[i]]@exprs)<-unlist(lapply(marker.names, function(x){sapply(str_split(x,"_"),head,1)})) 
  }
  return(copy_flowset)
}

##Extract the values within a flowset in a list of dataframe
extract_values<-function(flowset){
  return(lapply(as.list(flowset@frames),function(x){x=as.data.frame(x@exprs)}))
}

#Numerical functions

##Compute a tsne object from a flowset
tsne_flowset<-function(flowset, scale,colstoignore,dims,perplexity,theta,verbose,max_iter,eta){
  expr=lapply(as.list(flowset@frames),function(x){x=as.data.frame(x@exprs[,-colstoignore])})#get values of flowset
  binded_expr=list.rbind(expr)#bind
  if (scale==TRUE){#conditional : scaling ?
    binded_expr=as.data.frame(scale(as.matrix(binded_expr),scale=TRUE,center=TRUE))} #scaling
  return(Rtsne(binded_expr, dims = dims, perplexity=perplexity,theta=theta, verbose=verbose, max_iter = max_iter,eta=eta))#return tsne
}

##Compute a umap object from a flowset
umap_flowset<-function(flowset,scale,colstoignore){
  expr=lapply(as.list(flowset@frames),function(x){x=as.data.frame(x@exprs)})
  binded_expr=list.rbind(expr)
  binded_expr=as.data.frame(scale(as.matrix(binded_expr),scale=TRUE,center=TRUE))
  if (scale==TRUE){
    binded_expr=as.data.frame(scale(as.matrix(binded_expr),scale=TRUE,center=TRUE))}
  return(umap(binded_expr[,-colstoignore]))
}

##Compute phenograph for a flowset
Rphenograph_flowset=function(flowset, colstoignore, k, scale){
  expr=lapply(as.list(flowset@frames),function(x){x@exprs})
  binded_expr=list.rbind(expr)#bind
  if (scale==TRUE){ #Conditional scaling
    binded_expr=as.data.frame(scale(as.matrix(binded_expr),scale=TRUE,center=TRUE))}
  return(Rphenograph(binded_expr[,-colstoignore],k))
}

##Apply a minmax normalization on a list. This function is used in another function
minmax_normalize <- function(x)
{
  return((x- min(x)) /(max(x)-min(x)))
}


##This function compute the pearson correlation coefficient between a cell and the profiles given by the user. This function is called within another functions
compute_corr<-function(cell_profile,celltypes_profile){
  celltypes=celltypes_profile[,c(3:ncol(celltypes_profile))]
  res=as.vector(apply(celltypes,2,function(x){
    cor(cell_profile[which(!is.na(x))],x[which(!is.na(x))])}))
  names(res)<-colnames(celltypes)
  return(res)
}

##This function returns a list of dataframe by cluster containing all the pearson correlation coefficient values between cells and given profiles. It performs zscore and minmax normalizations
correlation_cluster_celltype<-function(flowset,cluster,celltable,colstoignore){
  flowset=rename_markers(flowset)
  facsdata=list.rbind(lapply(as.list(flowset@frames),function(x){x=as.data.frame(x@exprs[,-colstoignore])}))
  facsdata=as.data.frame(scale(as.matrix(facsdata),scale=TRUE,center=TRUE))
  facsdata=as.data.frame(apply(facsdata,2, minmax_normalize))
  rownames(celltable)=toupper(celltable$X)
  colnames(facsdata)=toupper(colnames(facsdata))
  
  intersect_markers=intersect(colnames(facsdata),rownames(celltable))
  intersect_markers=intersect_markers[order(intersect_markers)]
  celltable=celltable[intersect_markers,]
  facsdata=facsdata[,intersect_markers]
  
  cluster=paste0("cluster_",cluster)
  facsdata$Cluster=as.factor(cluster)
  correlation=vector(mode="list",length(unique(cluster)))
  names(correlation)=unique(cluster)
  
  for (i in unique(cluster)){
    currclust=facsdata[facsdata$Cluster==i,]
    currclust=currclust[,-ncol(currclust)]
    clustcorr=t(apply(currclust,1,compute_corr,celltable))
    correlation[[i]]=clustcorr
  }
  return(correlation)
}

##Get a table for repartition of the sample within the clusters
sample_per_cluster=function(pheno){
  totable=cbind(sample,pheno) #Associates the vector of sample and the results of the clustering
  colnames(totable)=c("sample","cluster") 
  totable=as.data.frame(totable) # Formatting
  return(table(totable$sample,totable$cluster)) #contingency table
}

##Compute the mean value of each phenograph cluster
mean_per_cluster=function(flowset,pheno,colstoignore){
  expr=as.data.frame(list.rbind(lapply(as.list(flowset@frames),function(x){x=as.data.frame(x@exprs[,-c(colstoignore)])})))
  tohm=cbind(expr,as.factor(pheno))#To Heatmap : binds the data and the corresponding clusters
  tres=data.frame(row.names =unique(tohm$`as.factor(pheno)`))#declate a dataframe with rows correponding to cluster
  #2 for loops for clarity between index / cluster number 
  for (i in colnames(tohm)){
    for (j in unique(tohm$`as.factor(pheno)`)){
      tres[rownames(tres)==j,i]=mean(tohm[tohm$`as.factor(pheno)`==j,i])#compute the mean value of each marker through clusters
    }
  }
  return(tres[,-ncol(tres)])
}


#Graphical functions

##Plot density of each flowframe within a flowset
plotdensity_flowset<-function(flowset){ ggplot(melt(lapply(as.list(flowset@frames),function(x){x=as.data.frame(x@exprs)})), aes(x=value,y=L1,fill=L1)) + geom_density_ridges(alpha=.4,verbose=FALSE) +facet_wrap(~variable)}

##Plot a tsne with a title
plot_tsne<-function(tsne,title){
  colnames(tsne$Y)<-c("X","Y")
  ggplot(as.data.frame(tsne$Y)) + geom_point(aes(x=X, y=Y,col=sample),size=0.5)+ guides(colour = guide_legend(override.aes = list(size=3)))+ scale_colour_manual(values  = c("deepskyblue1","orange","violetred2","maroon2","brown2","dodgerblue","cyan","yellow3","darkgoldenrod2","brown3","maroon2","violetred2"))+ggtitle(title)
}

##Plot a umap with a title
plot_umap<-function(umap,title){
  colnames(umap$layout)<-c("X","Y")
  ggplot(as.data.frame(umap$layout)) + geom_point(aes(x=X, y=Y,col=sample),size=0.5)+ guides(colour = guide_legend(override.aes = list(size=3)))+ scale_colour_manual(values  = c("deepskyblue1","orange","violetred2","maroon2","brown2","dodgerblue","cyan","yellow3","darkgoldenrod2","brown3","maroon2","violetred2"))+ggtitle(title)
}

##Plot a tsne colored by phenograph cluster
plot_tsne_clust<-function(tsne,pheno,title){
  cluster=as.factor(pheno)
  colnames(tsne$Y)<-c("X","Y")
  print(ggplot(as.data.frame(tsne$Y))+
          geom_point(aes(x=X,y=Y,col=cluster),size=0.5)+
          theme_minimal()+guides(colour = guide_legend(override.aes = list(size=3)))+ scale_colour_manual(values  = c("#36dee6", # bright blue
                                                                                                                      "#6d71d8", # blue
                                                                                                                      "#5c3788", # dark purple blue 
                                                                                                                      "#b1457b", # purple
                                                                                                                      "#ca73c6", # other purple
                                                                                                                      "#006400", # dark green
                                                                                                                      "#54b06c", # green
                                                                                                                      "#799e43", # olive green
                                                                                                                      "#c1a339", # yellow mustard
                                                                                                                      "#b86738", # orange
                                                                                                                      "Khaki3", 
                                                                                                                      "brown3",
                                                                                                                      "darkgrey",
                                                                                                                      "darkorchid3",
                                                                                                                      "chocolate",
                                                                                                                      "chartreuse3",
                                                                                                                      "bisque3",
                                                                                                                      "aquamarine2",
                                                                                                                      "mediumspringgreen",
                                                                                                                      "peru",
                                                                                                                      "peachpuff4","yellowgreen","thistle2","olivedrab1","powderblue","rosybrown1","rosybrown3","lawngreen","mistyrose","gold4"
          ))+ggtitle(title))
}


##This function is used to plot the marker intensity on one tsne. It is not used directly but in a loop in another function
plot_marker<-function(marker,tsne,lexpr){
  intensity=lexpr[,marker]
  colnames(tsne$Y)<-c("X","Y")
  print(ggplot(as.data.frame(tsne$Y))+
          geom_point(aes(x=X,y=Y,col=intensity),size=0.2)+scale_colour_gradient2(
            low = muted("blue"),
            mid = "green",
            high = muted("red"),
            midpoint = 0,
            space = "Lab",
            na.value = "grey50",
            guide = "colourbar",
            aesthetics = "colour"
          )+ggtitle(ggtitle(paste0("\nmarker: ",colnames(lexpr)[marker])))+
          theme_minimal())
}

plot_marker_tsne=function(flowset,tsne,colstoignore){
  flo=rename_markers(flowset)
  expr=as.data.frame(list.rbind(lapply(as.list(flo@frames),function(x){x=as.data.frame(x@exprs[,-c(colstoignore)])})))
  myplots<-lapply(c(1:length(colnames(expr))), plot_marker, tsne=tsne,lexpr=expr)
  n <- length(myplots) #Get the dimensions so the final figure looks like a square
  nCol <- floor(sqrt(n))
  do.call("grid.arrange", c(myplots, ncol=nCol))
}


heatmap_cluster<-function(flowset, pheno_clustering,colstoignore){
  expr=cbind(as.data.frame(list.rbind(lapply(as.list(flowset@frames),function(x){x=as.data.frame(x@exprs[,-c(colstoignore)])}))),as.factor(pheno_clustering))
  expr=expr[order(expr$`as.factor(pheno_clustering)`),]
  expr$cell=as.factor(c(1:nrow(expr)))
  to_hm=as.data.frame(melt(expr))
  colnames(to_hm)<-c("cluster","cell","marker","value")
  to_hm=to_hm[order(to_hm$cluster),]
  ggplot(to_hm,aes(x=cell,y=marker,fill=value))+geom_raster(size=0.5)+  scale_fill_gradient2(low="cyan2", high ="darkorchid3")+scale_x_discrete(labels=to_hm$cluster,position = "top")+theme(axis.text.x = element_blank())+facet_grid(~cluster,scales="free",drop=TRUE)
}
