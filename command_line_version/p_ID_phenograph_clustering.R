args=commandArgs(trailingOnly = TRUE)
if (args[1]=="-h"){
  cat("Help : \n 
      This script is used for clustering multiple fcs files with phenograph algorithm\n
      \t parameters: 
      1: Input absolute path of fcs files (preprocessing recommended)
      2: Output absolute path for the clustering result
      3: K nearest neigbors for determining a relative size of communities
      4: Boolean: TRUE/FALSE indicating wether or not scaling and centering (z-score normalization) should be applied\n")
  quit()
}

require(flowCore)
require(Rphenograph)
require(rlist)



data_path=args[1]
output_path=args[2]
given_k=as.integer(args[3])
scaling=args[4]

fs=read.flowSet(path=data_path,phenoData = "annotation.txt",transformation =FALSE ,emptyValue = FALSE)

Rphenograph_flowset=function(flowset, colstoignore, k, scale){
  expr=lapply(as.list(flowset@frames),function(x){x@exprs})
  binded_expr=list.rbind(expr)#bind
  if (scale==TRUE){ #Conditional scaling
    binded_expr=as.data.frame(scale(as.matrix(binded_expr),scale=TRUE,center=TRUE))}
  return(Rphenograph(binded_expr[,-colstoignore],k))
}

clustering=Rphenograph_flowset(flowset=fs,colstoignore=c(1,2,9),k=given_k,scale=FALSE)
write.table(clustering[[2]],paste0(output_path,Sys.Date(),"-phenograph.csv"))
