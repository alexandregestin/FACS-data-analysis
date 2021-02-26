#---PREPROCESSING---#

#Help message
args=commandArgs(trailingOnly = TRUE)
if (args[1]=="-h"){
  cat("\nHelp : \n 
      This script is used for subsetting and transform multiple fcs files (flowset).\n
      \t parameters:
      1: Input absolute path to folder containing multiple fcs files only
      2: Output absolute path for storing subset transformed aligned fcs files
      3: Desired size for subsetting ('min' for using the smallest sample, 'all' for not subsetting)
      4: Type of transformation : \n\t\t. 'biexp' for biexponential\n\t\t. 'logicle' for logicle\n\t\t. 'log' for logarithm 'none' for none)\n\n
      Call example: \n Rscript subsample_transform.R /home/user/FACS_DATA/Experience1/ /home/user/FACS_OUTPUT/PHENO_ID1/ 25000 biexp\n\n")
  quit()
}

#Load libraries
require("flowCore",quietly=TRUE, warn.conflicts = FALSE)
require("flowClust",quietly=TRUE,warn.conflicts = FALSE)
require("flowWorkspace",quietly=TRUE, warn.conflicts = FALSE)
require("flowTrans",quietly=TRUE, warn.conflicts = FALSE)
require("flowStats",quietly=TRUE, warn.conflicts = FALSE)
require("rlist",quietly=TRUE, warn.conflicts = FALSE)
require("stringr",quietly=TRUE, warn.conflicts = FALSE)
require("dplyr",quietly=TRUE, warn.conflicts = FALSE)
require("ggplot2",quietly=TRUE, warn.conflicts = FALSE)
require("ggridges",quietly=TRUE, warn.conflicts = FALSE)
require("reshape2",quietly=TRUE, warn.conflicts = FALSE)

#Get parameters
input_path=args[1]
output_path=args[2]
desired_size=args[3]
transformation_type=args[4]

#Check inputs

possible_transformations=c("logicle","log","biexp","none")

if (length(args)<4){
  stop("4 arguments needed (input path, output path, desired size and type of transformation).\n\n\t Call 'subsample_transform.R -h' for help. ",call. = FALSE)
}

if(dir.exists(input_path)==FALSE){
  stop(paste0("Input path not found : ",input_path),call. = FALSE)
}

if(transformation_type %in% possible_transformations==FALSE){
  stop(paste0("Transformation: ",transformation_type," not found."),call. = FALSE)
}

if(dir.exists(output_path)==FALSE){#Check if output path exists, if not it creates it.
  dir.create(output_path) #create directory
}

#Read in the data
fs=read.flowSet(path=input_path,transformation =FALSE ,emptyValue = FALSE,package="flowWorkspace")
fs=fsApply(fs,function(x){x[ ,grepl("^[{'FJComp'}|{'FCS'}|{'SSC'}].*A$",colnames(fs))]})
#Writes counts and details in csv files
write.table(fsApply(fs,function(x){dim(x@exprs)})[,1],file=paste0(output_path,"number_cells.tsv"),row.names = TRUE) #Counts the number of cells inside each flowframe and writes it into the output folder.
write.table(apply(as.data.frame(list.rbind(lapply(sampleNames(fs),function(name){fs[[name]]@description}))),2,as.character),file=paste0(output_path,"flowset_description.csv"),col.names = TRUE)#Save descriptions of each flowframe in a csv file

#Filter by size
if (desired_size=='min'){
  desired_size=min(fsApply(fs,function(x){nrow(x@exprs)}))
  sf=sampleFilter(filterId = "SizeFilter", size =as.numeric(desired_size))#establish a "filter" to subset
  fs=fsApply(fs,function(x){Subset(x,sf)})
} else if (desired_size!='all'){
  set.seed(42) #Set a seed for reproducibility 
  sf=sampleFilter(filterId = "SizeFilter", size =as.numeric(desired_size))#establish a "filter" to subset
  fs=fsApply(fs,function(x){Subset(x,sf)})#apply the filter
}

#Transformation
if (transformation_type!='none'){
  if (transformation_type=='biexp'){
    biexp <- biexponentialTransform("biexp transform")
    transformed_fs <- transform(fs, transformList(colnames(fs), biexp))
  }
  if (transformation_type=="logicle"){
    logicle<-logicleTransform("logicle transform")
    transformed_fs <- transform (fs, transformList(colnames(fs),logicle))
  }
  if (transformation_type=="log"){
    log<-logTransform("log transform")
    transformed_fs <- transform (fs, transformList(colnames(fs),log))
  }
  dir.create(paste0(output_path,"transformed_flowset",as.character(Sys.Date())))
  write.flowSet(transformed_fs,outdir=paste0(output_path,"transformed_flowset",as.character(Sys.Date())))
}

rename_markers<-function(flowset){#Defines a function to use marker names 
  copy_flowset=flowset[seq(along=flowset)]
  for (i in 1:length(copy_flowset)){
    marker.names=copy_flowset[[i]]@parameters@data$desc
    marker.names=lapply(marker.names,function(x){str_replace_all(x,"-","_")})
    colnames(copy_flowset[[i]]@exprs)<-unlist(lapply(marker.names, function(x){sapply(str_split(x,"_"),head,1)})) 
  }
  return(copy_flowset)
}

#Writes a csv and a Rdata file to save dataset
list_of_flowframes=fsApply(rename_markers(transformed_fs),function(x){as.data.frame(x@exprs)})#Makes a list of dataframes
list_names=names(list_of_flowframes)#extract names for not calling names() function at each loop
for (index in seq_along(list_of_flowframes)){ #Iterates along the index
  list_of_flowframes[[index]]=list_of_flowframes[[index]] %>% #Using tidyverse for adding features
    mutate(Batch=list_names[index])
}
ds=list.rbind(list_of_flowframes)
ds$cell=as.factor(unlist(lapply(as.list(c(1:length(transformed_fs))),function(x){c(1:nrow(transformed_fs[[x]]@exprs))})))#add cell IDs 
write.csv(ds,file=paste0(output_path,"dataset.csv"))
saveRDS(ds,file=paste0(output_path,"dataset.rds"))

sys=Sys.info()
fileConn<-file(paste0(output_path,"call_description.txt")) #Writes txt files that describe the work 
writeLines(c(paste0(as.character(Sys.time())," ",Sys.timezone()),
            paste0("Import ", "8"," fcs files :"),
            "",
            sampleNames(fs),
            "",
            paste0("data is subset to ",as.character(desired_size)," cells per fcs file"),
            paste0("transformation used : ", transformation_type),
            "",
            paste0("written in ", output_path),
            "",
            "system informations:",
            paste0("sysname : ",sys[1]),
            paste0("release : ",sys[2]),
            paste0("version : ",sys[3]),
            paste0("node name : ",sys[4]),
            paste0("machine : ",sys[5]),
            paste0("login : ",sys[6]),
            paste0("user : ",sys[7]),
            paste0("effective_user : ",sys[8])
            ),
          fileConn)
close(fileConn)
save.image(paste0(output_path,"workspace.RData"))
