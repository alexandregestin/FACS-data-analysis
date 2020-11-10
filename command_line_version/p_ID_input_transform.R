args=commandArgs(trailingOnly = TRUE)

if (args[1]=="-h"){
  cat("Help : \n 
      This script is used for subsetting, arcsinh transform and align multiple fcs files\n
      \t parameters:
      1: Input absolute path to folder containing multiple fcs files only
      2: Output path for storing subset transformed aligned fcs files
      3: Desired size for subsetting (min for using the smallest sample, -1 for not subsetting)\n")
  quit()
}

require("flowCore")
require("flowWorkspace")
require("flowTrans")
require("flowStats")

#Get parameters
data_path=args[1]
output_path=args[2]
desired_size=args[3]

#Read in the data
fs=read.flowSet(path=data_path,transformation =FALSE ,emptyValue = FALSE,package="flowWorkspace")
fs=fsApply(fs,function(x){x[ ,grepl("^[{'FJComp'}|{'FCS'}|{'SSC'}].*A$",colnames(fs))]})

#Filter by size
if (desired_size=='min'){
  desired_size=min(fsApply(fs,function(x){nrow(x@exprs)}))
  sf=sampleFilter(filterId = "SizeFilter", size =as.numeric(desired_size))#establish a "filter" to subset
  fs=fsApply(fs,function(x){Subset(x,sf)})
} else if (desired_size!=-1) {
  set.seed(42) #Set a seed for reproducibility 
  sf=sampleFilter(filterId = "SizeFilter", size =as.numeric(desired_size))#establish a "filter" to subset
  fs=fsApply(fs,function(x){Subset(x,sf)})#apply the filter
}

#biexp function with default parameters, equivalent to arcsinh function. More parametrizable, it could potentially use more parameters in further versions.
biexp <- biexponentialTransform("biexp transform")
transformed_fs <- transform(fs, transformList(colnames(fs), biexp))

#Alignment step, 2 and 1 peaks
normtr=gaussNorm(transformed_fs,colnames(transformed_fs)[c(3,5:6,9:length(colnames(fs)))],max.lms = 2,peak.density.thr = 0.01)
fs_trans_norm2=normtr$flowset

normtr=gaussNorm(fs_trans_norm2,colnames(fs_trans_norm2)[c(4,7:8)],max.lms = 1,peak.density.thr = 0.05)
fs_trans_al=normtr$flowset


log.text=cat(as.character(Sys.Date()),"\n\nImport ",length(fs)," fcs files : \n \t",as.vector(sampleNames(fs)),"\n\n Data is subset to",as.character(desired_size),"cells per fcs files, arcsinh transformed and written in",output_path)

write.flowSet(fs_trans_al,outdir=output_path)

#write.table(log.text,output_path,filename=paste0("logfile-",Sys.Date()),row.names=FALSE)


