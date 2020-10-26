require("flowCore")
require("flowWorkspace")
require("flowTrans")
require("scales")
args=commandArgs(trailingOnly = TRUE)

if (args[1]=="-h"){
  cat("Help : \n 
      \t parameters: \t
      This script is used for subsetting and arcsinh transform multiple fcs files
      1:Input path\n
      2:Output Path\n
      3:desired_size (min for using the smallest sample, -1 for not subsetting)\n
      ")
  quit()
}

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
} else if (desired_size!=-1) {
  set.seed(42) #Set a seed for reproducibility 
  sf=sampleFilter(filterId = "SizeFilter", size =as.numeric(desired_size))#establish a "filter" to subset
  fs=fsApply(fs,function(x){Subset(x,sf)})#apply the filter
}

#biexp function with default parameters, equivalent to arcsinh function. More parametrizable, it could potentially use more parameters in further versions.
biexp <- biexponentialTransform("biexp transform")
transformed_fs <- transform(fs, transformList(colnames(fs), biexp))



log.text=cat(as.character(Sys.Date()),"\n\nImport ",length(fs)," fcs files : \n \t",as.vector(sampleNames(fs)),"\n\n Data is subset to",as.character(desired_size),"cells per fcs files, arcsinh transformed and written in",output_path)

write.text(log.text,output_path,filename=paste0("logfile-",Sys.Date()))
write.flowSet(transformed_fcsfiles,outdir=output_path)


