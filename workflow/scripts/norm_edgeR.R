
### load libraries
library(edgeR)

### configs
set.seed(42)

# input
data_path <- snakemake@input[["filtered_counts"]]

# parameters
result_path <- snakemake@params[["result_path"]]
split <- snakemake@params[["split"]]

# config
norm_parameters <- snakemake@config[["edgeR_parameters"]]
feature_annotation_path <- snakemake@config[["feature_annotation"]]

### load data
data <- read.csv(file=file.path(data_path), row.names=1)

# load subset feature annotation, if provided
if(feature_annotation_path!=""){
    feature_annotation <- read.csv(file=file.path(feature_annotation_path), row.names=1)
    feature_annotation <- feature_annotation[rownames(data),]
}

# normalize data using CalcNormFactors and CPM/RPKM for each configured method
for(method in norm_parameters[["method"]]){
    dge_tmp <- DGEList(data)
    
    # calculate normalization factors (normalized library size factors)
    dge_tmp <- calcNormFactors(dge_tmp,
                                    method = method,
                                    refColumn = if(norm_parameters[["refColumn"]]!="NULL") norm_parameters[["refColumn"]] else NULL,
                                    logratioTrim = norm_parameters[["logratioTrim"]],
                                    sumTrim = norm_parameters[["sumTrim"]],
                                    doWeighting = if(norm_parameters[["doWeighting"]]=="TRUE") TRUE else FALSE,
                                    Acutoff = norm_parameters[["Acutoff"]],
                                    p = norm_parameters[["p"]]
                                   )

    # quantification
    if (norm_parameters[["quantification"]]=="RPKM" & feature_annotation_path!=""){
        # determine (log2) RPKM values
        norm_data <- rpkm(dge_tmp,
                          gene.length = feature_annotation[[norm_parameters[["gene.length"]]]],
                          normalized.lib.sizes = TRUE,
                          log = if(norm_parameters[["log"]]=="TRUE") TRUE else FALSE,
                          prior.count = norm_parameters[["prior.count"]]
                          )
        
        if(method=="none"){
            method <- "RPKM" #norm_parameters[["quantification"]]
        }
        
    } else{
        # if (norm_parameters[["quantification"]]=="CPM"){}
        # determine (log2) CPM values
        norm_data <- cpm(dge_tmp,
                         normalized.lib.sizes = TRUE,
                         log = if(norm_parameters[["log"]]=="TRUE") TRUE else FALSE,
                         prior.count = norm_parameters[["prior.count"]]
                        )
        
        if(method=="none"){
            method <- "CPM" #norm_parameters[["quantification"]]
        }
    }
    
    # save normalized data
    write.csv(norm_data, file=file.path(result_path,split,paste0("norm",method,".csv")), row.names=TRUE)
}
