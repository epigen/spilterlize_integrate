
### load libraries
library("cqn")
library("data.table")

### configs
set.seed(42)

# input
data_path <- snakemake@input[["filtered_counts"]]
annot_path <- snakemake@input[["annotation"]]

# output
result_path <- snakemake@output[["normalized_counts"]]
plot_path <- snakemake@output[["cqn_plot"]]

# parameters
split <- snakemake@params[["split"]]

# config
norm_parameters <- snakemake@config[["cqn_parameters"]]
feature_annotation_path <- snakemake@config[["feature_annotation"]]

### load data
# data <- read.csv(file=file.path(data_path), row.names=1)
# annot <- read.csv(file=file.path(annot_path), row.names=1)
# feature_annotation <- read.csv(file=file.path(feature_annotation_path), row.names=1)
data <- data.frame(fread(file.path(data_path), header=TRUE), row.names=1)
annot <- data.frame(fread(file.path(annotation_path), header=TRUE), row.names=1)
feature_annotation <- data.frame(fread(file.path(feature_annotation_path), header=TRUE), row.names=1)


# subset feature annotation
feature_annotation <- feature_annotation[rownames(data),]

### normalize data using conditional quantile normalization (cqn)
cqn_result <- cqn(counts = data,
                  x = feature_annotation[rownames(data), norm_parameters[["x"]]],
                  lengths = feature_annotation[rownames(data), norm_parameters[["lengths"]]],
                  sizeFactors = if(norm_parameters[["sizeFactors"]]=="NULL") NULL else annot[colnames(data), norm_parameters[["sizeFactors"]]],
                  subindex = NULL,
                  tau = norm_parameters[["tau"]],
                  sqn = if(norm_parameters[["sqn"]]=="TRUE") TRUE else FALSE,
                  lengthMethod = norm_parameters[["lengthMethod"]],
                  verbose = TRUE
                 )
# determine corrected values
norm_data <- cqn_result$y + cqn_result$offset
    
# save normalized data
# write.csv(norm_data, file=file.path(result_path), row.names=TRUE)
fwrite(as.data.frame(norm_data), file=file.path(result_path), row.names=TRUE)

### plot model fit
png(file=file.path(plot_path), height=4, width=8, units="in", res=300)

if(norm_parameters[["lengthMethod"]]=="smooth"){
    par(mfrow=c(1,2))
    cqnplot(cqn_result, n = 2, xlab = norm_parameters[["lengths"]])   
}

cqnplot(cqn_result, n = 1, xlab = norm_parameters[["x"]])

x <- dev.off()
