
### load libraries
library("edgeR")
library("data.table")

### configs

# input
data_path <- snakemake@input[["data"]]
annotation_path <- snakemake@input[["annotation"]]

# output
filtered_path <- snakemake@output[["filtered_counts"]]

# parameters
filter_parameters <- snakemake@params[["filter_parameters"]]

### load data
# data <- read.csv(file=file.path(data_path), row.names=1)
# annot <- read.csv(file=file.path(annotation_path), row.names=1)
data <- data.frame(fread(file.path(data_path), header=TRUE), row.names=1)
annot <- data.frame(fread(file.path(annotation_path), header=TRUE), row.names=1)

# set group variable
if(filter_parameters[["group"]]!="NULL"){
    group <- annot[[filter_parameters[["group"]]]]
}else{
    group <- NULL
}

# identify and filter features using edgeR
keep.exprs <- filterByExpr(data,
                           group = group,
                           min.count = filter_parameters[["min.count"]],
                           min.total.count = filter_parameters[["min.total.count"]],
                           large.n = filter_parameters[["large.n"]],
                           min.prop = filter_parameters[["min.prop"]]
                          )

filtered <- data[keep.exprs,]

# save data
# write.csv(filtered, file=file.path(filtered_path), row.names=TRUE)
fwrite(as.data.frame(filtered), file=file.path(filtered_path), row.names=TRUE)

# print stats
print(paste0("before: ",dim(data)[1]," x ",dim(data)[2]))
print(paste0("after: ",dim(filtered)[1]," x ",dim(filtered)[2]))