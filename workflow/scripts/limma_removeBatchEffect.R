
### load libraries
library("limma")
library("edgeR")
library("data.table")

### configs
set.seed(42)

# input
data_path <- snakemake@input[["normalized_data"]]
annot_path <- snakemake@input[["annotation"]]

# output
result_path <- snakemake@output[["integrated_data"]]

# parameters
desired_cols <- snakemake@params[["desired"]]
unwanted_categorical_cols <- snakemake@params[["unwanted_categorical"]]
unwanted_numerical_cols <- snakemake@params[["unwanted_numerical"]]

# load data
data <- data.frame(fread(file.path(data_path), header=TRUE), row.names=1)
annot <- data.frame(fread(file.path(annot_path), header=TRUE), row.names=1)

# create design based on desired effects
design <- model.matrix(reformulate(desired_cols), data = annot)

# extract unwanted categorical effects for batch variables
if(length(unwanted_categorical_cols) == 2){
    batch <- annot[[unwanted_categorical_cols[1]]]
    batch2 <- annot[[unwanted_categorical_cols[2]]]
} else if(length(unwanted_categorical_cols) == 1){
    batch <- annot[[unwanted_categorical_cols[1]]]
    batch2 <- NULL
} else {
    batch <- NULL
    batch2 <- NULL
}

# extract unwanted numerical effects for covariates variables
if(length(unwanted_numerical_cols) > 0){
    covariates <- as.matrix(annot[, unwanted_numerical_cols, drop=FALSE])
} else {
    covariates <- NULL
}

# integrate data / remove batch effect
integrated_data <- removeBatchEffect(as.matrix(data), 
                                     batch = batch, 
                                     batch2 = batch2, 
                                     covariates = covariates, 
                                     design = design
                                    )

# save integrated data
fwrite(as.data.frame(integrated_data), file = file.path(result_path), row.names = TRUE)