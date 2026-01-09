
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
data <- data.frame(fread(file.path(data_path), header=TRUE), row.names=1, check.names=FALSE)
annot <- data.frame(fread(file.path(annot_path), header=TRUE), row.names=1, check.names=FALSE)

# We can only integrate data for splits where the desired_cols have more than one level.
# This allows desired_cols to be over-specified in the config and impossible cols are dropped.

valid_desired_cols <- character(0)
if (length(desired_cols) > 0) {
    # Check number of levels for each desired column
    n_levels <- sapply(desired_cols, function(col) length(unique(annot[[col]])))

    # Identify and warn about columns with a single level that will be dropped
    cols_to_drop <- n_levels <= 1
    if(any(cols_to_drop)){
        warning(paste("The following desired_cols had only a single level and were dropped:",
                      paste(desired_cols[cols_to_drop], collapse = ", ")))
    }
    valid_desired_cols <- desired_cols[!cols_to_drop]
}

# create design based on desired effects
if (length(valid_desired_cols) > 0) {
    # If valid columns remain, create the design matrix from the formula
    design <- model.matrix(reformulate(valid_desired_cols), data = annot)
} else {
    # If no valid desired columns remain, write a warning the the output file.
    warning_message <- "No valid desired_cols with multiple levels were specified. Batch correction not possible for this split."
    warning(warning_message)
    writeLines(warning_message, con = file.path(result_path))
    exit()
}
                       
# extract unwanted categorical effects for batch variables
# if if(is.null(batch) && is.null(batch2) && is.null(covariates)) removeBatchEffect returns as.matrix(x)
                       
batch <- NULL
batch2 <- NULL
if (length(unwanted_categorical_cols) >= 1) {
    if (length(unique(annot[[unwanted_categorical_cols[1]]])) > 1) {
        batch <- annot[[unwanted_categorical_cols[1]]]
    } else {
        warning("Not all batch effects can be modelled for this split (not enough levels in batch).")
    }
}
if (length(unwanted_categorical_cols) == 2) {
    if (length(unique(annot[[unwanted_categorical_cols[2]]])) > 1) {
        batch2 <- annot[[unwanted_categorical_cols[2]]]
    } else {
        warning("Not all batch effects can be modelled for this split (not enough levels in batch2).")
    }
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