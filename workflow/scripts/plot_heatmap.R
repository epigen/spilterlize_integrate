#### load libraries
library("ComplexHeatmap")
library("circlize")
library("RColorBrewer")
library("data.table")
library("fastcluster")

### configurations
set.seed(42)
ht_opt(fast_hclust = TRUE)

# inputs
data_path <- snakemake@input[["data"]]
metadata_path <- snakemake@input[["metadata"]]

# output
plot_path <- snakemake@output[["heatmap"]]

# parameters
metadata_cols <- c(snakemake@config[["visualization_parameters"]][["annotate"]])
label <- snakemake@wildcards[["label"]]

### load data
data <- data.frame(fread(file.path(data_path), header=TRUE), row.names=1)
metadata <- data.frame(fread(file.path(metadata_path), header=TRUE), row.names=1)

# raw counts have to be log-normalized
if (label=="counts" | label=="filtered"){
    data <- log2(data + 1)
}

# determine sample wise correlation
data <- cor(data, method = "pearson")

# perform hierarchical clustering using fastcluster
hc <- fastcluster::hclust(dist(data, method = "euclidean"), method = "complete")

# prepare metadata
if(length(metadata_cols)==0){
    metadata_cols <- colnames(metadata)[1:2]
}

for (metadata_col in metadata_cols){
    # check if metadata column is only NA and switch to the first that is not
    if(all(is.na(metadata[[metadata_col]]))){
        for(col in colnames(metadata)){
            if(all(is.na(metadata[[col]]))){
                next
            }else{
                metadata_col <- col
                break
            }
        }
    }

    if (is.numeric(metadata[[metadata_col]]) & length(unique(metadata[[metadata_col]]))<=25){
        if(all(metadata[[metadata_col]] == round(metadata[[metadata_col]]))){
            metadata[metadata_col] <- as.factor(metadata[[metadata_col]])
        }
    }
    # if a metadata class is empty ("") fill with "unknown"
    if (!any(is.na(metadata[[metadata_col]]))){
        if (any(metadata[[metadata_col]]=="")){
            metadata[metadata[[metadata_col]]=="", metadata_col] <- "unknown"
        }
    }

    # replace NA values
    if (any(is.na(metadata[[metadata_col]]))){
        if (is.numeric(metadata[[metadata_col]])){
            metadata[is.na(metadata[[metadata_col]]),metadata_col] <- 0
        }else{
            metadata[is.na(metadata[[metadata_col]]),metadata_col] <- "NA"
        }
    }
}

# plot specifications
plot_dim <- min(0.2 * nrow(data) + 2,20) #maximum of 20in

# make colors for values
col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# create a color mapping for each metadata column
colors_list <- list()
for (col in metadata_cols) {
    if (!is.numeric(metadata[[col]])) {
        n_cat <- length(unique(metadata[[col]]))
        qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
        colors <- sample(unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))), n_cat, replace=TRUE)
        names(colors) <- unique(metadata[[col]])
        colors_list[[col]] <- colors
    }
}

# Create the row annotation
row_annot <- HeatmapAnnotation(df = metadata[, metadata_cols, drop = FALSE], which = "row", col = colors_list)

### make & save heatmap
# options(repr.plot.width=plot_dim, repr.plot.height=plot_dim)
png(filename=plot_path, width=plot_dim, height=plot_dim, units = "in", res=300)

Heatmap(data,
        name = "Pearson\nCorrelation",
        column_title = paste0("Hierarchically clustered heatmap of Pearson correlation matrix using method 'complete' with distance metric 'euclidean'."),
        col = col_fun,
        left_annotation = row_annot,
        show_column_dend = FALSE, 
        show_row_names = ifelse(nrow(data)>100, FALSE, TRUE),
        show_column_names = FALSE,
        cluster_rows = as.dendrogram(hc),
        cluster_columns = as.dendrogram(hc),
        row_dend_reorder = TRUE,
        column_dend_reorder = TRUE,
        use_raster = TRUE,
        raster_quality = 9
       )

dev.off()
