
### load libraries
library("ggplot2")
library("reshape2")
library("patchwork")
library("data.table")
library("hexbin")

### configs
set.seed(42)

# input
data_path <- snakemake@input[["data"]]
annot_path <- snakemake@input[["annotation"]]

# output
plot_path <- snakemake@output[["diag_plot"]]
cfa_plot_path <- snakemake@output[["cfa_plot"]]
cfa_results_path <- snakemake@output[["cfa_results"]]

# parameters
annot_vars <- snakemake@config[["visualization_parameters"]][["annotate"]]
split <- snakemake@wildcards[["split"]]
label <- snakemake@wildcards[["label"]]

### load data
data <- data.frame(fread(file.path(data_path), header=TRUE), row.names=1, check.names=FALSE)
annot <- data.frame(fread(file.path(annot_path), header=TRUE), row.names=1, check.names=FALSE)

# Convert numerical metadata with fewer than 25 unique values to factor AND ensure all categorical metadata are factors
annot <- as.data.frame(lapply(annot, function(x) {
  if ((is.numeric(x) && length(unique(x)) <= 25 && nrow(annot) > 25) || !is.numeric(x)) {
    return(factor(x))
  } else {
    return(x)
  }
}),
    row.names = rownames(annot),
    check.names = FALSE
)

# need to handle empty data, e.g. if no HVFs
if (nrow(data) == 0) {
    error_message <- paste("Input data is empty!",
                           "Most likely you have no HVFs, check your config under 'hvf_parameters'.",
                           sep = "\n")

    # Create a ggplot object that displays the error message
    error_plot <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = error_message, size = 5, color = "darkred", fontface="bold") +
      theme_void() +
      labs(title = "Data Processing Error") +
      theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))

    ggsave(plot_path, error_plot, width=8, height=12, dpi = 300)
    ggsave(cfa_plot_path, error_plot, width=8, height=6, dpi = 300)

    message(error_message)
    quit(save = "no", status = 0)
}

if (length(annot_vars)<1){
    annot_vars <- c(colnames(annot)[1], colnames(annot)[2])
    sample_col <- "sample"
} else if (length(annot_vars)<2){
    annot_vars <- c(annot_vars, colnames(annot)[1])
    sample_col <- annot_vars[1]
} else{
    sample_col <- annot_vars[1]
}

# raw counts have to be log-normalized
if (label=="counts" | label=="filtered"){
    data <- log2(data + 1)
}

# transform data into long format for plotting
data_long <- reshape2::melt(data = data, value.name = "counts", variable.name = "sample")

# add metadata for plotting
if(sample_col!="sample"){
    data_long$annot <- annot[data_long$sample, sample_col]
}else{
    data_long$annot <- data_long$sample
}

# Calculate mean and variance for each feature
data_mean <- rowMeans(data)
data_var <- apply(data, 1, sd)

#### DIAGNOSTIC PLOT ####

# Create mean-variance plot
mean_var_p <- ggplot(data.frame(data_mean, data_var), aes(x=data_mean, y=data_var)) +
#   geom_point(alpha=0.1, size=0.1) +
  geom_hex() +
  scale_fill_viridis_c(option="plasma", trans="log10", name=NULL, direction=-1) +
  labs(x="Mean", y="Standard Deviation", title="Mean-Variance Relationship") +
  theme_minimal() +
 theme(legend.position = c(0.9, 0.9), legend.key.size = unit(0.5, "cm"), legend.text = element_text(size=8))

# Create density plot of log-normalized counts per feature
density_p <- ggplot(data_long, aes(x=counts, color=annot, group = sample)) +
  geom_density() +
  labs(x="Log-normalized Counts", title=paste0("Density of Log-normalized Counts per Sample\ncolored by ", sample_col)) +
  theme_minimal() + 
  theme(plot.title = element_text(size = 10)) +
  guides(color="none")

# Create boxplots of (log normalized) counts of all samples
boxplots_p <- ggplot(data_long, aes(x=sample, y=counts, color=annot)) +
  geom_boxplot(outlier.size=1, outlier.stroke = 0) +
  labs(x="Sample", y="Log-normalized Counts", title=paste0("Boxplots of Log-normalized Counts per Sample colored by ", sample_col)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
       plot.title = element_text(size = 10)) + 
  guides(color="none")

# Create PCA plots colored by "batch" and "metadata"

# remove features with 0 variance to avoid error in prcomp()
data <- data[data_var!=0, ]

# Perform PCA
pca <- prcomp(t(data), center = TRUE, scale. = TRUE)
# Get variance explained by each PC
var_explained <- pca$sdev^2 / sum(pca$sdev^2)
# Prepare plotting dataframe
pca_data <- data.frame(sample=colnames(data), pca$x[colnames(data),c('PC1','PC2')], annot[colnames(data), annot_vars], check.names=FALSE)

pca_p <- list()

for (var in annot_vars){    
    p <- ggplot(pca_data, aes(x=PC1, y=PC2, color=.data[[var]])) +
        geom_point() +
        labs(x=paste0("Principal Component 1 (", round(var_explained[1] * 100, 2), "%)"),
             y=paste0("Principal Component 2 (", round(var_explained[2] * 100, 2), "%)"),
             title="Principal Component Analysis"
            ) +
        theme_minimal()

    if (length(unique(na.omit(pca_data[[var]]))) > 15) {
        p <- p + theme(legend.position = "none")
    }
    pca_p[[var]] <- p
}

# Combine all plots into one figure
figure <- (mean_var_p + density_p) / boxplots_p / (pca_p[[annot_vars[1]]] + pca_p[[annot_vars[2]]]) +
    plot_annotation(title = paste0(split,": ", label), subtitle = paste(dim(data)[2], "samples with", dim(data)[1], "features"))

# Set figure size
# options(repr.plot.width = 8, repr.plot.height = 12)
# figure

ggsave(plot_path, figure, width=8, height=12, dpi = 300)

#### CONFOUNDING FACTOR ANALYSIS PLOT ####

# don't save all PCs since very low var PCs are mostly noise and waste of storage space
# filesize like this is not a problem since we only have a metadata x PCs matrix -> kilobytes
pc_n <- min(max(which(cumsum(var_explained) >= 0.95)[1], 10), ncol(pca$x))
print(paste("Number of PCs to keep:", pc_n))
pc_data <- data.frame(pca$x[, 1:pc_n])

# Calculate percentage variance explained
var_explained_percent <- round(var_explained[1:pc_n] * 100, 1)

# Remove metadata without variation
annot <- annot[, apply(annot, 2, function(x) { length(unique(na.omit(x))) > 1 })]

# Split metadata into numeric and categorical
numeric_metadata <- annot[sapply(annot, is.numeric)]
categorical_metadata <- annot[sapply(annot, is.factor)]

# Calculate p-values for each PC and numeric metadata
if (ncol(numeric_metadata)>0){
    p_values_numeric <- sapply(pc_data, function(pc) {
        apply(numeric_metadata, 2, function(meta){
            cor.test(pc, meta, method="kendall")$p.value
        })
    })
    
    # ensure matrix
    p_values_numeric <- matrix(p_values_numeric, nrow = ncol(numeric_metadata), ncol = ncol(pc_data))
    rownames(p_values_numeric) <- colnames(numeric_metadata)
    colnames(p_values_numeric) <- colnames(pc_data)
}else{
    p_values_numeric <- matrix(nrow = 0, ncol = ncol(pc_data))
}


# Calculate p-values for each PC and categorical metadata
if (ncol(categorical_metadata)>0){
    p_values_categorical <- sapply(pc_data, function(pc) {
      apply(categorical_metadata, 2, function(meta){
          kruskal.test(pc ~ meta)$p.value
      })
    })
    
    # ensure martix
    p_values_categorical <- matrix(p_values_categorical, nrow = ncol(categorical_metadata), ncol = ncol(pc_data))
    rownames(p_values_categorical) <- colnames(categorical_metadata)
    colnames(p_values_categorical) <- colnames(pc_data)
}else{
    p_values_categorical <- matrix(nrow = 0, ncol = ncol(pc_data))
}

                       
# Combine p-values
p_values <- rbind(p_values_numeric, p_values_categorical)
min_nonzero <- min(p_values[p_values > 0])
p_values[p_values == 0] <- min_nonzero

# Adjust p-values for multiple testing
p_values_adjusted <- p.adjust(as.vector(p_values), method = "BH")
p_values_adjusted <- matrix(p_values_adjusted, nrow = nrow(p_values), ncol = ncol(p_values))
rownames(p_values_adjusted) <- rownames(p_values)
colnames(p_values_adjusted) <- colnames(p_values)

# Transform p-values to -log10(p-values)
log_p_values <- -log10(p_values_adjusted)

# save p-values and var_explained table
cfa_results <- rbind(var_explained, log_p_values)
write.csv(cfa_results, file=cfa_results_path, row.names=TRUE)
print(paste("Results saved to", cfa_results_path))

# Add variance explained to column names for plotting, and keep only first 10 PCs
pc_n_plot <- min(10, pc_n)
log_p_values_for_plot <- log_p_values[, 1:pc_n_plot]
colnames(log_p_values_for_plot) <- paste0("PC", 1:pc_n_plot, "\n(", var_explained_percent[1:pc_n_plot], "%)")

# Melt the data for ggplot
log_p_values_melted <- reshape2::melt(log_p_values_for_plot, varnames = c("Metadata", "PC"))
                       
# Perform hierarchical clustering on the rows (metadata)
hclust_rows <- hclust(dist(log_p_values_for_plot))
ordered_metadata <- rownames(log_p_values_for_plot)[hclust_rows$order]
log_p_values_melted$Metadata <- factor(log_p_values_melted$Metadata, levels = ordered_metadata)
                       
# Create a new column in log_p_values_melted to annotate metadata as numeric or categorical
log_p_values_melted$Type <- ifelse(log_p_values_melted$Metadata %in% colnames(numeric_metadata), "Numeric", "Categorical")

# Plot using ggplot2
cfa_plot <- ggplot(log_p_values_melted, aes(x = PC, y = Metadata, fill = value)) +
                       geom_tile(color="black") +
                       geom_text(aes(label = round(value, 0)), color = "black", size = 3) +
                       scale_fill_gradient2(low = "royalblue4", high = "firebrick2", mid = "white", midpoint = 0, name = "") +
                       facet_grid(Type ~ ., scales = "free", space = "free") +
                       theme_minimal() +
                       labs(title = "Statistical Association between PCs and Metadata as -log10 Adjusted P-values") +
                       theme_minimal(base_size = 10) +
                       theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
                             axis.title.x = element_blank(),
                             axis.title.y = element_blank(),
                             axis.ticks.x = element_blank(),
                             axis.ticks.y = element_blank(),
                             plot.title = element_text(size = 10),
                             legend.position = "none"
                            ) 

# determine plot size
heigth_hm <- dim(log_p_values_for_plot)[1] * 0.2 + 1
width_hm <- dim(log_p_values_for_plot)[2] * 0.75 + 1
                       
# save plot
ggsave(cfa_plot_path, cfa_plot, width=width_hm, height=heigth_hm, dpi = 300)