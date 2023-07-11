
### load libraries
library(ggplot2)
library(reshape2)
library(patchwork)

### configs
set.seed(42)

# input
data_path <- snakemake@input[["data"]]
annot_path <- snakemake@input[["annotation"]]

# output
plot_path <- snakemake@output[["plot"]]

# parameters
annot_vars <- snakemake@config[["visualization_parameters"]][["annotate"]]
split <- snakemake@wildcards[["split"]]
label <- snakemake@wildcards[["label"]]

### load data
data <- read.csv(file=file.path(data_path), row.names=1)
annot <- read.csv(file=file.path(annot_path), row.names=1)

if (length(annot_vars)<1){
    annot_vars <- c(colnames(annot)[1], colnames(annot)[2]) 
} else if (length(annot_vars)<2){
    annot_vars <- c(annot_vars, colnames(annot)[1]) 
}

# raw counts have to be log-normalized
if (label=="counts" | label=="filtered"){
    data <- log2(data + 1)
}


# transform data into long format for plotting
data_long <- melt(data = data, value.name = "counts", variable.name = "sample")

# Calculate mean and variance for each feature
data_mean <- rowMeans(data)
data_var <- apply(data, 1, sd)

# Create mean-variance plot
mean_var_p <- ggplot(data.frame(data_mean, data_var), aes(x=data_mean, y=data_var)) +
  geom_point(alpha=0.1, size=0.1) +
  labs(x="Mean", y="Standard Deviation", title="Mean-Variance Relationship") +
  theme_minimal()

# Create density plot of log-normalized counts per feature
density_p <- ggplot(data_long, aes(x=counts, color=sample)) +
  geom_density() +
  labs(x="Log-normalized Counts", title="Density of Log1p-Values per Sample") +
  theme_minimal() + guides(color="none")

# Create boxplots of (log normalized) counts of all samples
boxplots_p <- ggplot(data_long, aes(x=sample, y=counts, color=sample)) +
  geom_boxplot() +
  labs(x="Sample", y="Log-normalized Counts", title="Boxplots of Log1p-Values per Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + guides(color="none")

# Create PCA plots colored by "batch" and "metadata"

# remove features with 0 variance to avoid error in prcomp()
data <- data[data_var!=0, ]

# Perform PCA
pca <- prcomp(t(data), center = TRUE, scale. = TRUE)
# Get variance explained by each PC
var_explained <- pca$sdev^2 / sum(pca$sdev^2)
# Prepare plotting dataframe
pca_data <- data.frame(sample=colnames(data), pca$x[colnames(data),c('PC1','PC2')], annot[colnames(data), annot_vars])

pca_p <- list()

for (var in annot_vars){    
    pca_p[[var]] <- ggplot(pca_data, aes_string(x="PC1", y="PC2", color=var)) +
        geom_point() +
        labs(x=paste0("Principal Component 1 (", round(var_explained[1] * 100, 2), "%)"),
             y=paste0("Principal Component 2 (", round(var_explained[2] * 100, 2), "%)"),
             title="Principal Component Analysis"
            ) +
        theme_minimal()
}

# Combine all plots into one figure
figure <- (mean_var_p + density_p) / boxplots_p / (pca_p[[annot_vars[1]]] + pca_p[[annot_vars[2]]]) +
    plot_annotation(title = paste0(split,": ", label), subtitle = paste(dim(data)[2], "samples with", dim(data)[1], "features"))

# Set figure size
# options(repr.plot.width = 8, repr.plot.height = 12)
# figure

ggsave(plot_path, figure, width=8, height=12, dpi = 300)

