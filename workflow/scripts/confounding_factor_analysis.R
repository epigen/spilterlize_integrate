
### load libraries
library("ggplot2")
library("reshape2")
library("data.table")

### configs
set.seed(42)

# input
annot_path <- snakemake@input[["annotation"]]

# output
plot_path <- snakemake@output[["cfa_plot"]]
cfa_results_path <- snakemake@output[["cfa_results"]]

# load data
annot <- data.frame(fread(file.path(annot_path), header=TRUE), row.names=1)

# Remove metadata without variation
annot <- annot[, apply(annot, 2, function(x) { length(unique(na.omit(x))) > 1 })]

# Convert numerical metadata with fewer than 25 unique values to factor AND ensure all categorical metadata are factors
annot <- as.data.frame(lapply(annot, function(x) {
  if ((is.numeric(x) && length(unique(x)) <= 25 && nrow(annot) > 25) || !is.numeric(x)) {
    return(factor(x))
  } else {
    return(x)
  }
}))

# perform pairwise statistical association testing between all annotation variables
p_values <- data.frame(var1=character(), var2=character(), p_value=numeric(), stringsAsFactors=FALSE)
var_names <- colnames(annot)

for(i in 1:(length(var_names)-1)){
  for(j in (i+1):length(var_names)){
      x <- annot[[var_names[i]]]
      y <- annot[[var_names[j]]]
      type_x <- if(is.numeric(x)) "numeric" else "categorical"
      type_y <- if(is.numeric(y)) "numeric" else "categorical"
    
    if(type_x=="numeric" && type_y=="numeric"){
        pvalue <- cor.test(x, y, method="kendall")$p.value
    } else if(type_x=="categorical" && type_y=="categorical"){
        # Perform the Fisher's exact test on a contingency table from two categorical variables using Monte Carlo simulation to compute p-values in reasonable time
        pvalue <- fisher.test(table(x, y), simulate.p.value=TRUE, B=10000)$p.value
        # alternatively, but requires count of >=5 per cell, which is not always the case with NGS samples e.g., 3 replicates per condition
        # pvalue <- chisq.test(table(x, y))$p.value
    } else {
      if(type_x=="numeric"){
        pvalue <- kruskal.test(x ~ as.factor(y))$p.value
      } else {
        pvalue <- kruskal.test(y ~ as.factor(x))$p.value
      }
    }
    p_values <- rbind(p_values, data.frame(var1=var_names[i], var2=var_names[j], p_value=pvalue, stringsAsFactors=FALSE))
  }
}

# adjust p-values for multiple testing
p_values$p_values_adjusted <- p.adjust(as.vector(p_values$p_value), method = "BH")

# handle adjusted p-value 0 exception
min_nonzero <- min(p_values$p_values_adjusted[p_values$p_values_adjusted > 0])
p_values$p_values_adjusted[p_values$p_values_adjusted == 0] <- min_nonzero

# Transform p-values to -log10(p-values)
p_values$log_p_values <- -log10(p_values$p_values_adjusted)

# create symmetric matrix of -log10(adjusted p values)
var_names <- unique(c(p_values$var1, p_values$var2))
mat <- matrix(NA, nrow=length(var_names), ncol=length(var_names), dimnames=list(var_names, var_names))
for(i in seq_len(nrow(p_values))){
  mat[p_values$var1[i], p_values$var2[i]] <- p_values$log_p_values[i]
  mat[p_values$var2[i], p_values$var1[i]] <- p_values$log_p_values[i]
}

# set diagonal to the maximum of -log10(adjusted p values) in the data for clustering
diag(mat) <- max(p_values$log_p_values, na.rm=TRUE)

# hierarchical clustering to reorder matrix
hc <- hclust(dist(mat))
mat <- mat[hc$order, hc$order]

# melt the matrix into long format
df_long <- melt(mat, varnames=c("Var1", "Var2"), value.name="log_p")

# Overlay diagonal tiles with a blocking color
df_long$diag <- df_long$Var1 == df_long$Var2

# plot statistical association between metadata using ggplot2
cfa_plot <- ggplot(df_long, aes(x=Var1, y=Var2, fill=log_p)) +
    geom_tile(color="black") +
    geom_tile(data = subset(df_long, diag), aes(x=Var1, y=Var2), fill="grey90", color="black") +
    geom_text(data = subset(df_long, !diag), aes(label = round(log_p, 0)), color = "black", size = 3) +
    scale_fill_gradient2(low = "royalblue4", high = "firebrick2", mid = "white", midpoint = 0, name = "") +
    labs(title = "Pairwise statistical association between metadata as -log10 adjusted P-values", x="", y="") +
    theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
       plot.title = element_text(size = 8),
        legend.position = "none")

# determine plot size
heigth_hm <- nrow(mat) * 0.2 + 2
                       
# save plot
# options(repr.plot.width = heigth_hm, repr.plot.height = heigth_hm)
# cfa_plot

ggsave(plot_path, cfa_plot, width=heigth_hm, height=heigth_hm, dpi = 300)
write.csv(df_long, file=cfa_results_path, row.names=TRUE)