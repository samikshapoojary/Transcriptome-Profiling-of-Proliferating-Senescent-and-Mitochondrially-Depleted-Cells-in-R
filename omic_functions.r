# Installing necessary packages and commenting  out when installed
# install.packages("ggplot2")   # For creating plots
# install.packages("ggrepel")   # For adding labels without overlap on plots
# install.packages("eulerr")    # For venn diagram
# install.packages("reshape2")  # For reshaping data (wide to long formats)
# install.packages("amap")      # For advanced clustering methods
# install.packages("BiocManager")   # Bioconductor package manager
# BiocManager::install("clusterProfiler")  # Gene ontology and pathway enrichment analysis
# BiocManager::install("org.Hs.eg.db")    # Human gene annotation package
# BiocManager::install("STRINGdb")        # STRING database for protein-protein interactions

# Loading the installed libraries
library(ggplot2)       # For data visualization 
library(ggrepel)       # For repelling labels to avoid overlap in plots
library(eulerr)        # For creating venn diagrams 
library(reshape2)      # For reshaping data (useful when working with wide-format data)
library(amap)          # For clustering algorithms
library(clusterProfiler)  # For performing gene enrichment analysis
library(org.Hs.eg.db)     # Human gene annotation database
library(STRINGdb)         # STRING database for protein interactions

# Function to read CSV files into R as a table
csv_reader = function(file) {
  data = read.table(file, header = TRUE, row.names = 1, sep = "\t") # Read the file using read.table() with header and row names
  return(data)  # Return the data as a data frame
}

# Function to merge expression matrix, differential expression results, and annotations
create_master_table = function(em, de_table, annotations, exclude_cols = NULL) {
  master_temp = merge(em, annotations, by.x=0, by.y=0) # Merge the expression matrix (EM) and gene annotations 
  master_de = merge(de_table, master_temp, by.x=0, by.y=1) # Merge the differential expression results with the master temp
  row.names(master_temp) = master_temp[,1]
  row.names(master_de) = master_de[,1] # Set gene id as row names
  master_temp = master_temp[,-1]
  master_de = master_de[,-1] # Remove the gene id column as it is already row names
  master_de = na.omit(master_de)   # Remove rows with missing data or NA values
  if (!is.null(exclude_cols)) {
    master_de = master_de[, -exclude_cols]  # Remove specified columns
  }
  sorted_order = order(master_de[,"p.adj"], decreasing=FALSE) # Sort by adjusted p-value to prioritize significant results
  master_de = master_de[sorted_order, ]
  # Add additional useful columns: mean expression, -log10(p-value), and significance 
  master_de$mean = rowMeans(master_de[, c(4:9)])  # Calculate the mean expression
  master_de$mlog10p = -log10(master_de$p)         # Calculate -log10 of the p-values
  master_de$sig = as.factor(master_de$p.adj < 0.001 & abs(master_de$log2fold) > 2.0)  # Define significance
  return(master_de)   # Return the processed master table
}

# Function to create a table of significant genes based on adjusted p-value and fold change thresholds
create_sig_genes_table = function(master_table, p_adj_threshold, log2fc_threshold) {
  sig_genes_table = master_table[master_table$p.adj < p_adj_threshold & abs(master_table$log2fold) > log2fc_threshold, ] # Filter the genes that meet the significance criteria
  return(sig_genes_table)  # Return the subset of significant genes
}

# Function to create a table with gene symbols as rownames
create_symbols_table = function(master_table, select_cols = NULL, rownames_col = NULL) {
  symbols = master_table[, select_cols] # Extract the selected columns
  row.names(symbols) = symbols[, rownames_col] # Set gene symbols as row names
  return(symbols)  # Return the table of gene symbols
}

# Function to merge two differential expression tables 
merge_de_tables = function(de_table_1, de_table_2, suffix_1, suffix_2) {
  merged_table = merge(de_table_1, de_table_2, by.x=0, by.y=0, suffixes = c(suffix_1, suffix_2)) # Merge the tables and use custom suffix
  row.names(merged_table) = merged_table[, 1] # Set the row names as gene ids
  merged_table = merged_table[, -1]# Remove the first column containing gene ids
  return(merged_table)  # Return the merged table
}

# Function to merge an already merged DE table with an expression matrix
merge_merged_de_with_em = function(merged_de, em) {
  merged_table = merge(merged_de, em, by.x=0, by.y=0) # Merge the DE table with the expression matrix (EM)
  row.names(merged_table) = merged_table[, 1] #Set the row names as gene ids
  merged_table = merged_table[, -1] # Remove the first column containing gene ids
  
  return(merged_table)  # Return the merged table
}

# Function to merge a merged DE table + expression matrix with annotations
merge_with_annotations = function(merged_de_em, annotations) {
  master = merge(merged_de_em, annotations, by.x=0, by.y=0) # Merge the already merged DE table with annotations
  names(master)[1] = "Gene_ID"  # Rename the first column to "Gene_ID" for clarity
  row.names(master) = master[, "SYMBOL"] # Set row names based SYMBOL column
  return(master)  # Return the fully merged table with annotations
}

# Creating a custom theme function for consistent and attractive formatting
my_theme = theme(
  plot.title = element_text(size = 18, face = "bold", color = "midnightblue", vjust = 1),  # Customize plot title
  axis.text.x = element_text(size = 10, color = "darkslateblue", hjust = 0.5),  # Customize x-axis labels
  axis.text.y = element_text(size = 10, color = "darkslateblue"),  # Customize y-axis labels
  axis.title.x = element_text(size = 12, color = "darkslateblue"),  # Customize x-axis title
  axis.title.y = element_text(size = 12, color = "darkslateblue"),  # Customize y-axis title
  legend.title = element_text(size = 14, face = "bold", color = "darkorchid"),  # Customize legend title
  legend.text = element_text(size = 10, color = "purple4"),  # Customize legend text
  panel.background = element_rect(fill = "whitesmoke"),  # Set light background
  panel.grid.major = element_line(color = "lavender", size = 0.5),  # Major gridlines with subtle color
  panel.grid.minor = element_line(color = "lavender", size = 0.25),  # Minor gridlines with subtle color
  plot.margin = margin(20, 20, 20, 20),  # Add margins around the plot for spacing
  strip.background = element_rect(color = "lightsteelblue", fill = "transparent"),  # Facet strip background
  strip.text = element_text(size = 12, color = "midnightblue", face = "bold")  # Facet label styling
)

# Function to save plots as PNG files to provided file path
save_plot = function(ggp, file_path, height = 800, width = 1000) { #setting default measurements incase not provided
  png(file_path, height = height, width = width)  # Open a PNG device
  print(ggp)  # Print the plot to the device
  dev.off()  # Close the PNG device and save the file
}

# Function to generate a volcano plot from differential expression data
plot_volcano = function(de_table, log2fold = log2fold, p.adj = p.adj, p_threshold, fold_threshold, save_path, height = 600, width = 800, title = "Volcano Plot of Differential Gene Expression") {
  if (sum(de_table$p.adj == 0) > 0) { #if any p.adj value is 0, changing it to avoid Inf values while taking log
    de_table$p.adj[de_table$p.adj == 0] = 0.0000000000000001  # Replacing 0 with a very small value
  }
  # Adding a 'direction' column to categorize genes into upregulated, downregulated, or non-significant
  de_table$direction = "non-significant"  # Default to non-significant
  de_table$direction[de_table$log2fold > fold_threshold & de_table$p.adj < p_threshold] = "upregulated"
  de_table$direction[de_table$log2fold < -fold_threshold & de_table$p.adj < p_threshold] = "downregulated"
  
  # Subset the significant genes (top 10 upregulated and downregulated)
  de_table_sig_up = subset(de_table, direction == "upregulated") # Upregulated genes
  de_table_sig_down = subset(de_table, direction == "downregulated") # Downregulated genes
  de_table_sig_up_top10 = de_table_sig_up[1:10,]  # Top 10 upregulated genes
  de_table_sig_down_top10 = de_table_sig_down[1:10,]  # Top 10 downregulated genes
  
  # Create the volcano plot
  ggp = ggplot(de_table, aes(x = log2fold, y = -log10(p.adj), colour = direction)) + # Plot points with the 'direction' variable for coloring
    geom_point(na.rm=T) +  
    scale_colour_manual(values = c("non-significant" = "darkgrey", "upregulated" = "midnightblue", "downregulated" = "darkviolet"), # Custom color scale for the categories
                        labels = c("Downregulated", "Non-Significant", "Upregulated"), # Add index labels
                        name = "Type") +  # Add label title
    geom_vline(xintercept = -fold_threshold, linetype = "dashed", color = "orchid", size = 0.5) +  # Threshold for downregulated genes
    geom_vline(xintercept = fold_threshold, linetype = "dashed", color = "violet", size = 0.5) +  # Threshold for upregulated genes
    geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "plum", size = 0.5) +  # p-value threshold
    labs(title = title, x = "Log2 Fold Change", y = "-Log10 P-Value") + #adding plot labels
    geom_text_repel(data = de_table_sig_up_top10, aes(label = SYMBOL), colour = "navy", size = 3, vjust = 2, hjust= -1, nudge_x = 5, nudge_y = 75, max.overlaps = 10, force = 10, show.legend = FALSE) +  # Label upregulated genes and adjust their position
    geom_text_repel(data = de_table_sig_down_top10, aes(label = SYMBOL), colour = "purple4", size = 3, vjust = 2, hjust = 1,nudge_x = -5, nudge_y = 75, max.overlaps = 10, force = 10, show.legend = FALSE) +  # Label downregulated genes and adjust their position
    xlim(c(-15, 15)) +  # Set the limits for the x-axis 
    ylim(c(-5, 350)) +  # Set the limits for the y-axis 
    my_theme  # Apply the custom theme
  save_plot(ggp, save_path, height = height, width = width)  # Save the plot to file
  return(ggp)  # Return the plot object
}

# Function to create MA Plot for Differential Expression Analysis
plot_ma = function(de_table, log2fold = log2fold, p.adj = p.adj, mean = mean, p_threshold, fold_threshold, save_path, height = 600, width = 800, title = "MA Plot of Differential Gene Expression") {
  
  # Adding a 'direction' column to categorize genes into upregulated, downregulated, or non-significant
  de_table$direction = "non-significant"  # Default category is non-significant
  de_table$direction[de_table$log2fold > fold_threshold & de_table$p.adj < p_threshold] = "upregulated"  # Upregulated genes
  de_table$direction[de_table$log2fold < -fold_threshold & de_table$p.adj < p_threshold] = "downregulated"  # Downregulated genes
  
  # Subset the significant genes (top 10 upregulated and downregulated)
  de_table_sig_up = subset(de_table, direction == "upregulated")  # Upregulated genes
  de_table_sig_down = subset(de_table, direction == "downregulated")  # Downregulated genes
  de_table_sig_up_top10 = de_table_sig_up[1:10,]  # Top 10 upregulated genes
  de_table_sig_down_top10 = de_table_sig_down[1:10,]  # Top 10 downregulated genes
  
  # Create the MA plot using ggplot
  ggp = ggplot(de_table, aes(x = log10(mean), y = log2fold, colour = direction)) + # Plot all points, color by 'direction'
    geom_point(size = 2, alpha = 0.6, shape = 19) +   
    scale_colour_manual(values = c("non-significant" = "darkgrey", "upregulated" = "midnightblue", "downregulated" = "darkviolet"), # Custom color scale for the categories
                        labels = c("Downregulated", "Non-Significant", "Upregulated"), # Add index labels
                        name = "Type") +  # Add label title
    geom_vline(xintercept = 0, linetype = "dashed", color = "orchid", size = 0.5) +  # Dashed line at mean expression = 0
    geom_hline(yintercept = 0, linetype = "dashed", color = "violet", size = 0.5) +  # Dashed line at log2fold = 0
    labs(title = title,  # Plot title
      x = "Log10 Mean Expression",  # X-axis label
      y = "Log2 Fold Change"  # Y-axis label
    ) +
    geom_text_repel(data = de_table_sig_up_top10, aes(label = SYMBOL), colour = "navy", size = 3, vjust = -0.5, hjust= 0.5, nudge_x = 2, nudge_y = 1, max.overlaps = 10, force = 10, show.legend = FALSE) +    # Label for upregulated genes and adjust their position
    geom_text_repel(data = de_table_sig_down_top10, aes(label = SYMBOL), colour = "purple4", size = 3, vjust = 0.5, hjust = 0.5,nudge_x = 2, nudge_y = -1, max.overlaps = 10, force = 10, show.legend = FALSE) +  # Label for downregulated genes and adjust their position
    ylim(c(-12.5, 12.5)) + # Set limit for y-axis to accommodate labels
    my_theme  # Apply custom theme
  save_plot(ggp, save_path, height = height, width = width) # Save the plot to the specified path
  return(ggp)  # Return the plot object
}

# Function for Hypergeometric Test between two significant groups for overlap analysis
hypergeo_test = function(sig_group_1, sig_group_2, group_1_size, group_2_size, total_size) {
  overlap = length(intersect(sig_group_1, sig_group_2)) # Calculate the overlap between the two gene groups
  p_value = phyper(overlap - 1, group_1_size, total_size - group_1_size, group_2_size, lower.tail = FALSE) # Perform the hypergeometric test to calculate the p-value for overlap significance
  return(p_value)  # Return the p-value from the hypergeometric test
}

# Function to create a Fold vs Fold plot for comparison of two conditions
create_fvf_plot = function(fold_change_data, comparison_x, comparison_y, color, title, x_label, y_label, save_path, height = 600, width = 800) {
  ggp = ggplot(fold_change_data, aes_string(x = comparison_x, y = comparison_y)) +
    geom_point(alpha = 0.6, color = color) +  
    my_theme +  # Apply custom theme
    labs(
      title = title, #Plot title
      x = x_label,  # X-axis label
      y = y_label   # Y-axis label
    ) +
    xlim(-15,15) + # Set the limits for the x-axis 
    theme(
      plot.title = element_text(hjust = 0.5),  # Position the title
      axis.title = element_text(size = 12),    # Set axis title size
      axis.text = element_text(size = 10)      # Set axis text size
    )
  save_plot(ggp, save_path, height = height, width = width)  # Save the plot to the specified path
  return(ggp)  # Return the plot object
}

# Function to perform PCA and plot results
plot_pca = function(data, groups, group_colors, title, save_path, height = 600, width = 600) {
  data_matrix = as.matrix(sapply(data, as.numeric)) # Convert the data matrix to numeric and transpose as samples should be rows for PCA
  data_matrix_transposed = t(data_matrix) 
  pca_result = prcomp(data_matrix_transposed) # Perform PCA using prcomp
  pca_coordinates = data.frame(pca_result$x)   # Extract PCA coordinates
  pca_coordinates$group = groups # Add the sample group information to PCA coordinates
  
  # Calculate the variance explained by each principal component 
  vars = apply(pca_result$x, 2, var)
  prop_x = round(vars["PC1"] / sum(vars), 4) * 100  # Percentage variance explained by PC1
  prop_y = round(vars["PC2"] / sum(vars), 4) * 100  # Percentage variance explained by PC2
  
  # Custom axis labels showing variance explained by PC1 and PC2
  x_axis_label = paste("PC1 (", prop_x, "%)", sep = "")
  y_axis_label = paste("PC2 (", prop_y, "%)", sep = "")
  
  # Create the PCA plot
  pca_plot = ggplot(pca_coordinates, aes(x = PC1, y = PC2, colour = group)) +
    geom_point() +
    scale_color_manual(values = group_colors) +  # Apply custom colors as per group
    labs(
      title = title, #Plot title
      x = x_axis_label, y = y_axis_label, colour = "Sample Group") +
    my_theme +
    theme(
      plot.title = element_text(hjust = 0.5),  # Position the title
      axis.title = element_text(size = 12),    # Set axis title size
      axis.text = element_text(size = 10)      # Set axis text size
    )
  save_plot(pca_plot, save_path, height = height, width = width)  # Save the plot to the specified path
  return(pca_plot)  # Return the plot object
}

plot_exp_density = function(data, title = "Expression Density Plot", fill_color = "lavenderblush", line_color = "midnightblue", columns = 3, x_axis_title = "Log10(Expression)", y_axis_title = "Density", font_size = 8, save_path, height = 600, width = 800) {
  melted_data = melt(data) # Melt the expression matrix as data should be in a matrix format

  # Plot the density plot
  ggp = ggplot(melted_data, aes(x = log10(value))) +
    geom_density(size = 0.5, alpha = 0.6, fill = fill_color, color = line_color) +  # Add fill color and adjust line size
    facet_wrap(~variable, ncol = columns) + #create faceted plot
    labs(title = title, x = x_axis_title, y = y_axis_title) +
    my_theme +
    theme(
      # Remove gridlines and background
      panel.grid = element_line(color = "gray90", size = 0.2),  # Subtle gridlines
      panel.background = element_rect(fill = "white"),  # Light background color
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),  # Border around the panels
      strip.text.x = element_text(size = font_size + 2, family = "sans", face = "bold", color = line_color, vjust = 1),
      plot.title = element_text(size = font_size + 4, family = "sans", face = "bold", color = "black", hjust = 0.5),
      axis.title = element_text(size = font_size, family = "sans", face = "bold"),
      axis.text = element_text(size = font_size, family = "sans"),
      panel.spacing = unit(1.5, "lines"),  # Increased space between facets
      legend.position = "none" 
    ) 
  save_plot(ggp, save_path, height = height, width = width)  # Save the plot to the specified path
  return(ggp)# Return the plot object
}

# Define a function to create a boxplot for specific gene expression
create_boxplot_for_gene = function(gene_name, em, sample_groups, group_order, group_colors, save_path, height = 600, width = 800) {
  gene_data = em[gene_name, ]  # Extract the gene's expression data for all samples
  gene_data = data.frame(t(gene_data))  # # Transpose the data so each row corresponds to one sample
  gene_data$sample_group = sample_groups # Add sample group information to the data frame
  gene_data$sample_group = factor(gene_data$sample_group, levels = group_order) # Set the order of sample groups
  names(gene_data) = c("expression", "sample_group")  # Rename the columns to "expression" and "sample_group" for clarity
  
  # Create a boxplot to show gene expression across groups, including median and interquartile range
  ggp = ggplot(gene_data, aes(x = sample_group, y = expression, fill = sample_group)) +
    geom_boxplot(size = 0.5, outlier.size = 0, alpha = 0.5, colour = "black") +  
    labs(title = paste("Boxplot of", gene_name, "Expression"), x = "Sample Group", y = "Expression") +  # Add labels
    scale_fill_manual(values = group_colors) +  # Apply custom colors to the boxes for each group
    my_theme  # Apply a custom theme
  
  save_plot(ggp, save_path, height = height, width = width) # Save the plot to the specified file pat
  return(ggp) # Return the plot object
}

# Define a function to create a violin plot for gene expression
create_violin_plot = function(gene_name, em, sample_groups, group_order, group_colors, save_path, height = 600, width = 800) {
  gene_data = em[gene_name, ]  # Extract the gene's expression data for all samples
  gene_data = data.frame(t(gene_data))  # # Transpose the data so each row corresponds to one sample
  gene_data$sample_group = sample_groups # Add sample group information to the data frame
  gene_data$sample_group = factor(gene_data$sample_group, levels = group_order) # Set the order of sample groups
  names(gene_data) = c("expression", "sample_group")  # Rename the columns to "expression" and "sample_group" for clarity
   
  # Create a violin plot using ggplot2 to show the distribution of data with density estimates on both sides of the central axis
    ggp = ggplot(gene_data, aes(x = sample_group, y = expression, fill = sample_group)) +
    geom_violin(width = 0.1, alpha = 0.5, trim = TRUE, show.legend = FALSE) +  # 
    labs(title = paste("Violin Plot of", gene_name, "Expression"), x = "Sample Group", y = "Expression") +
    scale_fill_manual(values = group_colors) +  # Apply custom group colors
    my_theme  # Apply custom theme
  
    save_plot(ggp, save_path, height = height, width = width) # Save the plot to the specified file pat
    return(ggp) # Return the plot object
}

# Define a function to create a jitter plot for gene expression
create_jitter_plot = function(gene_name, em, sample_groups, group_order, group_colors, save_path, height = 600, width = 800) {
  gene_data = em[gene_name, ]  # Extract the gene's expression data for all samples
  gene_data = data.frame(t(gene_data))  # # Transpose the data so each row corresponds to one sample
  gene_data$sample_group = sample_groups # Add sample group information to the data frame
  gene_data$sample_group = factor(gene_data$sample_group, levels = group_order) # Set the order of sample groups
  names(gene_data) = c("expression", "sample_group")  # Rename the columns to "expression" and "sample_group" for clarity
  
  # Create a jitter plot to show the distribution of individual points
  ggp = ggplot(gene_data, aes(x = sample_group, y = expression, color = sample_group)) +
    geom_jitter(width = 0.1, colour = "black") +  # Jitter points with slight horizontal displacement for visibility
    labs(title = paste("Jitter Plot of", gene_name, "Expression"), x = "Sample Group", y = "Expression") +
    scale_color_manual(values = group_colors) +  # Apply custom colors to the points
    my_theme  # Apply custom theme
  
  save_plot(ggp, save_path, height = height, width = width) # Save the plot to the specified file pat
  return(ggp) # Return the plot object
}

# Define a function to create a combined violin and jitter plot for gene expression
create_violin_jitter_plot = function(gene_name, em, sample_groups, group_order, group_colors, save_path, height = 600, width = 800) {
  gene_data = em[gene_name, ]  # Extract the gene's expression data for all samples
  gene_data = data.frame(t(gene_data))  # # Transpose the data so each row corresponds to one sample
  gene_data$sample_group = sample_groups # Add sample group information to the data frame
  gene_data$sample_group = factor(gene_data$sample_group, levels = group_order) # Set the order of sample groups
  names(gene_data) = c("expression", "sample_group")  # Rename the columns to "expression" and "sample_group" for clarity
  
  # Create a combined violin and jitter plot to visualize both the distribution and individual points
  ggp = ggplot(gene_data, aes(x = sample_group, y = expression, fill = sample_group)) +
    geom_violin(width = 0.1, alpha = 0.5) +  # Add the violin plot 
    geom_jitter(width = 0.1, color = "black") +  # Add jitter 
    labs(title = paste("Violin + Jitter Plot of", gene_name, "Expression"), x = "Sample Group", y = "Expression") +
    scale_fill_manual(values = group_colors) +  # Apply custom colors to the violins
    scale_color_manual(values = group_colors) +  # Apply custom colors to the jitter points
    my_theme  # Apply custom theme
  
  save_plot(ggp, save_path, height = height, width = width) # Save the plot to the specified file pat
  return(ggp) # Return the plot object
}

# Function to extract top genes with the lowest p-values from a master table
get_top_genes = function(master_table, top_no) {
  candidate_genes = rownames(master_table)[1:top_no] # Assuming the table is ordered in increasing p-values order
  return(candidate_genes) # Return the list of top N candidate genes
}

# Function to create boxplots for a set of genes to visualize their expression across sample groups
create_boxplots_for_genes = function(gene_names, em, ss, annot, sample_group_col, group_order, fill_colors, ncol, save_path, height = 800, width = 1000) {
  em = data.frame(scale(em)) # Subset the expression matrix for the chosen genes
  data = em[gene_names, ]
  data = data.frame(t(data))  # Transpose the data so that each row represents a sample
  data$sample_group = ss[[sample_group_col]]  # Add sample group information from the provided sample sheet
  data$sample_group = factor(data$sample_group, levels = group_order) # Set the order of the sample groups based on the specified group order
  gene_data_melted = melt(data, id.vars = "sample_group") # Melt the data into a long format suitable for ggplot
  colnames(gene_data_melted) = c("sample_group", "gene", "expression") # Rename columns for clarity
  gene_data_melted$gene_symbol = annot$SYMBOL[match(gene_data_melted$gene, rownames(annot))] # Annotate gene names with the provided annotation table
  
  # Create the boxplot with custom fill colors for each sample group
  ggp = ggplot(gene_data_melted, aes(x = sample_group, y = expression, fill = sample_group)) +
    geom_boxplot() +
    facet_wrap(~gene_symbol, ncol = ncol) +  # Facet by gene symbols, showing them in multiple columns
    labs(title = "Boxplot of Top Genes", x = "Sample Group", y = "Expression") +
    scale_fill_manual(values = fill_colors) +  # Apply custom colors to each sample group
    my_theme +  # Apply custom theme 
    theme(
      axis.text.x = element_text(size = 8, angle = 90, hjust = 1),  # Rotate x-axis labels for readability
      axis.text.y = element_text(size = 8),  # Adjust size of y-axis labels
      axis.title.x = element_text(size = 10),  # Adjust x-axis title size
      axis.title.y = element_text(size = 10),  # Adjust y-axis title size
      plot.title = element_text(size = 12, face = "bold"),  # Increase plot title size
      strip.text = element_text(size = 8),  # Adjust the size of facet labels (gene names)
      legend.title = element_text(size = 10),  # Adjust legend title size
      legend.text = element_text(size = 8)  # Adjust legend text size
    )
  save_plot(ggp, save_path, height = height, width = width) # Save the plot to the specified file pat
  return(ggp) # Return the plot object
}
 

# Function to create a clustered heatmap for a set of significant genes
create_clustered_heatmap = function(em, sig_genes, colour_palette, cluster_axis = "both", save_path, title = "Clustered Heatmap of Significant Genes", height = 800, width = 800) {
  em_sig = em[sig_genes, ] #Subset the expression matrix for the selected genes
  em_sig_scaled = scale(em_sig) #Scale the data to make genes comparable across samples
  hm.matrix = as.matrix(em_sig_scaled)  # Convert scaled data to matrix
  y.dist = Dist(hm.matrix, method = "spearman")  # Calculate pairwise distance between genes using Spearman correlation
  
  if (cluster_axis == "both" | cluster_axis == "y") { # Perform hierarchical clustering for genes (rows) if specified
    y.cluster = hclust(y.dist, method = "average")  
    y.dd = as.dendrogram(y.cluster)
    y.dd.reorder = reorder(y.dd, 0, FUN = "average") # Reorder the genes according to the dendrogram
    y.order = order.dendrogram(y.dd.reorder) # Get the new order of genes (rows) after clustering
    hm.matrix_clustered = hm.matrix[y.order, ] # Reorder the matrix according to the clustered gene order
  } 
  else {
    hm.matrix_clustered = hm.matrix # If no clustering on the y-axis specified, keep original gene order
  }
  
  if (cluster_axis == "both" | cluster_axis == "x") { # Perform hierarchical clustering on the x-axis (columns) if specified
    x.dist = Dist(t(hm.matrix_clustered), method = "spearman")  
    x.cluster = hclust(x.dist, method = "average")  
    x.dd = as.dendrogram(x.cluster)
    x.dd.reorder = reorder(x.dd, 0, FUN = "average") # Reorder the columns according to the dendrogram
    x.order = order.dendrogram(x.dd.reorder) # Get the new sample order (column order) after clustering
    hm.matrix_clustered = hm.matrix_clustered[, x.order] # Reorder the matrix according to the clustered gene order
  }
  
  hm.matrix_clustered = melt(hm.matrix_clustered) # Melt the clustered matrix into a long format for ggplot
  
  # Create the heatmap plot using ggplot2
  ggp = ggplot(hm.matrix_clustered, aes(x = Var2, y = Var1, fill = value)) + 
    geom_tile() +  # Use tiles to visualize the heatmap
    scale_fill_gradientn(colours = colorRampPalette(colour_palette)(100)) +  # Apply custom color gradient
    my_theme + 
    labs(title = title, x = "Samples", y = "Genes") +  # Add labels for the heatmap
    ylab("") +  # Remove y-axis label for a cleaner plot
    xlab("") +  # Remove x-axis label for a cleaner plot
    theme(
      axis.text.y = element_blank(),  # Hide y-axis text for a cleaner plot
      axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis labels for readability
      axis.ticks = element_blank(),  # Remove axis ticks for a cleaner look
      legend.title = element_blank(),  # Remove legend title
      legend.spacing.x = unit(0.25, 'cm')  # Adjust spacing for the legend
    )
  
  save_plot(ggp, save_path, height = height, width = width) # Save the plot to the specified file pat
  return(ggp) # Return the plot object
}

# Function to generate a rug plot to visualize sample groups
create_rug_plot = function(sample_column, group_colors, save_path, height = 600, width = 800) {
  sample_column = factor(sample_column, levels = names(group_colors)) 
  groups_data = as.matrix(as.numeric(sample_column))  # Convert to numeric factors
  groups_data = melt(groups_data)  # Melt the data into a long format suitable for ggplot
  groups_data$value = factor(groups_data$value, labels = names(group_colors)) # Set the factor labels in the melted data to the group names
  
  # Create the rug plot
  ggp = ggplot(groups_data, aes(y = Var2, x = value, fill = value)) + 
    geom_tile(linetype = "blank") + # Avoid tiles for cleaner look
    scale_fill_manual(values = group_colors) +  # Use the custom group colors
    my_theme + 
    labs(title = "Rug Plot for Sample Groups", x = "Groups", y = "Samples") + # Add labels
    theme(
      legend.position = "none",  
      legend.title = element_blank(),  # Remove legend title
      axis.text.x = element_text(hjust = 0.5),  # Rotate x-axis labels
      axis.text.y = element_blank(),  # Remove y-axis labels for a cleaner plot
      axis.ticks = element_blank()  # Remove axis ticks for a cleaner plot
    )
  
  save_plot(ggp, save_path, height = height, width = width) # Save the plot to the specified file pat
  return(ggp) # Return the plot object
}

# Function to generate and plot the STRING network for a list of candidate genes
plot_string_network = function(candidate_genes) {
  candidate_genes_table = data.frame(candidate_genes) # Convert the candidate genes into a data frame with a single column "gene"
  names(candidate_genes_table) = "gene"
  string_db = STRINGdb$new(version = "11.5", species = 9606, score_threshold = 200, network_type = "full", input_directory = "") # Load the STRING database for human species = 9606
  string_mapped = string_db$map(candidate_genes_table, "gene", removeUnmappedRows = TRUE) # Map the genes to the STRING database and filter out unmapped genes
  plot_str = string_db$plot_network(string_mapped) # Plot the STRING network for the candidate genes
  return(plot_str) # Return the plot object 
}

# Function to perform Over-representation analysis of various ontology types
do_ora = function(sig_genes, organism_db, sample_group_col, group_order, fill_colors, ncol, em, ss, annot, ontology_type, colour_palette, box_path, heatmap_path, heatmap_title) {
  results_list = list() # Create a list to store all results
  
  # Run ORA (Enrichment analysis) using enrichGO based on the ontology_type parameter
  ora_results = enrichGO(gene = sig_genes, OrgDb = organism_db, ont = ontology_type,  readable = T, pvalueCutoff = 0.05,  qvalueCutoff = 0.10)  # BP, MF, CC OR ALL
  results_list$ora_results = ora_results # Store ORA results in the list

  ora_results_table = data.frame( 
    gene_sets = ora_results$geneID,
    p.adj = ora_results$p.adjust
  )
  rownames(ora_results_table) = ora_results$Description   # Extract the ORA results table and rename the row names
  results_list$ora_results_table = ora_results_table  # Store the results table
  
  # Extract and split the most enriched gene set from the first row
  enriched_gene_set = as.character(ora_results_table[1, "gene_sets"])
  candidate_genes = unlist(strsplit(enriched_gene_set, "/"))
  results_list$candidate_genes = candidate_genes # Store the candidate genes for further analysis
  
  boxplot = create_boxplots_for_genes(candidate_genes, em, ss, annot, sample_group_col, group_order, fill_colors, ncol, save_path = box_path) # Create boxplots for candidate genes
  heatmap = create_clustered_heatmap(em, candidate_genes, colour_palette, save_path = heatmap_path, title = heatmap_title ) # Create a clustered heatmap for the candidate genes
  string_network = plot_string_network(candidate_genes) # Generate the STRING network for the candidate genes
  
  # Plot other General Enrichment Plots (bar, dot, go, and cnet plots)
  ggp_bar = barplot(ora_results, showCategory = 10)
  ggp_dot = dotplot(ora_results, showCategory = 10)
  ggp_goplot = goplot(ora_results, showCategory = 10)
  ggp_cnet = cnetplot(ora_results, categorySize = "pvalue")
  
  # Store the plots
  results_list$boxplot = boxplot
  results_list$heatmap = heatmap
  results_list$string_network = string_network
  results_list$ggp_bar = ggp_bar
  results_list$ggp_dot = ggp_dot
  results_list$ggp_goplot = ggp_goplot
  results_list$ggp_cnet = ggp_cnet
  
  return(results_list)  # Return all results
}

# Function to perform Gene Set Enrichment Analysis for different conditions
do_gsea = function(organism_db, condition, master_all, sample_group_col, group_order, fill_colors, ncol, em, ss, annot, ontology_type, colour_palette,box_path, heatmap_path, heatmap_title) {
  results_list = list() # Create a list to store all results
  gsea_input = master_all[[paste0("log2fold.", condition)]] # Prepare the GSEA input 
  names(gsea_input) = row.names(master_all)  # Assign gene names 
  gsea_input = na.omit(gsea_input)  # Remove NA values
  gsea_input = sort(gsea_input, decreasing = TRUE)  # Sort the genes in decreasing order
  
  # Run GSEA for the specified ontology type
  gseGO_result = gseGO(
    geneList = gsea_input,  
    ont = ontology_type,  
    keyType = "ENSEMBL",  
    nPerm = 10000,  
    minGSSize = 3,  
    maxGSSize = 800,  
    pvalueCutoff = 0.05,  
    verbose = TRUE,  
    OrgDb = organism_db,  
    pAdjustMethod = "none"
  )
  
  results_list$gsea_results = gseGO_result # Store GSEA results in the results list
  
  ridgeplot_result <- ridgeplot(gseGO_result) # Generate a ridge plot for the selected ontology
  
  # Extract the most enriched ontology for the selected condition
  gse_results_sorted = gseGO_result[order(gseGO_result$p.adjust), ]
  gse_results_sorted = gse_results_sorted[order(-gse_results_sorted$NES), ]
  most_enriched_ontology = gse_results_sorted[1, ]
  core_genes = unlist(strsplit(most_enriched_ontology$core_enrichment, "/"))  # Get core enrichment genes for the enriched ontology
  results_list$core_enrichment_genes = core_genes # Store the core enrichment genes in the results list
  
  boxplot_result = create_boxplots_for_genes(core_genes, em, ss, annot, sample_group_col, group_order, fill_colors, ncol, save_path = box_path)   # Generate a boxplot for the core enrichment genes
  heatmap = create_clustered_heatmap(em, core_genes, colour_palette, save_path = heatmap_path, title = heatmap_title)   # Generate a clustered heatmap for the core enrichment genes
  string_network = plot_string_network(core_genes)  # Generate STRING network plot for the core enrichment genes 
  
  # Store the plots
  results_list$ridgeplot = ridgeplot_result
  results_list$boxplot = boxplot_result
  results_list$heatmap = heatmap
  results_list$string_network = string_network
  
  return(results_list) # Return the results list
}

# Function to generates a boxplot for metagene expression based on a given list of genes
make_metagene_boxplot = function(gene_list, em_scaled, sample_groups, fill_colors) {
  metagene = data.frame(colMeans(em_scaled[gene_list, ]))  # Select rows corresponding to the gene list
  colnames(metagene) = "MetageneExpression"  # Name the resulting column "MetageneExpression"
  metagene$Group = sample_groups  # Add group information (sample_groups) to the metagene expression data.

  # Create the boxplot 
  ggp = ggplot(metagene, aes(x = Group, y = MetageneExpression, fill = Group)) + # Fill colors as per sample groups
    geom_boxplot() +  
    labs(title = "Metagene Boxplot",  # Add a title 
         y = "Metagene Expression", x = "Group") +  # Label the axes
    scale_fill_manual(values = fill_colors) +  # Custom colors for the groups
    my_theme  # Apply custom theme 

  return(ggp) # Return the ggplot object
}

# Function to run a Signature Analysis (ORA) and generates various plots for the provided signatures
run_sig_analysis = function(sig_genes, organism_db, sample_group_col, group_order, fill_colors, ncol, em, ss, annot, ontology_type, colour_palette, group_colors, save_path_prefix, box_path, heatmap_path, heatmap_title) {
  
  ora_results = do_ora(sig_genes, organism_db, sample_group_col, group_order, fill_colors, ncol, em, ss, annot, ontology_type, colour_palette, box_path, heatmap_path, heatmap_title) # Run ORA using the custom function defined earlier
  
  # Extract the results from the ORA output
  ora_results_table = ora_results$ora_results_table  # Table of ORA results
  candidate_genes = ora_results$candidate_genes  # Candidate genes identified by ORA
  boxplot = ora_results$boxplot  # Boxplot visualizing gene expression across groups
  ggp_bar = ora_results$ggp_bar  # Bar plot of ORA results
  ggp_dot = ora_results$ggp_dot  # Dot plot for gene expression data
  ggp_goplot = ora_results$ggp_goplot  # GO plot for gene ontology enrichment
  ggp_cnet = ora_results$ggp_cnet  # Co-expression network plot
  string_network = ora_results$string_network  # STRING network showing protein-protein interactions
  
  metagene_boxplot = make_metagene_boxplot(candidate_genes, em, sample_group_col, group_colors) # Generate a boxplot for the candidate genes' metagene expression.
  
  # Dynamically generate file paths to save the plots that will include the provided prefix and different suffixes for each plot type
  save_path_barplot = paste0(save_path_prefix, "_barplot.png")
  save_path_dotplot = paste0(save_path_prefix, "_dotplot.png")
  save_path_goplot = paste0(save_path_prefix, "_goplot.png")
  save_path_cnetplot = paste0(save_path_prefix, "_cnetplot.png")
  save_path_metagene_boxplot = paste0(save_path_prefix, "_metagene_boxplot.png")
  save_path_stringnet = paste0(save_path_prefix, "_stringnet.png")
  
  # Save the generated plots using a custom function `save_plot()`.
  # Each plot is saved to the corresponding file path with specified dimensions.
  save_plot(metagene_boxplot, save_path_metagene_boxplot, height = 600, width = 800)
  save_plot(string_network, save_path_stringnet, height = 1000, width = 1000)
  save_plot(ggp_bar, save_path_barplot, height = 600, width = 800)
  save_plot(ggp_dot, save_path_dotplot, height = 600, width = 800)
  save_plot(ggp_goplot, save_path_goplot, height = 800, width = 800)
  save_plot(ggp_cnet, save_path_cnetplot, height = 800, width = 800)
  
  # Return a list containing all the results for
  results_list = list(
    ora_results_table = ora_results_table,  # Table with ORA results
    candidate_genes = candidate_genes,  # List of candidate genes identified
    heatmap = heatmap,  # Heatmap plot
    boxplot = boxplot,  # Boxplot plot
    metagene_boxplot = metagene_boxplot,  # Metagene expression boxplot
    ggp_bar = ggp_bar,  # Bar plot
    ggp_dot = ggp_dot,  # Dot plot
    ggp_goplot = ggp_goplot,  # GO enrichment plot
    ggp_cnet = ggp_cnet,  # Co-expression network plot
    string_network = string_network  # STRING network plot
  )
  
  return(results_list) # Return the results list containing all outputs
}


