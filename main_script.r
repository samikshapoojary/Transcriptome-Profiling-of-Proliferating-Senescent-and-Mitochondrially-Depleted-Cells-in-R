# Loading the created functions from the script for analysis
source("omic_functions.R") 

# Load data from various CSV files containing experimental data for different conditions using the custom csv_reader function
em = csv_reader("EM.csv")  # Expression matrix
de_m_vs_p = csv_reader("DE_Senes_MtD_vs_Prolif.csv")  # Differential expression between Senes_MtD (m) and Prolif (p)
de_m_vs_s = csv_reader("DE_Senes_MtD_vs_Senes.csv")  # Differential expression between Senes_MtD and Senes(s)
de_s_vs_p = csv_reader("DE_Senes_vs_Prolif.csv")  # Differential expression between Senes and Prolif
annotations = csv_reader("Human_Background_GRCh38.p13.csv")  # Annotations for gene symbols and other information
ss = csv_reader("sample_sheet.csv")  # Sample groups

# Creating master tables for each pairwise comparisons by combining the gene expression data with differential expression and annotations
master_m_vs_p = create_master_table(em, de_m_vs_p, annotations, exclude_cols = c(7, 8, 9))
master_m_vs_s = create_master_table(em, de_m_vs_s, annotations, exclude_cols = c(4, 5, 6))
master_s_vs_p = create_master_table(em, de_s_vs_p, annotations, exclude_cols = c(10, 11, 12))

# Merging the master tables from pairwise comparisons and cleaning it
master1 = merge_de_tables(master_m_vs_p, master_m_vs_s, ".mvp", ".mvs")
master_temp = merge_de_tables(master1, master_s_vs_p, "", ".svp")

# Cleaning the master table
cols_to_keep = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 15, 16, 17, 18, 19, 20, 21, 22, 23, 32, 33, 34, 35, 36, 37, 44, 45, 46, 47, 49, 50,51)
master_all = master_temp[, cols_to_keep] # Keeping only relevant columns
colnames(master_all)[c(7, 8, 9, 22, 23, 24, 29, 30,31)] = c("Senes_MtD_1", "Senes_MtD_2", "Senes_MtD_3","log2fold.svp", "p.svp", "p.adj.svp", "mean.svp","mlog10p.svp", "sig.svp") # Renaming columns for ease
master_all = master_all[, c(4, 5, 6, 7, 8, 9, 16, 17, 18, 1, 2, 3, 10, 11, 12, 13, 14, 15, 19, 20, 21, 22, 23, 24, 29, 30, 31, 25, 26, 27, 28)] # Reordering to group similar information columns

# Ordering the rows by the minimum p-value to priorite genes that are statistically significant across all conditions
master_all = master_all[order(pmin(master_all$p.adj.mvp, master_all$p.adj.mvs, master_all$p.adj.svp)), ]

# Creating a mean column of DGE results as there are 3 columns of 3 comparisons
#master_all$mean_log2fold = rowMeans(master_all[, c("log2fold.mvp", "log2fold.mvs", "log2fold.svp")])
#master_all$mean_p.adj = rowMeans(master_all[, c("p.adj.mvp", "p.adj.mvs", "p.adj.svp")])
#master_all$mean_mean =rowMeans(master_all[, c("mean.mvp", "mean.mvs", "mean.svp")])

# Creating a table that uses the gene symbol as row names
em_symbols = create_symbols_table(master_table = master_all,select_cols = c(1:9,28),  rownames_col = "SYMBOL")
em_symbols = em_symbols[,-10]

# Scaling the gene expression data to standardize values
em_symbols_scaled = data.frame(scale(em_symbols))
em_scaled = data.frame(scale(em))

# Creating tables of significantly differentially expressed genes in each comparison based on strict thresholds and obtaining their names as vectors
#sig_m_vs_p_table = create_sig_genes_table(master_m_vs_p, 0.001, 2)
#sig_m_vs_s_table = create_sig_genes_table(master_m_vs_s, 0.001, 2)
#sig_s_vs_p_table = create_sig_genes_table(master_s_vs_p, 0.001, 2)
#sig_m_vs_p = row.names(sig_m_vs_p_table)
#sig_m_vs_s = row.names(sig_m_vs_s_table)
#sig_s_vs_p = row.names(sig_s_vs_p_table)

# Subset genes that are significant in all three comparisons
sig_all_table = subset(master_all, abs(log2fold.svp) > 2 & p.adj.svp < 0.001 & abs(log2fold.mvs) > 2 & p.adj.mvs < 0.001 &abs(log2fold.mvp) > 2 & p.adj.mvp < 0.001)
sig_all = row.names(sig_all_table)
sig_all_symbols = sig_all_table$SYMBOL #getting symbols of genes for future usage

# Subset genes that are significant in any of the comparisons
#sig_in_either_table = subset(master_all, abs(log2fold.svp) > 2 & p.adj.svp < 0.001 |abs(log2fold.mvs) > 2 & p.adj.mvs < 0.001 | abs(log2fold.mvp) > 2 & p.adj.mvp < 0.001)
#sig_in_either = row.names(sig_in_either_table)
#sig_in_either_symbols = sig_in_either_table$SYMBOL #adding a symbols column for future usage

# Defining various parameters to be used in future plots, etc. 
group_colors = c("Prolif" = "royalblue", "Senes" = "darkviolet", "Senes_MtD" = "darkgreen") #colors for the experimental groups for visualization
sample_groups = ss$SAMPLE_GROUP #extracting sample group information from the sample sheet
group_order = c("Prolif", "Senes", "Senes_MtD") #order in which the experimental groups should appear in plots
sample_group_col = "SAMPLE_GROUP"  #column name for sample groups in the sample sheet
fill_colors = group_colors #fill colors for the group plots
ncol = 5  #number of columns for faceting in plots
colour_palette = c("plum3", "white", "royalblue", "turquoise", "black") #color palette for the heatmaps
organism_db = "org.Hs.eg.db"  # Defining the Human gene database database for Gene Ontology and GSEA as dealing with human fibroblasts

# Generating volcano plots for each pairwise comparison and all comparisons to visualize differentially expressed genes based on p-value and fold change
volcano_m_vs_p = plot_volcano(master_m_vs_p, log2fold, p.adj, 0.001, 2, "volcano_m_vs_p.png", title = "Volcano Plot of DGE between Senes_MtD v/s Prolif")
volcano_m_vs_s = plot_volcano(master_m_vs_s, log2fold, p.adj, 0.001, 2, "volcano_m_vs_s.png", title = "Volcano Plot of DGE between Senes_MtD v/s Senes")
volcano_s_vs_p = plot_volcano(master_s_vs_p, log2fold, p.adj, 0.001, 2, "volcano_s_vs_p.png", title = "Volcano Plot of DGE between Senes v/s Prolif")
#volcano_all = plot_volcano(master_all, log2fold = mean_log2fold, p.adj = mean_p.adj, 0.001, 2, "volcano_all.png", title = "Volcano Plot of mean DGE between Senes_MtD, Senes and Prolif")

# Generating MA plots for each pairwise comparison and all comparisons to show the relationship between mean expression and fold change
ma_m_vs_p = plot_ma(master_m_vs_p, log2fold, p.adj, mean, 0.001, 2, "ma_m_vs_p.png", title = "MA Plot of DGE between Senes_MtD v/s Prolif")
ma_m_vs_s = plot_ma(master_m_vs_s, log2fold, p.adj, mean, 0.001, 2, "ma_m_vs_s.png", title = "MA Plot of DGE between Senes_MtD v/s Senes")
ma_s_vs_p = plot_ma(master_s_vs_p, log2fold, p.adj, mean, 0.001, 2, "ma_s_vs_p.png", title = "MA Plot of DGE between Senes v/s Prolif")
#ma_all = plot_ma(master_all, log2fold = mean_log2fold, p.adj = mean_p.adj, mean = mean_mean, 0.001, 2, "volcano_all.png", title = "MA Plot of mean DGE between Senes_MtD, Senes and Prolif")

#Plotting a expression density plot to check the distribution across all samples
exp_density = plot_exp_density(em, columns = 3, save_path = "exp_density.png")

# Creating a list that contains the significant genes from each pairwise comparison to make the venn diagram
#venn_data = list("M_vs_P" = sig_m_vs_p, "M_vs_S" = sig_m_vs_s, "S_vs_P" = sig_s_vs_p)
#colors = c("M_vs_P" = "lightblue",  "M_vs_S" = "lavender", "S_vs_P" = "honeydew") #assigning colors for each comparison for better visualization
#venn_plot = plot(euler(venn_data, shape = "ellipse"), quantities = TRUE, fill = colors) #plotting the venn diagram using euler 
#save_plot(venn_plot, "venn_diagram.png", height = 600, width = 600) #saving the Venn diagram as a PNG file

# Noting the total numbers of significant genes in each of the pairwise comparisons and the total of entire dataset for calculating overlaps
#group_mvs = nrow(sig_m_vs_p_table)   #genes significant in M_vs_P
#group_mvp = nrow(sig_m_vs_s_table)   #genes significant in M_vs_S
#group_svp = nrow(sig_s_vs_p_table)   #genes significant in S_vs_P
#total = nrow(master_all) #total number of genes in the entire dataset

# Calculating the overlap using a hypergeometric test to check if the observed overlap is greater than expected by chance
#p_value_mvs_mvp = hypergeo_test(sig_m_vs_p, sig_m_vs_s, group_mvs, group_mvp, total) #p-value for genes that are in both M_vs_P and M_vs_S
#p_value_mvs_svp = hypergeo_test(sig_m_vs_p, sig_s_vs_p, group_mvs, group_svp, total) #p-value for genes that are in both M_vs_P and S_vs_P
#p_value_mvp_svp = hypergeo_test(sig_m_vs_s, sig_s_vs_p, group_mvp, group_svp, total) #p-value for genes that are in both M_vs_S and S_vs_P

# Combining the results of the overlap p-values into a data frame for easier display and interpretation
#p_values_overlap = data.frame(comparison = c("M_vs_P and M_vs_S", "M_vs_P and S_vs_P", "M_vs_S and S_vs_P"), p_value = c(p_value_mvs_mvp, p_value_mvs_svp, p_value_mvp_svp))
  
# Creating a data frame containing fold change values for each pairwise comparison
#fold_change_data = data.frame(m_vs_p = master_all$log2fold.mvp, m_vs_s = master_all$log2fold.mvs, s_vs_p = master_all$log2fold.svp)

# Creating fold vs fold (FVF) plots to compare fold changes between different comparisons to visualize how genes with significant fold changes in one comparison behave in another comparison
#fold change (M_vs_P) vs fold change (M_vs_S), colored in dark blue
#mvp_mvs_fvf = create_fvf_plot(fold_change_data, comparison_x = "m_vs_p", comparison_y = "m_vs_s", color = "darkblue", title = "Fold Change vs Fold Change (M vs P vs M vs S)", x_label = "Fold Change (M vs P)", y_label = "Fold Change (M vs S)", "mvp_mvs_fvf.png")
#fold change (M_vs_P) vs fold change (S_vs_P), colored in maroon
#mvp_svp_fvf = create_fvf_plot(fold_change_data, comparison_x = "m_vs_p", comparison_y = "s_vs_p", color = "maroon", title = "Fold Change vs Fold Change (M vs P vs S vs P)", x_label = "Fold Change (M vs P)", y_label = "Fold Change (S vs P)", "mvp_svp_fvf.png")
#fold change (M_vs_S) vs fold change (S_vs_P), colored in seagreen
#mvs_svp_fvf = create_fvf_plot(fold_change_data, comparison_x = "m_vs_s", comparison_y = "s_vs_p", color = "seagreen", title = "Fold Change (M vs S) vs Fold Change (S vs P)", x_label = "Fold Change (M vs S)", y_label = "Fold Change (S vs P)", "mvs_svp_fvf.png")

# Perform correlation tests to assess whether the fold changes in different comparisons are related 

cor_test_m_vs_p_vs_m_vs_s = cor.test(master_all$log2fold.mvp, master_all$log2fold.mvs) #correlation between fold change in M_vs_P and M_vs_S
cor_test_m_vs_p_vs_s_vs_p = cor.test(master_all$log2fold.mvp, master_all$log2fold.svp) #correlation between fold change in M_vs_P and S_vs_P
cor_test_m_vs_s_vs_p = cor.test(master_all$log2fold.mvs, master_all$log2fold.svp) #correlation between fold change in M_vs_S and S_vs_P

# Doing Principal Component Analysis (PCA) to identify patterns and clusters in the data corresponding to biological groups
# Generating PCA plot for original expression matrix to seehow the samples group based on their raw expression data
#pca_plot_all = plot_pca(em_symbols, ss$SAMPLE_GROUP, group_colors, title = "PCA Plot of Expressions", "pca_plot_all.png")

# Generating PCA plot for the scaled expression data i.e. normalized for better comparison
pca_plot_scaled = plot_pca(em_symbols_scaled, ss$SAMPLE_GROUP, group_colors, title = "PCA Plot of Scaled Expressions", "pca_plot_scaled.png")

# Creating boxplot, violin plot, jitter plot, combined violin + jitter plot for selected gene
#gene_name = "TP53" #important in cell cycle regulation
#boxplot_tp53 = create_boxplot_for_gene(gene_name, em_symbols, sample_groups=ss$SAMPLE_GROUP, group_order, group_colors, "boxplot_tp53.png")
#violin_plot_tp53 = create_violin_plot(gene_name, em_symbols, sample_groups, group_order, group_colors, "violin_plot_tp53.png")
#jitter_plot_tp53 = create_jitter_plot(gene_name,em_symbols, sample_groups, group_order, group_colors, "jitter_plot_tp53.png")
#violin_jitter_plot_tp53 = create_violin_jitter_plot(gene_name, em_symbols,sample_groups, group_order, group_colors, "violin_jitter_plot_tp53.png")

# Identifying the top 10 most significant genes in each comparison group to be further analyzed for patterns and biological relevance
candidate_genes_mvp = get_top_genes(master_m_vs_p, top_no=10)  # Top 10 genes in M_vs_P comparison
candidate_genes_mvs = get_top_genes(master_m_vs_s, top_no=10)  # Top 10 genes in M_vs_S comparison
candidate_genes_svp = get_top_genes(master_s_vs_p, top_no=10)  # Top 10 genes in S_vs_P comparison
candidate_genes_all = get_top_genes(master_all, top_no=10)      # Top 10 genes across all comparisons

# Creating faceted boxplots for the selected genes to visualize their expression levels across samples
top10_boxplot_mvp = create_boxplots_for_genes(candidate_genes_mvp, em, ss, annotations, sample_group_col, group_order, fill_colors, ncol, "top10_boxplot_mvp.png")
top10_boxplot_mvs = create_boxplots_for_genes(candidate_genes_mvs, em, ss, annotations, sample_group_col, group_order, fill_colors, ncol, "top10_boxplot_mvs.png")
top10_boxplot_svp = create_boxplots_for_genes(candidate_genes_svp, em, ss, annotations, sample_group_col, group_order, fill_colors, ncol, "top10_boxplot_svp.png")
top10_boxplot_all = create_boxplots_for_genes(candidate_genes_all, em, ss, annotations, sample_group_col, group_order, fill_colors, ncol, "top10_boxplot_all.png")

# Creating a clustered heatmap to visualize the expression patterns of significant genes across samples
heatmap_sig_all = create_clustered_heatmap(em, sig_all, colour_palette, save_path = "heatmap_sig_all.png")

# Creating a rug plot to identify sample groups in heatmap
groups_rugplot = create_rug_plot(ss$SAMPLE_GROUP, group_colors, "groups_rugplot.png")

# Subsetting the upregulated genes across all comparisons based on the log2 fold change and adjusted p-value thresholds
upregulated_genes = subset(master_all, log2fold.svp > 0 & p.adj.svp < 0.001 & log2fold.mvs > 0 & p.adj.mvs < 0.001 & log2fold.mvp > 0 & p.adj.mvp < 0.001)
upregulated_gene_names = row.names(upregulated_genes)  #extracting gene names for upregulated genes
upregulated_gene_names_symbols= upregulated_genes$SYMBOL #extracting symbols for upregulated genes

# Subsetting the downregulated genes across all comparisons based on the log2 fold change and adjusted p-value thresholds
downregulated_genes = subset(master_all, log2fold.svp < 0 & p.adj.svp < 0.001 & log2fold.mvs < 0 & p.adj.mvs < 0.001 & log2fold.mvp < 0 & p.adj.mvp < 0.001)
downregulated_gene_names = row.names(downregulated_genes) #extracting gene names for downregulated genes
downregulated_gene_names_symbols = downregulated_genes$SYMBOL #extracting gene names for downregulated genes

# Converting the significant genes from SYMBOL to ENTREZID using the 'bitr' function for performing Gene Ontology (GO) enrichment analysis (ORA) and Gene Set Enrichment Analysis (GSEA), which use ENTREZ IDs
# The bitr function retrieves the ENTREZ ID based on the gene SYMBOL, using the org.Hs.eg.db database for human genes
sig_genes_entrez = bitr(sig_all_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) 
upregulated_genes_entrez = bitr(upregulated_gene_names_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
downregulated_genes_entrez = bitr(downregulated_gene_names_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Performing Gene Ontology (GO) enrichment analysis for the significant genes for each ontology category
# ORA tests the overrepresentation of gene sets in specific biological processes (BP), molecular functions (MF), or cellular components (CC) to identify which biological functions or cellular processes are enriched among the significant genes
# Perform GO enrichment analysis for the upregulated and downregulated genes using the same ontology categories

# For BP ontology
ora_results_sig_genes_BP = do_ora(sig_genes = sig_genes_entrez$ENTREZID, organism_db, sample_group_col, group_order, fill_colors, ncol, em = em_symbols, ss, annotations, ontology_type = "BP",colour_palette, "ora_sig_BP_boxplot.png", "ora_sig_BP_heatmap.png", "Clustered Heatmap of Enriched BP of Significant Genes")
ora_results_upregulated_genes_BP = do_ora(sig_genes = upregulated_genes_entrez$ENTREZID, organism_db, sample_group_col, group_order, fill_colors, ncol, em = em_symbols, ss, annotations, ontology_type = "BP",colour_palette, "ora_upreg_BP_boxplot.png", "ora_upreg_BP_heatmap.png", "Clustered Heatmap of Enriched BP of Upregulated Genes")
ora_results_downregulated_genes_BP = do_ora(sig_genes = downregulated_genes_entrez$ENTREZID, organism_db, sample_group_col, group_order, fill_colors, ncol, em = em_symbols, ss, annotations, ontology_type = "BP",colour_palette, "ora_downreg_BP_boxplot.png", "ora_downreg_BP_heatmap.png", "Clustered Heatmap of Enriched BP of Downregulated Genes")

# Saving the barplot resuts
save_plot(ora_results_sig_genes_BP$ggp_bar, "ora_sig_BP_bar.png", height = 600, width = 800)
save_plot(ora_results_upregulated_genes_BP$ggp_bar, "ora_upregulated_BP_bar.png", height = 600, width = 800)
save_plot(ora_results_downregulated_genes_BP$ggp_bar, "ora_downregulated_BP_bar.png", height = 600, width = 800)

#ora_results_sig_genes_MF = do_ora(sig_genes = sig_genes_entrez$ENTREZID, organism_db, sample_group_col, group_order, fill_colors, ncol, em = em_symbols, ss, annotations, ontology_type = "MF",colour_palette, "ora_sig_MF_boxplot.png", "ora_sig_MF_heatmap.png")
#ora_results_sig_genes_CC = do_ora(sig_genes = sig_genes_entrez$ENTREZID, organism_db, sample_group_col, group_order, fill_colors, ncol, em = em_symbols, ss, annotations, ontology_type = "CC",colour_palette, "ora_sig_CC_boxplot.png", "ora_sig_CC_heatmap.png")
#ora_results_sig_genes_ALL = do_ora(sig_genes = sig_genes_entrez$ENTREZID, organism_db, sample_group_col, group_order, fill_colors, ncol, em = em_symbols, ss, annotations, ontology_type = "ALL",colour_palette, "ora_sig_ALL_boxplot.png", "ora_sig_ALL_heatmap.png")
#ora_results_upregulated_genes_MF = do_ora(sig_genes = upregulated_genes_entrez$ENTREZID, organism_db, sample_group_col, group_order, fill_colors, ncol, em = em_symbols, ss, annotations, ontology_type = "MF",colour_palette, "ora_upreg_MF_boxplot.png", "ora_upreg_MF_heatmap.png")
#ora_results_upregulated_genes_CC = do_ora(sig_genes = upregulated_genes_entrez$ENTREZID, organism_db, sample_group_col, group_order, fill_colors, ncol, em = em_symbols, ss, annotations, ontology_type = "CC",colour_palette, "ora_upreg_CC_boxplot.png", "ora_upreg_CC_heatmap.png")
#ora_results_upregulated_genes_ALL = do_ora(sig_genes = upregulated_genes_entrez$ENTREZID, organism_db, sample_group_col, group_order, fill_colors, ncol, em = em_symbols, ss, annotations, ontology_type = "ALL",colour_palette, "ora_upreg_ALL_boxplot.png", "ora_upreg_ALL_heatmap.png")
#ora_results_downregulated_genes_MF = do_ora(sig_genes = downregulated_genes_entrez$ENTREZID, organism_db, sample_group_col, group_order, fill_colors, ncol, em = em_symbols, ss, annotations, ontology_type = "MF",colour_palette, "ora_downreg_MF_boxplot.png", "ora_downreg_MF_heatmap.png")
#ora_results_downregulated_genes_CC = do_ora(sig_genes = downregulated_genes_entrez$ENTREZID, organism_db, sample_group_col, group_order, fill_colors, ncol, em = em_symbols, ss, annotations, ontology_type = "CC",colour_palette, "ora_downreg_CC_boxplot.png", "ora_downreg_CC_heatmap.png")
#ora_results_downregulated_genes_ALL = do_ora(sig_genes = downregulated_genes_entrez$ENTREZID, organism_db, sample_group_col, group_order, fill_colors, ncol, em = em_symbols, ss, annotations, ontology_type = "ALL",colour_palette, "ora_downreg_ALL_boxplot.png", "ora_downreg_ALL_heatmap.png")

# Performing Gene Set Enrichment Analysis (GSEA) for each condition (M_vs_P, M_vs_S, S_vs_P) using the same ontology categories
# GSEA assesses whether a predefined set of genes shows statistically significant differences between two biological states
#gsea_results_mvp_BP = do_gsea(organism_db = "org.Hs.eg.db", condition = "mvp", master_all, sample_group_col, group_order, fill_colors, ncol, em, ss, annotations, ontology_type = "BP",colour_palette, "gsea_mvp_BP_boxplot.png", "gsea_mvp_BP_heatmap.png", "Clustered Heatmap of Enriched BP in Senes_MtD and Prolif")
#gsea_results_mvs_BP = do_gsea(organism_db = "org.Hs.eg.db", condition = "mvs", master_all, sample_group_col, group_order, fill_colors, ncol, em, ss, annotations, ontology_type = "BP",colour_palette, "gsea_mvs_BP_boxplot.png", "gsea_mvs_BP_heatmap.png", "Clustered Heatmap of Enriched BP in Senes_MtD and Senes")
#gsea_results_svp_BP = do_gsea(organism_db = "org.Hs.eg.db", condition = "svp", master_all, sample_group_col, group_order, fill_colors, ncol, em, ss, annotations, ontology_type = "BP",colour_palette, "gsea_svp_BP_boxplot.png", "gsea_svp_BP_heatmap.png", "Clustered Heatmap of Enriched BP in Prolif and Senes")
#gsea_results_mvp_ALL = do_gsea(organism_db = "org.Hs.eg.db", condition = "mvp", master_all, sample_group_col, group_order, fill_colors, ncol, em, ss, annotations, ontology_type = "ALL",colour_palette,"gsea_mvp_ALL_boxplot.png", "gsea_mvp_ALL_heatmap.png")
#gsea_results_mvp_MF = do_gsea(organism_db = "org.Hs.eg.db", condition = "mvp", master_all, sample_group_col, group_order, fill_colors, ncol, em, ss, annotations, ontology_type = "MF",colour_palette, "gsea_mvp_MF_boxplot.png", "gsea_mvp_MF_heatmap.png")
#gsea_results_mvp_CC = do_gsea(organism_db = "org.Hs.eg.db", condition = "mvp", master_all, sample_group_col, group_order, fill_colors, ncol, em, ss, annotations, ontology_type = "CC",colour_palette, "gsea_mvp_CC_boxplot.png", "gsea_mvp_CC_heatmap.png")
#gsea_results_mvs_ALL = do_gsea(organism_db = "org.Hs.eg.db", condition = "mvs", master_all, sample_group_col, group_order, fill_colors, ncol, em, ss, annotations, ontology_type = "ALL",colour_palette, "gsea_mvs_ALL_boxplot.png", "gsea_mvs_ALL_heatmap.png")
#gsea_results_mvs_MF = do_gsea(organism_db = "org.Hs.eg.db", condition = "mvs", master_all, sample_group_col, group_order, fill_colors, ncol, em, ss, annotations, ontology_type = "MF",colour_palette, "gsea_mvs_MF_boxplot.png", "gsea_mvs_MF_heatmap.png")
#gsea_results_mvs_CC = do_gsea(organism_db = "org.Hs.eg.db", condition = "mvs", master_all, sample_group_col, group_order, fill_colors, ncol, em, ss, annotations, ontology_type = "CC",colour_palette,"gsea_mvs_CC_boxplot.png", "gsea_mvs_CC_heatmap.png")
#gsea_results_svp_ALL = do_gsea(organism_db = "org.Hs.eg.db", condition = "svp", master_all, sample_group_col, group_order, fill_colors, ncol, em, ss, annotations, ontology_type = "ALL",colour_palette, "gsea_svp_ALL_boxplot.png", "gsea_svp_ALL_heatmap.png")
#gsea_results_svp_MF = do_gsea(organism_db = "org.Hs.eg.db", condition = "svp", master_all, sample_group_col, group_order, fill_colors, ncol, em, ss, annotations, ontology_type = "MF",colour_palette, "gsea_svp_MF_boxplot.png", "gsea_svp_MF_heatmap.png")
#gsea_results_svp_CC = do_gsea(organism_db = "org.Hs.eg.db", condition = "svp", master_all, sample_group_col, group_order, fill_colors, ncol, em, ss, annotations, ontology_type = "CC",colour_palette, "gsea_svp_CC_boxplot.png", "gsea_svp_ALL_heatmap.png")

# Create a clustered heatmap of significant genes in either of the comparisons
#heatmap_sig_in_either = create_clustered_heatmap(em, sig_in_either, colour_palette, save_path = "heatmap_sig_in_either.png")

# Signature 1 - genes that are upregulated in 'svp' but downregulated in 'mvp
signature_1_table = subset(master_all, sig.svp == TRUE & log2fold.svp > 0 & sig.mvp == TRUE & log2fold.mvp < 0)
signature_1 = row.names(signature_1_table) #extracting gene names
signature_1_entrez = bitr(signature_1, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) #converting gene names to entrezid for ORA

# Signature 2 - genes that are downregulated in 'mvs' but upregulated in 'svp'
signature_2_table = subset(master_all, sig.mvs == TRUE & log2fold.mvs < 0 & sig.svp == TRUE & log2fold.svp > 0)
signature_2 = row.names(signature_2_table)  #extracting gene names
signature_2_entrez = bitr(signature_2, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) #converting gene names to entrezid for ORA

# Signature 3 - genes that are upregulated in 'mvs' but downregulated in 'svp'
signature_3_table = subset(master_all, sig.mvs == TRUE & log2fold.mvs > 0 & sig.svp == TRUE & log2fold.svp < 0)
signature_3 = row.names(signature_3_table) #extracting gene names
signature_3_entrez = bitr(signature_3, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) #converting gene names to entrezid for ORA

# Signature 4 - genes that are downregulated in 'mvs' but upregulated in 'mvp'
signature_4_table = subset(master_all, sig.mvs == TRUE & log2fold.mvs < 0 & sig.mvp == TRUE & log2fold.mvp > 0)
signature_4 = row.names(signature_4_table) #extracting gene names
signature_4_entrez = bitr(signature_4, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) #converting gene names to entrezid for ORA

#Biological Process (BP) Ontology enrichment analysis (ORA) for all Signatures using custom Signature analysis function
ora_results_sign_1_bp = run_sig_analysis(sig_genes = signature_1_entrez$ENTREZID, organism_db, sample_group_col,group_order, fill_colors, ncol, em, ss, annotations, ontology_type = "BP", colour_palette, group_colors, save_path_prefix = "sign1_bp", box_path="ora_sign1_BP_boxplot.png", heatmap_path="ora_sign1_BP_heatmap.png", heatmap_title="Clustered Heatmap of Signature 1")
ora_results_sign_2_bp = run_sig_analysis(sig_genes = signature_2_entrez$ENTREZID, organism_db, sample_group_col = "SAMPLE_GROUP", group_order, fill_colors, ncol = 2, em, ss, annotations, ontology_type = "BP", colour_palette, group_colors, save_path_prefix = "sign2_bp", box_path="ora_sign2_BP_boxplot.png", heatmap_path="ora_sign2_BP_heatmap.png", heatmap_title="Clustered Heatmap of Signature 2")
ora_results_sign_3_bp = run_sig_analysis(sig_genes = signature_3_entrez$ENTREZID, organism_db, sample_group_col = "SAMPLE_GROUP", group_order, fill_colors, ncol = 2, em, ss, annotations, ontology_type = "BP", colour_palette, group_colors, save_path_prefix = "sign3_bp", box_path="ora_sign3_BP_boxplot.png", heatmap_path="ora_sign3_BP_heatmap.png", heatmap_title="Clustered Heatmap of Signature 3")
ora_results_sign_4_bp = run_sig_analysis(sig_genes = signature_4_entrez$ENTREZID, organism_db, sample_group_col = "SAMPLE_GROUP", group_order, fill_colors, ncol = 2, em, ss, annotations, ontology_type = "BP", colour_palette, group_colors, save_path_prefix = "sign4_bp", box_path="ora_sign4_BP_boxplot.png", heatmap_path="ora_sign4_BP_heatmap.png", heatmap_title="Clustered Heatmap of Signature 4")

# For MF, CC, ALL ontology type of signature 1
#ora_results_sign_1_mf = run_sig_analysis(sig_genes = signature_1, organism_db, sample_group_col = "SAMPLE_GROUP", group_order, fill_colors, ncol = 2, em, ss, annotations, ontology_type = "MF", colour_palette, group_colors, save_path_prefix = "sign1_mf")
#ora_results_sign_1_cc = run_sig_analysis(sig_genes = signature_1, organism_db, sample_group_col = "SAMPLE_GROUP", group_order, fill_colors, ncol = 2, em, ss, annotations, ontology_type = "CC", colour_palette, group_colors, save_path_prefix = "sign1_cc")
#ora_results_sign_1_all = run_sig_analysis(sig_genes = signature_1, organism_db, sample_group_col = "SAMPLE_GROUP", group_order, fill_colors, ncol = 2, em, ss, annotations, ontology_type = "ALL", colour_palette, group_colors, save_path_prefix = "sign1_all")

# For MF, CC, ALL ontology type of signature 2
#ora_results_sign_2_mf = run_sig_analysis(sig_genes = signature_2, organism_db, sample_group_col = "SAMPLE_GROUP",group_order, fill_colors, ncol = 2, em, ss, annotations, ontology_type = "MF", colour_palette, group_colors, save_path_prefix = "sign2_mf")
#ora_results_sign_2_cc = run_sig_analysis(sig_genes = signature_2, organism_db, sample_group_col = "SAMPLE_GROUP",group_order, fill_colors, ncol = 2, em, ss, annotations, ontology_type = "CC", colour_palette, group_colors, save_path_prefix = "sign2_cc")
#ora_results_sign_2_all = run_sig_analysis(sig_genes = signature_2, organism_db, sample_group_col = "SAMPLE_GROUP",group_order, fill_colors, ncol = 2, em, ss, annotations, ontology_type = "ALL", colour_palette, group_colors, save_path_prefix = "sign2_all")

# For MF, CC, ALL ontology type of signature 3
#ora_results_sign_3_mf = run_sig_analysis(sig_genes = signature_3, organism_db, sample_group_col = "SAMPLE_GROUP",group_order, fill_colors, ncol = 2, em, ss, annotations, ontology_type = "MF", colour_palette, group_colors, save_path_prefix = "sign3_mf")
#ora_results_sign_3_cc = run_sig_analysis(sig_genes = signature_3, organism_db, sample_group_col = "SAMPLE_GROUP",group_order, fill_colors, ncol = 2, em, ss, annotations, ontology_type = "CC", colour_palette, group_colors, save_path_prefix = "sign3_cc")
#ora_results_sign_3_all = run_sig_analysis(sig_genes = signature_3, organism_db, sample_group_col = "SAMPLE_GROUP",group_order, fill_colors, ncol = 2, em, ss, annotations, ontology_type = "ALL", colour_palette, group_colors, save_path_prefix = "sign3_all")

# For MF, CC, ALL ontology type of signature 4
#ora_results_sign_4_mf = run_sig_analysis(sig_genes = signature_4, organism_db, sample_group_col = "SAMPLE_GROUP",group_order, fill_colors, ncol = 2, em, ss, annotations, ontology_type = "MF", colour_palette, group_colors, save_path_prefix = "sign4_mf")
#ora_results_sign_4_cc = run_sig_analysis(sig_genes = signature_4, organism_db, sample_group_col = "SAMPLE_GROUP",group_order, fill_colors, ncol = 2, em, ss, annotations, ontology_type = "CC", colour_palette, group_colors, save_path_prefix = "sign4_cc")
#ora_results_sign_4_all = run_sig_analysis(sig_genes = signature_4, organism_db, sample_group_col = "SAMPLE_GROUP",group_order, fill_colors, ncol = 2, em, ss, annotations, ontology_type = "ALL", colour_palette, group_colors, save_path_prefix = "sign4_all")

# For MF, CC, ALL ontology type of signature 5
#ora_results_sign_5_mf = run_sig_analysis(sig_genes = signature_5, organism_db, sample_group_col = "SAMPLE_GROUP",group_order, fill_colors, ncol = 2, em, ss, annotations, ontology_type = "MF", colour_palette, group_colors, save_path_prefix = "sign5_mf")
#ora_results_sign_5_cc = run_sig_analysis(sig_genes = signature_5, organism_db, sample_group_col = "SAMPLE_GROUP",group_order, fill_colors, ncol = 2, em, ss, annotations, ontology_type = "CC", colour_palette, group_colors, save_path_prefix = "sign5_cc")
#ora_results_sign_5_all = run_sig_analysis(sig_genes = signature_5, organism_db, sample_group_col = "SAMPLE_GROUP",group_order, fill_colors, ncol = 2, em, ss, annotations, ontology_type = "ALL", colour_palette, group_colors, save_path_prefix = "sign5_all")

# For MF, CC, ALL ontology type of signature 6
#ora_results_sign_6_mf = run_sig_analysis(sig_genes = signature_6, organism_db, sample_group_col = "SAMPLE_GROUP", group_order, fill_colors, ncol = 2, em, ss, annotations, ontology_type = "MF", colour_palette, group_colors, save_path_prefix = "sign6_mf")
#ora_results_sign_6_cc = run_sig_analysis(sig_genes = signature_6, organism_db, sample_group_col = "SAMPLE_GROUP", group_order, fill_colors, ncol = 2, em, ss, annotations, ontology_type = "CC", colour_palette, group_colors, save_path_prefix = "sign6_cc")
#ora_results_sign_6_all = run_sig_analysis(sig_genes = signature_6, organism_db, sample_group_col = "SAMPLE_GROUP",  group_order, fill_colors, ncol = 2, em, ss, annotations, ontology_type = "ALL", colour_palette, group_colors, save_path_prefix = "sign6_all")

