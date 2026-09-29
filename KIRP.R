setwd("D:/KIRP")
library(dplyr) 
library(tidyverse)
library(ggplot2) 
library(RColorBrewer) 

library(devtools)
library(ggplotify) 

library(ggcorrplot) 

library(corrplot)
install.packages("BiocManager")

library(org.Hs.eg.db)
library(tibble)

library(readr)
# Read in data
metadat <- read.csv("TCGA.KIRP.csv")
exp_data <- read.csv("HiSeqV2.csv")
surv_data <- read.delim("https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/survival%2FKIRP_survival.txt")
library(dplyr)
# Remove normal tissue samples
metadat_filtered <- metadat[metadat$sample_type != "Solid Tissue Normal", ]
# Keep only male patients
male_metadat <- metadat_filtered %>%
  filter(gender == "MALE")

#3. Compute ssGSEA LOY scores (for male samples only)
LOY_set <- c("RPS4Y1","ZFY","CDHL3","TBL1Y","USP9Y",
             "DDX3Y","UTY","TMSB4Y","NLGN4Y","HSFY2","KDM5D",
             "EIF1AY","RBMY1A1","PRY2")

# Custom ssGSEA function
ssgsea <- function(X, gene_sets, alpha = 0.25, scale = T, norm = F, single = T) {
  row_names <- rownames(X)
  num_genes <- nrow(X)
  gene_sets <- lapply(gene_sets, function(genes) {which(row_names %in% genes)})
  R <- matrixStats::colRanks(X, preserveShape = T, ties.method = 'average')
  es <- apply(R, 2, function(R_col) {
    gene_ranks <- order(R_col, decreasing = TRUE)
    es_sample <- sapply(gene_sets, function(gene_set_idx) {
      indicator_pos <- gene_ranks %in% gene_set_idx
      indicator_neg <- !indicator_pos
      rank_alpha <- (R_col[gene_ranks] * indicator_pos) ^ alpha
      step_cdf_pos <- cumsum(rank_alpha) / sum(rank_alpha)
      step_cdf_neg <- cumsum(indicator_neg) / sum(indicator_neg)
      step_cdf_diff <- step_cdf_pos - step_cdf_neg
      if (scale) step_cdf_diff <- step_cdf_diff / num_genes
      if (single) sum(step_cdf_diff) else step_cdf_diff[which.max(abs(step_cdf_diff))]
    })
    unlist(es_sample)
  })
  if (length(gene_sets) == 1) es <- matrix(es, nrow = 1)
  if (norm) es <- es / diff(range(es))
  rownames(es) <- names(gene_sets)
  colnames(es) <- colnames(X)
  return(es)
}
library(tibble)
# Prepare expression data
exp_data <- exp_data %>% column_to_rownames(var = "sample")

# Keep only male tumor samples in expression data
male_sample_ids <- male_metadat$sampleID
exp_data_male <- exp_data[, colnames(exp_data) %in% gsub("-", ".", male_sample_ids)]

# Run ssGSEA
LOY_ssgsea_male <- ssgsea(as.matrix(exp_data_male), LOY_set)
LOY_ssgsea_male <- na.omit(LOY_ssgsea_male)
LOY_means_male <- colMeans(LOY_ssgsea_male)

LOY_df_male <- data.frame(
  sample = names(LOY_means_male),
  LOY_score = as.numeric(LOY_means_male)
)
LOY_df_male$sample_dash <- gsub("\\.", "-", LOY_df_male$sample)
#4. Merge LOY scores with metadata
male_combined <- left_join(LOY_df_male, male_metadat, by = c("sample_dash" = "sampleID"))
#5. Classify into Y_high / Y_low (using male-specific median)
LOY_median_male <- median(male_combined$LOY_score, na.rm = TRUE)
male_combined <- male_combined %>%
  mutate(LOY_group = ifelse(LOY_score > LOY_median_male, "Y_high", "Y_low"))
# Wilcoxon test
p_value <- wilcox.test(LOY_score ~ LOY_group, data = male_combined)$p.value

# Adjust p-value (e.g., using FDR method)
adj_p_value <- p.adjust(p_value, method = "fdr")  # or "bonferroni"

library(ggplot2)
library(ggpubr)

ggplot(male_combined, aes(x = LOY_group, y = LOY_score, fill = LOY_group)) +
  geom_boxplot() +
  labs(
    title = "LOY Score by Group (Male Tumor Samples Only)",
    x = "LOY Group",
    y = "Mean LOY Score"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("Y_low" = "#1f77b4", "Y_high" = "#ff7f0e")) +
  stat_compare_means(
    method = "wilcox.test",     # Or "t.test" if you prefer
    label = "p.format",         # Shows the formatted p-value
    label.y = max(male_combined$LOY_score, na.rm = TRUE) * 1.05
  )
library(ggplot2)
library(ggpubr)

ggplot(male_combined, aes(x = LOY_group, y = LOY_score, fill = LOY_group)) +
  geom_boxplot() +
  labs(
    title = "LOY Score by Group (Male Tumor Samples Only)",
    x = "LOY Group",
    y = "Mean LOY Score"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("Y_low" = "#1f77b4", "Y_high" = "#ff7f0e")) +
  annotate(
    "text",
    x = 1.5,  # position between groups
    y = max(male_combined$LOY_score, na.rm = TRUE) * 1.05,
    label = paste0("Adjusted p = ", signif(adj_p_value, 3)),
    size = 5
  )
install.packages("survival")
library(survival)
#7. KM Survival Analysis (for male tumor samples only)
# Merge with survival data
male_surv <- left_join(male_combined, surv_data, by = c("sample_dash" = "sample"))

# Remove samples with missing survival info
male_surv_clean <- male_surv %>%
  filter(!is.na(OS.time), !is.na(OS), !is.na(LOY_group))


# Fit survival model
surv_fit_male <- survfit(Surv(OS.time, OS) ~ LOY_group, data = male_surv_clean)

install.packages("survminer")  # if not already installed
library(survminer)

ggsurvplot(
  surv_fit_male,
  data = male_surv_clean,
  pval = TRUE,
  risk.table = TRUE,
  legend.title = "LOY Group",
  legend.labs = c("Y_high", "Y_low"),
  palette = c("#ff7f0e", "#1f77b4")
)
ggsurvplot(
  surv_fit_male,
  data = male_surv_clean,
  pval = TRUE,
  risk.table = TRUE,
  legend.title = "LOY Group",
  legend.labs = c("Y_high", "Y_low"),
  palette = c("#ff7f0e", "#1f77b4"),
  title = "Overall Survival by LOY Group (Male Patients Only)"
)
# Make sure survminer and survival are loaded
library(survival)
library(survminer)

# Step 1: Prepare the survival dataframe (if not already)
male_surv <- left_join(LOY_df_male, male_metadat, by = c("sample_dash" = "sampleID"))
male_surv <- left_join(male_surv, surv_data, by = c("sample_dash" = "sample"))

# Remove missing survival data
male_surv_clean <- male_surv %>%
  filter(!is.na(OS.time), !is.na(OS), !is.na(LOY_score))

# Step 2: Determine optimal cutoff for LOY_score
cutpoint <- surv_cutpoint(
  male_surv_clean,
  time = "OS.time",
  event = "OS",
  variables = "LOY_score",
  minprop = 0.1  # ensures each group has at least 10% of samples
)

summary(cutpoint)  # Check the optimal cutoff value

# Step 3: Categorize samples based on optimal cutoff
male_surv_cut <- surv_categorize(cutpoint)



# Step 4: Fit survival model using optimal cutoff
surv_fit_optimal <- survfit(Surv(OS.time, OS) ~ LOY_score, data = male_surv_cut)

# Step 5: Plot KM curve
ggsurvplot(
  surv_fit_optimal,
  data = male_surv_cut,
  pval = TRUE,
  risk.table = TRUE,
  legend.title = "LOY Group (Optimal Cut)",
  legend.labs = c("Y_low", "Y_high"),
  palette = c("#1f77b4", "#ff7f0e"),
  title = "Overall Survival by Optimal LOY Cutoff"
)
# Check the available survival columns
colnames(surv_data)

# Make sure RFS data exists
# Example: using "PFI.time" and "PFI" if present

# Merge LOY score with metadata and survival
male_surv_rfs <- left_join(LOY_df_male, male_metadat, by = c("sample_dash" = "sampleID"))
male_surv_rfs <- left_join(male_surv_rfs, surv_data, by = c("sample_dash" = "sample"))

# Filter out missing RFS values
male_surv_rfs_clean <- male_surv_rfs %>%
  filter(!is.na(PFI.time), !is.na(PFI), !is.na(LOY_score))  # replace with actual column names if different

# Determine optimal LOY cutoff for RFS
cutpoint_rfs <- surv_cutpoint(
  male_surv_rfs_clean,
  time = "PFI.time",
  event = "PFI",
  variables = "LOY_score",
  minprop = 0.1
)

summary(cutpoint_rfs)

# Categorize using optimal cutoff
male_rfs_cut <- surv_categorize(cutpoint_rfs)

# Fit RFS model
rfs_fit <- survfit(Surv(PFI.time, PFI) ~ LOY_score, data = male_rfs_cut)

# Plot KM curve for RFS
ggsurvplot(
  rfs_fit,
  data = male_rfs_cut,
  pval = TRUE,
  risk.table = TRUE,
  legend.title = "LOY Group (RFS Optimal)",
  legend.labs = c("Y_low", "Y_high"),
  palette = c("#1f77b4", "#ff7f0e"),
  title = "Progression Free Interval by Optimal LOY Cutoff"
)


# Load required libraries
# Load libraries
library(limma)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(readr)

# Step 1: Read expression data (CSV format: genes in rows, samples in columns)
exp_data <- read.csv("HiSeqV2.csv", header = TRUE, stringsAsFactors = FALSE)
exp_data <- exp_data %>% column_to_rownames(var = "sample")

# Step 2: Load metadata and fix sample names
male_rfs_cut$sample_dash <- male_surv_rfs_clean$sample_dash
male_rfs_cut <- as_tibble(male_rfs_cut) %>%
  dplyr::mutate(
    sample = gsub("-", ".", sample_dash),
    LOY_score = factor(tolower(LOY_score), levels = c("high", "low"))
  ) %>%
  dplyr::select(sample, LOY_score) %>%
  dplyr::filter(sample %in% colnames(exp_data))

# Step 3: Subset and reorder expression matrix to match metadata
exp_matrix_filtered <- exp_data[, male_rfs_cut$sample]
stopifnot(all(colnames(exp_matrix_filtered) == male_rfs_cut$sample))

# Step 4: Create design matrix and run limma
design <- model.matrix(~ LOY_score, data = male_rfs_cut)
fit <- lmFit(exp_matrix_filtered, design)
fit <- eBayes(fit)

# Step 5: Extract DEG results
deg_results <- topTable(fit, coef = "LOY_scorelow", number = Inf, adjust.method = "fdr") %>%
  rownames_to_column(var = "Gene")

# Step 6: Annotate significance
deg_results <- deg_results %>%
  dplyr::mutate(Significance = case_when(
    adj.P.Val < 0.05 & logFC > 1 ~ "Upregulated",
    adj.P.Val < 0.05 & logFC < -1 ~ "Downregulated",
    TRUE ~ "Not Significant"
  ))

# Step 7: Select top genes for labeling
top_up <- deg_results %>%
  dplyr::filter(Significance == "Upregulated") %>%
  dplyr::arrange(adj.P.Val) %>%
  dplyr::slice(1:10)

top_down <- deg_results %>%
  dplyr::filter(Significance == "Downregulated") %>%
  dplyr::arrange(adj.P.Val) %>%
  dplyr::slice(1:10)

top_genes <- dplyr::bind_rows(top_up, top_down)

# Step 8: Volcano Plot
ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c(
    "Upregulated" = "#D62728",
    "Downregulated" = "#1F77B4",
    "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(data = top_genes, aes(label = Gene), size = 3.5, max.overlaps = 50) +
  labs(
    title = "Volcano Plot of Differentially Expressed Genes",
    x = "Log2 Fold Change (Low vs High LOY Score)",
    y = "-Log10 Adjusted P-Value",
    color = "Significance"
  ) +
  theme_minimal()

#For Fariha only
ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c(
    "Upregulated" = "#D62728",
    "Downregulated" = "#1F77B4",
    "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(data = top_genes, aes(label = Gene), size = 3.5, max.overlaps = 50) +
  labs(
    title = "Volcano Plot of Differentially Expressed Genes",
    subtitle = "Comparison: Y_low vs Y_high (baseline = Y_high)",
    x = "Log2 Fold Change (Y_low vs Y_high)",
    y = "-Log10 Adjusted P-Value",
    color = "Significance"
  ) +
  theme_minimal()


# Step 9: Save results
write.csv(deg_results, "deg_results_loy_group.csv", row.names = FALSE)



# Load libraries
library(ggplot2)
library(ggrepel)
library(dplyr)
library(readr)

# Step 1: Read DEG results
deg <- read_csv("deg_results loy group.csv")
#for fariha only
deg <- read_csv("/Users/farihahasan/Desktop/KIRP/deg_results_loy_group.csv")

# Step 2: Check and clean column names (optional, but helpful)
colnames(deg) <- make.names(colnames(deg))

# Expected columns: logFC, adj.P.Val, Gene or rownames
# If gene names are in rownames, convert them to a column
if (!"Gene" %in% colnames(deg)) {
  deg <- deg %>% rownames_to_column(var = "Gene")
}

# Step 3: Add Significance category
deg <- deg %>%
  mutate(Significance = case_when(
    adj.P.Val < 0.05 & logFC > 1  ~ "Upregulated",
    adj.P.Val < 0.05 & logFC < -1 ~ "Downregulated",
    TRUE                          ~ "Not Significant"
  ))

# Step 4: Pick top 10 up and down regulated genes for labeling
top_up <- deg %>%
  filter(Significance == "Upregulated") %>%
  arrange(adj.P.Val) %>%
  slice(1:10)

top_down <- deg %>%
  filter(Significance == "Downregulated") %>%
  arrange(adj.P.Val) %>%
  slice(1:10)

top_genes <- bind_rows(top_up, top_down)

# Step 5: Volcano plot
ggplot(deg, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c(
    "Upregulated" = "#D62728",
    "Downregulated" = "#1F77B4",
    "Not Significant" = "grey"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(
    data = top_genes,
    aes(label = Gene),
    size = 3.5,
    max.overlaps = 50
  ) +
  labs(
    title = "Volcano Plot of DEGs (LOY Group)",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-Value",
    color = "Significance"
  ) +
  theme_minimal()
# Install if not already installed
install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "ReactomePA", "msigdbr", "enrichplot", "DOSE"))

# Install if not already installed
install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "ReactomePA", "msigdbr", "enrichplot", "DOSE"))

# Load libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(msigdbr)
library(dplyr)
library(tibble)
library(enrichplot)
library(ggplot2)


# Filter significant genes (adj.P.Val < 0.05)
sig_genes <- deg %>%
  filter(adj.P.Val < 0.05)


# Map gene symbols to Entrez IDs
gene_df <- bitr(sig_genes$Gene, fromType = "SYMBOL",
                toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Merge with logFC
sig_genes_entrez <- inner_join(sig_genes, gene_df, by = c("Gene" = "SYMBOL"))

# Create ranked gene list for GSEA
gene_list <- sig_genes_entrez$logFC
names(gene_list) <- sig_genes_entrez$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

# Get vector of significant Entrez IDs (for over-representation)
entrez_ids <- unique(sig_genes_entrez$ENTREZID)
kegg_res <- enrichKEGG(gene = entrez_ids,
                       organism = "hsa",
                       pvalueCutoff = 0.05)

dotplot(kegg_res, showCategory = 20) + ggtitle("KEGG Pathways")
go_bp <- enrichGO(gene = entrez_ids,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pvalueCutoff = 0.05,
                  readable = TRUE)

dotplot(go_bp, showCategory = 20) + ggtitle("GO Biological Process")
go_mf <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, ont = "MF", pvalueCutoff = 0.05, readable = TRUE)
go_cc <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, ont = "CC", pvalueCutoff = 0.05, readable = TRUE)
dotplot(go_mf, showCategory = 20) + ggtitle("GO Metabolic pathway")
dotplot(go_cc, showCategory = 20) + ggtitle("GO CC")

reactome_res <- enrichPathway(gene = entrez_ids,
                              organism = "human",
                              pvalueCutoff = 0.05,
                              readable = TRUE)

dotplot(reactome_res, showCategory = 20) + ggtitle("Reactome Pathways")
# Get MSigDB Hallmark gene sets
msig_hallmark <- msigdbr(species = "Homo sapiens", category = "H")

# Convert to list format
hallmark_list <- split(msig_hallmark$entrez_gene, msig_hallmark$gs_name)
hallmark_list <- split(msig_hallmark$entrez_gene, msig_hallmark$gs_name)
term2gene_hallmark <- msig_hallmark %>%
  dplyr::select(gs_name, entrez_gene)
gsea_hallmark <- GSEA(geneList = gene_list,
                      TERM2GENE = term2gene_hallmark,
                      pvalueCutoff = 0.05)
dotplot(gsea_hallmark, showCategory = 20, title = "GSEA - Hallmark Pathways") +
  theme_minimal()
gseaplot2(gsea_hallmark, geneSetID = gsea_hallmark@result$ID[1], title = "Top Hallmark GSEA Pathway")

####

library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

# Step 1: Extract result table
gsea_df <- gsea_hallmark@result

# Step 2: Extract core enriched ENTREZ IDs for each pathway
gsea_df$GeneSymbol <- sapply(gsea_df$core_enrichment, function(x) {
  entrez <- unlist(strsplit(x, "/"))
  
  # Convert to symbols
  symbol_map <- bitr(entrez,
                     fromType = "ENTREZID",
                     toType = "SYMBOL",
                     OrgDb = org.Hs.eg.db)
  
  symbols <- symbol_map$SYMBOL[match(entrez, symbol_map$ENTREZID)]
  paste(symbols, collapse = ", ")
})

# Step 3: View genes in top enriched pathways
# You can view or export:
View(gsea_df[, c("ID", "Description", "GeneSymbol")])
write.csv(gsea_df[, c("ID", "Description", "GeneSymbol")], "hallmark_gsea_genes.csv", row.names = FALSE)

getwd()




library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

# Extract result table
gsea_df <- gsea_hallmark@result

# Convert ENTREZ to gene symbols safely
gsea_df$GeneSymbol <- sapply(gsea_df$core_enrichment, function(x) {
  if (is.na(x) || x == "") return(NA)
  
  entrez <- unlist(strsplit(x, "/"))
  
  symbol_map <- bitr(entrez,
                     fromType = "ENTREZID",
                     toType = "SYMBOL",
                     OrgDb = org.Hs.eg.db)
  
  # Preserve input order
  symbols <- symbol_map$SYMBOL[match(entrez, symbol_map$ENTREZID)]
  paste(symbols[!is.na(symbols)], collapse = ", ")
})


# Step 1: Extract result table from your hallmark GSEA object
# Define the enriched hallmark pathways you're interested in
selected_hallmarks <- c(
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", 
  "HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_MITOTIC_SPINDLE"
)

# Filter and extract unique gene symbols
selected_genes <- gsea_df %>%
  filter(ID %in% selected_hallmarks) %>%
  pull(GeneSymbol) %>%
  strsplit(",\\s*") %>%     # handles spaces after commas
  unlist() %>%
  unique()

  
  # View results
  head(selected_genes)
  length(selected_genes) # Number of unique genes
#prepare expression matrix for these genes
# Subset male expression matrix to selected genes
expr_selected <- exp_data_male[rownames(exp_data_male) %in% selected_genes, ]

# Match columns to sample order in clinical data
expr_selected <- expr_selected[, match(gsub("-", ".", male_surv_clean$sample_dash), colnames(expr_selected))]

#univar cox reg for each gene (pfi)
library(survival)
library(broom)

# Run univariate Cox for each gene
cox_results <- apply(expr_selected, 1, function(gene_expr) {
  fit <- coxph(Surv(PFI.time, PFI) ~ gene_expr, data = male_surv_clean)
  tidy(fit)
})

# Convert to data frame
cox_df <- do.call(rbind, cox_results)
cox_df$gene <- rownames(expr_selected)

# Filter significant genes
cox_sig <- cox_df %>% filter(p.value < 0.05)
print(cox_sig)

# Save significant Cox results to CSV
write.csv(cox_sig, "significant_univariate_cox_hallmark_genes.csv", row.names = FALSE)


#lasso cox on 24 sig genes
library(glmnet)
library(survival)

# STEP 1: Subset and transpose the expression matrix for LASSO input
selected_genes <- cox_sig$gene
X <- t(expr_selected[selected_genes, ])
X <- scale(X)  # z-score normalization across samples

# STEP 2: Create survival object for PFI
y <- Surv(male_surv_clean$PFI.time, male_surv_clean$PFI)

# Remove samples with non-positive survival time
valid_idx <- which(male_surv_clean$PFI.time > 0)

# Subset both X and survival data
X <- X[valid_idx, ]
y <- Surv(male_surv_clean$PFI.time[valid_idx], male_surv_clean$PFI[valid_idx])


# STEP 3: Fit LASSO Cox
set.seed(123)
cv_fit <- cv.glmnet(X, y, family = "cox", alpha = 1, nfolds = 10)

# STEP 4: Plot CV error curve
plot(cv_fit)
title("LASSO Cross-Validation Curve (PFI)", line = 2.5)

# STEP 5: Extract genes with non-zero coefficients at optimal lambda
lasso_coefs <- coef(cv_fit, s = "lambda.min")
selected_lasso_genes <- rownames(lasso_coefs)[lasso_coefs[, 1] != 0]
selected_lasso_genes <- selected_lasso_genes[selected_lasso_genes != "(Intercept)"]

cat("Selected genes by LASSO Cox:\n")
print(selected_lasso_genes)


#build multivar cox model with selected genes
# Subset and scale expression matrix
X_lasso <- t(expr_selected[selected_lasso_genes, ])
X_lasso <- scale(X_lasso)

# Ensure matching samples
X_lasso <- X_lasso[valid_idx, ]  # use same valid_idx from before

# Create final data frame
cox_data <- data.frame(
  PFI.time = male_surv_clean$PFI.time[valid_idx],
  PFI = male_surv_clean$PFI[valid_idx],
  X_lasso
)

# Fit multivariate Cox model
library(survival)
fit_multi <- coxph(Surv(PFI.time, PFI) ~ ., data = cox_data)
summary(fit_multi)
# Save multivariate Cox model summary
cox_summary <- summary(fit_multi)
write.csv(as.data.frame(cox_summary$coefficients), "lasso_multivariate_cox_results.csv")





# --- Required libraries ---
library(survival)
library(boot)
library(dplyr)
library(ggplot2)

# --- STEP 1: Ensure data is aligned ---
# Assuming you already have:
# - cox_data: the multivariate expression + survival matrix (n x genes + 2 survival columns)
# - male_surv_final: the matched survival dataframe (same rows as cox_data)
# Get sample IDs from expression matrix used in cox_data
samples_in_cox <- rownames(cox_data)
samples_in_cox_dash <- gsub("\\.", "-", samples_in_cox)

# Now filter and match
male_surv_selected <- male_surv_clean %>%
  filter(sample_dash %in% samples_in_cox_dash)

# Identify unmatched ones (optional safety)
extra_in_surv <- setdiff(male_surv_selected$sample_dash, samples_in_cox_dash)



male_surv_final <- male_surv_selected %>%
  filter(!sample_dash %in% extra_in_surv) %>%
  arrange(match(sample_dash, samples_in_cox_dash))


stopifnot(nrow(male_surv_final) == nrow(cox_data))
stopifnot(all(gsub("-", ".", male_surv_final$sample_dash) == rownames(cox_data)))


# Overwrite survival info with cleaned values
cox_data$PFI.time <- male_surv_final$PFI.time
cox_data$PFI <- male_surv_final$PFI

# --- STEP 2: Fit multivariate Cox model ---
fit_multi <- coxph(Surv(PFI.time, PFI) ~ ., data = cox_data)
summary(fit_multi)

# --- STEP 3: Bootstrap validation ---
# Function to extract Cox coefficients
boot_cox <- function(data, indices) {
  boot_data <- data[indices, ]
  fit <- coxph(Surv(PFI.time, PFI) ~ ., data = boot_data)
  return(coef(fit))
}

# Run bootstrap
set.seed(42)
boot_results <- boot(data = cox_data, statistic = boot_cox, R = 1000)

# Compute percentile confidence intervals
boot.ci.results <- lapply(1:length(boot_results$t0), function(i) {
  boot.ci(boot_results, type = "perc", index = i)
})

# Compute standard errors
boot_se <- apply(boot_results$t, 2, sd)

# Combine results into summary table
boot_summary <- data.frame(
  Gene = names(coef(fit_multi)),
  Original = boot_results$t0,
  SE = boot_se,
  CI_Lower = sapply(boot.ci.results, function(x) x$percent[4]),
  CI_Upper = sapply(boot.ci.results, function(x) x$percent[5])
)

# Save summary
write.csv(boot_summary, "lasso_bootstrap_ci_results.csv", row.names = FALSE)
print(boot_summary)

#As my bootstrap Result shows low not that good confidence interval for Capg
#I decided to do penalized cox with cross validation
library(survival)
library(glmnet)

# Prepare input matrices
# Remove survival columns and keep predictors
predictors <- cox_data[, !(names(cox_data) %in% c("PFI.time", "PFI"))]

# Make sure categorical variables are factors (glmnet requires numeric/dummy)
# Convert factors to dummy variables using model.matrix
x <- model.matrix(~ ., data = predictors)[, -1]  # remove intercept column

# Survival object
y <- Surv(cox_data$PFI.time, cox_data$PFI)

set.seed(42)  # for reproducibility

# Fit elastic net Cox model (alpha=0.5 as example; alpha=1 for LASSO, 0 for ridge)
fit_cv <- cv.glmnet(x, y, family = "cox", alpha = 0.5, nfolds = 5)

# Plot cross-validation curve
plot(fit_cv)

# Lambda minimizing cross-validation error
best_lambda <- fit_cv$lambda.min
cat("Best lambda:", best_lambda, "\n")

# Coefficients at best lambda (sparse solution)
coef_best <- coef(fit_cv, s = best_lambda)
selected_vars <- rownames(coef_best)[as.numeric(coef_best) != 0]

cat("Selected variables:\n")
print(selected_vars)

# Coefficients values
print(coef_best[coef_best != 0])

# Optional: Predict risk scores for your samples
risk_scores <- predict(fit_cv, newx = x, s = best_lambda, type = "link")

# Add risk score to your dataframe
cox_data$risk_score <- as.numeric(risk_scores)

library(tibble)
library(dplyr)  # make sure this is loaded

# Convert sparse matrix to data frame
lasso_coef_df <- as.matrix(coef_best) %>%
  as.data.frame()

# Name the coefficient column
colnames(lasso_coef_df) <- "Coefficient"

# Add gene names
lasso_coef_df$Gene <- rownames(lasso_coef_df)

# Use base R for safety (if dplyr select fails)
lasso_coef_df <- lasso_coef_df[, c("Gene", "Coefficient")]
lasso_coef_df <- lasso_coef_df[lasso_coef_df$Coefficient != 0, ]

print(lasso_coef_df)
write.csv(lasso_coef_df, "lasso_lambda_min_coefficients.csv", row.names = FALSE)


# You can now proceed to Step 2 permutation testing or survival analysis with risk groups
#permutation test with c index
library(survival)
library(dplyr)

# Your data and risk scores from Step 1:
# cox_data with columns: PFI.time, PFI, risk_score

# Calculate observed C-index
cox_obj <- Surv(cox_data$PFI.time, cox_data$PFI)

library(survcomp)  # for concordance index calculation
observed_cindex <- concordance.index(cox_data$risk_score, cox_data$PFI.time, cox_data$PFI)$c.index
cat("Observed C-index:", observed_cindex, "\n")

# Number of permutations
n_perm <- 1000
set.seed(42)

# Initialize vector to store permuted c-indexes
perm_cindex <- numeric(n_perm)

for (i in 1:n_perm) {
  # Permute risk scores randomly
  perm_risk <- sample(cox_data$risk_score)
  
  # Calculate C-index for permuted risk scores
  perm_cindex[i] <- concordance.index(perm_risk, cox_data$PFI.time, cox_data$PFI)$c.index
}

# Calculate p-value: proportion of permuted c-index >= observed
p_val_perm <- mean(perm_cindex >= observed_cindex)
cat("Permutation test p-value:", p_val_perm, "\n")

# Plot histogram of permuted C-index with observed value line
library(ggplot2)
perm_df <- data.frame(cindex = perm_cindex)

ggplot(perm_df, aes(x = cindex)) +
  geom_histogram(binwidth = 0.01, fill = "lightblue", color = "black") +
  geom_vline(xintercept = observed_cindex, color = "red", linewidth = 1.5)+
  labs(
    title = "Permutation Test for C-index",
    subtitle = paste0("Observed C-index (red line) = ", round(observed_cindex, 3),
                      "\nPermutation p-value = ", signif(p_val_perm, 3)),
    x = "C-index",
    y = "Frequency"
  ) +
  theme_minimal()


#clinical adjustment Cox code
# Step 1: Combine risk score with selected clinical variables
# Step 0: Ensure matching order by sample IDs
# Step 0: Ensure samples are matched and ordered the same
# Assume these are your inputs:
# cox_data: data.frame with rows = samples, columns = PFI.time, PFI, risk_score, rownames = sample IDs with dots
# male_metadat (or clinical_data): clinical metadata for male tumor samples, with sample IDs matching cox_data but possibly with dashes
# clinical_vars: vector of clinical variable names to adjust for, e.g.:
clinical_vars <- c("age_at_initial_pathologic_diagnosis", "clinical_stage")

# Step 0: Prepare sample IDs for matching (replace dashes with dots in clinical metadata)
male_metadat$sample_id_dots <- gsub("-", ".", male_metadat$sampleID)

# Step 1: Subset clinical metadata to only those samples in cox_data
clinical_sub <- male_metadat %>%
  filter(sample_id_dots %in% rownames(cox_data)) %>%
  arrange(match(sample_id_dots, rownames(cox_data)))

# Step 2: Check sample order matches perfectly
stopifnot(all(clinical_sub$sample_id_dots == rownames(cox_data)))

# Step 3: Combine risk score, survival, and clinical variables into one dataframe
final_df <- cbind(
  cox_data[, c("PFI.time", "PFI", "risk_score")],
  clinical_sub[, clinical_vars]
)

# Step 4: (Optional) Convert 'stage' or other categorical variables to factor with correct levels
if ("stage" %in% clinical_vars) {
  final_df$stage <- factor(final_df$stage,
                           levels = c("Stage I", "Stage II", "Stage III", "Stage IV"))
}

# Step 5: Fit adjusted Cox proportional hazards model
library(survival)
fit_adj <- coxph(Surv(PFI.time, PFI) ~ risk_score + age_at_initial_pathologic_diagnosis + clinical_stage, data = final_df)


# Step 6: View results
summary(fit_adj)

# Roc Curve

# Assuming your Cox model is stored in 'fit_adj'
# and your data is 'final_df'

# Get linear predictor (risk score) from the fitted model
final_df$risk_pred <- predict(fit_adj, newdata = final_df, type = "lp")
library(timeROC)

# Define time points for ROC (in days): e.g., 1, 3, 5 years
times <- c(365, 1095, 1825)
roc_obj <- timeROC(T = final_df$PFI.time,
                   delta = final_df$PFI,
                   marker = final_df$risk_pred,
                   cause = 1,
                   weighting = "marginal",
                   times = times,
                   iid = TRUE)
print(roc_obj$AUC)
pplot(roc_obj, time = 365, col = "red", lwd = 2, title = FALSE)
plot(roc_obj, time = 1095, col = "blue", lwd = 2, add = TRUE)
plot(roc_obj, time = 1825, col = "green", lwd = 2, add = TRUE)
title("Time-dependent ROC curve for multivariate Cox model")

# Add legend
legend("bottomright",
       legend = c(
         paste0("1 year AUC = ", round(roc_obj$AUC[1], 3)),
         paste0("3 years AUC = ", round(roc_obj$AUC[2], 3)),
         paste0("5 years AUC = ", round(roc_obj$AUC[3], 3))
       ),
       col = c("red", "blue", "green"),
       lwd = 2)

#Fit Cox models for CAPG alone and full risk score, compare performance
library(survival)
library(survminer)

# Assuming your data frame is final_df, with columns:
# PFI.time, PFI (event), risk_score, CAPG expression, and clinical covariates

# Check rownames are present in both datasets and match
head(rownames(final_df))
head(rownames(cox_data))

# Create a data frame with CAPG expression and rownames from cox_data
capg_df <- data.frame(CAPG = cox_data$CAPG, sample_id = rownames(cox_data))

# Add sample_id as a column in final_df
final_df$sample_id <- rownames(final_df)

# Merge final_df with capg_df by sample_id
final_merged <- merge(final_df, capg_df, by = "sample_id")

# Check merged data
head(final_merged)

# Now you can run Cox model with CAPG
library(survival)
fit_CAPG <- coxph(Surv(PFI.time, PFI) ~ CAPG + age_at_initial_pathologic_diagnosis + clinical_stage, data = final_merged)
summary(fit_CAPG)

# 2. Cox model with risk_score (full score)
fit_risk <- coxph(Surv(PFI.time, PFI) ~ risk_score + age_at_initial_pathologic_diagnosis + clinical_stage, data = final_merged)

# Compare concordance
cindex_CAPG <- summary(fit_CAPG)$concordance[1]
cindex_risk <- summary(fit_risk)$concordance[1]
print(paste("Concordance CAPG:", round(cindex_CAPG, 3)))
print(paste("Concordance risk score:", round(cindex_risk, 3)))

# 3. KM plots for CAPG groups and risk_score groups
# Define groups by median split
final_merged$CAPG_group <- ifelse(final_merged$CAPG > median(final_merged$CAPG), "High", "Low")
final_merged$risk_group <- ifelse(final_merged$risk_score > median(final_merged$risk_score), "High", "Low")

# KM plot for CAPG
km_CAPG <- survfit(Surv(PFI.time, PFI) ~ CAPG_group, data = final_merged)
ggsurvplot(km_CAPG, data = final_merged, pval = TRUE, risk.table = TRUE, title = "KM plot by CAPG expression")

# KM plot for risk score
km_risk <- survfit(Surv(PFI.time, PFI) ~ risk_group, data = final_merged)
ggsurvplot(km_risk, data = final_merged, pval = TRUE, risk.table = TRUE, title = "KM plot by risk score")

# 4. Time-dependent ROC AUC at 1, 3, 5 years (365, 1095, 1825 days)
times <- c(365, 1095, 1825)

roc_CAPG <- timeROC(T = final_merged$PFI.time, delta = final_merged$PFI,
                    marker = final_merged$CAPG,
                    cause = 1, weighting = "marginal",
                    times = times, iid = TRUE)

roc_risk <- timeROC(T = final_merged$PFI.time, delta = final_merged$PFI,
                    marker = final_merged$risk_score,
                    cause = 1, weighting = "marginal",
                    times = times, iid = TRUE)

# Plot ROC curves for CAPG and risk score at 3 years (for example)
plot(roc_CAPG, time = 1095, col = "blue", title = FALSE)
plot(roc_risk, time = 1095, col = "red", add = TRUE)
legend("bottomright", legend = c("CAPG", "Risk score"), col = c("blue", "red"), lwd = 2)
title("Time-dependent ROC at 3 years")

# Print AUC values
print("AUC CAPG:")
print(roc_CAPG$AUC)
print("AUC risk score:")
print(roc_risk$AUC)


# 5. Correlation between CAPG and risk score
cor_test <- cor.test(final_merged$CAPG, final_merged$risk_score, method = "spearman")


#test association of these genes with loss of y
genes_11 <- c("CAPG", "MTHFD2", "KPNA2", "EGF", "FBXO5", "PTTG1", "PRDX4", "UCK2", "CENPA", "SPAG5", "PLAU")

# Prepare expression matrix subset for 11 genes (rows = genes, columns = samples)
exp_11genes <- exp_data_male[genes_11, ]

# Transpose to samples x genes
exp_11genes_t <- as.data.frame(t(exp_11genes))
exp_11genes_t$sample <- rownames(exp_11genes_t)

# Merge LOY scores
exp_11genes_LOY <- merge(exp_11genes_t, LOY_df_male[, c("sample", "LOY_score")], by = "sample")

# Optional: add LOY group (optimal cutoff)
library(survminer)
# Merge survival data by sample ID (use same ID format)
exp_11genes_LOY_surv <- merge(exp_11genes_LOY, surv_data, by.x = "sample", by.y = "sample")

# Check columns
head(exp_11genes_LOY_surv)

# Then run surv_cutpoint with proper survival columns:
cutpoint <- surv_cutpoint(
  exp_11genes_LOY_surv,
  time = "PFI.time",      # or "PFI.time" if you want progression-free interval
  event = "PFI",          # or "PFI"
  variables = "LOY_score",
  minprop = 0.1
)

opt_cut <- cutpoint$cutpoint$cutpoint
exp_11genes_LOY$LOY_group_opt <- ifelse(exp_11genes_LOY$LOY_score > opt_cut, "Y_high", "Y_low")
cor_results <- data.frame(Gene = genes_11, Correlation = NA, p.value = NA)

for (gene in genes_11) {
  test <- cor.test(exp_11genes_LOY[[gene]], exp_11genes_LOY$LOY_score, method = "spearman")
  cor_results[cor_results$Gene == gene, c("Correlation", "p.value")] <- c(test$estimate, test$p.value)
}

print(cor_results)
wilcox_results <- data.frame(Gene = genes_11, p.value = NA)

for (gene in genes_11) {
  test <- wilcox.test(exp_11genes_LOY[[gene]] ~ exp_11genes_LOY$LOY_group_opt)
  wilcox_results[wilcox_results$Gene == gene, "p.value"] <- test$p.value
}

print(wilcox_results)

library(tidyr)
library(dplyr)
library(ggplot2)

# Assuming exp_11genes_LOY has columns: sample, genes..., LOY_score, LOY_group_opt

# Select only your 11 genes + LOY group
plot_df <- exp_11genes_LOY %>%
  select(all_of(genes_11), LOY_group_opt) %>%
  pivot_longer(cols = all_of(genes_11), names_to = "Gene", values_to = "Expression")

# Now plot
ggplot(plot_df, aes(x = LOY_group_opt, y = Expression, fill = LOY_group_opt)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~ Gene, scales = "free_y") +    # separate plot per gene, free y axis
  theme_minimal() +
  labs(title = "Expression of 11 Genes by LOY group (optimal cutoff)",
       x = "LOY Group", y = "Expression") +
  scale_fill_manual(values = c("Y_low" = "#1f77b4", "Y_high" = "#ff7f0e")) +
  theme(legend.position = "bottom")

library(broom)

pvals_df <- plot_df %>%
  group_by(Gene) %>%
  summarise(
    p.value = wilcox.test(Expression ~ LOY_group_opt)$p.value
  ) %>%
  mutate(
    p.adj = p.adjust(p.value, method = "fdr"),        # adjust for multiple testing
    label = ifelse(p.adj < 0.001, "***",
                   ifelse(p.adj < 0.01, "**",
                          ifelse(p.adj < 0.05, "*", "ns")))
  )
library(dplyr)

pvals_df <- plot_df %>%
  group_by(Gene) %>%
  summarise(
    p.value = wilcox.test(Expression ~ LOY_group_opt)$p.value
  ) %>%
  mutate(
    p.adj = p.adjust(p.value, method = "fdr"),             # adjust p-values
    label = paste0("p = ", signif(p.adj, 3))               # format p-values as text labels
  )
library(ggplot2)

ggplot(plot_df, aes(x = LOY_group_opt, y = Expression, fill = LOY_group_opt)) +
  geom_boxplot() +
  facet_wrap(~ Gene, scales = "free_y") +
  theme_minimal() +
  labs(title = "Expression of 11 Genes by LOY group (optimal cutoff)",
       x = "LOY Group", y = "Expression") +
  scale_fill_manual(values = c("Y_low" = "#1f77b4", "Y_high" = "#ff7f0e")) +
  geom_text(data = pvals_df, aes(x = 1.5, y = Inf, label = label),
            inherit.aes = FALSE, vjust = 1.5, size = 4)
###check combined association
# Standardize gene expression if needed
gene_matrix <- scale(exp_11genes_LOY[, genes_11])
# Regress LOY score on all 11 genes
model_loy <- lm(LOY_score ~ ., data = as.data.frame(cbind(LOY_score = exp_11genes_LOY$LOY_score, gene_matrix)))
summary(model_loy)

#PCA ofpca <- prcomp(gene_matrix, scale. = TRUE)
pca <- prcomp(gene_matrix, scale. = TRUE)
pc1 <- pca$x[, 1]
cor.test(pc1, exp_11genes_LOY$LOY_score, method = "spearman")

#build loy expression score from these 11 genes
exp_11genes_LOY$LOY_signature_score <- rowMeans(gene_matrix)
cor.test(exp_11genes_LOY$LOY_signature_score, exp_11genes_LOY$LOY_score, method = "spearman")

library(ggplot2)

ggplot(exp_11genes_LOY, aes(x = LOY_signature_score, y = LOY_score)) +
  geom_point(alpha = 0.7, color = "#1f77b4") +
  geom_smooth(method = "lm", se = TRUE, color = "#ff7f0e") +
  labs(
    title = "Correlation Between 11-Gene Signature and LOY Score",
    x = "11-Gene Signature Score",
    y = "LOY Score"
  ) +
  theme_minimal()

#Venn diagram of overlapping gene
# Install & load ggVennDiagram
library(ggVennDiagram)

# Prepare gene list
gene_list <- selected_gsea %>%
  filter(ID %in% selected_hallmarks) %>%
  arrange(match(ID, selected_hallmarks)) %>%
  mutate(GeneSymbol = strsplit(GeneSymbol, ",\\s*")) %>%
  pull(GeneSymbol)

# Short names for clarity
names(gene_list) <- c("G2M", "E2F", "EMT", "KRAS", "Mitotic")

# Simple and clean plot
ggVennDiagram(gene_list, 
              label_alpha = 0,     # remove shaded areas
              edge_size = 0.3,     # thin lines
              label = "count") +   # show only overlap counts
  scale_fill_gradient(low = "white", high = "lightblue") +  # subtle fill
  theme_void() +  # remove axis and background
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size = 14)) +
  ggtitle("Core Gene Overlap (5 Hallmark Pathways)")


