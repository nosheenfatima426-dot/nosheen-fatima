# ==============================================================================
# Pipeline: Nanoparticle Exposure Transcriptomics in Lung Cancer (A549)
# Dataset: NCBI GEO GSE53700
# Target: SM30 3-hour Recovery vs. Untreated Controls
# ==============================================================================

# 1. Load Required Libraries ---------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
packages <- c("GEOquery", "limma", "clusterProfiler", "org.Hs.eg.db", 
              "pheatmap", "EnhancedVolcano", "dplyr")
BiocManager::install(setdiff(packages, installed.packages()[,"Package"]))

library(GEOquery)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(EnhancedVolcano)
library(pheatmap)
library(dplyr)

# 2. Download and Prepare Expression Data ---------------------------------------
# Download GSE53700 from GEO
gse <- getGEO("GSE53700", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]

# Define target samples from thesis methodology
# Treated (SM30 3h): GSM1299149, GSM1299150, GSM1299151
# Control: GSM1299152, GSM1299153, GSM1299154
target_samples <- c("GSM1299149", "GSM1299150", "GSM1299151", 
                    "GSM1299152", "GSM1299153", "GSM1299154")
sub_gse <- gse[, target_samples]

# Log2 transformation check
ex <- exprs(sub_gse)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
LogC <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
if (LogC) {
  ex[which(ex <= 0)] <- NaN
  exprs(sub_gse) <- log2(ex)
}

# 3. Differential Expression Analysis (limma) ----------------------------------
group <- factor(c("SM30", "SM30", "SM30", "Control", "Control", "Control"),
                levels = c("Control", "SM30"))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# Fit linear model
fit <- lmFit(sub_gse, design)
contrast_matrix <- makeContrasts(SM30_vs_Control = SM30 - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Extract full results
deg_results <- topTable(fit2, coef = "SM30_vs_Control", number = Inf, adjust.method = "fdr")

# Clean gene annotations (extract Gene Symbol from metadata)
deg_results$GeneSymbol <- fData(sub_gse)$`Gene symbol`
# Remove probes without gene symbols or with ambiguous mappings
deg_clean <- deg_results %>% 
  filter(!is.na(GeneSymbol) & GeneSymbol != "")

# Apply thesis cutoffs: adj.P.Val < 0.05 and |logFC| >= 1
significant_degs <- deg_clean %>% 
  filter(adj.P.Val < 0.05 & abs(logFC) >= 1)

upregulated_genes <- significant_degs %>% filter(logFC >= 1)
downregulated_genes <- significant_degs %>% filter(logFC <= -1)

# Export lists
write.csv(significant_degs, "significant_degs_SM30.csv", row.names = FALSE)

# 4. Data Visualization --------------------------------------------------------
# Volcano Plot
EnhancedVolcano(deg_clean,
                lab = deg_clean$GeneSymbol,
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                pointSize = 2.0,
                labSize = 4.0,
                title = 'SM30 Exposure vs Control (A549)',
                subtitle = 'Cutoff: FDR < 0.05, |log2FC| >= 1')

# Heatmap of Top 50 DEGs (Z-score normalized)
top50_genes <- bind_rows(
  upregulated_genes %>% slice_max(order_by = abs(logFC), n = 25),
  downregulated_genes %>% slice_max(order_by = abs(logFC), n = 25)
)

heatmap_mat <- exprs(sub_gse)[rownames(top50_genes), ]
rownames(heatmap_mat) <- top50_genes$GeneSymbol
pheatmap(heatmap_mat,
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_method = "average",
         show_colnames = TRUE,
         main = "Top 50 DEGs Expression Profile")

# 5. Pathway Enrichment Analysis (clusterProfiler) -----------------------------
# Convert Symbols to Entrez IDs
gene_mapping <- bitr(significant_degs$GeneSymbol, 
                     fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db)

# KEGG Pathway Enrichment
kegg_enrich <- enrichKEGG(gene = gene_mapping$ENTREZID,
                          organism = "hsa",
                          pvalueCutoff = 0.05)

# Gene Ontology (Biological Process)
go_enrich <- enrichGO(gene = gene_mapping$ENTREZID,
                      OrgDb = org.Hs.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      readable = TRUE)

# Plot Pathway Results
dotplot(kegg_enrich, showCategory = 10, title = "Top Enriched KEGG Pathways")