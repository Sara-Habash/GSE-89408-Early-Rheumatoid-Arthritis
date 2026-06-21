
###############################################################
# BULK RNA-SEQ DIFFERENTIAL EXPRESSION ANALYSIS
# RA_Early vs Normal
# Corrected + Publication-Grade Workflow
###############################################################

# =========================
# 1. Load Libraries
# =========================

library(DESeq2)
library(apeglm)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)

library(AnnotationDbi)
library(org.Hs.eg.db)


# =========================
# 2. Load raw counts
# =========================
counts <- read.table("GSE89408_GEO_count_matrix_rename.txt.gz",
                     header = TRUE,      # Does the first row contain column names?
                     row.names = 1,      # Use the first column (e.g., Gene IDs) as row names
                     sep = "\t",         # Use "\t" for tabs, " " for spaces, or "," for CSV
                     check.names = FALSE) # Prevents R from changing symbols like '-' to '.' header = TRUE, 


counts
counts <- as.matrix(counts)   # Must be a matrix of integers

# =========================
# 2. Load metadata
# =========================
coldata <- read.csv("RA_metadata.csv", row.names = 1)
# coldata must have samples as rows, conditions as columns

# ======================================================================
# 3. Verify sample order matches between counts and coldataLoad metadata
# ======================================================================
rownames(coldata) = colnames(counts)  # Must return TRUE then type the next command

all(rownames(coldata) == colnames(counts)) # Must return TRUE
# =================================
# 4. #change name of diseases in R
# =================================
coldata$condition =as.character(coldata$condition)
view(coldata)
coldata$condition [coldata$condition == "RA _ established"] = "RA_Est" # copy name of disease from coldata file bs spaces matter
coldata$condition [coldata$condition == "RA _early"] = "RA_Early"
# Then back to factor
coldata$condition = factor(coldata$condition)

# ====================
# 4. Check Conditions
# ====================

table(coldata$condition)

# =====================
# 5. Select Conditions
# =====================
selected_conditions <- c(
       "Normal",
     "RA_Early" )
# ===================
# 6. Subset Metadata
# ===================

coldata_subset <- coldata[
  
  coldata$condition %in% selected_conditions,
  
  ,
  
  drop = FALSE
]
# ========================
# 7. Check Sample Numbers
# ========================

table(coldata_subset$condition)
# =========================
# 8. Subset Count Matrix
# =========================

counts_subset <- counts[
  
  ,
  
  rownames(coldata_subset)
] 
# =====================
# 9. Ensure Same Order
# =====================

all(
  
  colnames(counts_subset) ==
    
    rownames(coldata_subset)
)

# Must return TRUE

# =========================
# 10. Convert Metadata Types
# =========================

coldata_subset$condition <- factor(
  
  coldata_subset$condition
)

coldata_subset$sex <- factor(
  
  coldata_subset$sex
)

coldata_subset$AGE <- as.numeric(
  
  coldata_subset$AGE
)
#Scale age (IMPORTANT for model stability)
coldata_subset$AGE_scaled <- scale(coldata_subset$AGE)


# =========================
# 11. Set Reference Group
# =========================

coldata_subset$condition <- relevel(
  
  coldata_subset$condition,
  
  ref = "Normal"
)
# =========================
# 12. Check Factor Levels
# =========================

levels(coldata_subset$condition)
# =============================
# 13. Rounding of counts_subset
# =============================

counts_subset <- round(counts_subset)  # DESeq2 needs integer counts

# =========================
# 14. Create DESeq2 Object
# =========================

dds <- DESeqDataSetFromMatrix(
  
  countData = counts_subset,
  
  colData = coldata_subset,
  
  design = ~ sex + AGE_scaled + condition
)
# =========================
# 15. Filter Low Count Genes
# =========================

keep <- rowSums(
  
  counts(dds) >= 10
  
) >= 3

dds <- dds[keep, ]

# =========================
# 16. Run DESeq2
# =========================

dds <- DESeq(dds)

# =========================
# 17. Check Model Coefficients
# =========================

resultsNames(dds)

# =========================
# 18. Extract Results
# =========================

res <- results(
  
  dds,
  
  contrast = c(
    
    "condition",
    
    "RA_Early",
    
    "Normal"
  )
)
# =========================
# 19. Shrink Log Fold Changes
# =========================

resLFC <- lfcShrink(
  
  dds,
  
  coef = "condition_RA_Early_vs_Normal",
  
  type = "apeglm"
)
# =========================
# 20. Convert to DataFrame
# =========================

res_df <- as.data.frame(resLFC)

# ===================
# 21. Add Gene Names
# ===================

res_df$gene <- rownames(res_df)

# =====================
# 22. Remove NA Values
# =====================

res_df <- na.omit(res_df)
# =========================
# 23. Order by Adjusted P-value
# =========================

res_df <- res_df[
  
  order(res_df$padj),
  
]
# =========================
# 24. View Top Results
# =========================

head(res_df)
# =========================
# 25. Gene Annotation
# =========================

res_df$entrez <- mapIds(
  org.Hs.eg.db,
  keys = res_df$gene,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)

res_df$ensembl <- mapIds(
  org.Hs.eg.db,
  keys = res_df$gene,
  column = "ENSEMBL",
  keytype = "SYMBOL",
  multiVals = "first"
)

res_df$gene_name <- mapIds(
  org.Hs.eg.db,
  keys = res_df$gene,
  column = "GENENAME",
  keytype = "SYMBOL",
  multiVals = "first"
)

# Replace NA
res_df$gene_name[is.na(res_df$gene_name)] <- res_df$gene[is.na(res_df$gene_name)]

# =========================
# 26. Significant Genes
# =========================

sig_genes <- res_df %>%
  
  filter(
    
    padj < 0.05 &
      
      abs(log2FoldChange) > 1
  )

# =========================
# 27. Number of Significant Genes
# =========================

nrow(sig_genes)

# =========================
# 28. Save Results
# =========================

write.csv(res_df, "RA_Early_vs_Normal_ALL_DEGs.csv", row.names = FALSE)
write.csv(sig_genes, "RA_Early_vs_Normal_SIG_DEGs.csv", row.names = FALSE)

# =========================
# 29. Variance Stabilization
# =========================

vsd <- vst(dds)

# =========================
# 30. PCA Plot
# =========================

plotPCA(vsd, intgroup = "condition")

pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(
    title = "PCA Plot",
    x = paste0("PC1: ", percentVar[1], "%"),
    y = paste0("PC2: ", percentVar[2], "%")
  )
# =========================
# 31. Sample Distance Heatmap
# =========================

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

pheatmap(sampleDistMatrix)

# =========================
# 32. MA Plot
# =========================

plotMA(resLFC, ylim = c(-5, 5))

# =========================
# 24. Volcano Plot (Annotated)
# =========================

res_df$significance <- case_when(
  res_df$padj < 0.05 & res_df$log2FoldChange > 1 ~ "Up",
  res_df$padj < 0.05 & res_df$log2FoldChange < -1 ~ "Down",
  TRUE ~ "NS"
)
top_labeled <- res_df %>% arrange(padj) %>% head(10)

ggplot(res_df, aes(log2FoldChange, -log10(padj), color = significance)) +
  geom_point(alpha = 0.5) +
  geom_point(data = top_labeled, color = "red") +
  geom_text_repel(data = top_labeled, aes(label = gene_name)) +
  theme_classic()

# =========================
# 25. Heatmap of Top Genes
# =========================

normalized_counts <- counts(dds, normalized = TRUE)

top_genes <- sig_genes %>%
  arrange(padj) %>%
  head(30) %>%
  pull(gene)

top_genes <- intersect(top_genes, rownames(normalized_counts))
heatmap_data <- normalized_counts[top_genes, rownames(coldata_subset)]
heatmap_scaled <- t(scale(t(heatmap_data)))
annotation_col <- data.frame(condition = coldata_subset$condition)
rownames(annotation_col) <- rownames(coldata_subset)

pheatmap(
  heatmap_scaled,
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 8,
  main = "Top DEGs Heatmap"
)

# =========================
# 26. Dispersion Plot
# =========================

plotDispEsts(dds)

# _____ FUNCTIONAL ENRICHMENT ANALYSIS____________

library(clusterProfiler)
library(enrichplot)
BiocManager::install("pathview")
library(pathview)
BiocManager::install("ReactomePA")
library(ReactomePA)

# =========================
# 27. Prepare Gene Lists
# =========================
sig_genes_enrich <- sig_genes %>%
  filter(!is.na(entrez))

gene_list <- sig_genes_enrich$entrez

background_genes <- res_df %>%
  filter(!is.na(entrez)) %>%
  pull(entrez)
# =========================
# 28. GO Enrichment (BP)
# =========================
ego_bp <- enrichGO(
  gene          = gene_list,
  universe      = background_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  
readable      = TRUE
)

write.csv(as.data.frame(ego_bp), "GO_BP_RA_Early_vs_Normal.csv")

dotplot(ego_bp, showCategory = 15) + ggtitle("GO Biological Process")
barplot(ego_bp, showCategory = 15)
emapplot(pairwise_termsim(ego_bp), showCategory = 15)

#GO Enrichment (MF)
ego_MF <- enrichGO(
  gene          = gene_list,
  universe      = background_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  
  readable      = TRUE
)
write.csv(as.data.frame(ego_MF), "GO_MF_RA_Early_vs_Normal.csv")
dotplot(ego_MF, showCategory = 15) + ggtitle("GO Molecular Function")
barplot(ego_MF, showCategory = 15)
emapplot(pairwise_termsim(ego_MF), showCategory = 15)

#GO Enrichment (CC)
ego_CC <- enrichGO(
  gene          = gene_list,
  universe      = background_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  
  readable      = TRUE
)

write.csv(as.data.frame(ego_CC), "GO_CC_RA_Early_vs_Normal.csv")
dotplot(ego_CC, showCategory = 15) + ggtitle("GO Cellular Component")
barplot(ego_CC, showCategory = 15)
emapplot(pairwise_termsim(ego_CC), showCategory = 15)

# =========================
# 29. KEGG Enrichment
# =========================
ekegg <- enrichKEGG(
  gene         = gene_list,
  universe     = background_genes,
  organism     = "hsa",
  pvalueCutoff = 0.05
)

ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

write.csv(as.data.frame(ekegg), "KEGG_RA_Early_vs_Normal.csv")

dotplot(ekegg, showCategory = 15) +
  ggtitle("KEGG Pathways")


# =========================
# 30. Reactome Enrichment
# =========================

ereact <- enrichPathway(
  gene         = gene_list,
  organism     = "human",
  pvalueCutoff = 0.05,
  readable     = TRUE
)

write.csv(as.data.frame(ereact), "Reactome_RA_Early_vs_Normal.csv")

dotplot(ereact, showCategory = 12) +
  ggtitle("Reactome Pathways")

# =========================
# 31. GSEA-ready Ranking
# =========================

gene_rank <- res_df %>%
  filter(!is.na(entrez)) %>%
  select(entrez, log2FoldChange)

gene_rank_vector <- gene_rank$log2FoldChange
names(gene_rank_vector) <- gene_rank$entrez

# =========================
# 32. GSEA GO (Optional)
# =========================

gsea_go <- gseGO(
  geneList     = sort(gene_rank_vector, decreasing = TRUE),
  OrgDb        = org.Hs.eg.db,
  keyType      = "ENTREZID",
  ont          = "BP",
  nPerm        = 1000,
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05
)

ridgeplot(gsea_go)

# =========================
# 33. Save Session
# =========================

save(
  res_df, sig_genes,
  ego_bp,ego_MF,ego_CC, ekegg, ereact,
  file = "RA_Early_vs_Normal_full_analysis.RData"
)

# =========================
# DONE
# =========================





