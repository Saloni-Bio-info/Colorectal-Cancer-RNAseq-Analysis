library(cluster)
install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)
BiocManager::install("GenomeInfoDbData")
BiocManager::valid()
options(repos = BiocManager::repositories())
BiocManager::install("GenomeInfoDbData", force = TRUE)
library(DESeq2)
sessionInfo()
BiocManager::install("GEOquery")
library(GEOquery)
install.packages("pheatmap")
library(pheatmap)
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
gset
library(GEOquery)
gset <- getGEO("GSE45827", GSEMatrix = TRUE)
library(DESeq2)
example(DESeq)
head(res)
summary(res)
plotMA(res)
plotDispEsts(dds)
expr <- read.table(
  "~/Downloads/GSE104836_gene_exp.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)
head(expr)
colnames(expr)
count_cols <- grep("_COUNT$", colnames(expr))
counts <- expr[, count_cols]
count_cols <- grep("_COUNT$", colnames(expr))
counts <- expr[, count_cols]
dim(counts)
head(counts)
rownames(counts) <- expr$gene
sum(duplicated(expr$gene))
gene_names <- make.unique(as.character(expr$gene))
rownames(counts) <- gene_names
gene_names <- make.unique(as.character(expr$gene))
rownames(counts) <- gene_names
head(rownames(counts))
sample_info <- data.frame(
  condition = factor(c(
    "Cancer","Normal",
    "Cancer","Normal",
    "Cancer","Normal",
    "Cancer","Normal",
    "Cancer","Normal",
    "Cancer","Normal",
    "Cancer","Normal",
    "Cancer","Normal",
    "Cancer","Normal",
    "Cancer","Normal"
  ))
)
rownames(sample_info) <- colnames(counts)
library(DESeq2)
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = sample_info,
  design = ~ condition
)
dds
dds <- DESeq(dds)
res <- results(dds)
summary(res)
head(res)
sig <- subset(res,
              padj < 0.05 &
                abs(log2FoldChange) > 1)
nrow(sig)
sig <- sig[order(sig$padj), ]
head(sig)
write.csv(
  as.data.frame(sig),
  "Breast_Cancer_DEGs.csv"
)
up_genes <- subset(sig, log2FoldChange < -1)
down_genes <- subset(sig, log2FoldChange > 1)
nrow(up_genes)
nrow(down_genes)
res_df <- as.data.frame(res)
res_df <- na.omit(res_df)
library(EnhancedVolcano)
EnhancedVolcano(
  res_df,
  lab = rownames(res_df),
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1,
  title = "Breast Cancer vs Normal",
  subtitle = "GSE104836"
)
pdf("VolcanoPlot_GSE104836.pdf",
    width = 8,
    height = 6)
EnhancedVolcano(
  res_df,
  lab = rownames(res_df),
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1
)
dev.off()
library(pheatmap)
norm_counts <- counts(dds, normalized = TRUE)
top50 <- rownames(sig)[1:50]
heatmap_data <- norm_counts[top50, ]
heatmap_data <- log2(heatmap_data + 1)
annotation_col <- data.frame(
  Condition = sample_info$condition
)
rownames(annotation_col) <- colnames(heatmap_data)
pheatmap(
  heatmap_data,
  scale = "row",
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = TRUE
)
pdf("Heatmap_GSE104836.pdf",
    width = 10,
    height = 8)
pheatmap(
  heatmap_data,
  scale = "row",
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = TRUE
)
dev.off()
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
library(clusterProfiler)
library(org.Hs.eg.db)
gene_list <- rownames(sig)
length(gene_list)
head(gene_list)
gene_df <- bitr(
  gene_list,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
head(gene_df)
nrow(gene_df)
ego <- enrichGO(
  gene = gene_df$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
library(enrichplot)
dotplot(ego, showCategory = 20)
dotplot(ego, showCategory = 20)
p <- dotplot(ego, showCategory = 20)
p
library(ggplot2)
p <- dotplot(ego, showCategory = 20)
print(p)
dev.off()
p <- dotplot(ego, showCategory = 20)
print(p)
pdf("GO_Dotplot_GSE104836.pdf",
    width = 10,
    height = 8)
print(p)
dev.off()
png("GO_Dotplot_GSE104836.png",
    width = 2000,
    height = 1600,
    res = 300)
print(p)
dev.off()
getwd()
gene_df
kegg <- enrichKEGG(
  gene = gene_df$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)
p_kegg <- dotplot(kegg, showCategory = 20)
print(p_kegg)
pdf("KEGG_Dotplot_GSE104836.pdf",
    width = 10,
    height = 8)
print(p_kegg)
dev.off()
kegg <- enrichKEGG(
  gene = gene_df$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)
nrow(as.data.frame(kegg))
[1] 41
up_genes <- subset(
  res,
  padj < 0.05 &
    log2FoldChange > 1
)
up_genes <- up_genes[order(-up_genes$log2FoldChange), ]
head(up_genes, 20)
down_genes <- subset(
  res,
  padj < 0.05 &
    log2FoldChange < -1
)
down_genes <- down_genes[order(down_genes$log2FoldChange), ]
head(down_genes, 20)
write.csv(
  as.data.frame(up_genes),
  "Top_Upregulated_Genes.csv"
)
write.csv(
  as.data.frame(down_genes),
  "Top_Downregulated_Genes.csv"
)
head(
  data.frame(
    Gene = rownames(up_genes),
    log2FC = up_genes$log2FoldChange,
    padj = up_genes$padj
  ),
  10
)
head(
  data.frame(
    Gene = rownames(down_genes),
    log2FC = down_genes$log2FoldChange,
    padj = down_genes$padj
  ),
  10
)
vsd <- vst(dds)
plotPCA(vsd, intgroup = "condition")
write.csv(
  as.data.frame(sig),
  "Significant_DEGs_GSE104836.csv"
)
write.csv(
  as.data.frame(up_genes),
  "Upregulated_Genes.csv"
)
write.csv(
  as.data.frame(down_genes),
  "Downregulated_Genes.csv"
)
write.csv(
  as.data.frame(ego),
  "GO_Enrichment.csv"
)
write.csv(
  as.data.frame(kegg),
  "KEGG_Enrichment.csv"
)
savehistory("GSE104836_commands.Rhistory")
