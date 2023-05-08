library(DESeq2)
library(apeglm)
library(org.Gg.eg.db)
library(EnhancedVolcano)

dat <- read.csv("CountMatrix.csv", header = T, row.names = 1)
info <- read.table("colData.txt", header = T, sep = '\t')

dds <- DESeqDataSetFromMatrix(dat, info, ~condition)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

ddsDE <- DESeq(dds)

normCounts <- counts(ddsDE, normalized = T)

res <- results(ddsDE, alpha = 0.05)
res.df <- as.data.frame(res)
res.df$symbol <- mapIds(org.Gg.eg.db, keys = rownames(res.df),keytype = "ENSEMBL", column = "SYMBOL")
res.dfO <- res.df[order(res.df$padj),]
write.csv(res.dfO, "DESeq.csv")

normCounts.df <- as.data.frame(normCounts)
normCounts.df$symbol <- mapIds(org.Gg.eg.db, keys = rownames(normCounts.df),keytype = "ENSEMBL", column = "SYMBOL")
write.csv(normCounts.df, "NormCounts.csv")

EnhancedVolcano(res.df, x = "log2FoldChange", y = "padj", lab = res.df$symbol,
                pCutoff = 1e-4, FCcutoff = 1)