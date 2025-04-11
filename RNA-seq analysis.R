library(DESeq2)
library(dplyr)
library(apeglm)
library(pheatmap)
library(viridis)

cts = read.csv("GSE216444_subtypes_TOC_postQC.csv",sep=",",head=TRUE,row.names = 1)
cts = as.matrix(cts)
coldata = read.csv("GSE216444_subtypes_colData_postQC.csv",row.names=1)
coldata$Sex = as.factor(coldata$Sex)
coldata$Condition = as.factor(coldata$Condition)
coldata$Timepoint = as.factor(coldata$Timepoint)
coldata$Batch = as.factor(coldata$Batch)
coldata$Population = as.factor(coldata$Population)
coldata$Time_Pop_Cond = as.factor(coldata$Time_Pop_Cond)
nfkbgenes <- read.table("HALLMARK_TNFA_SIGNALING_VIA_NFKB_geneset.csv", quote="\"", comment.char="")

dds = DESeqDataSetFromMatrix(countData=cts, colData=coldata, design=~ Sex + Time_Pop_Cond)
smallestGroupSize = 7
keep = rowSums(counts(dds) >= 10) >= smallestGroupSize
dds = dds[keep,]
dds$Time_Pop_Cond = relevel(dds$Time_Pop_Cond, ref="3D_MRTD_Cont")
dds = DESeq(dds)
res = results(dds)

res_order  = res[order(-res$log2FoldChange),]
res_order_df = as.data.frame(res_order)
res_order_df$gene_id = rownames(res_order_df)
res_order_df = res_order_df %>% filter(gene_id %in% nfkbgenes$V1)
res_order_df = res_order_df[1:20,]

samples = c("MRTD_3D_A_c","MRTD_3D_B_c","MRTD_3D_C_c","MRTD_3D_D_c","MRTD_3D_E_c","MRTD_3D_F_c","MRTD_3D_G_c","MRTD_3D_H_c","MRTD_3D_A_i","MRTD_3D_B_i","MRTD_3D_C_i","MRTD_3D_D_i","MRTD_3D_E_i","MRTD_3D_F_i","MRTD_3D_G_i","MRTD_3D_H_i")
samples = factor(samples, level=samples)
geneset = res_order_df$gene_id
geneset = factor(geneset, level=geneset)

vsd = varianceStabilizingTransformation(dds)
df = as.data.frame(colData(dds)[,"Condition"], row.names=rownames(coldata), col.names="Side")
colnames(df) = "Side"
pheatmap(assay(vsd)[rownames(assay(vsd)) %in% geneset, colnames(assay(vsd)) %in% samples], cluster_rows=FALSE, show_rownames=TRUE, show_colnames=FALSE, cluster_cols=FALSE, annotation_col=df, border_color = NA, annotation_names_col = FALSE, color=viridis(1000))

#Order of rows determined by order of rows in counts file
#I can't work out how to change this post-running DESeq2 so work out 20 most differentially expressed genes and place these at top of dataset before re-running DESeq2
to_plot = res_order_df$gene_id
to_plot = as.data.frame(to_plot)
cts_working = read.csv("GSE216444_subtypes_TOC_postQC.csv",sep=",",head=TRUE,row.names = 1)
cts_working$gene_id = rownames(cts_working)
to_plot = left_join(to_plot, cts_working, by=c("to_plot"="gene_id"))
rownames(to_plot) = to_plot[,1]
to_plot = to_plot[,-1]
cts_working = cts_working[,-17]
cts_working = cts_working[!rownames(cts_working) %in% rownames(to_plot),]
cts_working = rbind(to_plot, cts_working)
cts_working = as.matrix(cts_working)
dds_working = DESeqDataSetFromMatrix(countData=cts_working, colData=coldata, design=~ Sex + Condition)
keep_working = rowSums(counts(dds_working) >= 10) >= smallestGroupSize
dds_working = dds_working[keep_working,]
dds_working = DESeq(dds_working)
 
vsd_working = varianceStabilizingTransformation(dds_working)
vsd_working = assay(vsd_working)
z = t(scale(t(vsd_working)))
pheatmap(z[rownames(z) %in% geneset, colnames(z) %in% samples], cluster_rows=FALSE, show_rownames=TRUE, show_colnames=FALSE, cluster_cols=FALSE, annotation_col=df, border_color = NA, annotation_names_col = FALSE, color=viridis(1000), scale="row")