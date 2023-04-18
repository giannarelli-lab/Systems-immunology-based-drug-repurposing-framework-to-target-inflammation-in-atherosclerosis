setwd("~/Desktop/athero")

rm(list=ls())

library(DESeq)
library(RColorBrewer)

setwd("~/Desktop/athero")

source("loadData.r")

# DESeq analysis of RNA-seq data
# --------------------------------------------------------------------------------
cds = newCountDataSet(counts, sample_info$condition)  # DESeq count data set

cds = estimateSizeFactors(cds)

rna_norm = counts(cds, normalized=TRUE)  # normalized by size factor
rownames(rna_norm) = gene_symbols

x_zscore = t(scale(t(rna_norm)))  # standardized data

x_log = log2(rna_norm + 1)  # log transform


# Average normalized RNA-seq counts over experimental repeats
# ----------------------------------------------------------------------
ave_rna_norm = matrix(NA, nrow=nrow(rna_norm), ncol=length(unique(sample_info$exp)))

col_order = match(1:max(sample_info$exp), sample_info$exp)
colnames(ave_rna_norm) = colnames(rna_norm)[col_order]
rownames(ave_rna_norm) = rownames(rna_norm)

# Calculate replicate means
for (exp in unique(sample_info$exp)) {
	ave_rna_norm[,exp] = apply(rna_norm[, sample_info$exp == exp], 1, mean, na.rm=TRUE)
}

x_ave_log = log2(ave_rna_norm + 1)
sample_info_ave = sample_info[col_order,]


# PCA on normalized counts (by size factor) averaged over technical repeats.
colors = brewer.pal(9, "Set1")

pca_rna = prcomp(t(x_ave_log))  # PCA 


# 2D PCA plot
pdf("plots/rna_seq_cytof_pca2.pdf", height=4, width=7)

par(mfrow=c(1, 2))

col_vector = colors[as.numeric(sample_info_ave$condition)]

plot(pca_rna$x[,1], pca_rna$x[,2],
	main="RNA-seq",
	pch=16,
	col=col_vector,
	xlab="PC1", ylab="PC2"
)
legend("topright",
	pch=16,
	col=colors[c(1, 2)],
	legend=c("Healthy", "Atherosclerosis"))

pca_cytof = prcomp(t(cytof))

col_vector = rep(NA, nrow(cytof_sample_info))
col_vector[cytof_sample_info$status == "HEALTHY"] = colors[1]
col_vector[cytof_sample_info$status == "DISEASE_UNSTABLE"] = colors[2]
col_vector[cytof_sample_info$status == "DISEASE_STABLE"] = colors[2]

plot(pca_cytof$x[,1], pca_cytof$x[,2],
	pch=16, col=col_vector,
	main="CyTOF",
	xlab="PC1",
	ylab=""
	)
dev.off()

library(scatterplot3d)
source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')  # for addgrids


pdf("plots/rna_seq_cytof_pca2_3d.pdf", height=5, width=9)
par(mfrow=c(1, 2))
col_vector = colors[as.numeric(sample_info_ave$condition)]

grid_col = rgb(0.9, 0.9, 0.9)

s3d = scatterplot3d(pca_rna$x[,1:3],
	main="RNA-seq",
	pch="",
	cex.axis=0.7,
	box=FALSE, grid=FALSE,
	# tick.marks=FALSE
)

addgrids3d(pca_rna$x[, 1:3],
	grid=c("xy", "xz", "yz"),
	col=grid_col
)

s3d$points3d(pca_rna$x[, 1:3],
	col=col_vector,
	pch=16,
	type="h"
	)

legend("topright",
	pch=16,
	col=colors[c(1, 2)],
	bg="white",
	legend=c("Healthy", "Atherosclerosis"))

col_vector = rep(NA, nrow(cytof_sample_info))
col_vector[cytof_sample_info$status == "HEALTHY"] = colors[1]
col_vector[cytof_sample_info$status == "DISEASE_UNSTABLE"] = colors[2]
col_vector[cytof_sample_info$status == "DISEASE_STABLE"] = colors[2]

s3d = scatterplot3d(pca_cytof$x[,1:3],
	main="CyTOF",
	pch="",
	cex.axis=0.7,
	box=FALSE, grid=FALSE,
	# tick.marks=FALSE
)

addgrids3d(pca_cytof$x[, 1:3],
	grid=c("xy", "xz", "yz"),
	col=grid_col
)

s3d$points3d(pca_cytof$x[, 1:3],
	col=col_vector,
	pch=16,
	type="h"
)

dev.off()


# library(R.matlab)  # for importing matlab files

## Figure 3B
library(DESeq)
library(RColorBrewer)

library(data.table)
library(gplots)

rm(list=ls())

# setwd("~/Desktop/athero")
setwd("~/GoogleDrive/projects/athero-chiara")

source("loadData.r")
source("func.r")


# DESeq analysis of RNA-seq data
# --------------------------------------------------------------------------------
cds = newCountDataSet(counts, sample_info$condition)  # DESeq count data set

# cds = newCountDataSet(counts, sample_info$phenotype)  # DESeq count data set

# Size factors
cds = estimateSizeFactors(cds)
sizeFactors(cds)

rna_norm = counts(cds, normalized=TRUE)  # normalized by size factor
rownames(rna_norm) = gene_symbols

x_zscore = t(scale(t(rna_norm)))  # standardized data

x_log = log2(rna_norm + 1)  # log transform


# PCA on normalized counts (by size factor)
colors = brewer.pal(9, "Set1")
pca = prcomp(t(x_log))


pdf("plots/rna_seq_cytof_pca.pdf", height=5, width=9)
par(mfrow=c(1, 2))
col_vector = colors[as.numeric(sample_info$phenotype)]
col_vector[is.na(col_vector)] = colors[3]

point_pch = c(17, 16)[as.numeric(sample_info$condition)]
plot(pca$x[,1], pca$x[,2],
     main="RNA-seq. DESeq size factor norm., log2(x + 1) transf.",
     pch=point_pch,
     col=col_vector,
     xlab="PC1", ylab="PC2"
)
# text(pca$x[,1], pca$x[,2] + 1400, labels=sample_info$id, cex=0.5)
# text(pca$x[,1], pca$x[,2] + 1, labels=sample_info$id, cex=0.5)
legend("topright",
       pch=c(17, 16, 16),
       col=colors[c(3, 1, 2)],
       legend=c("Healthy", "Stable atherosclerosis", "Unstable atherosclerosis"))

pca_cytof = prcomp(t(cytof))
# pca_cytof = prcomp(t(cytof_zscore))

point_pch = rep(16, nrow(cytof_sample_info))
point_pch[cytof_sample_info$status == "HEALTHY"] = 17

col_vector = rep(NA, nrow(cytof_sample_info))
col_vector[cytof_sample_info$status == "HEALTHY"] = colors[3]
col_vector[cytof_sample_info$status == "DISEASE_UNSTABLE"] = colors[1]
col_vector[cytof_sample_info$status == "DISEASE_STABLE"] = colors[2]

plot(pca_cytof$x[,1], pca_cytof$x[,2],
     pch=point_pch, col=col_vector,
     main="CyTOF",
     xlab="PC1", ylab="PC2")
# text(pca_cytof$x[,1], pca_cytof$x[,2] + 3, label=cytof_sample_info$status, cex=0.5)
dev.off()




# Estimate dispersion for each gene
cds = estimateDispersions(cds, sharingMode="maximum")

# Fit information
str(fitInfo(cds))
plotDispEsts(cds)

# Negative binomial test for control vs disease.
res = nbinomTest(cds, "control", "disease")
res$gene = gene_symbols


write.csv(res, file="differential_expression/diff_expr_PBMC.csv", row.names=FALSE)
# res = read.csv("differential_expression/diff_expr_PBMC.csv")

res[res$gene == "ZHX2",]

head(res[order(res$padj), c("gene", "log2FoldChange", "padj")], n=30)

fc_thresh = 1.2
mean_count_thresh = 4.0
# mean_count_thresh = 10.0

# sig_gene_high_fc = res$gene[which(
# 	res$padj < 0.05 & abs(res$log2FoldChange) > fc_thresh
# 	)
# ]

sig_gene_high_fc = res$gene[which(
  res$padj < 0.05 & abs(res$log2FoldChange) > fc_thresh & res$baseMean > mean_count_thresh
)
]


rids = match(sig_gene_high_fc, gene_symbols)

# Plot heatmap with 

color_scale = colorRampPalette(rev(brewer.pal(11, "RdBu")), interpolate="spline")(100)
pdf("plots/rna-seq-z-scores/fc_filtered_sig.pdf", width=5.5, height=10)
par(cex.main=0.8)
hmap = heatmap.2(
  x_zscore[rids,],
  # main=strsplit(file_names[i], "[.]")[[1]][1],
  main="RNA-seq, sig highest effect",
  margins=c(5, 7),
  key.title="z-score",
  trace="none",
  col=color_scale,
  ColSideColors=colors[-as.numeric(sample_info$condition) + 3],
  cexRow=0.7,
  labCol="",
  dendrogram="row",
  tracecol="black"
)
dev.off()

log_fc_vals = res[rids, "log2FoldChange"][hmap$rowInd]
names(log_fc_vals) = res[rids, "gene"][hmap$rowInd]

log_fc_vals[is.infinite(log_fc_vals)] = NA

light_green = rgb(178, 223, 138, maxColorValue=255)  # light green
light_blue = rgb(166, 206, 227, maxColorValue=255)

col_vec = rep(light_green, length(log_fc_vals))
col_vec[log_fc_vals < 0] = light_blue

pdf("plots/rna-seq-z-scores/fc_filtered_sig_row_barplot.pdf", width=1.7)
barplot(abs(log_fc_vals),
        las=1,
        cex.names=0.3,
        xlab="abs. log2 fold change",
        col=col_vec,
        border=NA, horiz=TRUE, space=0)
dev.off()

#### Figure 3E
heat<-read.csv("~/common.genes_GO0006954_deseq_transcripts.csv")

#scale data
scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

phmat<-scale_rows(heat)

anno_df = data.frame(group=meta$group,
                     stringsAsFactors = FALSE)
anno_df$group <- factor(meta$group, levels = c("Disease", "Healthy"))
head(anno_df)

annotation_colors = list(
  group=c("Disease"="#D51F26","Healthy"="#80b1d3"))

head(annotation_colors)
phmat <- t(scale(t(heat)))
phmat[phmat> 2] <- 2
phmat[phmat < -2] <- -2

rownames(anno_df) <- colnames(phmat)

library(RColorBrewer)
p<-pheatmap(phmat,annotation_col=anno_df, annotation_colors=annotation_colors,
            fontsize_row = 12, cellwidth = 12,cellheight=11,cluster_cols = T,cluster_rows = T,
            main ="GO:0006954(Inflammatory genes) differential expressed in our dataset",# labels_col=empty.cols,
            fontsize_col = 12)
pdf("~/GO006954_genes_heatmap.pdf",  width = 12, height = 20)
print(p)
dev.off()

### Pathway analysis
library(ggplot2)
library(cowplot)

#BP pathway
bp_up<-read.csv("~/selected_up_BP_pathways.csv",sep=",",header=T,row.names=NULL)

pdf("~/top12_bp_up_padj.pdf", height=8,width=14)
bp_up$star <- ""
bp_up$star[bp_up$Adjusted.P.value <= .05]  <- "*"
bp_up$star[bp_up$Adjusted.P.value <= .01]  <- "**"
bp_up$star[bp_up$Adjusted.P.value <= .001] <- "***"
ggplot(bp_up, aes(x=reorder(Term,-log10(Adjusted.P.value)),y=-log10(Adjusted.P.value))) + 
  geom_bar(stat='identity',position="dodge", width=0.5, fill="orange") +
  geom_text(aes(label=star), colour="black", vjust=1.2,hjust=0.5, size=4,angle=90) +
  ggtitle("Up-regulated, Biological Process") +
  xlab("") +
  ylab("-log10(padj)") +
  theme_cowplot()+
  #scale_y_continuous(limits = c(1, 3), oob=rescale_none) +
  #ylab(expression(paste("Enrichr combined score"))) +
  coord_flip()+
  theme(plot.margin = margin(2,0.8,2,0.8, "cm"))
dev.off()

#MF pathway
mf_up<-read.csv("~/selected_up_MF_pathways.csv",sep=",",header=T,row.names=NULL)
head(mf_up)

pdf("~/top_mf_up_padj.pdf", height=6,width=9)
mf_up$star <- ""
mf_up$star[mf_up$Adjusted.P.value <= .05]  <- "*"
mf_up$star[mf_up$Adjusted.P.value <= .01]  <- "**"
mf_up$star[mf_up$Adjusted.P.value <= .001] <- "***"
ggplot(mf_up, aes(x=reorder(Term,-log10(Adjusted.P.value)),y=-log10(Adjusted.P.value))) + 
  geom_bar(stat='identity',position="dodge", width=0.5, fill="orange") +
  geom_text(aes(label=star), colour="black", vjust=1.2,hjust=0.5, size=4,angle=90) +
  ggtitle("Up-regulated, Molecular Function") +
  xlab("") +
  ylab("-log10(padj)") +
  theme_cowplot()+
  #scale_y_continuous(limits = c(1, 3), oob=rescale_none) +
  #ylab(expression(paste("Enrichr combined score"))) +
  coord_flip()+
  theme(plot.margin = margin(2,0.8,2,0.8, "cm"))
dev.off()

#Bio Planet pathway
bio_up<-read.csv("~/selected_up_Bioplanet_pathways.csv",sep=",",header=T,row.names=NULL)

pdf("~/top_bio_up_padj.pdf", height=6,width=12)
bio_up$star <- ""
bio_up$star[bio_up$Adjusted.P.value <= .05]  <- "*"
bio_up$star[bio_up$Adjusted.P.value <= .01]  <- "**"
bio_up$star[bio_up$Adjusted.P.value <= .001] <- "***"
ggplot(bio_up, aes(x=reorder(Term,-log10(Adjusted.P.value)),y=-log10(Adjusted.P.value))) + 
  geom_bar(stat='identity',position="dodge", width=0.5, fill="darkblue") +
  geom_text(aes(label=star), colour="black", vjust=1.2,hjust=0.5, size=4,angle=90) +
  ggtitle("Up-regulated, BioPlanet 2019") +
  xlab("") +
  ylab("-log10(padj)") +
  theme_cowplot()+
  #scale_y_continuous(limits = c(1, 3), oob=rescale_none) +
  #ylab(expression(paste("Enrichr combined score"))) +
  coord_flip()+
  theme(plot.margin = margin(2,0.8,2,0.8, "cm"))
dev.off()

## KEGG
kegg_up<-read.csv("~/selected_up_Keg_pathways.csv",sep=",",header=T,row.names=NULL)

pdf("~/top_kegg_up_padj.pdf", height=6,width=12)
kegg_up$star <- ""
kegg_up$star[kegg_up$Adjusted.P.value <= .05]  <- "*"
kegg_up$star[kegg_up$Adjusted.P.value <= .01]  <- "**"
kegg_up$star[kegg_up$Adjusted.P.value <= .001] <- "***"
ggplot(kegg_up, aes(x=reorder(Term,-log10(Adjusted.P.value)),y=-log10(Adjusted.P.value))) + 
  geom_bar(stat='identity',position="dodge", width=0.5, fill="darkblue") +
  geom_text(aes(label=star), colour="black", vjust=1.2,hjust=0.5, size=6,angle=90) +
  ggtitle("Up-regulated, KEGG 2021") +
  xlab("") +
  ylab("-log10(padj)") +
  theme_cowplot()+
  #scale_y_continuous(limits = c(1, 3), oob=rescale_none) +
  #ylab(expression(paste("Enrichr combined score"))) +
  coord_flip()+
  theme(plot.margin = margin(2,0.8,2,0.8, "cm"))
dev.off()



