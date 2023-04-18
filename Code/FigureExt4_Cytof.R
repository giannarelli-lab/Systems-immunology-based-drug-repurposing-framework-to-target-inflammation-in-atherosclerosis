setwd("~/Desktop/athero")

rm(list=ls())

library(DESeq)
library(RColorBrewer)

setwd("~/DataProjects/athero-chiara")
source("/Users/sk/Google Drive/projects/athero-chiara/loadData.r")

setwd("~/Google Drive/projects/athero-chiara")
source("func.r")




# DESeq analysis of RNA-seq data
# --------------------------------------------------------------------------------
cds = newCountDataSet(counts, sample_info$condition)  # DESeq count data set

# cds = newCountDataSet(counts, sample_info$phenotype)  # DESeq count data set

# Size factors
cds = estimateSizeFactors(cds)

rna_norm = counts(cds, normalized=TRUE)  # normalized by size factor
rownames(rna_norm) = gene_symbols

x_zscore = t(scale(t(rna_norm)))  # standardized data

x_log = log2(rna_norm + 1)  # log transform

cds = estimateDispersions(cds, sharingMode="maximum")

# Fit information
# str(fitInfo(cds))
# plotDispEsts(cds)

# Negative binomial test for control vs disease.
res = nbinomTest(cds, "control", "disease")
res$gene = gene_symbols

head(res[order(res$padj), c("gene", "log2FoldChange", "padj")], n=30)

write.table(
	res[order(res$padj)[1:sum(res$padj < 0.05)],],
	file="results/athero_healthy_sig_transcripts.csv",
	sep=",",
	row.names=FALSE,
	quote=FALSE)



# CyTOF t-test
# ------------------------------------------------------
cytof_ttest_athero = apply(cytof, 1, function(row) {
	t.test(
		row[cytof_sample_info$status == "DISEASE_STABLE" | cytof_sample_info$status == "DISEASE_UNSTABLE"],
		row[cytof_sample_info$status == "HEALTHY"]
		)
})

cytof_athero_tstat = sapply(cytof_ttest_athero, function(ttest) {
	ttest$statistic
})

cytof_athero_pval = sapply(cytof_ttest_athero, function(ttest) {
	ttest$p.value
})

# cytof_athero_tstat[cytof_athero_pval < 0.01]


# Read gene list filter from GO analysis
# --------------------------------------------------

file_names = list.files("go-all-lists")
go_lib = "go-all-lists"
out_folder = "plots/rna-seq-z-scores/all"

color_scale = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)  # Color scale for heatmaps

# Read the gene filters from the GO analysis
gene_filters = list()
for (name in file_names) {
	gene_filters[[name]] = fread(paste0(go_lib, "/", name))
	gene_filters[[name]] = as.data.frame(gene_filters[[name]])
}


# Average over normalized RNA-seq counts over experimental repeats
# ----------------------------------------------------------------------
ave_rna_norm = matrix(NA, nrow=nrow(rna_norm), ncol=length(unique(sample_info$exp)))

col_order = match(1:max(sample_info$exp), sample_info$exp)
colnames(ave_rna_norm) = colnames(rna_norm)[col_order]
rownames(ave_rna_norm) = rownames(rna_norm)

# Calculate replicate means
for (exp in unique(sample_info$exp)) {
	ave_rna_norm[,exp] = apply(rna_norm[, sample_info$exp == exp], 1, mean, na.rm=TRUE)
}

sample_info_ave = sample_info[col_order,]


# Cross-correlations between RNA-seq and
# Assumes that the columns of CyTOF and RNA-seq matrices have been matched one-to-one.
cross_cor = cor(t(cytof), t(ave_rna_norm))


# get all rids associated with the gene filters from the GO analysis.
all_cids = list()
for (i in 1:length(gene_filters)) {
	genes = gene_filters[[i]][,2]
	all_cids[[i]] = match(genes, colnames(cross_cor))
}
all_cids_unique = unique(unlist(all_cids))

# Calculate FC matrix from GO analysis for significant
go_clab = matrix("white",
	ncol=length(all_cids_unique),
	nrow=length(gene_filters)
)
colnames(go_clab) = gene_symbols_no_iso[all_cids_unique]
rownames(go_clab) = sapply(names(gene_filters), function(file_name) {
	return(strsplit(file_name, "[.]")[[1]][1])
})

go_clab = rbind(go_clab, rep("white", length(all_cids_unique)))
rownames(go_clab)[nrow(go_clab)] = "FC"


# Fill in GO x gene color based on log2 FC values from DESeq analysis.
for (i in 1:length(gene_filters)) {
	genes = gene_filters[[i]][,2]

	found = match(genes, colnames(go_clab))

	global_ids = match(genes, gene_symbols_no_iso)

	go_clab[i, found] = "black"
}

go_clab[nrow(go_clab),] = colorGradient(
	x=res$log2FoldChange[all_cids_unique],
	gradlim=c(-1.0, 1.0),
	colors=brewer.pal(9, "PRGn"))


# Load heatmap.3
library(devtools)
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

cols = brewer.pal(12, "Paired")
# cols2 = brewer.pal(8, "Set2")
# cols = c(cols1, cols2)
rlab = rbind(
	"Phosphosite"=cols[as.integer(cytof_feature_info$phos)],
	"Cell type"=cols[as.integer(cytof_feature_info$cell)]
)


# phospho_ids = order(apply(abs(cross_cor), 1, mean), decreasing=TRUE)[1:16]

phospho_ids = which(cytof_athero_pval < 0.01)

rlab_sub = rlab[,phospho_ids]

pdf("plots/cross-cor3.pdf", width=12, height=10)

mat = cross_cor[phospho_ids, all_cids_unique]

heatmap.3(
	mat,
	trace="none",
	distfun=function(x) {
		# dmat = dist(abs(x))
		dmat = dist(x, method="minkowski", p=1.5)
		return(dmat)
	},
	main="RNA-seq CyTOF correlations",
	xlab="DEGs with GO enrichment",
	ylab="CyTOF",
	margins=c(6, 12),
	dendrogram="none",
	cexRow=1.0,
	# cexCol=0.5,
	RowSideColors=rlab_sub,
	ColSideColors=t(go_clab),
	ColSideColorsSize=12,
	labCol=FALSE,
	# rep(NULL, length(all_cids_unique)),
	col=colorRampPalette(rev(brewer.pal(9, "RdBu")), interpolate="spline")(100),
	KeyValueName="Pearson's correlation",
)
legend("bottomleft",
	legend=c(
		levels(cytof_feature_info$phos),
		"",
		levels(cytof_feature_info$cell)
	),
	fill=c(
		cols[1:length(levels(cytof_feature_info$phos))],
		"white",
		cols[1:length(levels(cytof_feature_info$cell))]
	),
	# legend=c("Basal","LumA","LumB","Her2","Claudin","Normal","","Positive","Negative","NA","","Targeted","Chemo","","Approved","Experimental"),
	# fill=c("red","blue","cyan","pink","yellow","green","white","black","white","grey","white","darkorchid","darkred","white","green","darkgreen"),
	border=FALSE,
	bty="n",
	y.intersp = 0.7,
	cex=0.7
)
dev.off()

phos_order = order(apply(abs(mat), 1, median), decreasing=TRUE)

pdf("plots/cross-cor_mean.pdf", width=6, height=5)
k = 10
par(mar=c(6, 10, 4, 4), bty="n")
boxplot(
	t(abs(mat[phos_order[k:1],])),
	las=1,
	horizontal=TRUE,
	col=cols[as.integer(cytof_feature_info$cell[phospho_ids][phos_order[k:1]])],
	xlab="Absolute CyTOF RNA-seq cross-correlation"
)
dev.off()
