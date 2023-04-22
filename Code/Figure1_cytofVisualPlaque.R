## Simon Koplev

setwd("~/Google Drive/projects/athero-chiara")
rm(list=ls())

library(data.table)
library(RColorBrewer)
library(reshape)
library(stringr)
library(qvalue)

# Load heatmap.3
library(devtools)
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

data_dir = "/Users/sk/DataProjects/athero-chiara"	


# Load and parse CyTOF gated median data
# --------------------------------------
cytof = fread(file.path(data_dir, "151120_CGI002b_PM.csv"))
cytof = as.data.frame(cytof)


# Exclude IFNa and PMA controls; and ungated event counts.
cytof = cytof[cytof$Conditions %in% c("healthy", "auto") , colnames(cytof) != "Ungated_EventCounts"]

# Split column names by (phos), otherwise return as is
colnames(cytof) = sapply(
	strsplit(colnames(cytof), "[()]"),
	function(x) tail(x, 1))  # returns last element in split vector

# Replace "_" with "/" in column names
colnames(cytof) = str_replace(colnames(cytof), "_", "/")

# Reshape CyTOF data into matrix of phos_cell x condition.
cytof_melt = melt(cytof)
colnames(cytof_melt)[colnames(cytof_melt) == "variable"] = "Phos"

cytof_cast = cast(cytof_melt, Conditions + Individuals ~ SampleType + Phos)


# Paired t-tests for effect of auto serum on plaque
# -------------------------------------------------

# Test if all individuals are paired
paired = cytof_cast$Individuals[cytof_cast$Conditions == "auto"]  == cytof_cast$Individuals[cytof_cast$Conditions == "healthy"]

if (!all(paired)) {
	stop("CyTOF data is not paired")
}

n_start = 3  # data starting column
t_tests = apply(cytof_cast[, n_start:ncol(cytof_cast)], 2, function(dcol) {
	x = dcol[cytof_cast$Conditions == "auto"]
	y = dcol[cytof_cast$Conditions == "healthy"]
	t_test = t.test(x, y, paired=TRUE)

	# log2 fold changes, auto vs healthy stimulation
	log2_fc = log2(x/y)
	# fc = x/y
	names(log2_fc) = cytof_cast$Individuals[cytof_cast$Conditions == "auto"]

	log2_fc[is.na(log2_fc)] = 0.0
	t_test$log2_fc = log2_fc
	return(t_test)
})
names(t_tests) = colnames(cytof_cast)[n_start:ncol(cytof_cast)]

pvals = sapply(t_tests, function(x) x$p.value)
pvals[is.na(pvals)] = 1.0
tvals = sapply(t_tests, function(x) x$statistic)
qvals = qvalue(pvals)$qvalue

fc = sapply(t_tests, function(x) x$log2_fc)
# fc = sapply(t_tests, function(x) x$fc)

results = data.frame(p=pvals, q=qvals, t=tvals, t(fc))

write.csv(results, "cytofPlaque/tables/paired_ttests_wiht_auto_healthy_FC.csv", quote=FALSE, row.names=TRUE)


# Include only cell types (not total)
include = sapply(strsplit(colnames(fc), "_"), function(x) x[1]) != "Total"

sig = qvals < 0.1 & include
# sig = qvals < 0.05 & include

# Subset fold change matrix
fc_sig = fc[, sig]
fc_sig[is.infinite(fc_sig)] = NA
colnames(fc_sig) = str_replace(colnames(fc_sig), "_", ":")

# fc_sig = (fc_sig - 1) * 100  # percentage increase

# fc_sig = fc_sig[, ]

# cell_palette = c("black", brewer.pal(8, "Dark2"))
# cell_palette = c(brewer.pal(8, "Dark2"), "black")
# cell_palette = c(brewer.pal(10, "Paired"), "black")

cell_palette = c(
	rgb(40, 119, 177, maxColorValue=255),
	rgb(252, 128, 37, maxColorValue=255),
	rgb(51, 158, 53, maxColorValue=255),
	rgb(212, 41, 49, maxColorValue=255),
	rgb(149, 104, 186, maxColorValue=255),
	rgb(140, 87, 74, maxColorValue=255),
	rgb(225, 122, 194, maxColorValue=255),
	rgb(128, 127, 127, maxColorValue=255),
	rgb(41, 191, 204, maxColorValue=255),
	rgb(188, 188, 53, maxColorValue=255)
)

cell_type = sapply(strsplit(colnames(fc), "_"), function(x) x[1])

cell_levels = c("CD1c DCs", "CD14 monocytes", "CD16 monocytes", "pDCs", "B cells", "CD4 T cells", "CD8 T cells", "NK cells", "NKT cells", "Basophils")

rlab = data.frame(
	"Cell type"=cell_palette[as.integer(
		factor(cell_type, levels=cell_levels)
	)]
)

rlab = as.matrix(rlab)
rlab[is.na(rlab)] = "grey"

rlab_sig = rlab[sig, , drop=FALSE]


	
pdf("cytofPlaque/plots/sig_cell_phosph_FDR01.pdf", width=3.0, height=3.5)
sig_cell_counts = table(sapply(strsplit(colnames(fc_sig), ":"), function(x) x[1]))
sig_cell_counts = rev(sort(sig_cell_counts))

bar_cols = cell_palette[as.integer(factor(names(sig_cell_counts), levels=cell_levels))]
par(mar=c(8, 4, 2, 2))
barplot(sig_cell_counts, las=2,
	ylab="phospshosites (n) (FDR < 10%)",
	col=bar_cols)
abline(h=0)
dev.off()

pdf("cytofPlaque/plots/sig_heatmap_fc.pdf", width=7)
# col_scale = colorRampPalette(rev(brewer.pal(9, "Spectral")), interpolate="spline")(100)
# col_scale = colorRampPalette(rev(brewer.pal(9, "RdBu")), interpolate="spline")(50)
col_scale = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)
heatmap.3(t(fc_sig),
	col=col_scale,
	mar=c(6, 26),
	RowSideColors=t(rlab_sig),
	# breaks=seq(-100, 100, length.out=101),
	# breaks=seq(-50, 50, length.out=101),
	breaks=seq(-1.0, 1.0, length.out=101),
	# ColSideColors=clab,
	cexRow=1.0, cexCol=1.0,
	KeyValueName=expression("log"[2] * " FC (auto/healthy)"),
	# KeyValueName="Auto serum induction (%)",
	main="",
	xlab="Patients", ylab="Cell type:phosphosite"
)

legend("topright", legend=cell_levels, col=cell_palette, pch=15)
dev.off()
