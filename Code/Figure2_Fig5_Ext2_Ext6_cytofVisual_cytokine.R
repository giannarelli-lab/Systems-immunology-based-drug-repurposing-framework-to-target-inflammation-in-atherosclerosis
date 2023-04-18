# Vizualization of CyTOF data

setwd("~/Google Drive/projects/athero-chiara")
rm(list=ls())

data_dir = "/Users/sk/DataProjects/athero-chiara"


library(data.table)
library(sva) 
library(RColorBrewer)

# Load heatmap.3
library(devtools)
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

# Load 2 CyTOF experiments  
setwd("~/Simon_study_updates/")

 #cytof1 = fread("CGI002_HD_CV_phospho_summary/151030_CGI002_Mar-04-2016.csv")
 cytof1 = fread("CGI002_HD_CV_phospho_summary/151030-Table 1.csv")
 cytof1 = as.data.frame(cytof1)

 #cytof2 = fread("/151120_CGI002b_PM.csv")
 cytof2 = fread("CGI002_HD_CV_phospho_summary/160217-Table 1.csv")
 cytof2 = as.data.frame(cytof2)
data_dir="~/Simon_study_updates/"
cytof = fread(file.path(data_dir, "CGI002_HD_CV_phospho_summary/Combined-Table 1.csv"))
cytof = as.data.frame(cytof)


# Format data
cytof$Experiment = as.character(cytof$Experiment)
cytof$SampleType[cytof$SampleType == "CD4 T Cells"] = "CD4 T cells"


cytof_melt = melt(cytof)  # Default sample ids works
colnames(cytof_melt)[colnames(cytof_melt) == "variable"] = "Phos"  # specifying the phosphosite column



# Test for monocyte signal, which was observed in other CyTOF data.
 row_ids = which(cytof1_row_info$SampleType == "CD14 monocytes")
 row_ids = which(cytof1_row_info$SampleType == "NK cells")
 x_sub = x1[row_ids,]

 rownames(x_sub) = cytof1_row_info$Individuals[row_ids]

# # Combined data set
 cytof_mat = cytof[,4:ncol(cytof)]

# # cell_types 
 cytof_row_info = cytof[,1:3]

# # Correct capitalization
 cytof_row_info$SampleType[cytof_row_info$SampleType == "CD4 T Cells"] = "CD4 T cells"

 phos = unique(c(colnames(x1), colnames(x2)))


# phos = unique(colnames(cytof_mat))
phos = c("pSTAT1", "pSTAT3", "pSTAT5", "pp38", "pMAPKAP2", "pERK1/2", "pPLCg2", "IkBa", "pCREB", "pS6")

# cell_types = unique(cytof_row_info$SampleType)
cell_types = c("CD1c DCs", "CD14 monocytes", "CD16 monocytes", "pDCs", "B cells", "CD4 T cells", "CD8 T cells", "NK cells", "NKT cells", "Basophils")



# id = which(cytof_row_info$SampleType == "CD14 monocytes")
# hist(cytof_mat[id,8])

# plot(cytof_mat[id, 8], cytof_mat[id, 10])
# text(cytof_mat[id, 8], cytof_mat[id, 10]+2, labels=cytof_row_info$Individuals[id], cex=0.5)

# Reshape CyTOF data matrix
# -------------------------------------------
# make row ids of reshaped feature matrix

# cytof_cast = cast(cytof_melt, Individuals ~ SampleType + Phos)  # SampleType is cell line
cytof_cast = cast(cytof_melt, Experiment + Individuals ~ SampleType + Phos)  # SampleType is cell line
numeric_cols = 3:ncol(cytof_cast)

x = cytof_cast[, 3:ncol(cytof_cast)]




# Standardize per phosphosite over all 
 cytof_cast[, numeric_columns] = scale(cytof_cast[, numeric_columns])



 cell_types = unique(c(cytof1_row_info$SampleType, cytof2_row_info$SampleType))
 # Remove cell types
 cell_types = cell_types[cell_types != "Total"]

# Phosphosite ids
 phos = unique(c(colnames(x1), colnames(x2)))

 row_info = expand.grid(phos, cell_types)  # make all pairwise combinations
 colnames(row_info) = c("phos", "cell")

# # Individuals/patient samples
 exp_ids = unique(c(cytof2_row_info$Individuals, cytof1_row_info$Individuals))

 col_info = data.frame(
 	sample=unique(cytof_row_info$Individuals),
 	experiment=factor(
 			cytof_row_info$Experiment[
 			match(unique(cytof_row_info$Individuals), cytof_row_info$Individuals)
 		]
 	)
 )


 x = matrix(NA, nrow=nrow(row_info), ncol=nrow(col_info))

 rownames(x) = paste(row_info$cell, row_info$phos, sep="+")
 colnames(x) = col_info$sample


 for (i in 1:nrow(cytof_mat)) {
 	for (j in 1:ncol(cytof_mat)) {
 		# find the conditions and location in combined matrix
 		sample = cytof_row_info$Individuals[i]
 		cell = cytof_row_info$SampleType[i]
 		phos = colnames(cytof_mat)[j]

# 		# Find target position
 		target_j = which(col_info$sample == sample)

 		target_i = which(
 			row_info$phos == phos &
 			row_info$cell == cell
 		)

# 		# Retrieve data value from experiment 1
 		value = cytof_mat[i, j]

# 		# insert data into target
 		x[target_i, target_j] = value
 	}
 }


# Empirical Bayes correction for batch effect
mat = cytof_cast[, numeric_cols] # Add a small amount of noise to zero variance features

epsilon = 0.0000001
zero_var = which(apply(mat, 2, var) < epsilon)
mat[, zero_var] = mat[, zero_var] + rnorm(nrow(mat), 0, epsilon)

x_batch = ComBat(x, batch=factor(col_info$experiment))
mat_batch = t(ComBat(t(mat), batch=factor(cytof_cast$Experiment)))
colnames(mat_batch) = colnames(mat)

cytof_batch_cast = cytof_cast
cytof_batch_cast[, numeric_cols] = mat_batch

# cytof_batch_melt = melt(cytof_batch_cast)

# Standardize each phosphosite across cell type
cytof_batch_cell = cast(melt(cytof_batch_cast), Individuals + SampleType  ~ Phos)
numeric_cols = 3:ncol(cytof_batch_cell)

cytof_batch_cell[, numeric_cols] = scale(cytof_batch_cell[, numeric_cols])

x_batch_phos_zscore = cast(melt(cytof_batch_cell), Individuals ~ SampleType + Phos)


# CyTOF measurements highlighted
n1 = 1
n2 = 13

which(apply(x, 1, mean) < 30 & apply(x, 1, mean) > 10)


n = ncol(x)


# Principal component analysis
# x_imp = x
# x_imp[is.na(x_imp)] = 0.0
# pca = prcomp(t(x_imp))
pca = prcomp(t(x))
pca_batch = prcomp(t(x_batch))

# x_imp = x_batch_cor
# x_imp[is.na(x_imp)] = 0.0
# pca = prcomp(t(x_imp))

col_vector = rep(colors[1], nrow(col_info))
col_vector[col_info$sample %in% paste0("HD", 1:10)] = colors[2]


pdf("~/CyTOF_pool_pca.pdf", height=4)
par(mfrow=c(1, 2))
plot(pca$x[,1], pca$x[,2], pch=16, col=col_vector,
	main="CyTOF",
	xlab="PC1", ylab="PC2")
# text(pca$x[,1], pca$x[,2] + 1.5, label=col_info$sample, cex=0.3)

plot(pca_batch$x[,1], pca_batch$x[,2], pch=16, col=col_vector,
	main="CyTOF batch correction",
	xlab="PC1", ylab="PC2")
# text(pca_batch$x[,1], pca_batch$x[,2] + 1.5, label=col_info$sample, cex=0.3)
dev.off()


x_zscore = t(scale(t(x)))

# Separate z-score within batch
x_zscore = cbind(
	t(scale(t(x[, col_info$experiment == "151030"]))),
	t(scale(t(x[, col_info$experiment == "160217"])))
)

# x_zscore = t(scale(t(x[, !col_info$sample %in% c("IFNa", "PMA")])))

# x_zscore_batch = t(scale(t(x_batch_cor[, !col_info$sample %in% c("IFNa", "PMA")])))
x_zscore_batch = t(scale(t(x_batch)))


x_phos_zscore_batch = matrix(NA, nrow=nrow(x_batch), ncol=ncol(x_batch))
rownames(x_phos_zscore_batch) = rownames(x_batch)
colnames(x_phos_zscore_batch) = colnames(x_batch)

for (phos in unique(row_info$phos)) {
	row_ids = which(row_info$phos == phos)

	x_phos_zscore_batch[row_ids,] = (x[row_ids,] - mean(x[row_ids,])) / sd(x[row_ids,])
}


# # Check that reshaped data is consistent
# i = 19
# j = 12

# x[i, j]

# col_info[j,]
# row_info[i,]


# dist(t(x_zscore))

# cols = brewer.pal(12, "Paired")
# cols2 = brewer.pal(8, "Set2")
# cols = c(cols1, cols2)



# rownames(cytof_mat) = paste(cytof_row_info$Individuals, cytof_row_info$SampleType, sep=":")
# cytof_zscore = t(scale(t(cytof_mat)))
# heatmap.3(
# 	# as.matrix(cytof_mat),
# 	cytof_zscore,
# 	col=colorRampPalette(rev(brewer.pal(9, "RdBu")), interpolate="spline")(100),
# 	cexRow=0.2)

# x_zscore_batch_sub 

# pdf("plots/CyTOF_pool_heatmap.pdf", width=5, height=10)


# Test if some data is excluded
if (!all(1:ncol(x_zscore_batch) %in% col_order)) {
	stop("Not all data columns included")
}

t_test = apply(x_batch, 1, function(row) {
	hd = col_info$sample %in% paste0("HD", 1:10)

	return(t.test(row[!hd], row[hd]))
})

pval = sapply(t_test, function(test) {
	return(test$p.value)
})


# mat = x_zscore_batch[, col_order]


# mat = x_batch[, col_order]
# mat[mat < 0] = 0

# mat = x_batch[, col_order]
# mat[mat < 1] = 1
# mat = log2(mat)


# mat = x_phos_zscore_batch[, col_order]

mat = x_batch_phos_zscore[, 2:ncol(x_batch_phos_zscore)]
mat = t(mat)
colnames(mat) = x_batch_phos_zscore$Individuals
rownames(mat) = colnames(x_batch_phos_zscore)[2:ncol(x_batch_phos_zscore)]

# REMOVE THIS 
col_order_names = c(paste0("HD", 1:10), c("CV7004", "CV7007", "CV7009", "CV7012", "CV7013", "CV7015", "CV7016", "CV7017", "CV7019", "CV7036", "CV7095", "CV7096", "CV7099", "CV7100", "CV7101", "CV7102", "CV7106", "CV7107", "CV7110", "CV7111"))
col_order = match(col_order_names, colnames(mat))

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

cell_types = sapply(strsplit(rownames(mat), "_"), function(x) x[1])
cell_types = factor(cell_types, levels=c(
	"CD1c DCs",
	"CD14 monocytes",
	"CD16 monocytes",
	"pDCs",
	"B cells",
	"CD4 T cells",
	"CD8 T cells",
	"NK cells",
	"NKT cells",
	"Basophils"
	)
)

rlab = rbind(
	# "Phosphosite"=cols[as.integer(row_info$phos)],
	# "Cell type"=cell_palette[as.integer(row_info$cell)]
	"Cell type"=cell_palette[as.integer(cell_types)]
)


# col_scale = colorRampPalette(rev(brewer.pal(9, "RdBu")), interpolate="spline")(100)

# col_scale = colorRampPalette(rev(brewer.pal(8, "Set1")), interpolate="spline")(100)
col_scale = colorRampPalette(rev(brewer.pal(11, "Spectral")), interpolate="spline")(100)
# col_scale = rev(col_scale)
# col_scale = rainbow(100)

pdf("plots/CyTOF_pool_batch_heatmap8.pdf", width=6, height=10)
heatmap.3(
	mat[, col_order],
	Rowv=FALSE, 
	Colv=FALSE,
	trace="none",
	main="CyTOF pooled",
	xlab="Sample",
	ylab="Aggregate cell type CyTOF",
	margins=c(6, 12),
	dendrogram="none",
	cexRow=0.5,
	cexCol=0.5,
	colsep=c(10),
	# rowsep=(1:11) * 11,
	# sepcolor=rgb(0.7, 0.7, 0.7),
	sepcolor="white",
	RowSideColors=rlab,
	na.color=rgb(0.85, 0.85, 0.85),
	# ColSideColors=t(go_clab),
	# ColSideColorsSize=12,
	# labCol=FALSE,
	# rep(NULL, length(all_cids_unique)),
	col=col_scale,
	# KeyValueName="Within-experiment z-score"
	KeyValueName="z-score"
)

legend(
	"topright",
	legend=c(
		# levels(row_info$phos),
		# "",
		levels(cell_types)
	),
	fill=c(
		# cols[1:length(levels(row_info$phos))],
		# "white",
		cell_palette[1:length(levels(row_info$cell))]
	),
	border=FALSE,
	bty="n",
	y.intersp = 0.7,
	cex=0.7
)
dev.off()


# -----------------------------------------------

image(x2)

x1_zscore = scale(x1)

image(x1_zscore)
heatmap(x1_zscore)


##### Cytokine release and plasma data. Figure 2E
rm(list=ls())

library(data.table)
library(RColorBrewer)
library(sva)  # For ComBat implementation

library(dendextend)

# Load heatmap.3
library(devtools)
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

setwd("~/Google Drive/projects/athero-chiara")
data_dir = "/Users/sk/DataProjects/athero-chiara"

# Test if numerical values are inside a given range
inRange = function(x, interval) {
  return(x >= interval[1] & x <= interval[2])
}

# inRange(c(1, 3, 11), range(0, 10))

# Test which range a vector of numbers belongs to.
# Returns first match
# ranges: List of range objects
whichRange = function(x, ranges) {
  sapply(x, function(num) {
    for (i in 1:length(ranges)) {
      if (inRange(num, ranges[[i]])) {
        return(i)
      }
    }
    return(NA)
  })
}


# Load data
release = fread(
  file.path(data_dir, "MultiplexSiimon/net release_reshape-Table 1.csv"))

release = as.data.frame(release)

# Exclude experiment 1 due to bad quality
release = release[release$experiment != 1,]

release_row_info = release[,1:2]

release_row_info$status = NA
release_row_info$status[1:10] = "Disease"
release_row_info$status[11:19] = "Healthy"


release = release[,3:ncol(release)]  # main data matrix
rownames(release) = release_row_info$sample

plasma = fread(
  file.path(data_dir, "MultiplexSiimon/Plasma_reshape-Table 1.csv"))
plasma = as.data.frame(plasma)

plasma_row_info = plasma[,1:2]  

plasma_row_info$status = NA
plasma_row_info$status[c(1:6, 17:26)] = "Healthy"
plasma_row_info$status[c(7:16, 27:36)] = "Disease"

plasma = plasma[,3:ncol(plasma)]
rownames(plasma) = plasma_row_info$sample

# Remove outliers
plasma_rids = which(!rownames(plasma) %in% c("HD10"))
plasma = plasma[plasma_rids,]

plasma_row_info = plasma_row_info[plasma_rids,]



# Analyze data
# Cytokine release data
release_zscore = scale(release)
release_zscore[is.na(release_zscore)] = 0.0

# release_norm = scale(release, scale=FALSE)  # centered data
# release_norm[is.na(release_norm)] = 0.0

# Cytokine plasma data
# Plasma batch correction
plasma_batch = t(ComBat(t(plasma), batch=factor(plasma_row_info$experiment)))

# Standardize with minimum variation
plasma_batch_zscore = scale(plasma_batch)


# Errors may not be critical
t_test_plasma = apply(plasma_batch, 2, function(col) {
  # hd = release_row_info$sample %in% paste0("HD", 1:10)
  neg_ctrl = plasma_row_info$status == "Healthy"
  
  try({
    return(t.test(col[!neg_ctrl], col[neg_ctrl]))
  })
})

pval_plasma = sapply(t_test_plasma, function(test) {
  if (! "p.value" %in% names(test)) {
    return(NA)
  } else if (is.null(test$p.value)) {
    return(NA)
  } else {
    return(test$p.value)
  }
})


t_test_plasma_nobatch = apply(plasma, 2, function(col) {
  # hd = release_row_info$sample %in% paste0("HD", 1:10)
  neg_ctrl = plasma_row_info$status == "Healthy"
  
  try({
    return(t.test(col[!neg_ctrl], col[neg_ctrl]))
  })
})

pval_plasma_nobatch = sapply(t_test_plasma_nobatch, function(test) {
  if (! "p.value" %in% names(test)) {
    return(NA)
  } else if (is.null(test$p.value)) {
    return(NA)
  } else {
    return(test$p.value)
  }
})



# Reset unstable z-scores due to small variation
min_sd = 0.1  # minimum variation to be considered
plasma_cids_all = attr(plasma_batch_zscore, "scaled:scale") >= min_sd
plasma_batch_zscore_all = plasma_batch_zscore[, plasma_cids_all]

# plasma_zscore = scale(plasma)
# plasma_zscore[is.na(plasma_zscore)] = 0.0

# Filter based on p-values
min_pval = 0.1
plasma_cids_sig = which(pval_plasma < min_pval)

plasma_batch_zscore_sig = plasma_batch_zscore[, plasma_cids_sig]

# PCA
# -------------------------------------
x = plasma
x[is.na(x)] = 0.0
pca = prcomp(x)

colors = brewer.pal(9, "Set1")

plot(pca$x[,1], pca$x[,2], col=colors[release_row_info$experiment])
text(pca$x[,1], pca$x[,2], label=rownames(x), cex=0.4)


# release

healthy_means = apply(
  release[release_row_info$status == "Healthy",],
  2,
  mean, na.rm=TRUE
)
healthy_means[is.na(healthy_means)] = 0.0

release_fc = release

release_fc[release_row_info$status == "Disease",] =
  release[release_row_info$status == "Healthy",] %*% diag(healthy_means)	


release[is.na(release)] = 0.0
t_test = apply(release, 2, function(col) {
  # hd = release_row_info$sample %in% paste0("HD", 1:10)
  neg_ctrl = release_row_info$status == "Healthy"
  
  return(t.test(col[!neg_ctrl], col[neg_ctrl]))
})

pval = sapply(t_test, function(test) {
  return(test$p.value)
})

pval[is.na(pval)] = 1.0


# Cytokine dotplot
# --------------------------
# Order based on t-test
mat = release[,order(pval)]

ordered_pval = pval[order(pval)]

# Exlude zero entries
col_means = apply(mat, 2, mean, na.rm=TRUE)

mat = mat[,abs(col_means) > 0.001]


# Reorder disease status

release_row_info$status[1:10] = "Disease"
release_row_info$status[11:19] = "Healthy"

row_order = c(11:19, 1:10)
mat = mat[row_order,]


# mat = plasma_batch[, colnames(plasma_batch) %in% cytokine_sel]
# mat = scale(mat)
# plot(0, 0,
# 	xlim=c(1, ncol(mat)),
# 	ylim=range(mat),
# 	type="n")

# for (i in 1:nrow(mat)) {
# 	lines(mat[i, ],
# 		col=colors[as.integer(factor(plasma_row_info$status))[i]]
# 	)
# }


# col_scale = colorRampPalette(rev(brewer.pal(9, "RdBu")), interpolate="spline")(100)

mat = release
mat[is.na(mat)] = 0.0
# mat = sign(mat) * abs(mat)^(1/3)


sample_cols = data.frame(Status=colors[as.integer(factor(release_row_info$status))])
sample_cols = as.matrix(sample_cols)

healthy_means = apply(release[release_row_info$status == "Healthy",], 2, mean, na.rm=TRUE)
disease_means = apply(release[release_row_info$status == "Disease",], 2, mean, na.rm=TRUE)


ranges = list(
  range(1000, 10000),
  range(100, 1000),
  range(10, 100),
  range(0, 10),
  range(-10, 0),
  range(-100, -10),
  range(-1000, -100),
  range(-10000, - 1000)
)

col_interval_groups = colorRampPalette(rev(brewer.pal(11, "PuOr")), interpolate="spline")(length(ranges))

row_cols = data.frame(
  Healthy=col_interval_groups[whichRange(healthy_means, ranges)],
  Disease=col_interval_groups[whichRange(disease_means, ranges)]
)

# release

pdf("plots/cytokines/release_heatmap.pdf", width=4.5, height=4.5)
col_scale = colorRampPalette(rev(brewer.pal(11, "Spectral")), interpolate="spline")(100)
heatmap.3(
  t(release_zscore),
  # t(release_norm),
  # t(mat),
  # plasma_zscore,
  # release,
  # x_zscore,
  # x_zscore_batch,
  # x_batch,
  # Rowv=FALSE, 
  # Colv=FALSE,
  trace="none",
  # distfun=function(x) {
  # 	dmat = dist(x)
  # 	dmat[is.na(dmat)] = max(dmat, na.rm=TRUE) + 1
  # 	# dmat = dist(x, method="minkowski", p=1.5)
  # 	return(dmat)
  # },
  main="Release",
  # cex.main=0.5,
  cex.lab=0.5,
  xlab="Patient sample",
  ylab="Cytokines",
  margins=c(6, 12),
  dendrogram="none",
  cexRow=0.5,
  cexCol=0.5,
  # colsep=c(10),
  # rowsep=(1:11) * 11,
  # sepcolor=rgb(0.7, 0.7, 0.7),
  sepcolor="white",
  # RowSideColors=rlab,
  na.color=rgb(0.85, 0.85, 0.85),
  ColSideColors=sample_cols,
  RowSideColors=t(row_cols),
  # ColSideColorsSize=12,
  # rep(NULL, length(all_cids_unique)),
  col=col_scale,
  # KeyValueName="Within-experiment z-score"
  # KeyValueName="z-score"
  KeyValueName="z-score"
)
legend(
  "bottomleft",
  # -1, -1,
  legend=c(
    sapply(ranges, paste, collapse=", ")
  ),
  fill=c(
    # cols[1:length(levels(row_info$phos))],
    # "white",
    # cols[1:length(levels(row_info$cell))]
    col_interval_groups
  ),
  # legend=c("Basal","LumA","LumB","Her2","Claudin","Normal","","Positive","Negative","NA","","Targeted","Chemo","","Approved","Experimental"),
  # fill=c("red","blue","cyan","pink","yellow","green","white","black","white","grey","white","darkorchid","darkred","white","green","darkgreen"),
  border=FALSE,
  bty="n",
  y.intersp = 0.7,
  cex=0.3
)
dev.off()

# --------------------------------------------------------------------------------------

ranges = list(
  range(1000, 30000),
  range(100, 1000),
  range(10, 100),
  range(0, 10),
  range(-10, 0),
  range(-100, -10),
  range(-1000, -100),
  range(-30000, - 1000)
)
# whichRange(c(3, 11, -5000, 0), ranges)

# Plot selection
pdf("plots/cytokines/plasma_heatmap_all.pdf", height=6.0)
plasma_cids = plasma_cids_all
mat = plasma_batch_zscore_all


# Calculate side colors
col_interval_groups = colorRampPalette(rev(brewer.pal(11, "PuOr")), interpolate="spline")(length(ranges))
colors = brewer.pal(9, "Set1")

sample_cols = data.frame(Status=colors[as.integer(factor(plasma_row_info$status))])
sample_cols = as.matrix(sample_cols)

healthy_means = apply(plasma[plasma_row_info$status == "Healthy", plasma_cids], 2, mean, na.rm=TRUE)
disease_means = apply(plasma[plasma_row_info$status == "Disease", plasma_cids], 2, mean, na.rm=TRUE)

row_cols = data.frame(
  Healthy=col_interval_groups[whichRange(healthy_means, ranges)],
  Disease=col_interval_groups[whichRange(disease_means, ranges)]
)

col_scale = colorRampPalette(rev(brewer.pal(11, "Spectral")), interpolate="spline")(100)
heatmap.3(
  t(mat),
  # Rowv=FALSE, 
  # Colv=FALSE,
  trace="none",
  distfun=function(x) {
    # dmat = dist(x)
    # dmat[is.na(dmat)] = max(dmat, na.rm=TRUE) + 1
    dmat = dist(x, method="minkowski", p=1.5)
    # dmat = dist(x, method="minkowski", p=1.0)
    # dmat = dist(x, method="minkowski", p=1.2)
    return(dmat)
  },
  main="Plasma",
  # cex.main=0.5,
  cex.lab=0.5,
  xlab="Patient sample",
  ylab="Cytokines",
  margins=c(6, 12),
  dendrogram="column",
  # cexRow=0.5,
  # cexCol=0.5,
  # colsep=c(10),
  # rowsep=(1:11) * 11,
  # sepcolor=rgb(0.7, 0.7, 0.7),
  sepcolor="white",
  # RowSideColors=rlab,
  na.color=rgb(0.85, 0.85, 0.85),
  ColSideColors=sample_cols,
  RowSideColors=t(row_cols),
  # ColSideColorsSize=12,
  # rep(NULL, length(all_cids_unique)),
  col=col_scale,
  # KeyValueName="Within-experiment z-score"
  # KeyValueName="z-score"
  KeyValueName="z-score"
)
legend(
  "bottomleft",
  # -1, -1,
  legend=c(
    sapply(ranges, paste, collapse=", ")
  ),
  fill=c(
    # cols[1:length(levels(row_info$phos))],
    # "white",
    # cols[1:length(levels(row_info$cell))]
    col_interval_groups
  ),
  # legend=c("Basal","LumA","LumB","Her2","Claudin","Normal","","Positive","Negative","NA","","Targeted","Chemo","","Approved","Experimental"),
  # fill=c("red","blue","cyan","pink","yellow","green","white","black","white","grey","white","darkorchid","darkred","white","green","darkgreen"),
  border=FALSE,
  bty="n",
  y.intersp = 0.7,
  cex=0.5
)
dev.off()


#### Figure 2F&G Cytokine expression 
library(dplyr)  
library(data.table)
library(ggplot2)
library(reshape2)
library(magrittr)
library(tidyverse)
library(ggplot2)
library(ggpubr)

# file with clinical information, Luminex data
data <- read.csv("~/Multiplex/counts_data.csv",header=T)#,fill=TRUE)

dim(data)

# Create vector containing cytokine names
name <- (names(data[,3:43]))

#log normalize data
norm_data <- data
norm_data[name] <- log(norm_data[name]+1)

#write.csv(norm_data,"~/Multiplex/norm_data_drug.csv")

norm_data$group<- rep(c("Athero","Healthy","Athero"), c(10,15,10))

#compare two groups
type1 <- "Athero"
type2 <- "Healthy"

df.test <- rbind(norm_data[which(norm_data$group == type1),], 
                 norm_data[which(norm_data$group == type2),])

matab_num <- ncol(df.test)

test_results <- data.frame(cytokines = colnames(df.test)[3:matab_num])

# T-test
test_results$ttest <- apply(df.test[, 3:matab_num], 2,
                            function(x) unlist(t.test(as.numeric(x) ~ df.test$group, data = df.test)[3]))

# adjust p value using BH method
test_results$ttest_BH <- p.adjust(test_results$ttest, method = "BH")

# Log2 (FC)
cn <- paste("FC_", type1, "/", type2, sep = "")
test_results[cn] <- apply(df.test[,3:matab_num], 2, 
                          function(x) 
                            mean(as.numeric(x[which(norm_data$group == type1)]))/
                            mean(as.numeric(x[which(norm_data$group == type2)])))
test_results$LOG2FC <- log2(test_results$`FC_Athero/Healthy`)

#point plot
cyto.bind<-test_results
cyto.bind$type <- ifelse(cyto.bind$ttest< 0.05, "Significant","Not Significant")
cyto.bind$type <- factor(cyto.bind$type, levels = c("Significant","Not Significant"))


# save the cytokines
write.csv(cyto.bind,"~/Multiplex/pointplot_athero_healthy.csv")


# point plot
library(ggdendro)
p<-ggplot(cyto.mhu.bind, aes(group, cytokines)) +
  geom_point(aes(fill = LOG2FC, size = -log10(ttest), color = type, shape = type)) + 
  scale_fill_gradientn(colours = colorRampPalette(c("#3d67a3", "white","#ce1020")),
                       limits = c(-5,5), breaks = c(-5, -3, 0, 3, 5)) +
  scale_shape_manual(values = c(7,23,21)) +
  scale_color_manual(values = c("red","black","black"))+#, "grey50")) +
  scale_size_continuous(range = c(0,5), breaks = c(0.5,1,2,2.5)) +
  theme_dendro() +
  theme(plot.margin = margin(5, 3, 5, 3, "cm")) +
  coord_fixed(ratio = 0.7) +
  labs(fill = "Log2 (Fold change)", size = "-Log10 (P value)", color = "", shape = "") +
  theme(axis.text.x = element_text(size = 10, face = "plain", colour = "black", angle = 60,vjust = 1.1, hjust=1.1), 
        axis.text.y = element_text(size = 10, face = "plain", colour = "black")) +
  theme(legend.text = element_text(size = 10, face = "plain", colour = "black"), 
        legend.title = element_text(size = 10, face = "plain", colour = "black"), 
        legend.key.height = unit(0.3, "cm"), legend.key.width = unit(0.3, "cm"))   
pdf(file="/Multiplex/pointplot_signif_athero_in_drug.pdf",height = 8,width = 8)
p
dev.off()

#barplot
pdf(paste('/Multiplex/','barplot_atherovshealthy_color.pdf', sep=''), width = 4, height = 4)
cols<-c("#ce1020","#3d67a3")
for(i in name){
  v<-ggbarplot(norm_data, x = "group", y = i,add = c("mean_se","jitter"),
               width = 0.4, color ="group", palette = c("#ce1020","#3d67a3"),
               position = position_dodge(0.8))+
    theme(aspect.ratio = 1.5)+
    labs(y = paste0(i,"(pg/mL)"))
  print(v)
}
dev.off()

