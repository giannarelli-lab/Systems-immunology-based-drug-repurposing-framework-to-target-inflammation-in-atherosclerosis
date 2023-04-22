# Perturbational experiment of healthy PBMC + drug + diseased plasma.
# CyTOF data, median signaling markers.

rm(list=ls())

library(data.table)
library(RColorBrewer)
# library(penalized)
# library(glmnet)
library(reshape)
library(stringi)


library(devtools)
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

svgText = function(file, width, height) {
  library(RSvgDevice)
  devSVG(file = file, width = width, height = height, bg = "white", fg = "black",
         onefile = TRUE, xmlHeader = TRUE)
}


setwd("~/DataProjects/athero-chiara")
# Load cyTOF data
cytof = fread("160720_CGI002_exported_populations_Exported_Stats_Jul-28-2016_07-07-PM.csv")
cytof = as.data.frame(cytof)

# Rename columns of cytof data frame. From->to colname remap
rename_cols = list(
	"Medians_Ungated_Sm147Di--(pSTAT5)"="pSTAT5",
	"Medians_Ungated_Eu153Di--(pSTAT1)"="pSTAT1",
	"Medians_Ungated_Gd156Di--(pp38)"="pp38",
	"Medians_Ungated_Gd158Di--(pSTAT3)"="pSTAT3",
	"Medians_Ungated_Tb159Di--(pMAPKAP2)"="pMAPKAP2",
	"Medians_Ungated_Dy162Di--(pPLCg2)"="pPLCg2",
	"Medians_Ungated_Dy164Di--(IkBa)"="IkBa",
	"Medians_Ungated_Ho165Di--(pCREB)"="pCREB",
	"Medians_Ungated_Er167Di--(pERK1_2)"="pERK1/2",
	"Medians_Ungated_Lu175Di--(pS6)"="pS6"
)

colnames(cytof)[
	match(names(rename_cols), colnames(cytof))  # indices of from names
] = unlist(rename_cols)


# Rename conditions
for (i in 1:8) {
	cytof$Conditions[cytof$Conditions == paste("DRUG", i)] = i
}

# Format data
cytof_melt = melt(cytof)  # Default sample ids works
colnames(cytof_melt)[colnames(cytof_melt) == "variable"] = "Phos"  # specifying the phosphosite column
cytof_cast = cast(cytof_melt, Conditions + Individuals ~ SampleType + Phos)  # SampleType is cell line

# Get matrix of numerical data in casted format
k = 3
cytof_cast_mat = as.matrix(cytof_cast[,k:ncol(cytof_cast)])
colnames(cytof_cast_mat) = colnames(cytof_cast)[k:ncol(cytof_cast)]
rownames(cytof_cast_mat) = cytof_cast$Conditions

# Groups of experiments
basal_ids = which(cytof_cast$Conditions == "NO PLASMA")  # no plasma, no drug
plasma_ids = which(cytof_cast$Conditions == "VEHICLE")  # plasma treated, no drug


# basal_means = apply(cytof_cast_mat[basal_ids,], 2, mean)
# basal_sd = apply(cytof_cast_mat[basal_ids,], 2, sd)

# plasma_means = apply(cytof_cast_mat[plasma_ids,], 2, mean)
# plasma_sd = apply(cytof_cast_mat[plasma_ids,], 2, sd)

# plot(plasma_means / basal_means, pch=16)

# data.frame(basal_means, plasma_means)


setwd("~/Google Drive/projects/athero-chiara")

# t-tests for plasma vs basal
t_tests = list()
for (i in 1:ncol(cytof_cast_mat)) {
	t_tests[[i]] = t.test(cytof_cast_mat[plasma_ids, i], cytof_cast_mat[basal_ids, i])
}

# t-tests for each drug vs plasma
t_tests_drugs = list()

for (condition in cytof_cast$Conditions) {
	if (condition == "NO PLASMA" | condition == "VEHICLE") next
	t_tests_drugs[[condition]] = list()

	for (i in 1:ncol(cytof_cast_mat)) {
		drug_ids = which(cytof_cast$Conditions == condition)
		t_tests_drugs[[condition]][[i]] = t.test(cytof_cast_mat[drug_ids, i], cytof_cast_mat[plasma_ids, i])
	}

}

tDotPlot = function(t_tests, cell_cols, t_sig, ...) {

	t_values = sapply(t_tests, function(test) test$statistic)
	t_values[is.na(t_values)] = 0.0

	plot(t_values,
		type="n",
		ylab="t-statistic",
		bty="n",  # no border
		xaxt="n",
		...
	)

	points(
		t_values,
		pch=16,
		col=cell_cols,
		xpd=TRUE
	)
	polygon(
		c(0, length(cells) + 1, length(cells) + 1, 0),
		c(-t_sig, -t_sig, t_sig, t_sig),
		col=rgb(0.9, 0.9, 0.9, 0.7), border=NA)

	abline(h=t_sig, col="grey", lty=3)
	abline(h=-t_sig, col="grey", lty=3)

}

# Get cell types from column strings
cells = sapply(strsplit(colnames(cytof_cast_mat), "_"), function(x) x[1])
phos = sapply(strsplit(colnames(cytof_cast_mat), "_"), function(x) x[2])

# phos_palette = c("black", brewer.pal(9, "Set1"))
phos_palette = brewer.pal(10, "Set3")
cell_palette = c("black", brewer.pal(8, "Dark2"))
# cell_palette = c("black", brewer.pal(9, "Paired"))
# cell_palette = c("black", brewer.pal(9, "Set1"))

cell_cols = cell_palette[as.numeric(factor(cells))]
phos_cols = phos_palette[as.numeric(factor(phos))]

t_sig = 3.143  # .99 significance level with df=6

svg("cytofPerturbFig/pert_tstat.svg", width=4)
par(mfrow=c(9, 1), mar=c(1, 4.1, 1, 2.1))
tDotPlot(t_tests, cell_cols, t_sig,
	main=paste0("Vehicle - No plasma (",
		paste(unique(phos), collapse=", "),  # order of phosphosites
		") (",
		paste(unique(cells), collapse=", "),  # order of cells
		")"
		),
	cex.main=0.3
)

drug_names = c(
	"Ro 31-8220 mesylate",
	"Alvocidib",
	"AZD8055",
	"Saracatinib",
	"PF-562271 HCl",
	"Mevastatin",
	"Dasatinib",
	"CGP 60474"
)

for (i in 1:8) {
	tDotPlot(t_tests_drugs[[i]], cell_cols, t_sig,
		main=drug_names[i], cex.main=0.7)
}

dev.off()



# Heatmap of standardized median CyTOF data
# ---------------------------------------------------------------
mat = scale(cytof_cast_mat)

rownames(mat) = stri_trans_totitle(rownames(mat))  # first letter upper case, rest lower case
colnames(mat) = sub("_", ":", colnames(mat))

col_scale = colorRampPalette(rev(brewer.pal(9, "RdBu")), interpolate="spline")(100)
# col_scale = colorRampPalette(rev(brewer.pal(11, "Spectral")), interpolate="spline")(100)

rlab = rbind(
	"Cell type"=cell_cols,
	"Phosphosite"=phos_cols
)

# conditions_palette = brewer.pal(9, "Set1")
conditions_palette = brewer.pal(8, "Accent")
conditions_col = rep(conditions_palette[1], nrow(cytof_cast))
conditions_col[cytof_cast$Conditions == "NO PLASMA"] = conditions_palette[2]
conditions_col[cytof_cast$Conditions == "VEHICLE"] = conditions_palette[3]

clab = cbind(
	"Experiment"=conditions_col
)

pdf("cytofPerturbFig/cytof_drug_pert_heatmap.pdf", width=6)
heatmap.3(
	t(mat),
	distfun=function(x) {
		x[is.na(x)] = 0.0
		dmat = dist(x)
		return(dmat)
	},
	# cytof_zscore,
	col=col_scale,
	RowSideColors=rlab,
	ColSideColors=clab,
	cexRow=0.4, cexCol=0.6,
	KeyValueName="z-score",
	xlab="Perturbation", ylab="Cell type:phosphosite"
)

legend(
	"topright",
	legend=c(
		unique(cells),
		"",
		levels(factor(phos))
	),
	fill=c(
		cell_palette,
		"white",
		phos_palette
	),
	# legend=c("Basal","LumA","LumB","Her2","Claudin","Normal","","Positive","Negative","NA","","Targeted","Chemo","","Approved","Experimental"),
	# fill=c("red","blue","cyan","pink","yellow","green","white","black","white","grey","white","darkorchid","darkred","white","green","darkgreen"),
	border=FALSE,
	bty="n",
	y.intersp = 0.7,
	cex=0.7
)

dev.off()


# Restricted heatmap and t-statistics dotplots...
# --------------------------------
cell_sel = c("CD14 monocytes", "CD16 monocytes", "CD1c DCs")
cells_include = cells %in% cell_sel
phos_sel = c("pp38", "pMAPKAP2", "pCREB", "pERK1/2", "pS6")
phos_include = phos %in% phos_sel
cytof_include = cells_include & phos_include



# Subselection of perturbation experiments
pert_sel =c("VEHICLE", "1", "3", "5", "6")
pert_include = rownames(cytof_cast_mat) %in% pert_sel

mat = scale(cytof_cast_mat[pert_include,])  # Standarize for pert selection only

# Format names for output
rownames(mat) = stri_trans_totitle(rownames(mat))  # first letter upper case, rest lower case
colnames(mat) = sub("_", ":", colnames(mat))

cytof_include_all = phos_include & cells %in% "All"
mat_subset_all = mat[,cytof_include_all]

mat_subset = mat[,cytof_include]


pdf("cytofPerturbFig/cytof_drug_pert_subset_heatmap2.pdf", width=6, height=4.5)
heatmap.3(
	t(mat_subset),
	distfun=function(x) {
		x[is.na(x)] = 0.0
		dmat = dist(x)
		return(dmat)
	},
	# cytof_zscore,
	col=col_scale,
	RowSideColors=rlab[,cytof_include],
	ColSideColors=clab[pert_include,, drop=FALSE],
	cexRow=0.6, cexCol=0.6,
	KeyValueName="z-score",
	mar=c(6, 8),
	xlab="Perturbation", ylab="Cell type:phosphosite"
)

legend(
	"topright",
	legend=c(
		unique(cells),
		"",
		levels(factor(phos))
	),
	fill=c(
		cell_palette,
		"white",
		phos_palette
	),
	# legend=c("Basal","LumA","LumB","Her2","Claudin","Normal","","Positive","Negative","NA","","Targeted","Chemo","","Approved","Experimental"),
	# fill=c("red","blue","cyan","pink","yellow","green","white","black","white","grey","white","darkorchid","darkred","white","green","darkgreen"),
	border=FALSE,
	bty="n",
	y.intersp = 0.7,
	cex=0.7
)
dev.off()


# All cell types (aggregate CyTOF statistics)
pdf("cytofPerturbFig/cytof_drug_pert_subset_all_heatmap2.pdf", width=6, height=4.5)
heatmap.3(
	t(mat_subset_all),
	distfun=function(x) {
		x[is.na(x)] = 0.0
		dmat = dist(x)
		return(dmat)
	},
	# cytof_zscore,
	col=col_scale,
	RowSideColors=rlab[,cytof_include_all],
	ColSideColors=clab[pert_include, , drop=FALSE],
	cexRow=0.6, cexCol=0.6,
	KeyValueName="z-score",
	mar=c(6, 8),
	xlab="Perturbation", ylab="Cell type:phosphosite"
)

legend(
	"topright",
	legend=c(
		unique(cells),
		"",
		levels(factor(phos))
	),
	fill=c(
		cell_palette,
		"white",
		phos_palette
	),
	# legend=c("Basal","LumA","LumB","Her2","Claudin","Normal","","Positive","Negative","NA","","Targeted","Chemo","","Approved","Experimental"),
	# fill=c("red","blue","cyan","pink","yellow","green","white","black","white","grey","white","darkorchid","darkred","white","green","darkgreen"),
	border=FALSE,
	bty="n",
	y.intersp = 0.7,
	cex=0.7
)
dev.off()



svg("cytofPerturbFig/pert_tstat_subset.svg", width=3)
t_sig = 1.943  # df=6, p=0.05
par(mfrow=c(9, 1), mar=c(1, 4.1, 1, 2.1))
tDotPlot(t_tests[cytof_include], cell_cols[cytof_include], t_sig,
	main=paste0("Vehicle - No plasma (",
		paste(unique(phos[cytof_include]), collapse=", "),  # order of phosphosites
		") (",
		paste(unique(cells[cytof_include]), collapse=", "),  # order of cells
		")"
		),
	cex.main=0.3
)

drug_names = c(
	"Ro 31-8220 mesylate",
	"Alvocidib",
	"AZD8055",
	"Saracatinib",
	"PF-562271 HCl",
	"Mevastatin",
	"Dasatinib",
	"CGP 60474"
)

for (i in 1:8) {
	tDotPlot(t_tests_drugs[[i]][cytof_include], cell_cols[cytof_include], t_sig,
		main=drug_names[i], cex.main=0.7)
}
dev.off()

svg("cytofPerturbFig/pert_tstat_subset_all.svg", width=3)
t_sig = 1.943  # df=6, p=0.05
par(mfrow=c(9, 1), mar=c(1, 4.1, 1, 2.1))
tDotPlot(t_tests[cytof_include_all], cell_cols[cytof_include_all], t_sig,
	main=paste0("Vehicle - No plasma (",
		paste(unique(phos[cytof_include_all]), collapse=", "),  # order of phosphosites
		") (",
		paste(unique(cells[cytof_include_all]), collapse=", "),  # order of cells
		")"
		),
	cex.main=0.3
)

drug_names = c(
	"Ro 31-8220 mesylate",
	"Alvocidib",
	"AZD8055",
	"Saracatinib",
	"PF-562271 HCl",
	"Mevastatin",
	"Dasatinib",
	"CGP 60474"
)

for (i in 1:8) {
	tDotPlot(t_tests_drugs[[i]][cytof_include_all], cell_cols[cytof_include_all], t_sig,
		main=drug_names[i], cex.main=0.7)
}
dev.off()




# Format data
# Get matrix of numerical data
cytof_mat =cytof[,5:ncol(cytof)]


# Regression model
# -------------------------------------------------------------------------
m = 1  # phosphosite

cytof_mat[,m]

drug_ = cytof$Conditions
drug_[cytof$Conditions == "NO PLASMA"] = "None"
drug_[cytof$Conditions == "VEHICLE"] = "None"

# drug_[cytof$Conditions == "NO PLASMA"] = NA
# drug_[cytof$Conditions == "VEHICLE"] = NA
drug_ = factor(drug_)

# Separate encoding for plasma experiment, particular to each patient
plasma_ = cytof$Individuals
plasma_[cytof$Conditions == "NO PLASMA"] = "None"
plasma_ = factor(plasma_)

# lin_fit = lm(cytof_mat[,m] ~ 0 + drug_ : cell_)
lin_fit = lm(cytof_mat[,m] ~ 0 + drug_ : cell_ + plasma_ : cell_)

# lin_fit = penalized(cytof_mat[,m] ~ 0 + drug_ : cell_ + plasma_ : cell_)
lin_fit = glmnet(cytof_mat[,m] ~ 0 + drug_ : cell_ + plasma_ : cell_)

confint(lin_fit)


# cast(cytof_mat, Conditions~Individuals)


# Test if lm uses NA factors correctly by incoorporating in 
lm(
	c(1, 2, 10, 9) ~
	factor(c("a1", "a2", NA, NA))
)






