library(dplyr)  
library(data.table)
library(ggplot2)
library(reshape2)
library(RColorBrewer)



############################# NEW DATA ##################################
athero <-read.csv("~/AtherovsHealthy/NEW_Athero_vs_healthy_noplasma.csv")
names <- make.unique(athero$Sample)
rownames(athero) <- names
athero <- athero[,-1]

# Create vector containing cytokine names
name <- (names(athero[,2:49]))

norm_data <- athero
norm_data[name] <- log(norm_data[name]+1)


norm_data$group<- rep(c("NO PLASMA","HEALTHY","ATHERO"), c(6,8,8))

#compare two groups
type1 <- "HEALTHY"
type2 <- "NO PLASMA"

df.test <- rbind(norm_data[which(norm_data$group == type1),], 
                 norm_data[which(norm_data$group == type2),])


matab_num <- ncol(df.test)

test_results <- data.frame(cytokines = colnames(df.test)[2:matab_num])

n_start = 2  # data starting column
t_tests= apply(norm_data[, n_start:ncol(norm_data)], 2, function(dcol) {
  x = dcol[norm_data$group == "HEALTHY"]
  y = dcol[norm_data$group == "NO PLASMA"]
  
  t_test = try(t.test(x, y),silent=TRUE)
  if (is(t_test, "try-error")) return(NA) 
  
  return(t_test)
})
names(t_tests) = colnames(norm_data)[n_start:ncol(norm_data)]

res = as.data.frame(do.call(rbind,t_tests))
res$adjP = p.adjust(res$p.value,"BH")
res$cytokines<-rownames(res)

# Log2 (FC)
cn <- paste("FC_", type1, "/", type2, sep = "")
test_results[cn] <- apply(df.test[,2:matab_num], 2, 
                          function(x) 
                            mean(as.numeric(x[which(df.test$group == type1)]))/
                            mean(as.numeric(x[which(df.test$group == type2)])))
test_results$LOG2FC <- log2(test_results$`FC_HEALTHY/NO PLASMA`)
data.frame(test_results)


#compare two groups
type1 <- "ATHERO"
type2 <- "NO PLASMA"

df.test <- rbind(norm_data[which(norm_data$group == type1),], 
                 norm_data[which(norm_data$group == type2),])


matab_num <- ncol(df.test)

test_results <- data.frame(cytokines = colnames(df.test)[2:matab_num])

n_start = 2
t_tests= apply(norm_data[, n_start:ncol(norm_data)], 2, function(dcol) {
  x = dcol[norm_data$group == "ATHERO"]
  y = dcol[norm_data$group == "NO PLASMA"]
  
  t_test = try(t.test(x, y),silent=TRUE)
  if (is(t_test, "try-error")) return(NA)

  return(t_test)
})
names(t_tests) = colnames(norm_data)[n_start:ncol(norm_data)]

res = as.data.frame(do.call(rbind,t_tests))
res$adjP = p.adjust(res$p.value,"BH")
res$cytokines<-rownames(res)

# Log2 (FC)
cn <- paste("FC_", type1, "/", type2, sep = "")
test_results[cn] <- apply(df.test[,2:matab_num], 2, 
                          function(x) 
                            mean(as.numeric(x[which(df.test$group == type1)]))/
                            mean(as.numeric(x[which(df.test$group == type2)])))
test_results$LOG2FC <- log2(test_results$`FC_ATHERO/NO PLASMA`)


#compare two groups
type1 <- "ATHERO"
type2 <- "HEALTHY"

df.test <- rbind(norm_data[which(norm_data$group == type1),], 
                 norm_data[which(norm_data$group == type2),])


matab_num <- ncol(df.test)

test_results <- data.frame(cytokines = colnames(df.test)[2:matab_num])

n_start = 2  
t_tests= apply(norm_data[, n_start:ncol(norm_data)], 2, function(dcol) {
  x = dcol[norm_data$group == "ATHERO"]
  y = dcol[norm_data$group == "HEALTHY"]
  
  t_test = try(t.test(x, y),silent=TRUE)
  if (is(t_test, "try-error")) return(NA) 

  return(t_test)
})
names(t_tests) = colnames(norm_data)[n_start:ncol(norm_data)]

res = as.data.frame(do.call(rbind,t_tests))
res$adjP = p.adjust(res$p.value,"BH")
res$cytokines<-rownames(res)

# Log2 (FC)
cn <- paste("FC_", type1, "/", type2, sep = "")
test_results[cn] <- apply(df.test[,2:matab_num], 2, 
                          function(x) 
                            mean(as.numeric(x[which(df.test$group == type1)]))/
                            mean(as.numeric(x[which(df.test$group == type2)])))
test_results$LOG2FC <- log2(test_results$`FC_ATHERO/HEALTHY`)


dot<-merge(res,test_results,by.x="cyto",by.y="cytokines")

dot$type_pvalue <- ifelse(dot$adjP < 0.05, "Significant","Not Significant")
dot$type_pvalue <- factor(dot$type_pvalue, levels = c("Significant","Not Significant"))

#raw pvalue
library(ggdendro)
p<-ggplot(dot, aes(group, cytokines)) +
  geom_point(aes(fill = LOG2FC, size = -log10(adjP), color = type_pvalue, shape = type_pvalue)) + 
  scale_fill_gradientn(colours = colorRampPalette(c("#3d67a3", "white","#ce1020"))(99), 
                       limits = c(-7.6,7.6), breaks = c(-7.6, -5, 0, 5, 7.6)) +
  scale_shape_manual(values = c(8,21,23)) +
  scale_color_manual(values = c("black","grey50","black"))+
  scale_size_continuous(breaks=c(1,2,3,5), range=c(1,5)) +
  theme_dendro() +
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  coord_fixed(ratio = 0.7) +
  labs(fill = "Log2 (Fold change)", size = "-Log10 (q value)", color = "", shape = "") +
  theme(axis.text.x = element_text(size = 10, face = "plain", colour = "black", angle = 45,vjust = 1.1, hjust=1.1), 
        axis.text.y = element_text(size = 10, face = "plain", colour = "black")) +
  theme(legend.text = element_text(size = 10, face = "plain", colour = "black"), 
        legend.title = element_text(size = 10, face = "plain", colour = "black"), 
        legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.3, "cm"))  
pdf(file="NEWDATA_pointplot_netrelease_athero_healthy_BH.pdf",height = 9,width = 8)
p
dev.off()


# Analyze data
# Cytokine release data
zscore<-norm_data[,-c(1)]
release_zscore = scale(zscore)
release_zscore[is.na(release_zscore)] = 0.0


col_scale = colorRampPalette(rev(brewer.pal(12, "Set3")), interpolate="spline")(100)

anno_df = data.frame(group=norm_data$group,
                     stringsAsFactors = FALSE)
anno_df$group <- factor(norm_data$group, levels = c("NO PLASMA", "HEALTHY", "ATHERO" ))


annotation_colors = list(
  group=c("NO PLASMA"="#272E6A","HEALTHY"="#80b1d3","ATHERO"="#D51F26"))

head(annotation_colors)
rownames(anno_df) <- rownames(release_zscore)


myBreaks <- c(seq(min(release_zscore), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(release_zscore)/paletteLength, max(release_zscore), length.out=floor(paletteLength/2)))
p<-pheatmap(t(release_zscore),annotation_col=anno_df, annotation_colors=annotation_colors,
            cluster_cols = F,cluster_rows = T,cellwidth = 10,cellheight=10,fontsize_col = 8,
            color = colorRampPalette(c("blue", "white", "red"))(paletteLength),breaks=myBreaks)
pdf("AtherovsHealthy/NEWDATA_zscored_heatmap_athero_noplasma.pdf",  width = 10, height = 10)
print(p)
dev.off()


## Read Phospho-cytof data
phospho<-read.csv("~/combat_corrected_without_drug.csv")

cd14<-phospho[which(phospho$celltype == "CD14_Mono"),]
cd16<-phospho[which(phospho$celltype == "CD16_Mono"),]
Dc<-phospho[which(phospho$celltype == "CD1c_DCs"),]
cd4<-phospho[which(phospho$celltype == "CD4_Tcells"),]
cd8<-phospho[which(phospho$celltype == "CD8_Tcells"),]

####### compare sys vs asym
cd14_nothealthy<-cd14[which(cd14$group != "Healthy"),]
dc_nothealthy<-Dc[which(Dc$group != "Healthy"),]
cd16_nothealthy<-cd16[which(cd16$group != "Healthy"),]
cd4_nothealthy<-cd4[which(cd4$group != "Healthy"),]
cd8_nothealthy<-cd8[which(cd8$group != "Healthy"),]



# Difference between sys vs asym
cd14_sys_aym_ttest = lapply(name, function(cyto) {
  
  idx1 = cd14_nothealthy$group == "Symptomatic"
  idx2 = cd14_nothealthy$group == "Asymptomatic"
  

  x = cd14_nothealthy[idx1, cyto]
  y = cd14_nothealthy[idx2, cyto]
  
  x[x < 0] = 0
  y[y < 0] = 0
  
  try({
    result = t.test(x, y, paired=FALSE)
    result$estimate = c(
      mean(x, na.rm=TRUE),
      mean(y, na.rm=TRUE)
    )
    return(result)
  })
  return("Failed")
})
names(cd14_sys_aym_ttest) = name


cd16_sys_aym_ttest = lapply(name, function(cyto) {
  
  idx1 = cd16_nothealthy$group == "Symptomatic"
  idx2 = cd16_nothealthy$group == "Asymptomatic"

  x = cd16_nothealthy[idx1, cyto]
  y = cd16_nothealthy[idx2, cyto]
  
  
  x[x < 0] = 0
  y[y < 0] = 0
  
  
  try({
    result = t.test(x, y, paired=FALSE)
    result$estimate = c(
      mean(x, na.rm=TRUE),
      mean(y, na.rm=TRUE)
    )
    return(result)
  })
  return("Failed")
})
names(cd16_sys_aym_ttest) = name

dc_sys_aym_ttest = lapply(name, function(cyto) {
  
  idx1 = dc_nothealthy$group == "Symptomatic"
  idx2 = dc_nothealthy$group == "Asymptomatic"

  x = dc_nothealthy[idx1, cyto]
  y = dc_nothealthy[idx2, cyto]
  
  
  x[x < 0] = 0
  y[y < 0] = 0
  
  
  try({
    result = t.test(x, y, paired=FALSE)
    result$estimate = c(
      mean(x, na.rm=TRUE),
      mean(y, na.rm=TRUE)
    )
    return(result)
  })
  return("Failed")
})
names(dc_sys_aym_ttest) = name


cd4_sys_aym_ttest = lapply(name, function(cyto) {
  
  idx1 = cd4_nothealthy$group == "Symptomatic"
  idx2 = cd4_nothealthy$group == "Asymptomatic"
  
  x = cd4_nothealthy[idx1, cyto]
  y = cd4_nothealthy[idx2, cyto]
  
  x[x < 0] = 0
  y[y < 0] = 0
  
  try({
    result = t.test(x, y, paired=FALSE)
    result$estimate = c(
      mean(x, na.rm=TRUE),
      mean(y, na.rm=TRUE)
    )
    return(result)
  })
  return("Failed")
})
names(cd4_sys_aym_ttest) = name

cd8_sys_aym_ttest = lapply(name, function(cyto) {
  
  idx1 = cd8_nothealthy$group == "Symptomatic"
  idx2 = cd8_nothealthy$group == "Asymptomatic"
  
  x = cd8_nothealthy[idx1, cyto]
  y = cd8_nothealthy[idx2, cyto]
  
  x[x < 0] = 0
  y[y < 0] = 0
  
  try({
    result = t.test(x, y, paired=FALSE)
    return(result)
  })
  return("Failed")
})
names(cd8_sys_aym_ttest) = name


pvals = data.frame(
  cd14_sys_vs_asym<-sapply(cd14_sys_aym_ttest,function(x) x$p.value),
  cd16_sys_vs_asym<-sapply(cd16_sys_aym_ttest, function(x) x$p.value),
  dc_sys_vs_asym<-sapply(dc_sys_aym_ttest,function(x) x$p.value),
  cd4_sys_vs_asym<-sapply(cd4_sys_aym_ttest, function(x) x$p.value),
  cd8_sys_vs_asym<-sapply(cd8_sys_aym_ttest,function(x) x$p.value)
)
pvals[is.na(pvals)] = 1.0  


fold_changes = data.frame(
  cd14_sys_vs_asym=sapply(cd14_sys_aym_ttest, function(x) log2(x$estimate[1] / x$estimate[2])),
  cd16_sys_vs_asym=sapply(cd16_sys_aym_ttest, function(x) log2(x$estimate[1] / x$estimate[2])),
  dc_sys_vs_asym=sapply(dc_sys_aym_ttest, function(x) log2(x$estimate[1] / x$estimate[2])),
  cd4_sys_vs_asym=sapply(cd4_sys_aym_ttest, function(x) log2(x$estimate[1] / x$estimate[2])),
  cd8_sys_vs_asym=sapply(cd8_sys_aym_ttest, function(x) log2(x$estimate[1] / x$estimate[2]))
)
rownames(fold_changes) =name
fold_changes[is.na(fold_changes)] = 0.0
fold_changes[fold_changes > 10] = 10
fold_changes[fold_changes < -10] = -10

mat = fold_changes
mat = data.matrix(mat)

notes = matrix("", nrow=nrow(mat), ncol=ncol(mat))
notes[pvals < 0.1] = ""
notes[pvals < 0.05] = "*"
notes[pvals < 0.01] = "**"
notes[pvals < 0.001] = "***"

## Foldchange heatmap
pdf(file="/heatmap_allcells_sys_asym.pdf",height = 6,width = 6)
heatmap.2(mat,
          Colv=FALSE,
          Rowv = FALSE,
          mar=c(16, 14),
          trace="none",
          # col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)
          breaks=seq(-1,1, length.out=101),
          col=colorRampPalette(rev(brewer.pal(9, "Spectral")))(100),
          key.xlab=expression("log"[2] * " FC"), key.ylab=NA,
          key.title="",
          cellnote=notes,
          notecol="black",
          cexCol=1
)
dev.off()

