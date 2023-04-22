library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(ggrepel)
library(pheatmap)


#load DESeq object
rawcounts<-read.csv("~/DESeq/raw_counts.csv")
sample_table<-read.csv("~/DESeq/metadata.csv")

names <- make.unique(rawcounts$gene_name)
rownames(rawcounts) <- names
rawcounts <- rawcounts[,-1] 


#Run DESEQ analysis of RNA-seq data

rawcounts[is.na(rawcounts)] <- 0
cds <- DESeqDataSetFromMatrix(countData = rawcounts,
                              colData = sample_table,
                              design = ~ group)

cds = DESeq(cds)

#normalization and differential expression steps are done here
resultsNames(cds)

#relevel the comparisons
cds$group <- relevel(cds$group, ref='Vehicle')
#cds <- DESeq(cds, parallel=T, betaPrior = F)
cds <- nbinomWaldTest(cds)

save(cds,file="~/cds.rda")


#Normalized by size factor
rna_norm = counts(cds, normalized=TRUE) 
gene_symbols<-rownames(rna_norm)

sara<-rna_norm

# standardized data
x_zscore = t(scale(t(sara)))
x_zscore[is.na(x_zscore)] <- 0

#Differential expression analysis for each group comparison 
auto.plasma_sara.25uM_vs_auto.plasma <- results(cds, name="group_Saracatinib._vs_Vehicle",alpha=0.05)
auto.plasma_sara.25uM_vs_auto.plasma <- as.data.frame(lfcShrink(cds_rpmi, coef = 2, res = auto.plasma_sara.25uM_vs_auto.plasma,type='apeglm'))

#Subset using Padj value
fc_thresh <- 0
log2cutoff<-1.2
qvaluecutoff <- 0.001
mean_counts_thresh = 4.0

sigGenes <- unique(c(
  rownames(subset(auto.plasma_sara.25uM_vs_auto.plasma, padj<=qvaluecutoff & abs(log2FoldChange)>=log2cutoff & baseMean > mean_counts_thresh)))

#Heatmap
#create matrix for averaged expression-values 
myBreaks <- c(seq(min(x_zscore), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(x_zscore)/paletteLength, max(x_zscore), length.out=floor(paletteLength/2)))
pheatmap(x_zscore[sigGenes,], fontsize_row = 6, 
            cellwidth = 8,cellheight=6,cluster_cols = T,
            annotation_colors = mycolors,
            main ="Unique Sig Genes across comparisons",
            color = colorRampPalette(c("#440154" ,"#21908C", "#FDE725")),breaks=myBreaks,
            fontsize_col = 7)


### Cytokine 

#Load luminex data for study comparison Saracatinib vs Vehicle
# file with clinical information, Luminex data

data<-fread("~/Luminex/luminex_Tissue.csv",header=T)#,sep=",")
data<-as.data.frame(data)

dim(data)

#Create vector containing cytokine names
name <- (names(data[,3:50]))

#Normalzie data
norm_data <- data
norm_data[name] <- log(norm_data[name]+1)


#Run Unpaired T Test
norm_data$group<- rep(c("Vehicle","Saracatinib","Media"), c(6,6,4))


#compare two groups
type1 <- "Saracatinib"
type2 <- "Vehicle"

df.test <- rbind(norm_data[which(norm_data$group == type1),], 
                 norm_data[which(norm_data$group == type2),])


matab_num <- ncol(df.test)

test_results <- data.frame(cytokines = colnames(df.test)[3:matab_num])

n_start = 3  
t_tests= apply(norm_data[, n_start:ncol(norm_data)], 2, function(dcol) {
  x = dcol[norm_data$group == "Saracatinib"]
  y = dcol[norm_data$group == "Vehicle"]
  
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
test_results[cn] <- apply(df.test[,3:matab_num], 2, 
                          function(x) 
                            mean(as.numeric(x[which(df.test$group == type1)]))/
                            mean(as.numeric(x[which(df.test$group == type2)])))
test_results$LOG2FC <- log2(test_results$`FC_Saracatinib/Vehicle`)

dot<-merge(res,test_results,by.x="cyto",by.y="cytokines")


dot$type_pvalue <- ifelse(dot$p.value < 0.05, "Significant","Not Significant")
dot$type_pvalue <- factor(dot$type_pvalue, levels = c("Significant","Not Significant"))

#raw pvalue
svg(file="/Luminex/point_plot.svg",height = 10,width = 8)
ggplot(dot, aes(group, cytokines)) +
  geom_point(aes(fill = LOG2FC, size = -log10(p.value), color = type_pvalue, shape = type_pvalue)) + 
  scale_fill_gradientn(colours = colorRampPalette(c("#3d67a3", "white","#ce1020"))(99),
                       limits = c(-1,1), breaks = c(-1, -0.5, 0, 0.5, 1)) +
  scale_shape_manual(values = c(21,23)) +
  scale_color_manual(values = c("black","grey50","black"))+
  scale_size_continuous(breaks=c(0.5,1.5,3,4), range=c(1,4)) +
  theme_dendro() +
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  coord_fixed(ratio = 0.7) +
  labs(fill = "Log2 (Fold change)", size = "-Log10 (p value)", color = "", shape = "") +
  theme(axis.text.x = element_text(size = 10, face = "plain", colour = "black", angle = 45),
        axis.text.y = element_text(size = 10, face = "plain", colour = "black")) +
  theme(legend.text = element_text(size = 10, face = "plain", colour = "black"), 
        legend.title = element_text(size = 10, face = "plain", colour = "black"), 
        legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.3, "cm"))   
dev.off()

