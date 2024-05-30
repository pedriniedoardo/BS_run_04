# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(GGally)
library(cowplot)
library(ComplexHeatmap)
library(scales)
library(circlize)
library(DESeq2)
library(RNAseqQC)
library(limma)
library(ashr)
library(magick)
library(UpSetR)

# read in the final object ------------------------------------------------
#
scobj <- readRDS("../../out/object/sobj_processed_donor.rds")
DimPlot(scobj,label = T,raster = T,group.by = "seurat_clusters")
DimPlot(scobj,label = T,raster = T,group.by = "expertAnno.l1")
DimPlot(scobj,label = T,raster = T,group.by = "treat_full")

# wrangling ---------------------------------------------------------------
# explore the full meta
scobj@meta.data
# aggregate the expression per sample per treatment, donor and cell type
cts_sample_all <- AggregateExpression(object = scobj,
                                      group.by = c("id_sample_short","harmonized_donor2","treat_full", "expertAnno.l1"),
                                      assays = 'RNA',
                                      slot = "counts",
                                      return.seurat = FALSE)

# MG processing -----------------------------------------------------------
# focus on the immune cells
# 1. Get counts matrix
counts_MG <- cts_sample_all$RNA %>% 
  data.frame() %>% 
  rownames_to_column("gene") %>% 
  # pull only the MG cells
  dplyr::select("gene",contains("MG")) %>% 
  # remove the unassigned donors
  dplyr::select(-contains("unassigned")) %>% 
  column_to_rownames("gene")

# 2. generate sample level metadata
colData_MG <- data.frame(samples = colnames(counts_MG)) %>%
  separate(samples,into = c("id_sample_short","harmonized_donor2","treat_full", "expertAnno.l1"),remove = F,sep = "_") %>%
  mutate(treat_full = factor(treat_full,levels = c("CSF.ctrl.24h","BASELINE","CSF.MS.24h","CSF.MS.48h","cytokine","Fe","myelin","TBHP")))

# perform DESeq2 ----------------------------------------------------------
# build the model
treat_full <- colData_MG$treat_full
design_MG <- model.matrix(~ treat_full)
colnames(design_MG)[1] <- c("intercept")
# save the disign
saveRDS(design_MG,"../../out/object/design_MG_refCSFctrl.rds")

# Create DESeq2 object   
dds_MG <- DESeqDataSetFromMatrix(countData = counts_MG,
                                 colData = colData_MG,
                                 design = design_MG)

# filter
# filter at least 10 read per gene per sample
keep_MG <- rowSums(counts(dds_MG)) >=310
# keep_MG <- rowSums(counts(dds_MG)) >=10
dds_MG_filter <- dds_MG[keep_MG,]

# scale the data
vds_MG_filter <- vst(dds_MG_filter, blind = F)

# clustering samples ------------------------------------------------------
# set seed to control random annotation colors
# pdf("../../out/image/heatmap_SampleCluster_pseudobulk_MG_refCSFctrl.pdf",width = 12,height = 6)
# set.seed(1)
# hm_MG <- plot_sample_clustering(vds_MG_filter,
#                                 anno_vars = c("treat_full","harmonized_donor2"),
#                                 distance = "euclidean")
# draw(hm_MG,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(2, 2, 2, 30), "mm"))
# dev.off()

# PCA ---------------------------------------------------------------------
# plot_vsd_MG <- plotPCA(vds_MG_filter,
#                        intgroup = c("treat_full","harmonized_donor2")) +
#   theme_bw()
# 
# plot_vsd_MG$data %>%
#   # ggplot(aes(x=PC1,y=PC2,col=condition,label=clone,shape=gender_update)) +
#   ggplot(aes(x=PC1,y=PC2,col=treat_full,shape=harmonized_donor2)) +
#   geom_point(size =3,alpha=0.6) +
#   # ggrepel::geom_text_repel(show.legend = F)+
#   scale_x_continuous(expand = expansion(mult = 0.1))+
#   theme_bw() + ylab(plot_vsd_MG$labels[1]) + xlab(plot_vsd_MG$labels[2])
# ggsave("../../out/image/PCA_pseudobulk_MG_refCSFctrl.pdf",width = 6,height = 4)

# run DE ------------------------------------------------------------------
# run DESeq2
ddsHTSeq_MG_filter <- DESeq(dds_MG_filter)

# Check the coefficients for the comparison
resultsNames(ddsHTSeq_MG_filter)

# save the filtered object
saveRDS(ddsHTSeq_MG_filter,"../../out/object/ddsHTSeq_pseudobulk_MG_refCSFctrl.rds")

# print the contrast
resultsNames(ddsHTSeq_MG_filter)

# save the resut of the contrast of interest. notice that is also possible to set the alpha value for the significance (adj-p value)
contrast_MG <- makeContrasts(BASELINE_vs_CSF.ctrl.24h = treat_fullBASELINE,
                             CSF.MS.24h_vs_CSF.ctrl.24h = treat_fullCSF.MS.24h,
                             CSF.MS.48h_vs_CSF.ctrl.24h = treat_fullCSF.MS.48h,
                             cytokine_vs_CSF.ctrl.24h = treat_fullcytokine,
                             Fe_vs_CSF.ctrl.24h = treat_fullFe,
                             myelin_vs_CSF.ctrl.24h = treat_fullmyelin,
                             TBHP_vs_CSF.ctrl.24h = treat_fullTBHP,
                             levels = design_MG)

# pull the results table
res_BASELINE <- results(ddsHTSeq_MG_filter, contrast=contrast_MG[,"BASELINE_vs_CSF.ctrl.24h"],alpha = 0.05)
summary(res_BASELINE)
res_CSF.MS.24h <- results(ddsHTSeq_MG_filter, contrast=contrast_MG[,"CSF.MS.24h_vs_CSF.ctrl.24h"],alpha = 0.05)
summary(res_CSF.MS.24h)
res_CSF.MS.48h <- results(ddsHTSeq_MG_filter, contrast=contrast_MG[,"CSF.MS.48h_vs_CSF.ctrl.24h"],alpha = 0.05)
summary(res_CSF.MS.48h)
res_cytokine <- results(ddsHTSeq_MG_filter, contrast=contrast_MG[,"cytokine_vs_CSF.ctrl.24h"],alpha = 0.05)
summary(res_cytokine)
res_Fe <- results(ddsHTSeq_MG_filter, contrast=contrast_MG[,"Fe_vs_CSF.ctrl.24h"],alpha = 0.05)
summary(res_Fe)
res_myelin <- results(ddsHTSeq_MG_filter, contrast=contrast_MG[,"myelin_vs_CSF.ctrl.24h"],alpha = 0.05)
summary(res_myelin)
res_TBHP <- results(ddsHTSeq_MG_filter, contrast=contrast_MG[,"TBHP_vs_CSF.ctrl.24h"],alpha = 0.05)
summary(res_TBHP)

# add the gene symbols
df_res <- 
  list(BASELINE = res_BASELINE,
       CSF.MS.24h = res_CSF.MS.24h,
       CSF.MS.48h = res_CSF.MS.48h,
       cytokine = res_cytokine,
       Fe = res_Fe,
       myelin = res_myelin,
       TBHP = res_TBHP) %>%
  lapply(function(x){
    df <- x %>%
      data.frame() %>%
      rownames_to_column(var = "symbol") %>%
      arrange(pvalue)
    # left_join(LUT_gene,by = c("symbol"="gene_name"))
    df
  }) %>%
  bind_rows(.id = "conditionVsCSFctrl")
# save the table
df_res %>%
  write_tsv("../../out/table/DE_treatvsCSFctrl_pseudobulk_MG.tsv")

# shrink ------------------------------------------------------------------  
res_BASELINE_shr <- lfcShrink(ddsHTSeq_MG_filter, res = res_BASELINE, type = "ashr")
res_CSF.MS.24h_shr <- lfcShrink(ddsHTSeq_MG_filter, res = res_CSF.MS.24h, type = "ashr")
res_CSF.MS.48h_shr <- lfcShrink(ddsHTSeq_MG_filter, res = res_CSF.MS.48h, type = "ashr")
res_cytokine_shr <- lfcShrink(ddsHTSeq_MG_filter, res = res_cytokine, type = "ashr")
res_Fe_shr <- lfcShrink(ddsHTSeq_MG_filter, res = res_Fe, type = "ashr")
res_myelin_shr <- lfcShrink(ddsHTSeq_MG_filter, res = res_myelin, type = "ashr")
res_TBHP_shr <- lfcShrink(ddsHTSeq_MG_filter, res = res_TBHP, type = "ashr")

# save the table of DEGs
df_res_shr <- 
  list(BASELINE_shr = res_BASELINE_shr,
       CSF.MS.24h_shr = res_CSF.MS.24h_shr,
       CSF.MS.48h_shr = res_CSF.MS.48h_shr,
       cytokine_shr = res_cytokine_shr,
       Fe_shr = res_Fe_shr,
       myelin_shr = res_myelin_shr,
       TBHP_shr = res_TBHP_shr) %>%
  lapply(function(x){
    df <- x %>%
      data.frame() %>%
      rownames_to_column(var = "symbol") %>%
      arrange(pvalue)
    # left_join(LUT_gene,by = c("symbol"="gene_name"))
    df
  }) %>%
  bind_rows(.id = "conditionVsCSFctrl")
# save the table
df_res_shr %>%
  write_tsv("../../out/table/DE_treatvsCSFctrl_pseudobulk_MG_shr.tsv")

# Another useful diagnostic plot is the histogram of the p values (figure below). This plot is best formed by excluding genes with very small counts, which otherwise generate spikes in the histogram.
df_res %>%
  data.frame()%>%
  dplyr::filter(baseMean>1) %>%
  ggplot(aes(x=pvalue))+geom_histogram(breaks = 0:20/20) +
  facet_wrap(~conditionVsCSFctrl)+
  theme_bw()+
  theme(strip.background = element_blank())
ggsave("../../out/image/histogram_pvalue_pseudobulk_MG_refCSFctrl.pdf",width = 6,height = 6)

# volcano -----------------------------------------------------------------
# add the info of the genename
plot_volcano <- df_res %>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>1&!is.na(symbol)),yes = 1,no = 0)) %>% 
  filter(!is.na(col))

plot_volcano %>%
  ggplot(aes(x=log2FoldChange,y=-log(padj)))+
  # geom_point()
  geom_point(data = plot_volcano[plot_volcano$col==0,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.05)+
  geom_point(data = plot_volcano[plot_volcano$col==1,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.5)+
  geom_vline(xintercept = c(-1,1),col="red",linetype="dashed")+
  geom_hline(yintercept = (-log(0.05)),col="red",linetype="dashed")+
  scale_color_manual(values = c("black","red"))+theme(legend.position = "none")+
  # ggrepel::geom_text_repel(
  #   data = plot_volcano[plot_volcano$col==1,][1:1000,],
  #   aes(label = symbol),max.overlaps = 1,segment.alpha=0.4,
  #   size = 2,
  #   box.padding = unit(0.35, "lines"),
  #   point.padding = unit(0.3, "lines")) +
  ggrepel::geom_text_repel(
    data = plot_volcano %>% group_by(conditionVsCSFctrl) %>% arrange(padj) %>% dplyr::slice(1:10),
    aes(label = symbol),segment.alpha=0.4) +
  facet_wrap(~conditionVsCSFctrl)+
  theme_bw()+
  theme(strip.background = element_blank())+
  theme(legend.position = "none")
ggsave("../../out/image/vulcano_plot_pseudobulk_MG_refCSFctrl.pdf",width = 20,height = 20)

#
plot_volcano_shr <- df_res_shr %>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>1&!is.na(symbol)),yes = 1,no = 0)) %>% 
  filter(!is.na(col))

plot_volcano_shr %>%
  ggplot(aes(x=log2FoldChange,y=-log(padj)))+
  # geom_point()
  geom_point(data = plot_volcano_shr[plot_volcano_shr$col==0,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.05)+
  geom_point(data = plot_volcano_shr[plot_volcano_shr$col==1,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.5)+
  geom_vline(xintercept = c(-1,1),col="red",linetype="dashed")+
  geom_hline(yintercept = (-log(0.05)),col="red",linetype="dashed")+
  scale_color_manual(values = c("black","red"))+theme(legend.position = "none")+
  ggrepel::geom_text_repel(
    data = plot_volcano_shr %>% group_by(conditionVsCSFctrl) %>% arrange(padj) %>% dplyr::slice(1:10),
    aes(label = symbol),segment.alpha=0.4) +
  facet_wrap(~conditionVsCSFctrl)+
  theme_bw()+
  theme(strip.background = element_blank())+
  theme(legend.position = "none")
ggsave("../../out/image/vulcano_plot_pseudobulk_MG_refCSFctrl_shr.pdf",width = 20,height = 20)

# MA plot -----------------------------------------------------------------
df_res %>%
  filter(!is.na(padj)) %>% 
  data.frame()%>%
  mutate(color = ifelse(padj<0.05,yes = 1,no = 0))%>%
  mutate(color = factor(ifelse(is.na(color),yes = 0,no = color)))%>%
  arrange(color) %>% 
  ggplot(aes(baseMean,y = log2FoldChange,col=color)) + geom_point(alpha=0.2) + 
  scale_x_log10() + scale_color_manual(values = c("gray","red")) +
  facet_wrap(~conditionVsCSFctrl)+
  theme_bw()+
  theme(strip.background = element_blank()) +
  theme(legend.position = "none")
ggsave("../../out/image/MA_plot_pseudobulk_MG_refCSFctrl.pdf",width = 20,height = 20)

# shr
df_res_shr %>%
  filter(!is.na(padj)) %>% 
  data.frame()%>%
  mutate(color = ifelse(padj<0.05,yes = 1,no = 0))%>%
  mutate(color = factor(ifelse(is.na(color),yes = 0,no = color)))%>%
  arrange(color) %>% 
  ggplot(aes(baseMean,y = log2FoldChange,col=color)) + geom_point(alpha=0.2) + 
  scale_x_log10() + scale_color_manual(values = c("gray","red")) +
  facet_wrap(~conditionVsCSFctrl)+
  theme_bw()+
  theme(strip.background = element_blank()) +
  theme(legend.position = "none")
ggsave("../../out/image/MA_plot_pseudobulk_MG_refCSFctrl_shr.pdf",width = 20,height = 20)

# heatmaps ----------------------------------------------------------------
# pull the scaled values
mat_filter_MG <- assay(vds_MG_filter) %>%
  data.frame() %>%
  # dplyr::select(contains(c("_0_","_6_"))) %>%
  as.matrix()

# the DEGs plot
DEG_2 <- df_res_shr %>%
  data.frame()%>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>1),yes = 1,no = 0)) %>%
  dplyr::filter(col==1) %>%
  pull(symbol) %>%
  unique()

# mat_filter <- assay(vds_filter) %>%
#   data.frame() %>%
#   # dplyr::select(contains(c("_0_","_6_"))) %>%
#   as.matrix()

mat_shr_MG <- mat_filter_MG[rownames(vds_MG_filter) %in% DEG_2, ]
mat2_shr_MG <- (mat_shr_MG - rowMeans(mat_shr_MG))/rowSds(mat_shr_MG,useNames = TRUE)
#
meta_sample_MG <- data.frame(colname = colnames(mat2_shr_MG)) %>% 
  left_join(colData_MG,by=c("colname"="samples"))

# make the column of the matrix more readable
colnames(mat2_shr_MG) <- meta_sample_MG$colname

column_ha_shr_MG <- HeatmapAnnotation(treat = meta_sample_MG$treat_full,
                                      col = list(treat = c("BASELINE" = "blue",
                                                           "CSF.ctrl.24h" = "gray90",
                                                           "CSF.MS.24h" = "gray40",
                                                           "CSF.MS.48h" = "gray10",
                                                           "cytokine" = "purple",
                                                           "Fe" = "orange",
                                                           "myelin" = "cyan",
                                                           "TBHP" = "yellow"))) 

ht2_shr_MG <- Heatmap(mat2_shr_MG, show_column_names = T,raster_by_magick = T,show_row_dend = F, use_raster = T,
                      name = "exp", 
                      column_title = "MG_shr",
                      # row_names_gp = gpar(fontsize = 3),
                      top_annotation = column_ha_shr_MG,show_row_names = F,
                      # cluster_rows = F, 
                      # right_annotation = row_ha, 
                      # row_split = rep(c(1,2,3,4),c(2,3,4,7))
                      
) 
pdf("../../out/image/heatmap_DEG_plot_pseudobulk_MG_refCSFctrl_shr.pdf",width = 7,height = 14) 
draw(ht2_shr_MG,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(30,2,2, 2), "mm")) 
dev.off()

# upset plot --------------------------------------------------------------
# read in the table of DEGs
df_res_shr <- read_tsv("../../out/table/DE_treatvsCSFctrl_pseudobulk_MG_shr.tsv")

# build a list of common elements belonging to each set of fegs  
list_DE_up <- df_res_shr %>%
  split(f = .$conditionVsCSFctrl) %>%
  lapply(function(x){
    x %>%
      filter(padj < 0.05,log2FoldChange>1) %>% 
      pull(symbol) %>%
      unique()  
  })

glimpse(list_DE_up)  

list_DE_down <- df_res_shr %>%
  split(f = .$conditionVsCSFctrl) %>%
  lapply(function(x){
    x %>%
      filter(padj < 0.05,log2FoldChange<(-1)) %>% 
      pull(symbol) %>%
      unique()  
  })
glimpse(list_DE_down)  

# try the upset plot version 
# library(UpSetR) 
pdf("../../out/image/upset_DEG_UP_plot_pseudobulk_MG_refCSFctrl_shr.pdf",width = 14,height = 7) 
upset(fromList(list_DE_up), order.by = "freq",nsets = 7) 
dev.off()

pdf("../../out/image/upset_DEG_DOWN_plot_pseudobulk_MG_refCSFctrl_shr.pdf",width = 14,height = 7) 
upset(fromList(list_DE_down), order.by = "freq",nsets = 7) 
dev.off()

# pull the intersections
df1_UP <- lapply(list_DE_up,function(x){
  data.frame(gene = x)
}) %>% 
  bind_rows(.id = "path")

df1_DOWN <- lapply(list_DE_down,function(x){
  data.frame(gene = x)
}) %>% 
  bind_rows(.id = "path")

head(df1_UP)
head(df1_DOWN)

df2_UP <- data.frame(gene=unique(unlist(list_DE_up)))
df2_DOWN <- data.frame(gene=unique(unlist(list_DE_down)))

head(df2_UP)
head(df2_DOWN)

df_int_UP <- lapply(df2_UP$gene,function(x){
  # pull the name of the intersections
  intersection <- df1_UP %>% 
    dplyr::filter(gene==x) %>% 
    arrange(path) %>% 
    pull("path") %>% 
    paste0(collapse = "|")
  
  # build the dataframe
  data.frame(gene = x,int = intersection)
}) %>% 
  bind_rows()

df_int_DOWN <- lapply(df2_DOWN$gene,function(x){
  # pull the name of the intersections
  intersection <- df1_DOWN %>% 
    dplyr::filter(gene==x) %>% 
    arrange(path) %>% 
    pull("path") %>% 
    paste0(collapse = "|")
  
  # build the dataframe
  data.frame(gene = x,int = intersection)
}) %>% 
  bind_rows()

df_int_UP %>%
  write_tsv("../../out/table/upset_DEG_UP_plot_pseudobulk_MG_refCSFctrl_shr.tsv")

df_int_DOWN %>%
  write_tsv("../../out/table/upset_DEG_DOWN_plot_pseudobulk_MG_refCSFctrl_shr.tsv")

head(df_int_UP,n=20)
head(df_int_DOWN,n=20)

df_int_UP %>% 
  group_by(int) %>% 
  summarise(n=n()) %>% 
  arrange(desc(n))

df_int_DOWN %>% 
  group_by(int) %>% 
  summarise(n=n()) %>% 
  arrange(desc(n))

# PLOT DISPERSION ---------------------------------------------------------
pdf("../../out/image/dispersion_plot_pseudobulk_MG_refCSFctrl.pdf",width = 5,height = 5) 
plotDispEsts(ddsHTSeq_MG_filter)
dev.off()