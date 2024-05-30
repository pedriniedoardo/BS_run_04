# AIM ---------------------------------------------------------------------
# explore the dataset to answer specific queries

# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# read in the data --------------------------------------------------------
sobj_run04 <- readRDS("../../out/object/sobj_processed_donor.rds")
sobj_run03 <- readRDS("/beegfs/scratch/ric.cosr/pedrini.edoardo/scRNAseq_Brainsphere_Absinta/BS_batch_2023/out/object/sobj_processed_donor.rds")

(DimPlot(sobj_run04,raster = T) + ggtitle("Run 04"))+
  (DimPlot(sobj_run03,raster = T) + ggtitle("Run 03"))

# explore SOX10 expression ------------------------------------------------
# One option suggested by Francesca is that the cells are drifting and loosing the cassette of SOX10.
# GOI <- c("SOX10","GAPDH","ACTB")
GOI <- list(
  IMMUNE = c("AIF1","TYROBP","FTL","HLA-DRA","TREM2","CX3CR1","C3","CSF1R", "CD74", "RUNX1","C1QB","PTPRC"),
  OLIGOLINEAGE = c("CSPG4","OLIG1","OLIG2", "PDGFRA", "SOX6", "PLP1","SOX10","BCAS1","MBP","MAG"),
  ASTRO = c("AQP4", "GFAP", "VIM","SLC1A2","S100B"),
  # Neu = c("SYT1")
  NEURONS = c("GAD2", "TLE4", "CUX2","SYP", "NEFL","SYT1"),
  NPC = c("NES", "PAX6", "SOX1"),
  CYCLING = c("TOP2A", "CDK1", "CENPF")
) %>%
  unlist()

# save the lut of the genes
LUT_df_genes <- list(
  IMMUNE = c("AIF1","TYROBP","FTL","HLA-DRA","TREM2","CX3CR1","C3","CSF1R", "CD74", "RUNX1","C1QB","PTPRC"),
  OLIGOLINEAGE = c("CSPG4","OLIG1","OLIG2", "PDGFRA", "SOX6", "PLP1","SOX10","BCAS1","MBP","MAG"),
  ASTRO = c("AQP4", "GFAP", "VIM","SLC1A2","S100B"),
  # Neu = c("SYT1")
  NEURONS = c("GAD2", "TLE4", "CUX2","SYP", "NEFL","SYT1"),
  NPC = c("NES", "PAX6", "SOX1"),
  CYCLING = c("TOP2A", "CDK1", "CENPF")
) %>%
  lapply(function(x){
    data.frame(gene = x)
  }) %>%
  bind_rows(.id = "category")

Idents(sobj_run04) <- "orig.ident"
df_SOX10_run04 <- AverageExpression(sobj_run04,features = GOI)$RNA %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(names_to = "sample",values_to = "avg_exp",-gene) %>%
  mutate(dataset = "run04")

Idents(sobj_run03) <- "orig.ident"
df_SOX10_run03 <- AverageExpression(sobj_run03,features = GOI)$RNA %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(names_to = "sample",values_to = "avg_exp",-gene) %>%
  mutate(dataset = "run03")

# plot the average expression per gene per run
bind_rows(df_SOX10_run03,df_SOX10_run04) %>%
  ggplot(aes(x=dataset,y=avg_exp)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1))+theme_bw()+facet_wrap(~gene,scales = "free")+theme(strip.background = element_blank())
ggsave("../../out/plot/boxplot_panel_marker_run04vsrun03.pdf",width = 10,height = 10)

# plot the same as an heatmap
# to make the heatmap more readable, scale the expression values genewise
df_wide <- bind_rows(df_SOX10_run03,df_SOX10_run04) %>%
  mutate(dataset_sample = paste0(dataset,"|",sample)) %>%
  group_by(gene) %>%
  mutate(scale_exp = (avg_exp-mean(avg_exp))/sd(avg_exp)) %>%
  dplyr::select(gene,dataset_sample,scale_exp) %>%
  pivot_wider(names_from = dataset_sample,values_from = scale_exp) %>%
  column_to_rownames("gene")

# notice that genes are already ordered with the grouping of the marker set
sum(!(rownames(df_wide) == GOI))

# plot the data as heatmap
meta_df_wide <- data.frame(colname = colnames(df_wide)) %>% 
  separate(col = colname,into = c("dataset","sample"),sep = "\\|",remove = T)

meta_gene <- data.frame(gene = rownames(df_wide)) %>%
  left_join(LUT_df_genes,by = "gene")

column_meta_df_wide <- HeatmapAnnotation(run = meta_df_wide$dataset,
                                         col = list(run = c("run03" = "blue",
                                                            "run04" = "orange")))

row_meta_df_wide <- rowAnnotation(category = meta_gene$category,
                                  col = list(category = c("ASTRO" = "red",
                                                          "CYCLING" = "purple",
                                                          "IMMUNE" = "yellow",
                                                          "NEURONS" = "gray",
                                                          "NPC" = "cyan",
                                                          "OLIGOLINEAGE" = "blue")))

# ht2 <- Heatmap(mat2[match(GOI,rownames(mat2)),], 
#                name = "exp", 
#                column_title = "AMD RPE", 
#                top_annotation = column_ha, 
#                cluster_rows = F, 
#                right_annotation = row_ha, 
#                row_split = rep(c(1,2,3,4),c(2,3,4,7))) 
# pdf("heatmap_test_reinhold.pdf",width = 5,height = 4.5) 
# draw(ht2,heatmap_legend_side = "left",annotation_legend_side = "left") 
# dev.off()

ht2_shr_MG22 <- Heatmap(df_wide,
                        show_column_names = T,
                        raster_by_magick = T,
                        show_row_dend = F, use_raster = T,
                        name = "avg exp",
                        column_title = "BS panel markers",
                        # col = viridis::viridis(option = "turbo",n = 10),
                        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                        right_annotation = row_meta_df_wide,
                        # row_names_gp = gpar(fontsize = 3),
                        top_annotation = column_meta_df_wide,show_row_names = T
)

ht2_shr_MG23 <- Heatmap(df_wide,
                        show_column_names = T,
                        raster_by_magick = T,
                        show_row_dend = F, use_raster = T,
                        name = "avg exp",
                        column_title = "BS panel markers",
                        cluster_rows = F,
                        # col = viridis::viridis(option = "turbo",n = 10),
                        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                        right_annotation = row_meta_df_wide,
                        # row_names_gp = gpar(fontsize = 3),
                        top_annotation = column_meta_df_wide,show_row_names = T
)

pdf("../../out/plot/heatmap_panel_marker_run04vsrun03_clusterRows.pdf",width = 7,height = 12)
draw(ht2_shr_MG22,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(55,2,2, 2), "mm"))
dev.off()

pdf("../../out/plot/heatmap_panel_marker_run04vsrun03_NoclusterRows.pdf",width = 7,height = 12)
draw(ht2_shr_MG23,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(55,2,2, 2), "mm"))
dev.off()

# explore MBP expression in MG --------------------------------------------



# explore MBP expression correlation --------------------------------------



# subclustering of MG -----------------------------------------------------



