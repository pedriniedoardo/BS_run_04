# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(scales)
library(ComplexHeatmap)
library(limma)
library(multtest)
library(metap)
library(ggbreak)
library(patchwork)
# library(lemon)
library(future)
library(ggnewscale)

#  read in the data -------------------------------------------------------
test <- readRDS("../../out/object/data.combined_harmonySkipIntegration_AllSoupX_01000_06000_15.rds")
DefaultAssay(test)
DimPlot(test,raster = T,label = T)

test@meta.data

# get the coordinates
UMAP1_df <- test@reductions$umap@cell.embeddings %>%
  data.frame() %>%
  rownames_to_column(var = "barcodes")

# read in the GOIs --------------------------------------------------------
# df_GOI <- read_csv("data/shortlist_genes.csv")
GOI <- c("GPR17") %>% unique()

# test comparison average expression --------------------------------------
# # check the average expression per cluster cell type
# # pull also the percentage of expressing cells
# # define the grouping variable
# exp_cluster_scale <- DotPlot(test, features = GOI,cluster.idents = T,group.by = c("clusterCellType"))
# exp_cluster <- AverageExpression(test,features = GOI,group.by = c("clusterCellType"))
# # put the data togeher
# exp_cluster_fix <- exp_cluster$RNA %>% 
#   data.frame() %>% 
#   rownames_to_column("features.plot") %>% 
#   pivot_longer(names_to = "id",values_to = "avg.exp.function",-features.plot)
# 
# # join the two dataset
# test_join <- left_join(exp_cluster_scale$data,
#                        exp_cluster_fix,
#                        by = c("features.plot","id"))
# 
# # check if the average expression provided by the average function is the same privided by the dotplot function
# test_join %>% 
#   mutate(delta = avg.exp - avg.exp.function) %>% 
#   ggplot(aes(x=delta)) + geom_histogram()
# # yes the estimate are the same. therefore there is no need to calculate both. I can keep the exp_cluster_scale alone

# plot UMAP ---------------------------------------------------------------
# plot also in the UMAP
FeaturePlot(test,features = GOI,order = T,raster = T)
ggsave("out/image/UMAP_exp_cluster_shortlist_Martina.pdf",width = 8,height = 4)

FeaturePlot(test,features = GOI,order = T,split.by = "treat")
ggsave("out/image/UMAP_exp_cluster_shortlist2_Martina.pdf",width = 6,height = 8)

# calculate average expression for different metadata ---------------------
exp_cluster <- DotPlot(test, features = GOI,cluster.idents = T,group.by = c("seurat_clusters"))
# exp_cluster_sample <- DotPlot(test, features = GOI,cluster.idents = T,group.by = c("clusterCellType_sample"))

# # check the average expression per cluster and pathology
exp_pathology <- DotPlot(test, features = GOI,cluster.idents = T,group.by = c("treat"))
# exp_pathology_sample <- DotPlot(test, features = GOI,cluster.idents = T,group.by = c("pathology_sample"))
# 
# # check the average expression per cluster and pathology
exp_cluster_pathology <- DotPlot(test, features = GOI,cluster.idents = T,group.by = c("cluster_pathology"))
# exp_cluster_pathology_sample <- DotPlot(test, features = GOI,cluster.idents = T,group.by = c("clusterCellType_pathology_sample"))

test@meta.data$group <- test@meta.data |> 
  mutate(group = paste0(treat,".",seurat_clusters)) |> 
  pull(group)

tt <- AverageExpression(test,group.by = "group",features = GOI) %>%
  .$RNA |> 
  data.frame() |> 
  mutate(gene = GOI) |> 
  pivot_longer(names_to = "sample",values_to = "avg_exp",-gene) |> 
  mutate(seurat_clusters = str_extract(sample,pattern = "\\d+")) |> 
  mutate(treat = str_remove_all(sample,pattern = "\\.\\d+"))

# add the annotation for the number fo cells per cell type
lut_cells <- test@meta.data |> 
  group_by(treat,seurat_clusters) |> 
  summarise(n = n()) |>
  ungroup() |> 
  group_by(treat) |> 
  mutate(tot = sum(n)) |> 
  mutate(prop = n/tot)

rank_clusters <- tt |> 
  left_join(lut_cells,by = c("treat","seurat_clusters")) |> 
  group_by(seurat_clusters) |> 
  summarise(avg_prop = mean(prop)) |> 
  mutate(rank = rank(avg_prop))

tt |> 
  left_join(lut_cells,by = c("treat","seurat_clusters")) |> 
  left_join(rank_clusters,by = c("seurat_clusters")) |> 
  mutate(seurat_clusters = fct_reorder(seurat_clusters,rank,.desc = T)) |> 
  ggplot(aes(x=treat,y=avg_exp,col=prop))+geom_point()+facet_wrap(~seurat_clusters)+theme_bw()+theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 90))+scale_color_viridis_c(option = "turbo")

# plotting clusterCellType ------------------------------------------------
# plot the averages expressions

# col_pal <- c("#E6E6E6","#ffff09","#c1ce08","#446d05","#053c03","#4D4D4D","#06a8ce","#033b6d","#ff0ed7","#9a0404")
# show_col(col_pal)

# cluster the genes
# calculate the genes that are above the detection level avg.exp.scaled > 0 & `% Expressing` > 1
GOI_present <- exp_cluster$data %>%
  mutate(`% Expressing` = (pct.exp)) %>%
  # filter out the genes for wich no sample and per gene is crossing the threshold for the  plotting
  group_by(id,features.plot) %>% 
  # mutate(test = avg.exp.scaled > 0 & `% Expressing` > 1) %>%
  mutate(test = `% Expressing` > 1) %>%
  ungroup() %>% 
  group_by(features.plot) %>% 
  summarise(tot = sum(test)) %>% 
  filter(tot != 0)

# make data square to calculate euclidean distance only on the GOI filtered
exp_cluster_wide <- exp_cluster$data %>%
  # filter out the genes for wich no sample and per gene is crossing the threshold for the  plotting
  dplyr::filter(features.plot %in% GOI_present$features.plot) %>% 
  dplyr::select(features.plot,id,avg.exp.scaled) %>%  # drop unused columns to faciliate widening
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  data.frame() %>% 
  column_to_rownames("features.plot")

dim(exp_cluster_wide)

# calculate the distance
clust_exp_cluster_gene <- hclust(dist(exp_cluster_wide %>% as.matrix())) # hclust with distance matrix
clust_exp_cluster_sample <- hclust(dist(exp_cluster_wide %>% as.matrix() %>% t())) # hclust with distance matrix

# generate the dendrogram
ddgram_exp_cluster_gene <- as.dendrogram(clust_exp_cluster_gene) # create dendrogram
ddgram_exp_cluster_sample <- as.dendrogram(clust_exp_cluster_sample) # create dendrogram
ggtree_plot_exp_cluster_gene <- ggtree::ggtree(ddgram_exp_cluster_gene,hang = 0)
ggtree_plot_exp_cluster_sample <- ggtree::ggtree(ddgram_exp_cluster_sample,hang = 0) + layout_dendrogram()
ggtree_plot_exp_cluster_gene
ggtree_plot_exp_cluster_sample

# plot the data
dotplot_exp_cluster <- exp_cluster$data %>%
  dplyr::filter(features.plot %in% GOI_present$features.plot) %>% 
  # order the features and the samples according to the clustering
  mutate(`% Expressing` = (pct.exp),
         features.plot = factor(features.plot, levels = clust_exp_cluster_gene$labels[clust_exp_cluster_gene$order]),
         id = factor(id, levels = clust_exp_cluster_sample$labels[clust_exp_cluster_sample$order])) %>% 
  # filter(avg.exp.scaled > 0, `% Expressing` > 1) %>% 
  filter(`% Expressing` > 1) %>% 
  ggplot(aes(x=id, y = features.plot, color = avg.exp.scaled, size = `% Expressing`)) + 
  geom_point() +
  scale_size(range = c(0, 6))+
  scale_color_viridis_c(name = 'scale counts',option = "turbo") + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('')+
  scale_y_discrete(position = "right")
# theme(axis.ticks = element_blank(),
#         plot.margin = margin(0, 1, 0, 0, "cm"))

plot_grid(ggtree_plot_exp_cluster_gene, NULL, dotplot_exp_cluster, nrow = 1, rel_widths = c(0.5,-0.1, 2), align = 'h')
ggsave("out/image/heatmap_exp_cluster_shortlist_martina.pdf",width = 8,height = 5)

# try to add also the clustering of the samples
# column dendrogram
ggtree_plot_sample <- ggtree_plot_exp_cluster_sample + xlim2(dotplot_exp_cluster)
# row dendrogram
ggtree_plot_gene <- ggtree_plot_exp_cluster_gene + ylim2(dotplot_exp_cluster)

# define the layout
layout <- "
A
B
B
B
B
B
B
"
G1 <- (plot_spacer() / ggtree_plot_gene) + plot_layout(design = layout)
G2 <- (ggtree_plot_sample / dotplot_exp_cluster) + plot_layout(design = layout)

G1 | G2
ggsave("out/image/heatmap_exp_cluster_shortlist_martina2.pdf",width = 8,height = 6)
# 
# see the effect of not scaling counts ------------------------------------
# cluster the genes
# calculate the genes that are above the detection level avg.exp.scaled > 0 & `% Expressing` > 1

# make data square to calculate euclidean distance only on the GOI filtered
exp_cluster_wide_test <- exp_cluster$data %>%
  # filter out the genes for wich no sample and per gene is crossing the threshold for the  plotting
  dplyr::filter(features.plot %in% GOI_present$features.plot) %>%
  dplyr::select(features.plot,id,avg.exp) %>%  # drop unused columns to faciliate widening
  pivot_wider(names_from = id, values_from = avg.exp) %>%
  data.frame() %>%
  column_to_rownames("features.plot")

dim(exp_cluster_wide_test)

# calculate the distance
clust_exp_cluster_gene_test <- hclust(dist(exp_cluster_wide_test %>% as.matrix())) # hclust with distance matrix
clust_exp_cluster_sample_test <- hclust(dist(exp_cluster_wide_test %>% as.matrix() %>% t())) # hclust with distance matrix

# generate the dendrogram
ddgram_exp_cluster_gene_test <- as.dendrogram(clust_exp_cluster_gene_test) # create dendrogram
ddgram_exp_cluster_sample_test <- as.dendrogram(clust_exp_cluster_sample_test) # create dendrogram
ggtree_plot_exp_cluster_gene_test <- ggtree::ggtree(ddgram_exp_cluster_gene_test,hang = 0)
ggtree_plot_exp_cluster_sample_test <- ggtree::ggtree(ddgram_exp_cluster_sample_test,hang = 0) + layout_dendrogram()
ggtree_plot_exp_cluster_gene_test
ggtree_plot_exp_cluster_sample_test

# plot the data
dotplot_exp_cluster_test <- exp_cluster$data %>%
  dplyr::filter(features.plot %in% GOI_present$features.plot) %>%
  # order the features and the samples according to the clustering
  mutate(`% Expressing` = (pct.exp),
         features.plot = factor(features.plot, levels = clust_exp_cluster_gene_test$labels[clust_exp_cluster_gene_test$order]),
         id = factor(id, levels = clust_exp_cluster_sample_test$labels[clust_exp_cluster_sample_test$order])) %>%
  # filter(avg.exp.scaled > 0, `% Expressing` > 1) %>%
  filter(`% Expressing` > 1) %>%
  ggplot(aes(x=id, y = features.plot, color = avg.exp, size = `% Expressing`)) +
  geom_point() +
  scale_size(range = c(0, 6))+
  scale_color_viridis_c(name = 'avg counts',option = "turbo") +
  cowplot::theme_cowplot() +
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('')+
  scale_y_discrete(position = "right")
# theme(axis.ticks = element_blank(),
#         plot.margin = margin(0, 1, 0, 0, "cm"))

plot_grid(ggtree_plot_exp_cluster_gene_test, NULL, dotplot_exp_cluster_test, nrow = 1, rel_widths = c(0.5,-0.1, 2), align = 'h')
ggsave("out/image/heatmap_exp_cluster_shortlist_testCounts_martina.pdf",width = 6,height = 5)

# plotting clusterCellType pathology --------------------------------------
# plot the averages expressions

# col_pal <- c("#E6E6E6","#ffff09","#c1ce08","#446d05","#053c03","#4D4D4D","#06a8ce","#033b6d","#ff0ed7","#9a0404")
# show_col(col_pal)

# cluster the genes
# calculate the genes that are above the detection level avg.exp.scaled > 0 & `% Expressing` > 1
GOI_present2 <- exp_cluster_pathology$data %>%
  mutate(`% Expressing` = (pct.exp)) %>%
  # filter out the genes for wich no sample and per gene is crossing the threshold for the  plotting
  group_by(id,features.plot) %>%
  # mutate(test = avg.exp.scaled > 0 & `% Expressing` > 1) %>%
  mutate(test = `% Expressing` > 1) %>%
  ungroup() %>%
  group_by(features.plot) %>%
  summarise(tot = sum(test)) %>%
  filter(tot != 0)

# make data square to calculate euclidean distance only on the GOI filtered
exp_cluster_pathology_wide <- exp_cluster_pathology$data %>%
  # filter out the genes for wich no sample and per gene is crossing the threshold for the  plotting
  dplyr::filter(features.plot %in% GOI_present2$features.plot) %>%
  dplyr::select(features.plot,id,avg.exp.scaled) %>%  # drop unused columns to faciliate widening
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>%
  data.frame() %>%
  column_to_rownames("features.plot")

dim(exp_cluster_pathology_wide)

# calculate the distance
clust_exp_cluster_pathology_gene <- hclust(dist(exp_cluster_pathology_wide %>% as.matrix())) # hclust with distance matrix
clust_exp_cluster_pathology_sample <- hclust(dist(exp_cluster_pathology_wide %>% as.matrix() %>% t())) # hclust with distance matrix

# generate the dendrogram
ddgram_exp_cluster_pathology_gene <- as.dendrogram(clust_exp_cluster_pathology_gene) # create dendrogram
ddgram_exp_cluster_pathology_sample <- as.dendrogram(clust_exp_cluster_pathology_sample) # create dendrogram
ggtree_plot_exp_cluster_pathology_gene <- ggtree::ggtree(ddgram_exp_cluster_pathology_gene,hang = 0)
ggtree_plot_exp_cluster_pathology_sample <- ggtree::ggtree(ddgram_exp_cluster_pathology_sample,hang = 0) + layout_dendrogram()
ggtree_plot_exp_cluster_pathology_gene
ggtree_plot_exp_cluster_pathology_sample

# plot the data
dotplot_exp_cluster_pathology <- exp_cluster_pathology$data %>%
  # fix the name to match the downstream clustrering
  mutate(id = str_replace_all(id,pattern = "-",replacement = "\\.")) %>%
  dplyr::filter(features.plot %in% GOI_present2$features.plot) %>%
  # order the features and the samples according to the clustering
  mutate(`% Expressing` = (pct.exp),
         features.plot = factor(features.plot, levels = clust_exp_cluster_pathology_gene$labels[clust_exp_cluster_pathology_gene$order]),
         id = factor(id, levels = clust_exp_cluster_pathology_sample$labels[clust_exp_cluster_pathology_sample$order])) %>%
  # filter(avg.exp.scaled > 0, `% Expressing` > 1) %>%
  filter(`% Expressing` > 1) %>%
  ggplot(aes(x=id, y = features.plot, color = avg.exp.scaled, size = `% Expressing`)) +
  geom_point() +
  scale_size(range = c(0, 6))+
  scale_color_viridis_c(name = 'scale counts',option = "turbo") +
  cowplot::theme_cowplot() +
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('')+
  scale_y_discrete(position = "right")
# theme(axis.ticks = element_blank(),
#         plot.margin = margin(0, 1, 0, 0, "cm"))

plot_grid(ggtree_plot_exp_cluster_pathology_gene, NULL, dotplot_exp_cluster_pathology, nrow = 1, rel_widths = c(0.5,-0.05, 2), align = 'h')
ggsave("out/image/heatmap_exp_cluster_pathology_shortlist_clusterSample_martina.pdf",width = 15,height = 5)

# try the same as above but discard the sample ordering
# plot the data
dotplot_exp_cluster_pathology2 <- exp_cluster_pathology$data %>%
  # fix the name to match the downstream clustrering
  mutate(id = str_replace_all(id,pattern = "-",replacement = "\\.")) %>%
  dplyr::filter(features.plot %in% GOI_present2$features.plot) %>%
  # order the features and the samples according to the clustering
  mutate(`% Expressing` = (pct.exp),
         features.plot = factor(features.plot, levels = clust_exp_cluster_pathology_gene$labels[clust_exp_cluster_pathology_gene$order])) %>%
  # filter(avg.exp.scaled > 0, `% Expressing` > 1) %>%
  filter(`% Expressing` > 1) %>%
  ggplot(aes(x=id, y = features.plot, color = avg.exp.scaled, size = `% Expressing`)) +
  geom_point() +
  scale_size(range = c(0, 6))+
  scale_color_viridis_c(name = 'scale counts',option = "turbo") +
  cowplot::theme_cowplot() +
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('')+
  scale_y_discrete(position = "right")
# theme(axis.ticks = element_blank(),
#         plot.margin = margin(0, 1, 0, 0, "cm"))

plot_grid(ggtree_plot_exp_cluster_pathology_gene, NULL, dotplot_exp_cluster_pathology2, nrow = 1, rel_widths = c(0.5,-0.05, 2), align = 'h')
ggsave("out/image/heatmap_exp_cluster_pathology_shortlist_martina.pdf",width = 15,height = 5)

# try to add also the clustering of the samples
# column dendrogram
ggtree_plot_sample2 <- ggtree_plot_exp_cluster_pathology_sample + xlim2(dotplot_exp_cluster_pathology)
# row dendrogram
ggtree_plot_gene2 <- ggtree_plot_exp_cluster_pathology_gene + ylim2(dotplot_exp_cluster_pathology)

# define the layout
layout <- "
A
B
B
B
B
B
B
"
G12 <- (plot_spacer() / ggtree_plot_gene2) + plot_layout(design = layout)
G22 <- (ggtree_plot_sample2 / dotplot_exp_cluster_pathology) + plot_layout(design = layout)

layout2 <- "
ABBBB
"

(G12 | G22) + plot_layout(design = layout2)
ggsave("out/image/heatmap_exp_cluster_pathology_shortlist_clusterSample2_martina.pdf",width = 10,height = 6)

# plotting for pathology --------------------------------------------------
# plot the averages expressions

# col_pal <- c("#E6E6E6","#ffff09","#c1ce08","#446d05","#053c03","#4D4D4D","#06a8ce","#033b6d","#ff0ed7","#9a0404")
# show_col(col_pal)

# cluster the genes
# calculate the genes that are above the detection level avg.exp.scaled > 0 & `% Expressing` > 1
GOI_present3 <- exp_pathology$data %>%
  mutate(`% Expressing` = (pct.exp)) %>%
  # filter out the genes for wich no sample and per gene is crossing the threshold for the  plotting
  group_by(id,features.plot) %>%
  # mutate(test = avg.exp.scaled > 0 & `% Expressing` > 1) %>%
  mutate(test = `% Expressing` > 1) %>%
  ungroup() %>%
  group_by(features.plot) %>%
  summarise(tot = sum(test)) %>%
  filter(tot != 0)

# make data square to calculate euclidean distance only on the GOI filtered
exp_pathology_wide <- exp_pathology$data %>%
  # filter out the genes for wich no sample and per gene is crossing the threshold for the  plotting
  dplyr::filter(features.plot %in% GOI_present3$features.plot) %>%
  dplyr::select(features.plot,id,avg.exp.scaled) %>%  # drop unused columns to faciliate widening
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>%
  data.frame() %>%
  column_to_rownames("features.plot")

dim(exp_pathology_wide)

# calculate the distance
clust_exp_pathology_gene <- hclust(dist(exp_pathology_wide %>% as.matrix())) # hclust with distance matrix
clust_exp_pathology_sample <- hclust(dist(exp_pathology_wide %>% as.matrix() %>% t())) # hclust with distance matrix

# generate the dendrogram
ddgram_exp_pathology_gene <- as.dendrogram(clust_exp_pathology_gene) # create dendrogram
ddgram_exp_pathology_sample <- as.dendrogram(clust_exp_pathology_sample) # create dendrogram
ggtree_plot_exp_pathology_gene <- ggtree::ggtree(ddgram_exp_pathology_gene,hang = 0)
ggtree_plot_exp_pathology_sample <- ggtree::ggtree(ddgram_exp_pathology_sample,hang = 0) + layout_dendrogram()
ggtree_plot_exp_pathology_gene
ggtree_plot_exp_pathology_sample

# plot the data
dotplot_exp_pathology <- exp_pathology$data %>%
  # fix the name to match the downstream clustrering
  mutate(id = str_replace_all(id,pattern = "-",replacement = "\\.")) %>%
  dplyr::filter(features.plot %in% GOI_present3$features.plot) %>%
  # order the features and the samples according to the clustering
  mutate(`% Expressing` = (pct.exp),
         features.plot = factor(features.plot, levels = clust_exp_pathology_gene$labels[clust_exp_pathology_gene$order]),
         id = factor(id, levels = clust_exp_pathology_sample$labels[clust_exp_pathology_sample$order])) %>%
  # filter(avg.exp.scaled > 0, `% Expressing` > 1) %>%
  filter(`% Expressing` > 1) %>%
  ggplot(aes(x=id, y = features.plot, color = avg.exp.scaled, size = `% Expressing`)) +
  geom_point() +
  scale_size(range = c(0, 6))+
  scale_color_viridis_c(name = 'scale counts',option = "turbo") +
  cowplot::theme_cowplot() +
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('')+
  scale_y_discrete(position = "right")
# theme(axis.ticks = element_blank(),
#         plot.margin = margin(0, 1, 0, 0, "cm"))

plot_grid(ggtree_plot_exp_pathology_gene, NULL, dotplot_exp_pathology, nrow = 1, rel_widths = c(0.5,-0.1, 2), align = 'h')
ggsave("out/image/heatmap_exp_pathology_shortlist_martina.pdf",width = 6,height = 5)

# try to add also the clustering of the samples
# column dendrogram
ggtree_plot_sample3 <- ggtree_plot_exp_pathology_sample + xlim2(dotplot_exp_pathology)
# row dendrogram
ggtree_plot_gene3 <- ggtree_plot_exp_pathology_gene + ylim2(dotplot_exp_pathology)

# define the layout
layout <- "
A
B
B
B
B
B
B
"
G13 <- (plot_spacer() / ggtree_plot_gene3) + plot_layout(design = layout)
G23 <- (ggtree_plot_sample3 / dotplot_exp_pathology) + plot_layout(design = layout)

G13 | G23
ggsave("out/image/heatmap_exp_pathology_shortlist2_martina.pdf",width = 6,height = 6)
