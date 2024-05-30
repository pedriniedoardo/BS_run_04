# AIM ---------------------------------------------------------------------
# the aim is to plot the data after integration.
# this integration was run by skipping the seurat integration (by merging the matrices) and running Harmony

# libraries ---------------------------------------------------------------
library(harmony)
library(Seurat)
library(dplyr)
library(cowplot)
library(tidyverse)
library(ggrepel)
library(scales)
library(RColorBrewer)
library(SeuratWrappers)
library(dittoSeq)

# # read in the data --------------------------------------------------------
data.combined <- readRDS("../../out/object/data.combined_harmonySkipIntegration_AllSoupX_01000_06000_15.rds")
Idents(data.combined)<-"seurat_clusters"

# plots -------------------------------------------------------------------
# main umap
plot03 <- DimPlot(data.combined, reduction = "umap", group.by = "seurat_clusters",label = T,raster = T)
ggsave(plot=plot03,"../../out/plot/UMAPCluster_harmonySkipIntegration_AllSoupX_01000_06000_15.pdf",width = 6,height = 5)

plot03a <- DimPlot(data.combined, reduction = "umap", group.by = "seurat_clusters",label = T,raster = T,split.by = "ID",ncol = 3)
ggsave(plot=plot03a,"../../out/plot/UMAPClusterSplit_harmonySkipIntegration_AllSoupX_01000_06000_15.pdf",width = 13,height = 12)

# main umap cell cycle
plot04 <- DimPlot(data.combined, reduction = "umap", group.by = "Phase",raster = T)
ggsave(plot=plot04,"../../out/plot/UMAPPhase_harmonySkipIntegration_AllSoupX_01000_06000_15.pdf",width = 5,height = 4)

# split by sample
df_meta <- data.combined@meta.data %>%
  rownames_to_column("barcode")
df_UMAP <- data.combined@reductions$umap@cell.embeddings %>%
  data.frame() %>%
  rownames_to_column("barcode")

data2 <- left_join(df_UMAP,df_meta,"barcode")
data2_avg <- data2 %>% group_by(seurat_clusters) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)

# single plot
plot05 <- data2 %>%
  mutate(orig.ident = factor(orig.ident)) %>%
  # arrange(signature_score1) %>%
  # mutate(gene = "Ptx3") %>%
  ggplot() + 
  geom_point(aes(x=UMAP_1,y=UMAP_2,col=seurat_clusters),alpha=0.5,size=0.1)+
  geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = seurat_clusters)) +
  guides(colour = guide_legend(override.aes = list(size=5,alpha=1),ncol=2))+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))+
  theme_bw()
ggsave(plot=plot05,"../../out/plot/UMAPClusterggplot_harmonySkipIntegration_AllSoupX_01000_06000_15.pdf",width = 6,height = 4)
ggsave(plot=plot05,"../../out/plot/UMAPClusterggplot_harmonySkipIntegration_AllSoupX_01000_06000_15.png",width = 6,height = 4)

# save the same with theme cowplot
plot05_alt <- data2 %>%
  mutate(orig.ident = factor(orig.ident)) %>%
  # arrange(signature_score1) %>%
  # mutate(gene = "Ptx3") %>%
  ggplot() + 
  geom_point(aes(x=UMAP_1,y=UMAP_2,col=seurat_clusters),alpha=0.5,size=0.1)+
  geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = seurat_clusters)) +
  guides(colour = guide_legend(override.aes = list(size=5,alpha=1),ncol=2))+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))+
  theme_cowplot()
ggsave(plot=plot05_alt,"../../out/plot/UMAPClusterggplot2_harmonySkipIntegration_AllSoupX_01000_06000_15.pdf",width = 6,height = 4)
ggsave(plot=plot05_alt,"../../out/plot/UMAPClusterggplot2_harmonySkipIntegration_AllSoupX_01000_06000_15.png",width = 6,height = 4)

# single plot
plot051 <- data2 %>%
  mutate(orig.ident = factor(orig.ident)) %>%
  arrange(percent.ribo) %>%
  # mutate(gene = "Ptx3") %>%
  ggplot() + 
  geom_point(aes(x=UMAP_1,y=UMAP_2,col=percent.ribo),alpha=0.5,size=0.1)+
  geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = seurat_clusters)) +
  scale_color_viridis_c(option = "turbo")+
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))
ggsave(plot = plot051,"../../out/plot/UMAPRibo_harmonySkipIntegration_AllSoupX_01000_06000_15.pdf",width = 6,height = 4)
ggsave(plot = plot051,"../../out/plot/UMAPRibo_harmonySkipIntegration_AllSoupX_01000_06000_15.png",width = 6,height = 4)

# single plot
plot052 <- data2 %>%
  mutate(orig.ident = factor(orig.ident)) %>%
  arrange(percent.mt) %>%
  # mutate(gene = "Ptx3") %>%
  ggplot() + 
  geom_point(aes(x=UMAP_1,y=UMAP_2,col=percent.mt),alpha=0.5,size=0.1)+
  geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = seurat_clusters)) +
  scale_color_viridis_c(option = "turbo")+
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))
ggsave(plot = plot052,"../../out/plot/UMAPMito_harmonySkipIntegration_AllSoupX_01000_06000_15.pdf",width = 6,height = 4)
ggsave(plot = plot052,"../../out/plot/UMAPMito_harmonySkipIntegration_AllSoupX_01000_06000_15.png",width = 6,height = 4)

# single plot
plot053 <- data2 %>%
  mutate(orig.ident = factor(orig.ident)) %>%
  arrange(nCount_RNA) %>%
  # mutate(gene = "Ptx3") %>%
  ggplot() + 
  geom_point(aes(x=UMAP_1,y=UMAP_2,col=nCount_RNA),alpha=0.5,size=0.1)+
  geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = seurat_clusters)) +
  scale_color_viridis_c(option = "turbo")+
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))
ggsave(plot = plot053,"../../out/plot/UMAPnCount_harmonySkipIntegration_AllSoupX_01000_06000_15.pdf",width = 6,height = 4)
ggsave(plot = plot053,"../../out/plot/UMAPnCount_harmonySkipIntegration_AllSoupX_01000_06000_15.png",width = 6,height = 4)

# split by sample
plot06 <- data2 %>%
  mutate(orig.ident = factor(orig.ident)) %>%
  # arrange(signature_score1) %>%
  # mutate(gene = "Ptx3") %>%
  ggplot() + 
  geom_point(aes(x=UMAP_1,y=UMAP_2,col=seurat_clusters),alpha=0.5,size=0.1)+
  geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = seurat_clusters)) +
  facet_wrap(~ID,ncol=3)+
  guides(colour = guide_legend(override.aes = list(size=5,alpha=1),ncol=2))+
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))
ggsave(plot=plot06,"../../out/plot/UMAPClusterSplit_harmonySkipIntegration_AllSoupX_01000_06000_15.pdf",width = 13,height = 12)
ggsave(plot=plot06,"../../out/plot/UMAPClusterSplit_harmonySkipIntegration_AllSoupX_01000_06000_15.png",width = 13,height = 12)

# proportion of cell per cluster
df_summary <- df_meta %>%
  group_by(ID,seurat_clusters) %>%
  summarise(n = n()) %>%
  group_by(ID) %>%
  mutate(tot=sum(n)) %>%
  ungroup() %>%
  mutate(prop=n/tot)
write_tsv(df_summary,"../../out/table/summary_harmonySkipIntegration_AllSoupX_01000_06000_15.tsv")

# define a convenient palette of colors
show_col(hue_pal()(9))
RColorBrewer::display.brewer.all()
col_pal <- RColorBrewer::brewer.pal(name = "Paired",n = 9)
# col_pal <- c("#E6E6E6","#ffff09","#c1ce08","#446d05","#053c03","#4D4D4D","#06a8ce","#033b6d","#ff0ed7","#9a0404")
show_col(col_pal)

plot07 <- df_summary %>%
  mutate(ID = factor(ID)) %>%
  # mutate(treat = factor(treat,levels = c("Young","Adult","Old"))) %>%
  ggplot(aes(x=seurat_clusters,y=prop,fill=ID))+geom_col()+theme_bw()+
  scale_fill_manual(values = col_pal)
ggsave(plot=plot07,"../../out/plot/summary_harmonySkipIntegration_AllSoupX_01000_06000_15.pdf",width = 7,height = 4)

plot08 <- df_summary %>%
  mutate(ID = factor(ID)) %>%
  ggplot(aes(x=seurat_clusters,y=prop,fill=ID))+geom_col(position = "dodge")+theme_bw()+scale_fill_manual(values = col_pal)
ggsave(plot=plot08,"../../out/plot/summaryDodge_harmonySkipIntegration_AllSoupX_01000_06000_15.pdf",width = 8,height = 4)

# # add the grouping variable
# read_csv("LUT_sample.csv") %>%
#   left_join(df_summary,by = c("sample"="orig.sample")) %>%
#   mutate(treat = factor(treat,levels = c("Young","Adult","Old"))) %>%
#   ggplot(aes(x=seurat_clusters,y=prop,fill=treat))+
#   geom_boxplot(outlier.shape = NA)+
#   geom_point(position=position_jitterdodge(dodge.width=0.9,jitter.width = 0.2))+
#   theme_bw()+scale_y_log10()
# ggsave("images/fixed_threshold_doubletSinglet//summary_group_data.combined_NOT_annotated_norm_fix_regress_CC_resolution_log_DoubletSinglet.pdf",width = 8,height = 4)

# read_csv("LUT_sample.csv") %>%
#   left_join(df_summary,by = c("sample"="orig.sample")) %>%
#   mutate(treat = factor(treat,levels = c("Young","Adult","Old"))) %>%
#   ggplot(aes(x=seurat_clusters,y=prop,fill=treat))+
#   geom_boxplot(outlier.shape = NA)+
#   geom_point(position=position_jitterdodge(dodge.width=0.9,jitter.width = 0.2))+
#   theme_bw()
# ggsave("images/fixed_threshold_doubletSinglet/summary_group_data.combined_NOT_annotated_norm_fix_regress_CC_resolution_lin_DoubletSinglet.pdf",width = 8,height = 4)
shortlist_features_list <- list(
  IMMUNE = c("AIF1","TYROBP","FTL","HLA-DRA","TREM2","CX3CR1","C3","CSF1R", "CD74", "RUNX1","C1QB","PTPRC"),
  OLIGOLINEAGE = c("CSPG4","OLIG1","OLIG2", "PDGFRA", "SOX6", "PLP1","SOX10","BCAS1","MBP","MAG"),
  ASTRO = c("AQP4", "GFAP", "VIM","SLC1A2","S100B"),
  # Neu = c("SYT1")
  NEURONS = c("GAD2", "TLE4", "CUX2","SYP", "NEFL","SYT1"),
  NPC = c("NES", "PAX6", "SOX1"),
  CYCLING = c("TOP2A", "CDK1", "CENPF")
)

# martina asked to remove FTL
shortlist_features_list2 <- list(
  IMMUNE = c("AIF1","TYROBP","HLA-DRA","TREM2","CX3CR1","C3","CSF1R", "CD74", "RUNX1","C1QB","PTPRC"),
  OLIGOLINEAGE = c("CSPG4","OLIG1","OLIG2", "PDGFRA", "SOX6", "PLP1","SOX10","BCAS1","MBP","MAG"),
  ASTRO = c("AQP4", "GFAP", "VIM","SLC1A2","S100B"),
  # Neu = c("SYT1")
  NEURONS = c("GAD2", "TLE4", "CUX2","SYP", "NEFL","SYT1"),
  NPC = c("NES", "PAX6", "SOX1"),
  CYCLING = c("TOP2A", "CDK1", "CENPF")
)

shortlist_features_list_long <- list(
  IMMUNE = c("IGHG1","CD38","CD8A","CD2","SKAP1","LYVE1","CD163","MRC1","LINGO1","HSPA1A","MOBP","CD22","CD83","HIF1A","VEGFA","SOD1","TREM2","CX3CR1","P2RY12","C3","CSF1R", "CD74", "RUNX1","C1QB","PTPRC"),
  OLIGOLINEAGE = c("PLP1","MOG","PPP1R16B","TNS3","HMGB1","CD81","B2M","C1QL1","HLA-A","HLA-C","NLGN4X","OLIG1","OLIG2","CSPG4", "PDGFRA", "SOX6", "SOX10","BCAS1","MBP","MAG"),
  ASTRO = c("AQP4", "GFAP", "CD44", "AQP1", "VIM","APOE", "VCAN", "STAT3", "ABCA1", "TNC", "SDC4","SLC1A2","S100B"),
  NEURONS = c("GAD2", "PVALB", "SV2C", "VIP", "TLE4", "CUX2", "THY1", "SLC17A7", "NRGN", "SATB2", "RORB", "SST", "STX1A", "STX1B", "SYP", "TH", "NEFL","SYT1"),
  NPC = c("NES", "PAX6", "SOX1"),
  CYCLING = c("TOP2A", "CDK1", "CENPF"),
  ENDO = c("VWF","CDH5","TEK","PECAM1","FLT1","KDR","NOS3","MCAM","MMRN1","CLDN5","BMX","ANGPT2","GJA4","TIE1","ROBO4","ECSCR"),
  PERICYTE = c("PDGFRB","DES","ACTA2","ANPEP","RGS5","ABCC9","KCNJ8","CD248","DLK1","NT5E","ANGPT1"))

# plot the shortlisted feature per cluster
# notice that this is done only on the subset of the young (control) cells
test_short <- DotPlot(data.combined, features = shortlist_features_list, dot.scale = 8,cluster.idents = T) +
  RotatedAxis()
ggsave(plot = test_short,"../../out/plot/Dotplot_harmonySkipIntegration_AllSoupX_01000_06000_15.pdf",width = 13,height = 5)

test_short2 <- DotPlot(data.combined, features = shortlist_features_list2, dot.scale = 8,cluster.idents = T) +
  RotatedAxis()
ggsave(plot = test_short2,"../../out/plot/Dotplot2_harmonySkipIntegration_AllSoupX_01000_06000_15.pdf",width = 13,height = 5)

test_long <- DotPlot(data.combined, features = shortlist_features_list_long, dot.scale = 8,cluster.idents = T) +
  RotatedAxis()
ggsave(plot=test_long,"../../out/plot/DotplotLong_harmonySkipIntegration_AllSoupX_01000_06000_15.pdf",width = 25,height = 5)

df_meta %>%
  group_by(ID) %>%
  summarise(n=n())

# Identify conserved cell type markers ------------------------------------
# data
DefaultAssay(data.combined) <- "RNA"
# notice that in this case the data object of the RNA slot is already filled with the normalzied data, therefore in this case (following the Normalize workfloe for the integration) there is no need to run the NormalizeData on the RNA slot of the integrated object
# sobj_total_h@assays$RNA@data[20:50,1:10]
# dim(sobj_total_h@assays$RNA@data)
# 
# # scale the data see the note on evernote on why this can be done also at this point. this is needed because the scale.data is empty
# sobj_total_h@assays$RNA@scale.data
# all.genes <- rownames(sobj_total_h)
# sobj_total_h <- ScaleData(sobj_total_h,vars.to.regress = c("percent.mt","nCount_RNA","S.Score","G2M.Score"))
# # confirm now the slot is loaded
# sobj_total_h@assays$RNA@scale.data[1:10,1:10]

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
sobj_total_h.markers <- RunPrestoAll(data.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# save the table of all markers
sobj_total_h.markers %>%
  write_tsv("../../out/table/FindAllMarkers_harmonySkipIntegration_AllSoupX_01000_06000_15.tsv")

# sobj_total_h.markers <- read_tsv("../../out/table/FindAllMarkers_harmonySkipIntegration_AllSoupX_01000_06000_15.tsv")

# try plotting the top markers
top_specific_markers <- sobj_total_h.markers %>%
  # filter ribosomal and mt genes
  filter(str_detect(gene,pattern = "^MT-",negate=T)) %>%
  filter(str_detect(gene,pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA",negate=T)) %>%
  filter(str_detect(gene,pattern = "^HB[^(P)]",negate=T)) %>%
  group_by(cluster) %>%
  top_n(3, avg_log2FC)

# And generate e.g. a dotplot:
dittoSeq::dittoDotPlot(data.combined,
                       vars = unique(top_specific_markers$gene), 
                       group.by = "seurat_clusters")+scale_color_viridis_c(option = "turbo",name="relative \nexpression")
