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

# # read in the data --------------------------------------------------------
data.combined <- readRDS("../../out/object/sobj_processed_donor.rds")
Idents(data.combined)<-"expertAnno.l1"

# plots -------------------------------------------------------------------
# main umap
plot03 <- DimPlot(data.combined, reduction = "umap", group.by = "expertAnno.l1",label = T,raster = T)
ggsave(plot=plot03,"../../out/image/UMAPCellID_sobj_processed_donor.pdf",width = 5,height = 4)

# # main umap cell cycle
# plot04 <- DimPlot(data.combined, reduction = "umap", group.by = "Phase",raster = T)
# ggsave(plot=plot04,"../../out/image/UMAPPhase_sobj_processed_donor.pdf",width = 5,height = 4)

# split by sample
df_meta <- data.combined@meta.data %>%
  rownames_to_column("barcode")
df_UMAP <- data.combined@reductions$umap@cell.embeddings %>%
  data.frame() %>%
  rownames_to_column("barcode")

data2 <- left_join(df_UMAP,df_meta,"barcode")
data2_avg <- data2 %>% group_by(expertAnno.l1) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)

# single plot
plot05 <- data2 %>%
  mutate(orig.ident = factor(orig.ident)) %>%
  # arrange(signature_score1) %>%
  # mutate(gene = "Ptx3") %>%
  ggplot() + 
  geom_point(aes(x=UMAP_1,y=UMAP_2,col=expertAnno.l1),alpha=0.5,size=0.1)+
  geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = expertAnno.l1)) +
  guides(colour = guide_legend(override.aes = list(size=5,alpha=1),ncol=2))+
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))

ggsave(plot=plot05,"../../out/image/UMAPCellIDggplot_sobj_processed_donor.pdf",width = 5,height = 3)

# save the same with theme cowplot
plot05_alt <- data2 %>%
  mutate(orig.ident = factor(orig.ident)) %>%
  # arrange(signature_score1) %>%
  # mutate(gene = "Ptx3") %>%
  ggplot() + 
  geom_point(aes(x=UMAP_1,y=UMAP_2,col=expertAnno.l1),alpha=0.5,size=0.1)+
  geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = expertAnno.l1)) +
  guides(colour = guide_legend(override.aes = list(size=5,alpha=1),ncol=2))+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))+
  theme_cowplot()

ggsave(plot=plot05_alt,"../../out/image/UMAPCellIDggplot2_sobj_processed_donor.pdf",width = 5,height = 3)

# save also the smaller version from Seurat
plot05_alt <- DimPlot(data.combined,group.by = "expertAnno.l1",raster = T,label = T)
ggsave(plot=plot05_alt,"../../out/image/UMAPSeurat_sobj_processed_donor.pdf",width = 5,height = 4)

# single plot
# plot051 <- data2 %>%
#   mutate(orig.ident = factor(orig.ident)) %>%
#   arrange(percent.ribo) %>%
#   # mutate(gene = "Ptx3") %>%
#   ggplot() + 
#   geom_point(aes(x=UMAP_1,y=UMAP_2,col=percent.ribo),alpha=0.5,size=0.1)+
#   geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = seurat_clusters)) +
#   scale_color_viridis_c(option = "turbo")+
#   theme_bw()+
#   theme(strip.background = element_blank(), 
#         panel.border = element_rect(colour = "black", fill = NA))
# 
# ggsave(plot = plot051,"../../out/image/UMAPRibo_sobj_processed_donor.pdf",width = 5,height = 3)
# 
# # single plot
# plot052 <- data2 %>%
#   mutate(orig.ident = factor(orig.ident)) %>%
#   arrange(percent.mt) %>%
#   # mutate(gene = "Ptx3") %>%
#   ggplot() + 
#   geom_point(aes(x=UMAP_1,y=UMAP_2,col=percent.mt),alpha=0.5,size=0.1)+
#   geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = seurat_clusters)) +
#   scale_color_viridis_c(option = "turbo")+
#   theme_bw()+
#   theme(strip.background = element_blank(), 
#         panel.border = element_rect(colour = "black", fill = NA))
# 
# ggsave(plot = plot052,"../../out/image/UMAPMito_sobj_processed_donor.pdf",width = 5,height = 3)
# 
# # single plot
# plot053 <- data2 %>%
#   mutate(orig.ident = factor(orig.ident)) %>%
#   arrange(nCount_RNA) %>%
#   # mutate(gene = "Ptx3") %>%
#   ggplot() + 
#   geom_point(aes(x=UMAP_1,y=UMAP_2,col=nCount_RNA),alpha=0.5,size=0.1)+
#   geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = seurat_clusters)) +
#   scale_color_viridis_c(option = "turbo")+
#   theme_bw()+
#   theme(strip.background = element_blank(), 
#         panel.border = element_rect(colour = "black", fill = NA))
# 
# ggsave(plot = plot053,"../../out/image/UMAPnCount_sobj_processed_donor.pdf",width = 5,height = 3)

# split the sample by treatment
plot07_2 <- DimPlot(data.combined,group.by = "expertAnno.l1",split.by = "treat_full",raster = T,label = T,ncol=3)
ggsave(plot=plot07_2,"../../out/image/UMAPSeuratCellIDSplitTreat_sobj_processed_donor.pdf",width = 13,height = 9)

# split by sample
plot06 <- data2 %>%
  mutate(orig.ident = factor(orig.ident)) %>%
  # arrange(signature_score1) %>%
  # mutate(gene = "Ptx3") %>%
  ggplot() + 
  geom_point(aes(x=UMAP_1,y=UMAP_2,col=expertAnno.l1),alpha=0.5,size=0.1)+
  geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = expertAnno.l1)) +
  facet_wrap(~ID,ncol=4)+
  guides(colour = guide_legend(override.aes = list(size=5,alpha=1),ncol=2))+
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))
ggsave(plot=plot06,"../../out/image/UMAPCellIDSplit_sobj_processed_donor.pdf",width = 13,height = 9)

# seurat version
plot06_2 <- DimPlot(data.combined,group.by = "expertAnno.l1",split.by = "orig.ident",raster = T,label = T,ncol=3)
ggsave(plot=plot06_2,"../../out/image/UMAPSeuratCellIDSplitSample_sobj_processed_donor.pdf",width = 13,height = 12)

# plot the donor deconvolution
plot_08 <- DimPlot(data.combined,split.by = "harmonized_donor2",raster = T,group.by = "treat_full")+facet_wrap(~harmonized_donor2,nrow=2)
ggsave(plot=plot_08,"../../out/image/UMAPSeuratDonorSplit_sobj_processed_donor.pdf",width = 10,height = 6)

# plot as seurat
plot08_2 <- data2 %>%
  mutate(orig.ident = factor(orig.ident)) %>%
  # scramble the cells to intermix the donors
  sample_n(replace = F,size = nrow(.)) |> 
  # arrange(signature_score1) %>%
  # mutate(gene = "Ptx3") %>%
  ggplot() + 
  geom_point(aes(x=UMAP_1,y=UMAP_2,col=treat_full),alpha=0.5,size=0.02)+
  # geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = expertAnno.l1)) +
  facet_wrap(~harmonized_donor2,ncol=3) +
  theme_cowplot()+
  guides(colour = guide_legend(override.aes = list(size=5,alpha=1),ncol=1))+
  theme(strip.background = element_blank())
ggsave(plot=plot08_2,"../../out/image/UMAPggplotDonorSplit_sobj_processed_donor.pdf",width = 10,height = 6)

# proportions -------------------------------------------------------------
# proportion of cell per cluster
df_summary <- df_meta %>%
  group_by(ID,expertAnno.l1) %>%
  summarise(n = n()) %>%
  group_by(ID) %>%
  mutate(tot=sum(n)) %>%
  ungroup() %>%
  mutate(prop=n/tot)
write_tsv(df_summary,"../../out/table/summary_sobj_processed_donor.tsv")

# make the summary per cell type
df_summary2 <- df_meta %>%
  group_by(ID,harmonized_donor2,expertAnno.l1) %>%
  summarise(n = n()) %>%
  group_by(ID,harmonized_donor2) %>%
  mutate(tot=sum(n)) %>%
  ungroup() %>%
  mutate(prop=n/tot)
write_tsv(df_summary2,"../../out/table/summary_sobj_processed_donor2.tsv")

df_summary2 |> 
  group_by(harmonized_donor2,expertAnno.l1) |> 
  summarise(n = n()) |> 
  print(n=40)

df_summary2 |> 
  group_by(harmonized_donor2,expertAnno.l1) |> 
  summarise(avg = mean(prop)) |> 
  filter(expertAnno.l1=="OLIGO")

# define the color palette
# col_pal <- colorRampPalette(brewer.pal(9, "Set1"))
# show_col()

plot082 <- df_summary2 %>%
  mutate(harmonized_donor2 = factor(harmonized_donor2)) %>%
  ggplot(aes(x=expertAnno.l1,y=prop,col=harmonized_donor2))+
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.8)) +
  theme_bw()+
  scale_y_sqrt(breaks = c(0,0.1,0.2,0.4,0.6,0.8,1))
ggsave(plot=plot082,"../../out/image/summaryDodge_sobj_processed_donor2.pdf",width = 8,height = 4)

# focus on the cell type per treat 
# make the summary per cell type
df_summary3 <- df_meta %>%
  group_by(ID,treat,treat_full,harmonized_donor2,expertAnno.l1) %>%
  summarise(n = n()) %>%
  group_by(ID,harmonized_donor2) %>%
  mutate(tot=sum(n)) %>%
  ungroup() %>%
  mutate(prop=n/tot)
write_tsv(df_summary3,"../../out/table/summary_sobj_processed_donor3.tsv")

plot083 <- df_summary3 %>%
  filter(!str_detect(harmonized_donor2,pattern="doublet|unassigned")) |> 
  mutate(harmonized_donor2 = factor(harmonized_donor2)) %>%
  ggplot(aes(x=expertAnno.l1,y=prop,col=treat_full))+
  geom_boxplot(outlier.shape = NA,position = position_dodge(width = 0.8),width=0.5)+
  geom_point(position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.8),alpha=0.5) +
  # geom_boxplot(outlier.shape = NA) +
  # geom_point(position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.8),alpha=0.5) +
  theme_bw()+
  scale_y_sqrt(breaks = c(0,0.05,0.1,0.2,0.3,0.4,0.8,1))
ggsave(plot=plot083,"../../out/image/summaryDodge_sobj_processed_donor3.pdf",width = 10,height = 4)

# Dotplot -----------------------------------------------------------------
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
ggsave(plot = test_short,"../../out/image/Dotplot_sobj_processed_donor.pdf",width = 14,height = 4)

test_short2 <- DotPlot(data.combined, features = shortlist_features_list2, dot.scale = 8,cluster.idents = T) +
  RotatedAxis()
ggsave(plot = test_short2,"../../out/image/Dotplot2_sobj_processed_donor.pdf",width = 14,height = 4)

test_long <- DotPlot(data.combined, features = shortlist_features_list_long, dot.scale = 8,cluster.idents = T) +
  RotatedAxis()
ggsave(plot=test_long,"../../out/image/DotplotLong_sobj_processed_donor.pdf",width = 25,height = 4)

df_meta %>%
  group_by(orig.ident) %>%
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
  write_tsv("../../out/table/FindAllMarkers_sobj_processed_donor.tsv")
