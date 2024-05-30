# AIM ---------------------------------------------------------------------
# the aim of the script is to count the senescent cells according to the SIT algorithm

# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(scales)
library(ggrepel)
library(cowplot)

# read the data -----------------------------------------------------------
# read in the dataset
data.combined <- readRDS("../../out/object/sobj_processed_donor.rds")
DimPlot(data.combined,label = T,raster = T,group.by = "seurat_clusters")
DimPlot(data.combined,label = T,raster = T,group.by = "expertAnno.l1")
DimPlot(data.combined,label = T,raster = T,group.by = "treat_full")

# pull the annotation derived from the SIT algorithm (generated from 00_test_SIT.R)
df_SIT <- read_tsv("../../out/table/meta_sobj_processed_donor_SIT.tsv")

# merge the data
# add also the coordinates for the UMAP
meta_full <- data.combined@meta.data %>%
  rownames_to_column("rowname") %>%
  left_join(df_SIT %>%
              dplyr::select("rowname",
                            "KEGG_CELL_CYCLE","Cell_cycle_arrest_signatureK","REACTOME_CELL_CYCLE","Cell_cycle_arrest_signatureR","WP_CELL_CYCLE","Cell_cycle_arrest_signatureWP",
                            "Senescence_scoreK",
                            "Senescence_scoreWP",
                            "Senescence_scoreR",
                            "Senescence_scoreKQUANTILE",
                            "Senescence_scoreRQUANTILE",
                            "Senescence_scoreWPQUANTILE",
                            "SENEQUANTILE"),by = "rowname") %>%
  left_join(data.combined@reductions$umap@cell.embeddings %>%
              data.frame() %>%
              rownames_to_column("rowname"),by = "rowname")

# summarise the full metadata
dim(meta_full)
dim(data.combined)

meta_full %>%
  group_by(SENEQUANTILE) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot)

# per treatment
meta_full %>%
  group_by(treat_full,SENEQUANTILE) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot)

# use the sema criterion used with the senmayo approach
# test_summary_cellID <- df_modules_cellID %>% 
#   filter(signature %in% x) %>% 
#   left_join(df_90_cellID2,c("expertAnno.l1","signature")) %>% 
#   mutate(sen = case_when(signature_score1>thr~1,
#                          T~0)) %>% 
#   group_by(signature,orig.ident,origin,disease,expertAnno.l1,sen) %>% 
#   summarise(n = n()) %>% 
#   ungroup() %>% 
#   group_by(signature,orig.ident,origin,disease,expertAnno.l1) %>% 
#   mutate(tot = sum(n)) %>% 
#   ungroup() %>% 
#   mutate(prop = n/tot)

# fill with 0 for the combination without senescent cells
enumeration_cellID_ref <- crossing(treat_full = unique(meta_full$treat_full),
         harmonized_donor2 = unique(meta_full$harmonized_donor2),
         seurat_clusters = unique(meta_full$seurat_clusters),
         SENEQUANTILE = unique(meta_full$SENEQUANTILE))

# run the same but for imputed cell type
enumeration_cellID_ref2 <- crossing(treat_full = unique(meta_full$treat_full),
                                    harmonized_donor2 = unique(meta_full$harmonized_donor2),
                                    expertAnno.l1 = unique(meta_full$expertAnno.l1),
                                    SENEQUANTILE = unique(meta_full$SENEQUANTILE))

enumeration_cellID <- meta_full %>%
  group_by(treat_full,harmonized_donor2,seurat_clusters,SENEQUANTILE) %>%
  summarise(n = n())

enumeration_cellID2 <- meta_full %>%
  group_by(treat_full,harmonized_donor2,expertAnno.l1,SENEQUANTILE) %>%
  summarise(n = n())

enumeration_cellID_full <- left_join(enumeration_cellID_ref,enumeration_cellID,by = c("treat_full","harmonized_donor2","seurat_clusters","SENEQUANTILE")) %>%
  # fill the NA with 0
  mutate(n = case_when(is.na(n)~0,
                       T~n)) |> 
  group_by(treat_full,harmonized_donor2,seurat_clusters) %>%
  mutate(tot = sum(n)) %>%
  ungroup() %>%
  mutate(prop = n/tot) %>%
  # mutate(pathology_class = factor(pathology_class,levels = c("CX_Ctrl","CX_Demye","CX_Mye","WM_Ctrl","WM_NAWM","WM_CI","WM_CA","WM_Core"))) %>%
  mutate(signature ="SIT")

enumeration_cellID_full2 <- left_join(enumeration_cellID_ref2,
                                      enumeration_cellID2,by = c("treat_full","harmonized_donor2","expertAnno.l1","SENEQUANTILE")) %>%
  # fill the NA with 0
  mutate(n = case_when(is.na(n)~0,
                       T~n)) |> 
  group_by(treat_full,harmonized_donor2,expertAnno.l1) %>%
  mutate(tot = sum(n)) %>%
  ungroup() %>%
  mutate(prop = n/tot) %>%
  # mutate(pathology_class = factor(pathology_class,levels = c("CX_Ctrl","CX_Demye","CX_Mye","WM_Ctrl","WM_NAWM","WM_CI","WM_CA","WM_Core"))) %>%
  mutate(signature ="SIT")

enumeration_cellID_full %>% 
  write_tsv("../../out/table/enumeration_SIT_donor_full_seuratClusters.tsv")

enumeration_cellID_full2 %>% 
  write_tsv("../../out/table/enumeration_SIT_donor_full_expertAnno.tsv")

# plot the data to show the proportions of senescent cells per condition
# remove the uassigned
enumeration_cellID_full %>%
  # filter(expertAnno.l1 %in% c(3,5,11)) %>%
  filter(SENEQUANTILE == "YES",
         !harmonized_donor2 %in% c("unassigned")) %>% 
  # filter(signature %in% c("CLASSICAL_SASP","Induces","Inhibits","SAEPHIA_CURATED_SASP","senmayo")) %>% 
  # ggplot(aes(x=pathology_class,y=prop,col=origin,label=orig.ident)) +
  ggplot(aes(x=treat_full,y=prop)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1)) +
  # geom_text_repel()+
  # facet_grid(signature~expertAnno.l1,scales = "free")+
  facet_wrap(signature~seurat_clusters,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text.y.right = element_text(angle = 0))+
  # geom_hline(yintercept = 0.1,linetype="dotted",col="black")+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))
ggsave("../../out/plot/enumeration_SIT_donor_full_seuratClusters.pdf",width = 15,height =15)

enumeration_cellID_full2 %>%
  # filter(expertAnno.l1 %in% c(3,5,11)) %>%
  filter(SENEQUANTILE == "YES",
         !harmonized_donor2 %in% c("unassigned")) %>% 
  # filter(signature %in% c("CLASSICAL_SASP","Induces","Inhibits","SAEPHIA_CURATED_SASP","senmayo")) %>% 
  # ggplot(aes(x=pathology_class,y=prop,col=origin,label=orig.ident)) +
  ggplot(aes(x=treat_full,y=prop)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1)) +
  # geom_text_repel()+
  # facet_grid(signature~expertAnno.l1,scales = "free")+
  facet_wrap(signature~expertAnno.l1,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text.y.right = element_text(angle = 0))+
  # geom_hline(yintercept = 0.1,linetype="dotted",col="black")+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))
ggsave("../../out/plot/enumeration_SIT_donor_full_expertAnno.pdf",width = 10,height =6)

# color the donors
enumeration_cellID_full %>%
  # filter(expertAnno.l1 %in% c(3,5,11)) %>%
  filter(SENEQUANTILE == "YES",
         !harmonized_donor2 %in% c("unassigned")) %>% 
  # filter(signature %in% c("CLASSICAL_SASP","Induces","Inhibits","SAEPHIA_CURATED_SASP","senmayo")) %>% 
  # ggplot(aes(x=pathology_class,y=prop,col=origin,label=orig.ident)) +
  ggplot(aes(x=treat_full,y=prop)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),aes(col=harmonized_donor2)) +
  # geom_text_repel()+
  # facet_grid(signature~expertAnno.l1,scales = "free")+
  facet_wrap(signature~seurat_clusters,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text.y.right = element_text(angle = 0))+
  # geom_hline(yintercept = 0.1,linetype="dotted",col="black")+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))
ggsave("../../out/plot/enumeration_SIT_donor_full2_seuratClusters.pdf",width = 15,height =15)

enumeration_cellID_full2 %>%
  # filter(expertAnno.l1 %in% c(3,5,11)) %>%
  filter(SENEQUANTILE == "YES",
         !harmonized_donor2 %in% c("unassigned")) %>% 
  # filter(signature %in% c("CLASSICAL_SASP","Induces","Inhibits","SAEPHIA_CURATED_SASP","senmayo")) %>% 
  # ggplot(aes(x=pathology_class,y=prop,col=origin,label=orig.ident)) +
  ggplot(aes(x=treat_full,y=prop)) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),aes(col=harmonized_donor2)) +
  # geom_text_repel()+
  # facet_grid(signature~expertAnno.l1,scales = "free")+
  facet_wrap(signature~expertAnno.l1,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text.y.right = element_text(angle = 0))+
  # geom_hline(yintercept = 0.1,linetype="dotted",col="black")+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))
ggsave("../../out/plot/enumeration_SIT_donor_full2_expertAnno.pdf",width = 10,height =6)

# try plot UMAPs for SIT --------------------------------------------------
meta_full %>%
  arrange(SENEQUANTILE) %>%
  # mutate(pathology_class = factor(pathology_class,levels = c("CX_Ctrl","CX_Demye","CX_Mye","WM_Ctrl","WM_NAWM","WM_CI","WM_CA","WM_Core"))) %>%
  ggplot(aes(x=UMAP_1,y=UMAP_2,col=SENEQUANTILE))+
  geom_point(alpha=0.1,size=0.1)+scale_color_manual(values = c("gray","red"))+
  facet_grid(SENEQUANTILE~treat_full)+
  theme_void() +
  guides(color = guide_legend(override.aes = list(size=1,alpha=1)))+
  theme(strip.background = element_blank())
ggsave(paste0("../../out/plot/UMAP_SIT_donor_full_seuratCluster.pdf"),width = 16,height = 4)

# balance the number of cells per sample
# what is the minimum nuber of cells per pathology
meta_full %>% group_by(treat_full) %>% summarise(n = n())

set.seed(2144)
meta_full %>%
  group_by(treat_full) %>%
  sample_n(9000) %>%
  arrange(SENEQUANTILE) %>%
  # mutate(pathology_class = factor(pathology_class,levels = c("CX_Ctrl","CX_Demye","CX_Mye","WM_Ctrl","WM_NAWM","WM_CI","WM_CA","WM_Core"))) %>%
  ggplot(aes(x=UMAP_1,y=UMAP_2,col=SENEQUANTILE))+
  geom_point(alpha=0.2,size=0.1)+scale_color_manual(values = c("gray","red"))+
  facet_grid(SENEQUANTILE~treat_full)+
  theme_void() +
  guides(color = guide_legend(override.aes = list(size=1,alpha=1)))+
  theme(strip.background = element_blank())
ggsave("../../out/plot/UMAP_SIT_donor_full_capped_seuratCluster.pdf",width = 16,height = 4)


















# # check the correlation between SIT and senmayo ---------------------------
# # join the two datasets pull only the senescence cells from both
# df_senmayo <- read_tsv("../../out/table/ManualClean/modules_SENESCENCE/HVG4000_expertAnno/enumeration_ALL_090_threshold_MSStatus_refCX_cellID_sampleWise.tsv") %>%
#   filter(signature == "senmayo") %>%
#   filter(sen == 1)
# 
# 
# df_SIT <- read_tsv("../../out/table/ManualClean/enumeration_SIT_WM_CX_harmonySkipIntegAllSoupX_expertAnno.tsv") %>%
#   filter(SENEQUANTILE == "YES")
# 
# # check the reason for the discrepancies in the two dataset. the call for senescence is not the same in both, therefore it might be absent in one estimate compard to the other
# full_join(df_senmayo %>%
#             group_by(orig.ident,origin,pathology_class) %>%
#             summarise(n = n()),
#           df_SIT %>%
#             group_by(orig.ident,origin,pathology_class) %>%
#             summarise(n = n()),by = c("orig.ident","origin","pathology_class"),suffix = c(".senmayo",".SIT")) %>%
#   mutate(delta = n.senmayo - n.SIT) %>%
#   filter(delta != 0)
# 
# # s20        wm     WM_NAWM                15    18    -3
# df_senmayo %>%
#   filter(orig.ident == "s20")
# 
# df_SIT %>%
#   filter(orig.ident == "s20")
# 
# # s27        wm     WM_Ctrl                17    16     1
# df_senmayo %>%
#   filter(orig.ident == "s27")
# 
# # confirm that the total number of cell is the same for both 
# full_join(read_tsv("../../out/table/ManualClean/modules_SENESCENCE/HVG4000_expertAnno/enumeration_ALL_090_threshold_MSStatus_refCX_cellID_sampleWise.tsv") %>%
#             filter(signature == "senmayo") %>%
#             group_by(orig.ident,expertAnno.l1,pathology_class) %>%
#             summarise(tot = sum(n)),
#           read_tsv("../../out/table/ManualClean/enumeration_SIT_WM_CX_harmonySkipIntegAllSoupX_expertAnno.tsv") %>%
#             group_by(orig.ident,expertAnno.l1,pathology_class) %>%
#             summarise(tot = sum(n)),by = c("orig.ident","expertAnno.l1","pathology_class"),suffix = c(".senmayo",".SIT")) %>%
#   mutate(delta = tot.senmayo - tot.SIT) %>%
#   filter(delta != 0)
# 
# 
# #
# df_plot <- full_join(df_senmayo %>%
#                        select(orig.ident,origin,pathology_class,expertAnno.l1,n,tot,prop),
#                      df_SIT %>%
#                        select(orig.ident,origin,pathology_class,expertAnno.l1,n,tot,prop),by = c("orig.ident","origin","pathology_class","expertAnno.l1"),suffix = c(".senmayo",".SIT"))
# 
# # it is safe to assume that if one is missing it means it is a 0%
# df_plot2 <- df_plot %>%
#   mutate(prop.SIT = case_when(is.na(prop.SIT)~0,
#                               T~prop.SIT),
#          prop.senmayo = case_when(is.na(prop.senmayo)~0,
#                                   T~prop.senmayo))
# 
# 
# # plot the correlation
# # notice that in this case I am discarding the samples for which there is a estimate of 0 senescence
# df_plot %>%
#   ggplot(aes(x=prop.senmayo,
#              y=prop.SIT))+
#   geom_smooth(method = "lm") +
#   geom_point(alpha=0.5)+theme_bw()+
#   facet_wrap(~expertAnno.l1,scales = "free")+theme(strip.background = element_blank())
# ggsave(paste0("../../out/image/ManualClean/modules_SENESCENCE/correlation_SITSenmayo_WM_CX_harmonySkipIntegAllSoupX_expertAnno_split.pdf"),width = 9,height =8)
# 
# # notice that in this case I am discarding the samples for which there is a estimate of 0 senescence
# df_plot2 %>%
#   ggplot(aes(x=prop.senmayo,
#              y=prop.SIT))+
#   geom_smooth(method = "lm") +
#   geom_point(alpha=0.5)+theme_bw()+
#   facet_wrap(~expertAnno.l1,scales = "free")+theme(strip.background = element_blank())
# ggsave(paste0("../../out/image/ManualClean/modules_SENESCENCE/correlation_SITSenmayo_WM_CX_harmonySkipIntegAllSoupX_expertAnno_split_keepZero.pdf"),width = 9,height =8)
# 
# # plot a globa correlation
# df_plot %>%
#   filter(!is.na(prop.SIT),
#          !is.na(prop.senmayo)) %>%
#   ggplot(aes(x=prop.senmayo,
#              y=prop.SIT))+
#   geom_smooth(method = "lm") +
#   geom_point(alpha=0.5)+theme_bw()+
#   theme(strip.background = element_blank())
# ggsave(paste0("../../out/image/ManualClean/modules_SENESCENCE/correlation_SITSenmayo_WM_CX_harmonySkipIntegAllSoupX_expertAnno_global.pdf"),width = 5,height =5)
# 
# df_plot %>%
#   filter(!is.na(prop.SIT),
#          !is.na(prop.senmayo)) %>%
#   ggplot(aes(x=prop.senmayo,
#              y=prop.SIT))+
#   geom_smooth(method = "lm") +
#   geom_point(alpha=0.5)+theme_bw()+
#   theme(strip.background = element_blank())+
#   scale_x_log10()+
#   scale_y_log10()
# 
# df_plot %>%
#   filter(!is.na(prop.SIT),
#          !is.na(prop.senmayo)) %>%
#   ggplot(aes(x=prop.senmayo,
#              y=prop.SIT))+
#   geom_smooth(method = "lm") +
#   geom_point(alpha=0.5)+theme_bw()+
#   theme(strip.background = element_blank())+
#   scale_x_sqrt()+
#   scale_y_sqrt()
# 
# df_plot
# cor.test(df_plot$prop.senmayo,df_plot$prop.SIT)
# 
# df_plot2 %>%
#   filter(!is.na(prop.SIT),
#          !is.na(prop.senmayo)) %>%
#   ggplot(aes(x=prop.senmayo,
#              y=prop.SIT))+
#   geom_smooth(method = "lm") +
#   geom_point(alpha=0.5)+theme_bw()+
#   theme(strip.background = element_blank())
# ggsave(paste0("../../out/image/ManualClean/modules_SENESCENCE/correlation_SITSenmayo_WM_CX_harmonySkipIntegAllSoupX_expertAnno_global_keepZero.pdf"),width = 5,height =5)
# 
# df_plot2
# cor.test(df_plot2$prop.senmayo,df_plot2$prop.SIT)
# 
# df_plot %>%
#   ggplot(aes(x=prop.senmayo,
#              y=prop.SIT))+
#   geom_smooth(method = "lm") +
#   geom_point()+theme_bw()+
#   facet_wrap(~pathology_class,scales = "free")+theme(strip.background = element_blank())
