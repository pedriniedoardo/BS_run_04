# # libraries ---------------------------------------------------------------
# library(tidyverse)
# library(Seurat)
# # library(SeuratData)
# # library(ggridges)
# # library(ComplexHeatmap)
# 
# # read in the data --------------------------------------------------------
# list_out <- readRDS("../../out/object/list_out_label_transfer_02_60_IMM_WM.rds")
# 
# # read in the reference metadata for the run 01 and run 02
# df_ref <- read_tsv(file = "../../out/table/ref_WM_meta.tsv")
# 
# # wrangling ---------------------------------------------------------------
# # compare all the filtering options
# df_all <- list_out %>% 
#   bind_rows(.id = "sample") %>% 
#   mutate(feature_filter = str_extract(sample,pattern = "02_60|06_60|10_60"))
# 
# # plot all the samples over the ref plot umap
# plot <- ggplot(label= TRUE) +
#   # reference layer
#   geom_point(data = df_ref %>% dplyr::select(UMAP_1,UMAP_2),aes(x=UMAP_1,y=UMAP_2),col="gray",alpha=0.3,size=0.1)+
#   # transfer layer
#   geom_point(data = df_all,aes(x = refUMAP_1.WM,y = refUMAP_2.WM, col = predicted.id),size=0.3,alpha=0.8) +
#   # labs(color= "Clusters") +
#   theme_bw() +
#   facet_wrap(feature_filter~ID)+
#   guides(colour = guide_legend(override.aes = list(size=5),ncol = 1))+theme(strip.background = element_blank())
# ggsave(plot = plot,"../../out/image/UMAP_label_transfer_WM_panel_woConfidence.pdf",width = 12,height = 9)
# 
# # plot the confident annotations
# plot2 <- ggplot(label= TRUE) +
#   # reference layer
#   geom_point(data = df_ref %>% dplyr::select(UMAP_1,UMAP_2),aes(x=UMAP_1,y=UMAP_2),col="gray",alpha=0.3,size=0.1)+
#   # transfer layer
#   geom_point(data = df_all,aes(x = refUMAP_1.WM,y = refUMAP_2.WM, col = robust_score),size=0.3,alpha=0.8) +
#   # labs(color= "Clusters") +
#   theme_bw() +
#   facet_wrap(feature_filter~ID)+
#   guides(colour = guide_legend(override.aes = list(size=5),ncol = 1))+theme(strip.background = element_blank())
# ggsave(plot = plot2,"../../out/image/UMAP_label_transfer_WM_panel_wConfidence.pdf",width = 12,height = 9)
# 
# # plot the counts with or without confidence
# # df_max <- df_all %>% 
# #   group_by(ID,feature_filter) %>% 
# #   summarise(tot = n()) %>% 
# #   ungroup() %>% 
# #   filter(feature_filter == "02_60") %>% 
# #   dplyr::select(ID,tot)
# 
# df_summary1 <- df_all %>% 
#   group_by(ID,feature_filter,predicted.id) %>% 
#   summarise(n = n()) %>% 
#   ungroup() %>% 
#   group_by(predicted.id) %>% 
#   mutate(avg_count = mean(n)) %>% 
#   ungroup()
# 
# df_summary2 <- df_all %>% 
#   group_by(ID,feature_filter,robust_score) %>% 
#   summarise(n = n()) %>% 
#   ungroup() %>% 
#   group_by(robust_score) %>% 
#   mutate(avg_count = mean(n)) %>% 
#   ungroup()
# 
# # check the proportions
# df_summary1 %>% 
#   mutate(predicted.id = fct_reorder(predicted.id, avg_count,.desc = T)) %>% 
#   ggplot(aes(x=predicted.id,y=n,fill=feature_filter))+geom_col(position = "dodge")+facet_wrap(~ID,scales = "free")+theme_bw()+theme(strip.background = element_blank())
# ggsave("../../out/image/barplot_label_transfer_WM_panel_woConfidence.pdf",width = 15,height = 10)
# 
# df_summary2 %>% 
#   mutate(robust_score = fct_reorder(robust_score, avg_count,.desc = T)) %>% 
#   ggplot(aes(x=robust_score,y=n,fill=feature_filter))+geom_col(position = "dodge")+facet_wrap(~ID,scales = "free")+theme_bw()+theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
# ggsave("../../out/image/barplot_label_transfer_WM_panel_wConfidence.pdf",width = 15,height = 10)
