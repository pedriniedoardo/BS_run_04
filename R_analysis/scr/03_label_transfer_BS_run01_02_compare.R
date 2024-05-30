# # libraries ---------------------------------------------------------------
# library(tidyverse)
# library(Seurat)
# # library(SeuratData)
# # library(ggridges)
# # library(ComplexHeatmap)
# 
# # read in the data --------------------------------------------------------
# list_out <- readRDS("../../out/object/list_out_label_transfer_BS0102.rds")
# 
# # read in the reference metadata for the run 01 and run 02
# df_ref <- read_tsv(file = "../../out/table/ref_BS_meta.tsv")
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
#   geom_point(data = df_all,aes(x = refUMAP_1,y = refUMAP_2, col = predicted.id),size=0.3,alpha=0.8) +
#   # labs(color= "Clusters") +
#   theme_bw() +
#   facet_grid(feature_filter~ID)+
#   guides(colour = guide_legend(override.aes = list(size=5),ncol = 1))+theme(strip.background = element_blank())
# ggsave(plot = plot,"../../out/image/UMAP_label_transfer_BS_panel.pdf",width = 40,height = 9)
# 
# # plot the ralative proportions across filtering stages
# df_max <- df_all %>% 
#   group_by(ID,feature_filter) %>% 
#   summarise(tot = n()) %>% 
#   ungroup() %>% 
#   filter(feature_filter == "02_60") %>% 
#   dplyr::select(ID,tot)
#   
# df_summary <- df_all %>% 
#   group_by(ID,feature_filter,predicted.id) %>% 
#   summarise(n = n()) %>% 
#   ungroup() %>% 
#   left_join(df_max,by = "ID") %>% 
#   mutate(prop = n/tot) %>% 
#   group_by(predicted.id) %>% 
#   mutate(avg_prop = mean(prop)) %>% 
#   ungroup()
# 
# # check the proportions
# df_summary %>% 
#   mutate(predicted.id = fct_reorder(predicted.id, avg_prop,.desc = T)) %>% 
#   ggplot(aes(x=predicted.id,y=prop,fill=feature_filter))+geom_col(position = "dodge")+facet_wrap(~ID,scales = "free")+theme_bw()+theme(strip.background = element_blank())
# ggsave("../../out/image/barplot_label_transfer_BS_panel.pdf",width = 15,height = 10)
