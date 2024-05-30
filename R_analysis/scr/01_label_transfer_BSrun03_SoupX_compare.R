# AIM ---------------------------------------------------------------------
# compare the results of the label transfer on the BS run 03 reference

# libraries ---------------------------------------------------------------
library(tidyverse)
library(Seurat)
# library(SeuratData)
# library(ggridges)
# library(ComplexHeatmap)

# read in the data --------------------------------------------------------
# list_out1 <- readRDS("../../out/object/list_out_label_transfer_BSrun0102_SoupX_00200_06000_15.rds") %>%
#   bind_rows(.id = "sample") %>% 
#   mutate(filter_param = "00200_06000_15")
# list_out2 <- readRDS("../../out/object/list_out_label_transfer_BSrun0102_SoupX_00600_06000_15.rds") %>% 
#   bind_rows(.id = "sample") %>% 
#   mutate(filter_param = "00600_06000_15")
list_out3 <- readRDS("../../out/object/list_out_label_transfer_BSrun0102_SoupX_01000_06000_15.rds") %>% 
  bind_rows(.id = "sample") %>% 
  mutate(filter_param = "01000_06000_15")

# read in the reference metadata for the run 03
df_ref <- read_tsv(file = "../../out/table/ref_BSrun03_meta.tsv")


# plot the label trasfer --------------------------------------------------
# compare all the filtering options
# df_all <- list(list_out1,
#                list_out2,
#                list_out3) %>%
#   bind_rows()

df_all <- list(list_out3) %>%
  bind_rows()

# plot all the samples over the ref plot umap unfiltered
plot <- ggplot(label= TRUE) +
  # reference layer
  geom_point(data = df_ref %>% dplyr::select(UMAP_1,UMAP_2),aes(x=UMAP_1,y=UMAP_2),col="gray",alpha=0.1,size=0.1)+
  # transfer layer
  geom_point(data = df_all,aes(x = refUMAP_1,y = refUMAP_2, col = predicted.id),size=0.1,alpha=0.8) +
  # labs(color= "Clusters") +
  theme_bw() +
  facet_grid(filter_param~ID)+
  guides(colour = guide_legend(override.aes = list(size=5),ncol = 1))+theme(strip.background = element_blank())
# ggsave(plot = plot,"../../out/image/UMAP_label_transfer_BSrun03_panel.pdf",width = 40,height = 9)
ggsave(plot = plot,"../../out/plot/UMAP_label_transfer_BSrun03_panel_raw.png",width = 29,height = 4)

# plot all the samples over the ref plot umap unfiltered
plot2 <- ggplot(label= TRUE) +
  # reference layer
  geom_point(data = df_ref %>% dplyr::select(UMAP_1,UMAP_2),aes(x=UMAP_1,y=UMAP_2),col="gray",alpha=0.1,size=0.1)+
  # transfer layer
  geom_point(data = df_all,aes(x = refUMAP_1,y = refUMAP_2, col = robust_score),size=0.1,alpha=0.8) +
  # labs(color= "Clusters") +
  theme_bw() +
  facet_grid(filter_param~ID)+
  guides(colour = guide_legend(override.aes = list(size=5),ncol = 1))+theme(strip.background = element_blank())
# ggsave(plot = plot,"../../out/image/UMAP_label_transfer_BSrun03_panel.pdf",width = 40,height = 9)
ggsave(plot = plot2,"../../out/plot/UMAP_label_transfer_BSrun03_panel_robust.png",width = 29,height = 4)

# summarise global dataset ------------------------------------------------
# global distribution of the transfer score
df_all %>%
  ggplot(aes(x=predicted.celltype.score))+geom_histogram()+geom_vline(xintercept = 0.7,linetype="dashed",col="red")+theme_bw()

# summarise the global trend in the dataset
df_summary_run_04_raw <- bind_rows(df_all %>%
            group_by(predicted.id) %>%
            summarise(n = n()) %>%
            mutate(tot = sum(n)) %>%
            mutate(prop = n/tot) %>%
            mutate(dataset = "BSrun04_raw") %>%
            dplyr::rename(ID = predicted.id),
          df_ref %>%
            group_by(expertAnno.l1) %>%
            summarise(n = n()) %>%
            mutate(tot = sum(n)) %>%
            mutate(prop = n/tot) %>%
            mutate(dataset = "BSrun03") %>%
            dplyr::rename(ID = expertAnno.l1))

df_summary_run_04_robust <- bind_rows(df_all %>%
                                     group_by(robust_score) %>%
                                     summarise(n = n()) %>%
                                     mutate(tot = sum(n)) %>%
                                     mutate(prop = n/tot) %>%
                                     mutate(dataset = "BSrun04_robust") %>%
                                     dplyr::rename(ID = robust_score),
                                   df_ref %>%
                                     group_by(expertAnno.l1) %>%
                                     summarise(n = n()) %>%
                                     mutate(tot = sum(n)) %>%
                                     mutate(prop = n/tot) %>%
                                     mutate(dataset = "BSrun03") %>%
                                     dplyr::rename(ID = expertAnno.l1))

df_summary_run_04_raw %>%
  mutate(ID = fct_reorder(ID, prop,.desc = T)) %>%
  ggplot(aes(x=ID,y=prop,fill=dataset))+
  geom_col(position = "dodge") +
  theme_bw() +
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))+scale_y_continuous(breaks = seq(from=0,to=0.6,by=0.05))

df_summary_run_04_robust %>%
  mutate(ID = fct_reorder(ID, prop,.desc = T)) %>%
  ggplot(aes(x=ID,y=prop,fill=dataset))+
  geom_col(position = "dodge") +
  theme_bw() +
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))+scale_y_continuous(breaks = seq(from=0,to=0.6,by=0.05))

# summarise split dataset -------------------------------------------------
# plot the ralative proportions across filtering stages
df_max_raw <- df_all %>%
  group_by(ID,filter_param) %>%
  summarise(tot = n()) %>%
  ungroup() %>%
  # filter(filter_param == "00200_06000_15") %>%
  dplyr::select(ID,tot)

df_summary_raw <- df_all %>%
  group_by(ID,filter_param,predicted.id) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  left_join(df_max_raw,by = "ID") %>%
  mutate(prop = n/tot) %>%
  group_by(predicted.id) %>%
  mutate(avg_prop = mean(prop)) %>%
  ungroup()

# check the proportions
df_summary_raw %>%
  mutate(predicted.id = fct_reorder(predicted.id, avg_prop,.desc = T)) %>%
  ggplot(aes(x=predicted.id,y=prop,fill=filter_param))+
  geom_col(position = "dodge")+facet_wrap(~ID,scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("../../out/plot/barplot_label_transfer_BSrun03_panel_raw.pdf",width = 15,height = 10)

# use a different scale
df_summary_raw %>%
  mutate(predicted.id = fct_reorder(predicted.id, avg_prop,.desc = T)) %>%
  ggplot(aes(x=predicted.id,y=prop,fill=filter_param))+
  geom_col(position = "dodge")+facet_wrap(~ID,scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))+
  scale_y_sqrt(breaks = c(c(0,0.01,0.05,0.1,0.2,0.4,0.6)))

# the filter for high confidence is score > 0.7
df_max_robust <- df_all %>%
  group_by(ID,filter_param) %>%
  summarise(tot = n()) %>%
  ungroup() %>%
  # filter(filter_param == "00200_06000_15") %>%
  dplyr::select(ID,tot)

df_summary_robust <- df_all %>%
  group_by(ID,filter_param,robust_score) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  left_join(df_max_raw,by = "ID") %>%
  mutate(prop = n/tot) %>%
  group_by(robust_score) %>%
  mutate(avg_prop = mean(prop)) %>%
  ungroup()

# check the proportions
df_summary_robust %>%
  mutate(robust_score = fct_reorder(robust_score, avg_prop,.desc = T)) %>%
  ggplot(aes(x=robust_score,y=prop,fill=filter_param))+
  geom_col(position = "dodge")+facet_wrap(~ID,scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("../../out/plot/barplot_label_transfer_BSrun03_panel_robust.pdf",width = 15,height = 10)


# use a different scale
df_summary_robust %>%
  mutate(robust_score = fct_reorder(robust_score, avg_prop,.desc = T)) %>%
  ggplot(aes(x=robust_score,y=prop,fill=filter_param))+
  geom_col(position = "dodge")+facet_wrap(~ID,scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))+
  scale_y_sqrt(breaks = c(c(0,0.01,0.05,0.1,0.2,0.4,0.6)))

# check if it is losing or gaining clusters in the comparison
# df_summary2 <- df_summary %>% 
#   # mutate(f = paste0(ID,"_",predicted.id)) %>% 
#   split(f = paste0(.$ID,"_",.$predicted.id)) %>% 
#   lapply(function(x){
#     x %>%
#       arrange(filter_param) %>% 
#       mutate(FC_prop = prop/prop[1]) %>% 
#       mutate(logFC_prop = log(FC_prop))
#   }) %>% 
#   bind_rows() %>% 
#   ungroup() %>% 
#   group_by(predicted.id) %>% 
#   mutate(avg_logFC_prop = mean(logFC_prop)) %>% 
#   ungroup()

# plot the general FC per cluster
# df_summary2 %>%
#   mutate(predicted.id = fct_reorder(predicted.id, avg_logFC_prop,.desc = T)) %>%
#   ggplot(aes(x=predicted.id,y=logFC_prop,fill=filter_param))+
#   geom_col(position = "dodge")+facet_wrap(~ID,scales = "free") +
#   theme_bw() +
#   theme(strip.background = element_blank())
# ggsave("../../out/image/barplot_label_transfer_BS_panel2.pdf",width = 15,height = 10)

# plot a summary trend
# df_summary2 %>%
#   mutate(predicted.id = fct_reorder(predicted.id, avg_logFC_prop,.desc = T)) %>%
#   ggplot(aes(x=predicted.id,y=logFC_prop,col=filter_param))+
#   geom_boxplot(outlier.shape = NA)+
#   geom_point(position = position_dodge(width = 0.7))+
#   theme_bw() +
#   theme(strip.background = element_blank())
# ggsave("../../out/image/barplot_label_transfer_BS_panel3.pdf",width = 6,height = 4)
