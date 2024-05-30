# extract the metadata to be more flexible with the plotting
# save the current meta add also the coordinates of the UMAP
df_umap <- Seurat.object@reductions$umap@cell.embeddings %>% 
  data.frame() %>% 
  rownames_to_column()

df_meta <- Seurat.object@meta.data %>% 
  data.frame() %>% 
  rownames_to_column()

df_meta_full <- left_join(df_umap,df_meta,"rowname")

df_meta_full %>%
  filter(treat %in% c("BASELINE","CSF.ctrl","CSF.MS")) %>% 
  arrange(SENEQUANTILE) %>% 
  ggplot(aes(x=UMAP_1,y=UMAP_2,col=seurat_clusters))+
  geom_point(size = 0.1,alpha=0.5)+theme_cowplot()+facet_wrap(~treat,ncol = 1)+
  # scale_color_manual(values = c("NO"="#F15A2B" , "YES"="#1D75BC"))+
  theme(strip.background = element_blank())+
  guides(colour = guide_legend(override.aes = list(size=5,alpha=1)))

# full map
df_meta_full %>%
  filter(treat %in% c("BASELINE","CSF.ctrl","CSF.MS")) %>% 
  arrange(SENEQUANTILE) %>% 
  ggplot(aes(x=UMAP_1,y=UMAP_2,col=SENEQUANTILE))+geom_point(size = 0.1,alpha=0.5)+theme_cowplot()+
  scale_color_manual(values = c("NO"="#F15A2B" , "YES"="#1D75BC"))+
  guides(colour = guide_legend(override.aes = list(size=5,alpha=1)))+
  facet_grid(SENEQUANTILE~treat) + theme(strip.background = element_blank())

# count the senescent cells per sample
df_summary_sample <- df_meta_full %>% 
  filter(treat %in% c("BASELINE","CSF.ctrl","CSF.MS")) %>% 
  group_by(ID,treat,clone,doxy,exposure,SENEQUANTILE) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  group_by(ID,treat,clone,doxy,exposure) %>% 
  mutate(tot = sum(n)) %>% 
  ungroup() %>% 
  mutate(prop_senescence = n/tot)

df_summary_sample %>% 
  ggplot(aes(x=treat,y=prop_senescence))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.2))+
  facet_wrap(~SENEQUANTILE,scales = "free")+theme_bw()+theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))

# do the same per cluster
df_summary_cluster <- df_meta_full %>% 
  filter(treat %in% c("BASELINE","CSF.ctrl","CSF.MS")) %>% 
  mutate(treat = as.factor(treat),
         ID= as.factor(ID),
         clone= as.factor(clone),
         doxy= as.factor(doxy),
         exposure= as.factor(exposure),
         seurat_clusters= as.factor(seurat_clusters),
         SENEQUANTILE= as.factor(SENEQUANTILE)) %>% 
  group_by(ID,treat,clone,doxy,exposure,seurat_clusters,SENEQUANTILE,.drop = FALSE) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  group_by(ID,treat,clone,doxy,exposure,seurat_clusters,.drop = FALSE) %>% 
  mutate(tot = sum(n)) %>% 
  ungroup() %>% 
  mutate(prop_senescence = n/tot)

df_summary_cluster %>%
  filter(SENEQUANTILE=="YES") %>% 
  ggplot(aes(x=treat,y=prop_senescence))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.2))+
  facet_wrap(~seurat_clusters,scales = "free")+theme_bw()+theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))

# try a cluster imputation
df_summary_cellID <- df_meta_full %>% 
  filter(treat %in% c("BASELINE","CSF.ctrl","CSF.MS")) %>% 
  mutate(cell_id = case_when(seurat_clusters%in%c(0,4,12)~"ASTRO",
                             seurat_clusters%in%c(3,11,14)~"IMMATURE",
                             seurat_clusters%in%c(8,15,17)~"PROG",
                             seurat_clusters%in%c(1,2,9,10)~"NEU",
                             seurat_clusters%in%c(5)~"OPC",
                             seurat_clusters%in%c(7)~"OLIGO",
                             seurat_clusters%in%c(6,13,16)~"MG")) %>% 
  mutate(treat = as.factor(treat),
         ID= as.factor(ID),
         clone= as.factor(clone),
         doxy= as.factor(doxy),
         exposure= as.factor(exposure),
         cell_id= as.factor(cell_id),
         SENEQUANTILE= as.factor(SENEQUANTILE)) %>% 
  group_by(ID,treat,clone,doxy,exposure,cell_id,SENEQUANTILE,.drop = F) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  group_by(ID,treat,clone,doxy,exposure,cell_id,.drop = F) %>% 
  mutate(tot = sum(n)) %>% 
  ungroup() %>% 
  mutate(prop_senescence = n/tot)

df_summary_cellID %>%
  filter(SENEQUANTILE=="YES") %>% 
  ggplot(aes(x=treat,y=prop_senescence))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.2))+
  facet_wrap(~cell_id,scales = "free")+theme_bw()+theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
s