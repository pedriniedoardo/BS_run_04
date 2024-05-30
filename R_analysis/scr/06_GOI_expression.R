# libraries ---------------------------------------------------------------
library(Seurat)
library(SeuratWrappers)
library(tidyverse)
library(scales)
library(ComplexHeatmap)
library(lemon)
library(finalfit)

# read in the dataset -----------------------------------------------------
data.combined <- readRDS("../../out/object/sobj_processed_donor.rds")
DimPlot(data.combined,label = T,raster = T,group.by = "seurat_clusters")

# define the gene of interest GOI
# GOI <- c("Irf7","Ddx58")
GOI <- c("DRD2","DRD3","5HT1A")

table(data.combined@meta.data$seurat_clusters)

# add a broader cell grouping
# data.combined@meta.data$cell_type2 <- data.combined@meta.data |> 
#   mutate(cell_type2 = case_when(# according to label transfer cluster 12 and 6 are Bipolar cells
#     cell_type %in% c("BP_0","BP_10","BP_11","BP_13","BP_16","BP_17","BP_5","BP_8","12","6")~"BP",
#     T~cell_type)) |> 
#   pull(cell_type2)

# generate the table for the plots ----------------------------------------
# get the metadata from the other object
meta <- data.combined@meta.data %>%
  rownames_to_column(var = "barcodes")

# extrac the expression value
df_exp <- FetchData(data.combined, vars = GOI,slot = "data") |> 
  rownames_to_column("barcodes") |> 
  pivot_longer(names_to = "gene",values_to = "exp",-barcodes) |> 
  # try to min/max normalize the count varaible per gene in order to rescale the difference in terms of expression
  group_by(gene) %>%
  # threshold of positiveness is based on the distriubtion of the expression of the signal in tihs case
  mutate(exp_min_max = ((exp - min(exp))/(max(exp)-min(exp))),
         exp_cat = case_when(exp > 0~"pos",
                             T~"neg")) %>%
  ungroup() %>%
  mutate(exp_fix = exp + rnorm(nrow(.))/100000)

# get the coordinates
UMAP1_df <- data.combined@reductions$umap@cell.embeddings %>%
  data.frame() %>%
  rownames_to_column(var = "barcodes")

# generate the dataset for mapping the data in the umamp
dim(UMAP1_df)
dim(df_exp)
dim(meta)

df_tot <- reduce(list(meta,UMAP1_df,df_exp),left_join, by="barcodes")
df_tot_avg <- df_tot %>% group_by(seurat_clusters) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)

dim(df_tot)

# plot the average expression per sample use the variable cell tyep per donor as grouping
# data.combined$group <- paste0(data.combined$orig.ident,".",data.combined$cell_type2)
data.combined$group <- paste0(data.combined$treat,"-",data.combined$seurat_clusters,"-",data.combined$harmonized_donor2)
# data.combined$group2 <- paste0(data.combined$orig.ident,".",data.combined$treat,".",data.combined$cell_type2)
Idents(data.combined) <- "group"
DefaultAssay(data.combined) <- "RNA"

average_GOI <- AverageExpression(data.combined,features = GOI,group.by = c("group"))

# plot general UMAP -------------------------------------------------------
# build the plot using both info
ggplot(label= TRUE) +
  # geom_point(data = data2,aes(x = UMAP_1,y = UMAP_2,col=cell_type),size=0.3,alpha=0.1) +
  geom_point(data = df_tot,aes(x = UMAP_1,y = UMAP_2,col=seurat_clusters),size=0.3) +
  # geom_point(data = data2_unc,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
  # geom_point(data = data2_defined,aes(x = UMAP_1,y = UMAP_2, col = robust_score),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = df_tot_avg,aes(x = UMAP_1,y = UMAP_2,label = seurat_clusters),col="black")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme_bw()
# facet_wrap(~infection)
# ggsave("../../out/image/ManualClean/UMAP_38_annotationConfident.pdf",width = 7,height = 5)

# no lab
ggplot(label= TRUE) +
  # geom_point(data = data2,aes(x = UMAP_1,y = UMAP_2,col=cell_type),size=0.3,alpha=0.1) +
  geom_point(data = df_tot,aes(x = UMAP_1,y = UMAP_2,col=seurat_clusters),size=0.3) +
  # geom_point(data = data2_unc,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
  # geom_point(data = data2_defined,aes(x = UMAP_1,y = UMAP_2, col = robust_score),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  # ggrepel::geom_text_repel(data = df_tot_avg,aes(x = UMAP_1,y = UMAP_2,label = cell_type2),col="black")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme_bw()
# facet_wrap(~infection)
# ggsave("../../out/image/ManualClean/UMAP_38_annotationConfident_noLab.pdf",width = 7,height = 5)

# split the sample
ggplot(label= TRUE) +
  # geom_point(data = data2,aes(x = UMAP_1,y = UMAP_2,col=cell_type),size=0.3,alpha=0.1) +
  geom_point(data = df_tot,aes(x = UMAP_1,y = UMAP_2,col=seurat_clusters),size=0.3) +
  # geom_point(data = data2_unc,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
  # geom_point(data = data2_defined,aes(x = UMAP_1,y = UMAP_2, col = robust_score),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = df_tot_avg,aes(x = UMAP_1,y = UMAP_2,label = seurat_clusters),col="black")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme_bw() +
  facet_wrap(~treat)+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA),axis.text.x = element_text(hjust = 1,angle = 45)
  )
# ggsave("../../out/image/ManualClean/UMAP_38_annotationConfident_split.pdf",width = 17,height = 15)
# DimPlot(data.combined,group.by = "seurat_clusters",split.by = "treat",raster = F)+facet_wrap(~treat,ncol=3)

# nolab
ggplot(label= TRUE) +
  # geom_point(data = data2,aes(x = UMAP_1,y = UMAP_2,col=cell_type),size=0.3,alpha=0.1) +
  geom_point(data = df_tot,aes(x = UMAP_1,y = UMAP_2,col=seurat_clusters),size=0.3) +
  # geom_point(data = data2_unc,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
  # geom_point(data = data2_defined,aes(x = UMAP_1,y = UMAP_2, col = robust_score),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  # ggrepel::geom_text_repel(data = df_tot_avg,aes(x = UMAP_1,y = UMAP_2,label = cell_type2),col="black")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme_bw() +
  facet_wrap(~treat)+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA),axis.text.x = element_text(hjust = 1,angle = 45)
  )
# ggsave("../../out/image/ManualClean/UMAP_38_annotationConfident_split_noLab.pdf",width = 17,height = 15)

# split by donors
ggplot(label= TRUE) +
  # geom_point(data = data2,aes(x = UMAP_1,y = UMAP_2,col=cell_type),size=0.3,alpha=0.1) +
  geom_point(data = df_tot,aes(x = UMAP_1,y = UMAP_2,col=treat),size=0.1) +
  # geom_point(data = data2_unc,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
  # geom_point(data = data2_defined,aes(x = UMAP_1,y = UMAP_2, col = robust_score),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  # ggrepel::geom_text_repel(data = df_tot_avg,aes(x = UMAP_1,y = UMAP_2,label = cell_type2),col="black")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme_bw() +
  facet_wrap(~harmonized_donor2)+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA),axis.text.x = element_text(hjust = 1,angle = 45)
  )

# expression distribution -------------------------------------------------
# crop the 0 expressing cells
df_exp %>%
  ggplot(aes(x=exp))+geom_histogram()+facet_grid(~gene)+theme_bw()+scale_x_log10()+geom_vline(xintercept = 1,col="red",linetype="dotted")

# keep the 0 expressing cells
df_exp %>%
  ggplot(aes(x=exp))+geom_histogram()+facet_wrap(~gene)+theme_bw()+
  # scale_x_log10()+
  geom_vline(xintercept = 2.5,col="red",linetype="dotted")

# library(scales)
# show_col(c("#4662D7FF","#FABA39FF","#7A0403FF"))

# plotting expression -----------------------------------------------------
# by counts
df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24")) %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  arrange(exp_cat) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2,col = exp)) + geom_point(alpha = 0.5,size = 0.2) +
  facet_wrap(gene~treat) +
  theme_bw() +
  scale_color_gradient(low = "gray",high = "blue") +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24")) %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  arrange(exp_cat) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2,col = exp)) + geom_point(alpha = 0.5,size = 0.05) +
  facet_wrap(gene~treat) +
  theme_bw() +
  # scale_color_gradient(low = "gray",high = "blue") +
  scale_color_viridis_c(option = "turbo") +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

# by min max normalized counts
df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24")) %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  arrange(exp_cat) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2,col = exp_min_max)) + geom_point(alpha = 0.5,size = 0.2) +
  facet_wrap(gene~treat) + theme_bw() + scale_color_gradient(low = "gray",high = "blue") +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

# plot the category. being 0 or non zero per cell
df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24","NMDA_36")) %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  arrange(exp_cat) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2,col = exp_cat)) + geom_point(alpha = 0.5,size = 0.05) +
  # facet_wrap(gene~NMDA_time,nrow = 2) +
  facet_rep_wrap(gene~treat,repeat.tick.labels = "all",nrow=3)+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme_bw() + scale_color_manual(values = c("gray","blue")) +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
# ggsave("../../out/image/ManualClean/UMAP_38_annotationConfident_TSPO_propPos.pdf",width = 16,height = 15)

# violin plot fro PTX3 expression use macro categories
df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24","NMDA_36")) %>%
  # this is the processing shown in the violinplot function
  # mutate(exp_fix = exp + rnorm(nrow(.))/100000) %>%
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=treat,y=exp_fix)) + 
  geom_violin(scale = "width")+
  geom_point(position=position_jitter(width = 0.2),alpha=0.1) +
  facet_wrap(~seurat_clusters) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1,angle = 90)) +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

# try to depict the average expression there is roughly one sample per condition
df_avg <- average_GOI$RNA %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(names_to = "group",values_to = "avg_exp",-gene) %>%
  filter(!str_detect(group,pattern="doublet|unassigned")) |> 
  mutate(treat = str_extract(group,pattern = c("BASELINE|CSF.ctrl|CSF.MS|cytokine|Fe|myelin|TBHP"))) |> 
  mutate(donor = str_extract(group,pattern = c("don\\w+"))) |> 
  mutate(seurat_clusters = str_remove_all(group,pattern = c("BASELINE.|CSF.ctrl.|CSF.MS.|cytokine.|Fe.|myelin.|TBHP.|.don\\w+")))

df_avg |>
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=treat,y=avg_exp))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),alpha = 0.6)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+
  facet_wrap(~seurat_clusters,scales = "free")
# scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")

# check is it is a marker in the cell type markes
# library(SeuratWrappers)
# data.combined$group2 <- paste0(data.combined$pathology_class,"-",data.combined$annotation_confident)
# DefaultAssay(data.combined) <- "RNA"

# make the stats using a pairwise test on the average expression
# in this case the origin and the sample are exclusives, therefore there is no need to generate another aggregation
df_avg |> 
  finalfit(formula = avg_exp~seurat_clusters+treat)

# cell_id <- "0"
id_test <- names(table(df_avg$seurat_clusters)) |> 
  # remove cluster 14 and 16 since there is no relicate
  str_subset(pattern = c("14|16"),negate = T)

list_pairwise <- lapply(id_test,function(cell_id){
  # check the progress
  print(cell_id)
  # subset the dataset to the single cell tyep
  test <- df_avg |> 
    filter(seurat_clusters %in% cell_id)
  
  # run the pairwise test
  test_out <- pairwise.t.test(test$avg_exp, test$treat,p.adjust.method = "fdr")
  
  # wrangle the table
  test2 <- test_out$p.value |> 
    data.frame() |> 
    rownames_to_column("first") |> 
    pivot_longer(names_to = "second",values_to = "padj",-first) |> 
    dplyr::filter(!is.na(padj)) |> 
    mutate(cell_id = cell_id)
  return(test2)
})

# build the facet for the interept
df_intercept <- list_pairwise |> 
  bind_rows() |>
  arrange(padj) |>
  group_by(first,second) |> 
  summarise() |> 
  mutate(intercept = -log(0.05))

list_pairwise |> 
  bind_rows() |>
  arrange(padj) |> 
  # mutate(comparison = paste(first,second,sep = " vs ")) |> 
  ggplot(aes(x=cell_id,y=-log(padj)))+
  facet_grid(first~second,drop = T)+
  geom_point()+
  geom_hline(data = df_intercept,aes(yintercept = intercept),col="red",linetype="dashed")+theme_bw()+theme(strip.background = element_blank())

# # use the crossing option
# crossing(first = data.combined$annotation_confident,second = data.combined$annotation_confident) |> 
#   # remove the equal ones
#   filter(first != second)

# barplot over time of the positive cells for Ptx3 per cell_type
# use the full panel of cells
# df_tot %>%
#   # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24","NMDA_36")) %>%
#   group_by(treat,seurat_clusters) %>%
#   summarise(cells = n(),
#             pos = sum(exp_cat=="pos")) %>%
#   mutate(prop_pos = pos/cells) %>%
#   mutate(log_number_cells = log10(cells)) |> 
#   # mutate(pathology_class = factor(pathology_class,levels = c("CX_Ctrl","CX_Mye","CX_Demye","WM_Ctrl","WM_NAWM","WM_CI","WM_CA","WM_Core"))) |> 
#   ggplot(aes(x=treat,y=prop_pos,fill=log_number_cells))+geom_col()+facet_wrap(~seurat_clusters)+ theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 90))+
#   theme(strip.background = element_blank(),
#         panel.border = element_rect(colour = "black", fill = NA)) +
#   scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
# ggsave("../../out/image/ManualClean/barplot_38_annotationConfident_TSPO_propPosCells.pdf",width = 9,height = 9)

# proportion on total sample
# use the full panel of cells
# df_tot %>%
#   # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24","NMDA_36")) %>%
#   group_by(pathology_class,expertAnno.l1) %>%
#   summarise(cells = n(),
#             pos = sum(exp_cat=="pos")) %>%
#   mutate(tot_cells = sum(cells)) %>%
#   mutate(prop_pos = pos/tot_cells) %>%
#   mutate(log_number_cells = log10(cells)) |> 
#   ggplot(aes(x=pathology_class,y=prop_pos,fill=log_number_cells))+geom_col()+facet_wrap(~expertAnno.l1)+ theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 90)) +
#   theme(strip.background = element_blank(),
#         panel.border = element_rect(colour = "black", fill = NA))+
#   scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
# ggsave("../../out/image/ManualClean/barplot_38_annotationConfident_TSPO_SamplepropPosCells.pdf",width = 9,height = 9)

# explore expression ------------------------------------------------------
data.combined$group_cellType <- paste0(data.combined$harmonized_donor2,"-",data.combined$seurat_clusters)
# data.combined$group2 <- paste0(data.combined$orig.ident,".",data.combined$treat,".",data.combined$cell_type2)
# Idents(data.combined) <- "group_cellType"
DefaultAssay(data.combined) <- "RNA"

average_group_cellType <- AverageExpression(data.combined,features = GOI,group.by = c("group_cellType"))

df_avg_cellType <- average_group_cellType$RNA %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  mutate(gene = "DRD2") |> 
  pivot_longer(names_to = "group",values_to = "avg_exp",-gene) |> 
  filter(!str_detect(group,pattern="doublet|unassigned")) |> 
  mutate(donor = str_extract(group,pattern = c("don\\w+"))) |> 
  mutate(seurat_clusters = str_remove_all(group,pattern = c("BASELINE.|CSF.ctrl.|CSF.MS.|cytokine.|Fe.|myelin.|TBHP.|don\\w+.")))

df_avg_cellType |>
  mutate(log1p = log1p(avg_exp)) |> 
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=seurat_clusters,y=avg_exp))+
  # ggplot(aes(x=annotation_confident,y=log1p,col=origin))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),alpha = 0.2)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 45))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
