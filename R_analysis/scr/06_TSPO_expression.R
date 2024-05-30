# libraries ---------------------------------------------------------------
library(Seurat)
library(SeuratWrappers)
library(tidyverse)
library(scales)
library(ComplexHeatmap)
library(lemon)
library(finalfit)
library(cowplot)

# read in the dataset -----------------------------------------------------
data.combined <- readRDS("../../out/object/sobj_processed_donor.rds")
DimPlot(data.combined,label = T,raster = T,group.by = "seurat_clusters")
DimPlot(data.combined,label = T,raster = T,group.by = "expertAnno.l1")
DimPlot(data.combined,label = T,raster = T,group.by = "treat_full")

# define the gene of interest GOI
# GOI <- c("Irf7","Ddx58")
GOI <- c("TSPO")

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
df_tot_avg <- df_tot %>% group_by(expertAnno.l1) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)

dim(df_tot)


head(data.combined@meta.data)
table(data.combined$ID)

# plot the average expression per sample use the variable cell tyep per donor as grouping
# data.combined$group <- paste0(data.combined$orig.ident,".",data.combined$cell_type2)
data.combined$group <- paste0(data.combined$treat_full,"-",data.combined$expertAnno.l1,"-",data.combined$harmonized_donor2)
# data.combined$group2 <- paste0(data.combined$orig.ident,".",data.combined$treat,".",data.combined$cell_type2)
Idents(data.combined) <- "group"
DefaultAssay(data.combined) <- "RNA"

average_GOI <- AverageExpression(data.combined,features = GOI,group.by = c("group"))

# plot general UMAP -------------------------------------------------------
# build the plot using both info
ggplot(label= TRUE) +
  # geom_point(data = data2,aes(x = UMAP_1,y = UMAP_2,col=cell_type),size=0.3,alpha=0.1) +
  geom_point(data = df_tot,aes(x = UMAP_1,y = UMAP_2,col=expertAnno.l1),size=0.3) +
  # geom_point(data = data2_unc,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
  # geom_point(data = data2_defined,aes(x = UMAP_1,y = UMAP_2, col = robust_score),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = df_tot_avg,aes(x = UMAP_1,y = UMAP_2,label = expertAnno.l1),col="black")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme_bw()
# facet_wrap(~infection)
# ggsave("../../out/image/ManualClean/UMAP_38_annotationConfident.pdf",width = 7,height = 5)

# no lab
ggplot(label= TRUE) +
  # geom_point(data = data2,aes(x = UMAP_1,y = UMAP_2,col=cell_type),size=0.3,alpha=0.1) +
  geom_point(data = df_tot,aes(x = UMAP_1,y = UMAP_2,col=expertAnno.l1),size=0.3) +
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
  geom_point(data = df_tot,aes(x = UMAP_1,y = UMAP_2,col=expertAnno.l1),size=0.3) +
  # geom_point(data = data2_unc,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
  # geom_point(data = data2_defined,aes(x = UMAP_1,y = UMAP_2, col = robust_score),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = df_tot_avg,aes(x = UMAP_1,y = UMAP_2,label = expertAnno.l1),col="black")+
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
  geom_point(data = df_tot,aes(x = UMAP_1,y = UMAP_2,col=expertAnno.l1),size=0.3) +
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
  facet_wrap(gene~treat_full) +
  theme_cowplot() +
  scale_color_gradient(low = "gray",high = "blue") +
  # theme(strip.background = element_blank(),
  #       panel.border = element_rect(colour = "black", fill = NA))
  theme(strip.background = element_blank())
ggsave("../../out/image/UMAPggplot_annotationConfident_TSPO_count.pdf",width = 13,height = 12)

# do the same using Seurat
FeaturePlot(data.combined,features = GOI,split.by = "treat_full",raster = T,order = T)
ggsave("../../out/image/UMAPSeurat_annotationConfident_TSPO_count.pdf",width = 25,height = 3)

df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24")) %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  arrange(exp_cat) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2,col = exp)) + geom_point(alpha = 0.5,size = 0.05) +
  facet_wrap(gene~treat_full) +
  theme_cowplot() +
  # scale_color_gradient(low = "gray",high = "blue") +
  scale_color_viridis_c(option = "turbo") +
  # theme(strip.background = element_blank(),
  #       panel.border = element_rect(colour = "black", fill = NA))
  theme(strip.background = element_blank())
ggsave("../../out/image/UMAPggplot_annotationConfident_TSPO_count2.pdf",width = 13,height = 12)

# by min max normalized counts
df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24")) %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  arrange(exp_cat) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2,col = exp_min_max)) + geom_point(alpha = 0.5,size = 0.2) +
  facet_wrap(gene~treat_full) +
  theme_cowplot() +
  scale_color_gradient(low = "gray",high = "blue") +
  # theme(strip.background = element_blank(),
  #       panel.border = element_rect(colour = "black", fill = NA))
  theme(strip.background = element_blank())
ggsave("../../out/image/UMAPggplot_annotationConfident_TSPO_minmax.pdf",width = 13,height = 12)

# plot the category. being 0 or non zero per cell
df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24","NMDA_36")) %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  arrange(exp_cat) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2,col = exp_cat)) + geom_point(alpha = 0.5,size = 0.05) +
  # facet_wrap(gene~NMDA_time,nrow = 2) +
  # facet_rep_wrap(gene~treat_full,repeat.tick.labels = "all",nrow=3)+
  facet_wrap(gene~treat_full)+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme_cowplot() +
  scale_color_manual(values = c("gray","blue")) +
  # theme(strip.background = element_blank(),
  #       panel.border = element_rect(colour = "black", fill = NA))
  theme(strip.background = element_blank())
ggsave("../../out/image/UMAPggplot_annotationConfident_TSPO_proppos.pdf",width = 13,height = 12)

# violin plot for GOI expression use macro categories
df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24","NMDA_36")) %>%
  # this is the processing shown in the violinplot function
  # mutate(exp_fix = exp + rnorm(nrow(.))/100000) %>%
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=treat_full,y=exp_fix)) + 
  geom_violin(scale = "width")+
  geom_point(position=position_jitter(width = 0.2),alpha=0.01) +
  facet_wrap(~expertAnno.l1) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1,angle = 90)) +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
ggsave("../../out/image/violin_annotationConfident_TSPO.pdf",width = 15,height = 10)

# try to depict the average expression there is roughly one sample per condition
df_avg <- average_GOI$RNA %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  mutate(gene = GOI) |> 
  pivot_longer(names_to = "group",values_to = "avg_exp",-gene) %>%
  filter(!str_detect(group,pattern="doublet|unassigned")) |> 
  mutate(treat_full = str_extract(group,pattern = c("BASELINE|CSF.ctrl.24h|CSF.MS.24h|CSF.MS.48h|cytokine|Fe|myelin|TBHP"))) |> 
  mutate(donor = str_extract(group,pattern = c("don\\w+"))) |> 
  mutate(expertAnno.l1 = str_remove_all(group,pattern = c("BASELINE.|CSF.ctrl.24h.|CSF.MS.24h.|CSF.MS.48h.|cytokine.|Fe.|myelin.|TBHP.|.don\\w+"))) |> 
  # add a more general treat annotation
  mutate(treat = case_when(treat_full%in%c("CSF.MS.24h","CSF.MS.48h")~"CSF.MS",
                           treat_full%in%c("CSF.ctrl.24h")~"CSF.ctrl",
                           T~treat_full))

# plot splitting by treat full
df_avg |>
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=treat_full,y=avg_exp))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),alpha = 0.6)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+
  facet_wrap(~expertAnno.l1,scales = "free")
# scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
ggsave("../../out/image/dotplot_annotationConfident_TSPO_expressionAvg_treatFull.pdf",width = 9,height = 9)

# plot splitting by treat
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
  facet_wrap(~expertAnno.l1,scales = "free")
# scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
ggsave("../../out/image/dotplot_annotationConfident_TSPO_expressionAvg_treat.pdf",width = 9,height = 9)

# try to keep the same scale
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
  facet_wrap(~expertAnno.l1)
# scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
ggsave("../../out/image/dotplot_annotationConfident_TSPO_expressionAvg_treat_scale.pdf",width = 9,height = 9)

# add also the doubtles as donor
df_avg2 <- average_GOI$RNA %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  mutate(gene = GOI) |> 
  pivot_longer(names_to = "group",values_to = "avg_exp",-gene) %>%
  filter(!str_detect(group,pattern="unassigned")) |> 
  mutate(treat_full = str_extract(group,pattern = c("BASELINE|CSF.ctrl.24h|CSF.MS.24h|CSF.MS.48h|cytokine|Fe|myelin|TBHP"))) |> 
  mutate(donor = str_extract(group,pattern = c("don\\w+|doublet"))) |> 
  mutate(expertAnno.l1 = str_remove_all(group,pattern = c("BASELINE.|CSF.ctrl.24h.|CSF.MS.24h.|CSF.MS.48h.|cytokine.|Fe.|myelin.|TBHP.|.don\\w+|.doublet"))) |> 
  # add a more general treat annotation
  mutate(treat = case_when(treat_full%in%c("CSF.MS.24h","CSF.MS.48h")~"CSF.MS",
                           treat_full%in%c("CSF.ctrl.24h")~"CSF.ctrl",
                           T~treat_full))

# save the table
df_avg2 |> 
  write_tsv("../../out/table/df_avg_cellTypeTreat_TSPO_BS.tsv")

# plot splitting by treat full
df_avg2 |>
  mutate(donor_id = case_when(donor =="doublet"~"doublet",
                              T~"donor")) |> 
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=treat_full,y=avg_exp))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=donor_id),position = position_jitter(width = 0.1),alpha = 0.6)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 45))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+scale_color_manual(values = c("black","red"))+
  # theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  # theme(strip.background = element_blank(),
  #       panel.border = element_rect(colour = "black", fill = NA))+
  facet_wrap(~expertAnno.l1,scales = "free")
# scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
ggsave("../../out/image/dotplot_annotationConfident_TSPO_expressionAvg_treatFull2.pdf",width = 9,height = 9)

# check is it is a marker in the cell type markes
# library(SeuratWrappers)
# data.combined$group2 <- paste0(data.combined$pathology_class,"-",data.combined$annotation_confident)
# DefaultAssay(data.combined) <- "RNA"

Idents(data.combined) <- "expertAnno.l1"
# subset the object to include only the one of interest
# cell_id <- "MG"

list_test <- lapply(names(table(data.combined$expertAnno.l1)),function(cell_id){
  # check the progress
  print(cell_id)
  test_obj <- subset(data.combined,subset = expertAnno.l1 == cell_id)
  Idents(test_obj) <- "treat_full"
  # return all the genes to show the realtive position of TSPO treatment
  test <- RunPrestoAll(test_obj, min.pct = 0.05, logfc.threshold = 0,return.thresh=1)
  
  test |> 
    mutate(cell_id = cell_id)
})

# save a table
df_test_subset <- list_test |> 
  bind_rows()

df_test_subset |> 
  write_tsv("../../out/table/FindAllMarkers_cellID_treat_minpct5_logfcthr01.tsv")

df_test_subset <- read_tsv("../../out/table/FindAllMarkers_cellID_treat_minpct5_logfcthr01.tsv")

# check if TSPO is in there
df_test_subset |> 
  filter(gene %in% GOI) |> 
  ggplot(aes(x = avg_log2FC,y = -log(p_val_adj),col=cluster,label=cluster))+geom_point()+facet_wrap(~cell_id,scales = "free")+theme_bw()+theme(strip.background = element_blank())+
  geom_vline(xintercept = c(0),linetype="dashed",col="gray")+
  geom_vline(xintercept = c(-0.5,0.5),linetype="dashed",col="red")+
  geom_hline(yintercept = c(-log(0.05)),linetype="dashed",col="red")+
  ggrepel::geom_text_repel(force = 100,segment.alpha=0.5,max.overlaps = 5)
ggsave("../../out/image/DE_annotationConfident_TSPO_sc.pdf",width = 12,height = 9)

# make the stats using a pairwise test on the average expression
# in this case the origin and the sample are exclusives, therefore there is no need to generate another aggregation
df_avg |> 
  finalfit(formula = avg_exp~expertAnno.l1+treat_full)

df_avg |> 
  finalfit(formula = avg_exp~expertAnno.l1+treat)

# cell_id <- "0"
id_test <- names(table(df_avg$expertAnno.l1)) 
  # remove cluster 14 and 16 since there is no relicate
  # str_subset(pattern = c("14|16"),negate = T)

list_pairwise <- lapply(id_test,function(cell_id){
  # check the progress
  print(cell_id)
  # subset the dataset to the single cell tyep
  test <- df_avg |> 
    filter(expertAnno.l1 %in% cell_id)
  
  # run the pairwise test
  test_out <- pairwise.t.test(test$avg_exp, test$treat_full,p.adjust.method = "fdr")
  
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
  geom_hline(data = df_intercept,aes(yintercept = intercept),col="red",linetype="dashed")+theme_bw()+theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,vjust = 1,angle = 45))
ggsave("../../out/image/DE_annotationConfident_TSPO_bulk_PairWiseTest.pdf",width = 9,height = 9)

# collapse the MS treatments
list_pairwise2 <- lapply(id_test,function(cell_id){
  # check the progress
  print(cell_id)
  # subset the dataset to the single cell tyep
  test <- df_avg |> 
    filter(expertAnno.l1 %in% cell_id)
  
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
df_intercept2 <- list_pairwise2 |> 
  bind_rows() |>
  arrange(padj) |>
  group_by(first,second) |> 
  summarise() |> 
  mutate(intercept = -log(0.05))

list_pairwise2 |> 
  bind_rows() |>
  arrange(padj) |> 
  # mutate(comparison = paste(first,second,sep = " vs ")) |> 
  ggplot(aes(x=cell_id,y=-log(padj)))+
  facet_grid(first~second,drop = T)+
  geom_point()+
  geom_hline(data = df_intercept2,aes(yintercept = intercept),col="red",linetype="dashed")+theme_bw()+theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,vjust = 1,angle = 45))
ggsave("../../out/image/DE_annotationConfident_TSPO_bulk_PairWiseTest2.pdf",width = 9,height = 9)

# explore expression ------------------------------------------------------
data.combined$group_cellType <- paste0(data.combined$harmonized_donor2,"-",data.combined$expertAnno.l1)
# data.combined$group2 <- paste0(data.combined$orig.ident,".",data.combined$treat,".",data.combined$cell_type2)
# Idents(data.combined) <- "group_cellType"
DefaultAssay(data.combined) <- "RNA"

average_group_cellType <- AverageExpression(data.combined,features = GOI,group.by = c("group_cellType"))

df_avg_cellType <- average_group_cellType$RNA %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  mutate(gene = GOI) |> 
  pivot_longer(names_to = "group",values_to = "avg_exp",-gene) |> 
  filter(!str_detect(group,pattern="doublet|unassigned")) |> 
  mutate(donor = str_extract(group,pattern = c("don\\w+"))) |> 
  mutate(expertAnno.l1 = str_remove_all(group,pattern = c("don\\w+.")))

df_avg_cellType |>
  mutate(log1p = log1p(avg_exp)) |> 
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=expertAnno.l1,y=avg_exp))+
  # ggplot(aes(x=annotation_confident,y=log1p,col=origin))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),alpha = 0.2)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 45))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
# scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
ggsave("../../out/image/dotplot_annotationConfident_TSPO_expressionAvg_cellType.pdf",width = 6,height = 5)

# martina asked to plot also the doublets
df_avg_cellType2 <- average_group_cellType$RNA %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  mutate(gene = GOI) |> 
  pivot_longer(names_to = "group",values_to = "avg_exp",-gene) |> 
  filter(!str_detect(group,pattern="unassigned")) |> 
  mutate(donor = str_extract(group,pattern = c("don\\w+|doublet"))) |> 
  mutate(expertAnno.l1 = str_remove_all(group,pattern = c("don\\w+.|doublet.")))

df_avg_cellType2 |>
  # label the doublets
  mutate(donor_id = case_when(donor =="doublet"~"doublet",
                           T~"donor")) |> 
  mutate(donor_id = factor(donor_id)) |> 
  mutate(log1p = log1p(avg_exp)) |> 
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=expertAnno.l1,y=avg_exp))+
  # ggplot(aes(x=annotation_confident,y=log1p,col=origin))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=donor_id),position = position_jitter(width = 0.1),alpha = 0.8)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 45))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+scale_color_manual(values = c("black","red"))
# scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
ggsave("../../out/image/dotplot_annotationConfident_TSPO_expressionAvg_cellType2.pdf",width = 6,height = 5)

# save the table
df_avg_cellType2 |> 
  write_tsv("../../out/table/df_avg_cellType_TSPO_BS.tsv")

# check is it is a marker in the cell type markes
# library(SeuratWrappers)
Idents(data.combined) <- "expertAnno.l1"
sobj_total_h.markers <- RunPrestoAll(data.combined, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.1)

sobj_total_h.markers |>
  write_tsv("../../out/table/FindAllMarkers_annotation_confident_PosOnly_minpct5_logfcthr01.tsv")

sobj_total_h.markers <- read_tsv("../../out/table/FindAllMarkers_annotation_confident_PosOnly_minpct5_logfcthr01.tsv")

sobj_total_h.markers |> 
  filter(gene %in% GOI)

# make the stats using a linear model on the average expression
# in this case the origin and the sample are exclusives, therefore there is no need to generate another aggregation
# use progenitors as reference
names(table(df_avg$expertAnno.l1)) 
df_avg_cellType |> 
  # mutate(expertAnno.l1 = factor(expertAnno.l1,levels = c("PROG","ASTRO","GLIA_IMM","MG","NEU","OLIGO","OPC"))) |> 
  mutate(expertAnno.l1 = factor(expertAnno.l1,levels = c("NEU","PROG","ASTRO","GLIA_IMM","MG","OLIGO","OPC"))) |> 
  finalfit(formula = avg_exp~expertAnno.l1)

df_avg_cellType |>
  # mutate(expertAnno.l1 = factor(expertAnno.l1,levels = c("PROG","ASTRO","GLIA_IMM","MG","NEU","OLIGO","OPC"))) |> 
  mutate(expertAnno.l1 = factor(expertAnno.l1,levels = c("NEU","PROG","ASTRO","GLIA_IMM","MG","OLIGO","OPC"))) |> 
  data.frame() |> 
  ff_plot(dependent = "avg_exp",explanatory = "expertAnno.l1")

# disease stage -----------------------------------------------------------
data.combined$group_disease <- paste0(data.combined$treat_full,"-",data.combined$harmonized_donor2)
# data.combined$group2 <- paste0(data.combined$orig.ident,".",data.combined$treat,".",data.combined$cell_type2)
# Idents(data.combined) <- "group_cellType"
DefaultAssay(data.combined) <- "RNA"

average_group_disease <- AverageExpression(data.combined,features = GOI,group.by = c("group_disease"))

df_avg_disease <- average_group_disease$RNA %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  mutate(gene = GOI) |> 
  pivot_longer(names_to = "group",values_to = "avg_exp",-gene) %>%
  filter(!str_detect(group,pattern="doublet|unassigned")) |> 
  mutate(donor = str_extract(group,pattern = c("don\\w+"))) |> 
  mutate(treat_full = str_extract(group,pattern = c("BASELINE|CSF.ctrl.24h|CSF.MS.24h|CSF.MS.48h|cytokine|Fe|myelin|TBHP"))) |> 
  mutate(treat = case_when(treat_full%in%c("CSF.MS.24h","CSF.MS.48h")~"CSF.MS",
                           treat_full%in%c("CSF.ctrl.24h")~"CSF.ctrl",
                           T~treat_full))

# save the table of average expression
df_avg_disease |> 
  write_tsv("../../out/table/df_avg_disease.tsv")

# keep the doublets
df_avg_disease2 <- average_group_disease$RNA %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  mutate(gene = GOI) |> 
  pivot_longer(names_to = "group",values_to = "avg_exp",-gene) %>%
  filter(!str_detect(group,pattern="unassigned")) |> 
  mutate(donor = str_extract(group,pattern = c("don\\w+|doublet"))) |> 
  mutate(treat_full = str_extract(group,pattern = c("BASELINE|CSF.ctrl.24h|CSF.MS.24h|CSF.MS.48h|cytokine|Fe|myelin|TBHP"))) |> 
  mutate(treat = case_when(treat_full%in%c("CSF.MS.24h","CSF.MS.48h")~"CSF.MS",
                           treat_full%in%c("CSF.ctrl.24h")~"CSF.ctrl",
                           T~treat_full))

# save the table
df_avg_disease2 |> 
  write_tsv("../../out/table/df_avg_disease_TSPO_BS.tsv")

df_avg_disease |>
  mutate(log1p = log1p(avg_exp)) |> 
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=treat_full,y=avg_exp))+
  # ggplot(aes(x=annotation_confident,y=log1p,col=origin))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),alpha = 0.2)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 45))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
# scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
ggsave("../../out/image/dotplot_annotationConfident_TSPO_expressionAvg_disease.pdf",width = 6,height = 5)

df_avg_disease2 |>
  # label the doublets
  mutate(donor_id = case_when(donor =="doublet"~"doublet",
                              T~"donor")) |> 
  mutate(donor_id = factor(donor_id)) |>
  mutate(log1p = log1p(avg_exp)) |> 
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=treat,y=avg_exp))+
  # ggplot(aes(x=annotation_confident,y=log1p,col=origin))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(col=donor_id),position = position_jitter(width = 0.1),alpha = 0.2)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 45))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+scale_color_manual(values = c("black","red"))
  # theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 45))+
  # theme(strip.background = element_blank(),
  #       panel.border = element_rect(colour = "black", fill = NA))
# scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
ggsave("../../out/image/dotplot_annotationConfident_TSPO_expressionAvg_disease2.pdf",width = 6,height = 5)

# check is it is a marker in the cell type markes
# library(SeuratWrappers)
Idents(data.combined) <- "treat_full"
sobj_total_h.markers2 <- RunPrestoAll(data.combined, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.1)

sobj_total_h.markers2 |> 
  filter(gene %in% GOI)

# make the stats using a linear model on the average expression
# in this case the origin and the sample are exclusives, therefore there is no need to generate another aggregation
df_avg_disease |> 
  finalfit(formula = avg_exp~treat_full)

df_avg_disease |>
  data.frame() |> 
  ff_plot(dependent = "avg_exp",explanatory = "treat_full")

df_avg_disease |>
  data.frame() |> 
  ff_plot(dependent = "avg_exp",explanatory = "treat")

# test pos TSPO vs neg TSPO -----------------------------------------------
# subset only the MG cells. define the positive and negative cells per TSPO and run a de analysis to compare pull the genes
# subset onlyt he MG cells
sobj_MG <- subset(data.combined,subset = expertAnno.l1 == 'MG')

sobj_MG@meta.data

# for each cell fetch TSPO expression
sobj_MG_TSPO <- FetchData(sobj_MG,"TSPO",slot = "data") |> 
  rownames_to_column("barcodes")

# add the TSPO expression to the original meta
sobj_MG$TSPO <- sobj_MG@meta.data |> 
  rownames_to_column("barcodes") |> 
  left_join(sobj_MG_TSPO) |> 
  pull("TSPO")

# build a category form the expression value
sobj_MG$TSPO_cat <- case_when(sobj_MG$TSPO == 0~"neg",
                              T~"pos")

# summarise the number of cells
table(sobj_MG$TSPO_cat)

# ggplot(aes(x=TSPO))+geom_histogram()

# run the DE
Idents(sobj_MG) <- "TSPO_cat"

# avg_logFC: log fold-chage of the average expression between the two groups. Positive values indicate that the gene is more highly expressed in the first group 
res_test_TSPO_MG <- RunPresto(object = sobj_MG,ident.1 = "pos",ident.2 = "neg")

# save the table of top markers
res_test_TSPO_MG |> 
  rownames_to_column("gene") |> 
  mutate(cell_id = "MG") |> 
  write_tsv("../../out/table/res_test_TSPO_MG.tsv")

# plot volcano
volcano_tot <- read_tsv("../../out/table/res_test_TSPO_MG.tsv") %>%
  mutate(DE_cat = case_when(avg_log2FC > 0.5 & p_val_adj < 0.01~"up",
                            avg_log2FC < (-0.5) & p_val_adj < 0.01~"down",
                            T~"no"))

ggplot() +
  geom_point(data = volcano_tot%>%filter(DE_cat=="no"),aes(x=avg_log2FC,y=-log10(p_val_adj)),size=1,alpha=0.2)+theme_bw()+
  geom_point(data = volcano_tot%>%filter(DE_cat!="no"),aes(x=avg_log2FC,y=-log10(p_val_adj)),size=1,alpha=0.2,col="red")+theme_bw()+
  ggrepel::geom_text_repel(data = volcano_tot%>%filter(DE_cat!="no"),aes(x=avg_log2FC,y=-log10(p_val_adj),label=gene))+theme_bw()+theme(strip.background = element_blank())
ggsave("../../out/image/volcano_test_TSPO_MG.pdf",width = 12,height = 9)

# explore the expression only inside the positive cells -------------------
# from the MG cell extract only the positive cells
sobj_MG_TSPOpos <- subset(sobj_MG,subset = TSPO_cat == 'pos')

sobj_MG_TSPOpos@meta.data

# run the DE
Idents(sobj_MG_TSPOpos) <- "treat_full"

# how many cells are in the different categories
table(Idents(sobj_MG_TSPOpos))

# avg_logFC: log fold-chage of the average expression between the two groups. Positive values indicate that the gene is more highly expressed in the first group
# vs baseline
df_res_TSPOpos <- lapply(c("myelin","CSF.MS.24h","CSF.ctrl.24h","cytokine","CSF.MS.48h","TBHP","Fe"),function(x){
  # track the progress
  print(x)
  res_MG_TSPOpos <- RunPresto(object = sobj_MG_TSPOpos,ident.1 = x,ident.2 = "BASELINE", logfc.threshold = 0) %>%
    rownames_to_column("gene") %>%
    mutate(comparison = paste0(x,"_vs_BASELINE"))
  return(res_MG_TSPOpos)
}) %>%
  bind_rows()

# extract the expression of BTK
df_res_TSPOpos %>%
  filter(gene %in% c("TSPO"))

# vs CSF.ctrl.24h
df_res_TSPOpos2 <- lapply(c("myelin","CSF.MS.24h","BASELINE","cytokine","CSF.MS.48h","TBHP","Fe"),function(x){
  # track the progress
  print(x)
  res_MG_TSPOpos <- RunPresto(object = sobj_MG_TSPOpos,ident.1 = x,ident.2 = "CSF.ctrl.24h",logfc.threshold = 0) %>%
    rownames_to_column("gene") %>%
    mutate(comparison = paste0(x,"_vs_CSF.ctrl.24h"))
  return(res_MG_TSPOpos)
}) %>%
  bind_rows()

# extract the expression of BTK
df_res_TSPOpos2 %>%
  filter(gene %in% c("TSPO"))

# # try a quick and dirty enrichR on the posirive and negatives
# library(tidyverse)
# library(enrichR)
# library(scales)
# library(patchwork)
# 
# # run enrichr with the list of genes in the module
# # DB selection ------------------------------------------------------------
# dbs <- listEnrichrDbs()
# #filter fo the db of interest
# dbs %>%
#   filter(str_detect(libraryName,pattern = "Atlas"))
# 
# dbs %>%
#   filter(str_detect(libraryName,pattern = "Cell"))
# 
# dbs_db <- c("KEGG_2021_Human","MSigDB_Hallmark_2020","Reactome_2016","HDSigDB_Human_2021","Azimuth_Cell_Types_2021","GO_Biological_Process_2023","Descartes_Cell_Types_and_Tissue_2021","CellMarker_Augmented_2021")
# 
# # query -------------------------------------------------------------------
# # seelct only the clusters with more than 10 genes as degs
# 
# # pull the gene names dividing the up regulated from the downregulated
# list_genes <- list(list_UP = volcano_tot %>% filter(DE_cat=="up") %>% pull(gene),
#                    list_DOWN = volcano_tot %>% filter(DE_cat=="down") %>% pull(gene))
# 
# # define the background
# # background <- df_modules$feature
# 
# # x <- list_res_tot_UP_filter$`DeMye_vs_Norm|clust_5`
# list_enrichr <- lapply(list_genes,function(x){
#   genes <- x
#   # out_enrich <- enrichr(genes, dbs_db,background = background,include_overlap=T)
#   out_enrich <- enrichr(genes, dbs_db)
#   
#   # filter out the annotations without an output
#   filter_out <- lapply(out_enrich,function(x){dim(x)[1]}) %>% 
#     unlist()
#   
#   out_enrich[filter_out>0] %>%
#     bind_rows(.id = "annotation")
# }) %>%
#   bind_rows(.id = "comparison")
# 
# list_enrichr %>%
#   write_tsv("../../out/table/enrichR_test_TSPO_MG.tsv")
# 
# # list  <- read_tsv("out/table/enrichR_out_all_filtered_clusters_CTRL_vs_CSF_harmonyMartina.tsv")
# 
# plot_list_UP <- list_enrichr %>%
#   split(f = .$comparison)
# 
# # library(scales)
# list_plot_UP <- pmap(list(plot_list_UP,names(plot_list_UP)), function(x,y){
#   x %>%
#     group_by(annotation) %>%
#     arrange(P.value) %>%
#     dplyr::slice(1:10) %>%
#     mutate(Term = str_sub(Term,start = 1,end = 30)) %>%
#     mutate(Term = fct_reorder(Term, Combined.Score,.desc = F)) %>%
#     # ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,scales = "free",ncol = 1)+theme_bw() +
#     ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,scales = "free",ncol = 1)+theme_bw() +
#     scale_color_gradientn(colors = c("red","blue"),
#                           values = rescale(c(0,1)),
#                           limits = c(0,0.2))+
#     theme(strip.background = element_blank(),
#           panel.border = element_rect(colour = "black", fill = NA))+
#     ggtitle(y)
#   # scale_color_gradient(low = "red",high = "blue")
#   
#   #ggsave(paste0("image/enrichR_out_",y,".pdf"),width = 7,height = 15)
# })
# 
# wrap_plots(list_plot_UP)
# ggsave("../../out/image/enrichR_test_TSPO_MG.pdf",width = 13,height = 25,limitsize = FALSE)
# 
