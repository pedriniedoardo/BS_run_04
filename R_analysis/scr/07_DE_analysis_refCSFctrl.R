# libraries ---------------------------------------------------------------
library(Seurat)
library(SeuratWrappers)
library(tidyverse)
library(scales)
library(ComplexHeatmap)
library(lemon)
library(finalfit)
library(enrichR)
library(patchwork)

# read in the dataset -----------------------------------------------------
data.combined <- readRDS("../../out/object/sobj_processed_donor.rds")
DimPlot(data.combined,label = T,raster = T,group.by = "seurat_clusters")
DimPlot(data.combined,label = T,raster = T,group.by = "expertAnno.l1")
DimPlot(data.combined,label = T,raster = T,group.by = "treat_full")

# run DE analysis ---------------------------------------------------------
# DefaultAssay(data.combined) <- "RNA"

Idents(data.combined) <- "expertAnno.l1"
# subset the object to include only the one of interest
# cell_id <- "MG"

# split the data per cell id wise and run DE comparing each conditions vs BASELINE.
list_test <- lapply(names(table(data.combined$expertAnno.l1)),function(cell_id){
  # check the progress
  print(cell_id)
  test_obj <- subset(data.combined,subset = expertAnno.l1 == cell_id)
  Idents(test_obj) <- "treat_full"

  # run the test using as reference CSF.ctrl.24h vs all the other treatments one by one
  # loop the testing
  list_DE <- lapply(c("myelin","CSF.MS.24h","BASELINE","cytokine","CSF.MS.48h","TBHP","Fe"),function(condition){
    print(paste(cell_id,condition))
    test <- RunPresto(test_obj, ident.1 = condition,ident.2 = "CSF.ctrl.24h",min.pct = 0.05, logfc.threshold = 0,return.thresh=1,only.pos = F) %>%
      mutate(cell_id = cell_id) %>%
      mutate(ident.1 = condition,
             ident.2 = "CSF.ctrl.24h") |>
      rownames_to_column("gene")
  })

  # pull all the tables in one
  list_DE |>
    bind_rows()
})

# save the list in a single table
df_test_subset <- list_test |>
  bind_rows()

df_test_subset |>
  write_tsv("../../out/table/DE_treatvsCSFctrl_cellIDwise_minpct5.tsv")

# plot volcano ------------------------------------------------------------
# plot the volcano for the cluster 12 DE
volcano_tot <- read_tsv("../../out/table/DE_treatvsCSFctrl_cellIDwise_minpct5.tsv") %>%
  mutate(DE_cat = case_when(avg_log2FC > 1 & p_val_adj < 0.01~"up",
                            avg_log2FC < (-1) & p_val_adj < 0.01~"down",
                            T~"no")) %>%
  mutate(test = paste0(ident.1,"_",ident.2))
# filter(cluster %in% c("IMM"))

ggplot() +
  geom_point(data = volcano_tot%>%filter(DE_cat=="no"),aes(x=avg_log2FC,y=-log10(p_val_adj)),size=1,alpha=0.2)+theme_bw()+
  geom_point(data = volcano_tot%>%filter(DE_cat!="no"),aes(x=avg_log2FC,y=-log10(p_val_adj)),size=1,alpha=0.2,col="red")+
  theme_bw()+
  ggrepel::geom_text_repel(data = volcano_tot%>%filter(DE_cat!="no"),aes(x=avg_log2FC,y=-log10(p_val_adj),label=gene))+
  theme_bw()+
  facet_grid(test~cell_id)+
  theme(strip.background = element_blank())
ggsave("../../out/image/volcano_DE_treatvsCSFctrl_cellIDwise_minpct5.pdf",width = 25,height = 20)

ggplot() +
  geom_point(data = volcano_tot%>%filter(DE_cat=="no"),aes(x=avg_log2FC,y=-log10(p_val_adj)),size=1,alpha=0.2)+theme_bw()+
  geom_point(data = volcano_tot%>%filter(DE_cat!="no"),aes(x=avg_log2FC,y=-log10(p_val_adj)),size=1,alpha=0.2,col="red")+
  theme_bw()+
  # ggrepel::geom_text_repel(data = volcano_tot%>%filter(DE_cat!="no"),aes(x=avg_log2FC,y=-log10(p_val_adj),label=gene))+
  theme_bw()+
  facet_grid(test~cell_id)+
  theme(strip.background = element_blank())
ggsave("../../out/image/volcano_DE_treatvsCSFctrl_cellIDwise_minpct5_noLab.pdf",width = 25,height = 20)

# check if TSPO is in any comparision
volcano_tot |>
  filter(gene %in% c("TSPO")) |>
  arrange(desc(avg_log2FC))

# run enrichR on each comparison ------------------------------------------
dbs <- listEnrichrDbs()

# filter fo the db of interest
# dbs %>%
#   filter(str_detect(libraryName,pattern = "Atlas"))
#
# dbs %>%
#   filter(str_detect(libraryName,pattern = "Cell"))

dbs_db <- c("KEGG_2021_Human","MSigDB_Hallmark_2020","Reactome_2016","HDSigDB_Human_2021","GO_Biological_Process_2023")

# query -------------------------------------------------------------------
# seelct only the clusters with more than 10 genes as degs

# pull the gene names dividing the up regulated from the downregulated by cell_id and group
list_genes_UP <- volcano_tot %>%
  mutate(test2 = paste0(.$cell_id,"|",.$test,"|UP")) %>%
  filter(DE_cat=="up") %>%
  split(f = .$test2) %>%
  map(function(x){
    x %>%
      pull(gene)
  })
# pull the downregulated
list_genes_DOWN <- volcano_tot %>%
  mutate(test2 = paste0(.$cell_id,"|",.$test,"|DOWN")) %>%
  filter(DE_cat=="down") %>%
  split(f = .$test2) %>%
  map(function(x){
    x %>%
      pull(gene)
  })

# define the background
# background <- df_modules$feature

# x <- list_res_tot_UP_filter$`DeMye_vs_Norm|clust_5`
list_enrichr_UP <- lapply(list_genes_UP,function(x){
  genes <- x
  # out_enrich <- enrichr(genes, dbs_db,background = background,include_overlap=T)
  out_enrich <- enrichr(genes, dbs_db)

  # filter out the annotations without an output
  filter_out <- lapply(out_enrich,function(x){dim(x)[1]}) %>%
    unlist()

  out_enrich[filter_out>0] %>%
    bind_rows(.id = "annotation")
}) %>%
  bind_rows(.id = "comparison")

# run the same for the downregulated
list_enrichr_DOWN <- lapply(list_genes_DOWN,function(x){
  genes <- x
  # out_enrich <- enrichr(genes, dbs_db,background = background,include_overlap=T)
  out_enrich <- enrichr(genes, dbs_db)

  # filter out the annotations without an output
  filter_out <- lapply(out_enrich,function(x){dim(x)[1]}) %>%
    unlist()

  out_enrich[filter_out>0] %>%
    bind_rows(.id = "annotation")
}) %>%
  bind_rows(.id = "comparison")

# save the tables
list_enrichr_UP %>%
  write_tsv("../../out/table/enrichR_DE_treatvsCSFctrl_cellIDwise_minpct5_UP.tsv")

list_enrichr_DOWN %>%
  write_tsv("../../out/table/enrichR_DE_treatvsCSFctrl_cellIDwise_minpct5_DOWN.tsv")

# plot the results --------------------------------------------------------
plot_list_UP <- list_enrichr_UP %>%
  split(f = .$comparison)

# library(scales)
list_plot_UP <- pmap(list(plot_list_UP,names(plot_list_UP)), function(x,y){
  x %>%
    group_by(annotation) %>%
    arrange(P.value) %>%
    dplyr::slice(1:10) %>%
    mutate(Term = str_sub(Term,start = 1,end = 30)) %>%
    mutate(Term = fct_reorder(Term, Combined.Score,.desc = F)) %>%
    # ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,scales = "free",ncol = 1)+theme_bw() +
    ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,scales = "free",ncol = 1)+theme_bw() +
    scale_color_gradientn(colors = c("red","blue"),
                          values = rescale(c(0,1)),
                          limits = c(0,0.2))+
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))+
    ggtitle(y)
  # scale_color_gradient(low = "red",high = "blue")

  #ggsave(paste0("image/enrichR_out_",y,".pdf"),width = 7,height = 15)
})

wrap_plots(list_plot_UP)
ggsave("../../out/image/enrichR_DE_treatvsCSFctrl_cellIDwise_minpct5_UP.pdf",width = 40,height = 60,limitsize = FALSE)

# down genes
plot_list_DOWN <- list_enrichr_DOWN %>%
  split(f = .$comparison)

# library(scales)
list_plot_DOWN <- pmap(list(plot_list_DOWN,names(plot_list_DOWN)), function(x,y){
  x %>%
    group_by(annotation) %>%
    arrange(P.value) %>%
    dplyr::slice(1:10) %>%
    mutate(Term = str_sub(Term,start = 1,end = 30)) %>%
    mutate(Term = fct_reorder(Term, Combined.Score,.desc = F)) %>%
    # ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,scales = "free",ncol = 1)+theme_bw() +
    ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,scales = "free",ncol = 1)+theme_bw() +
    scale_color_gradientn(colors = c("red","blue"),
                          values = rescale(c(0,1)),
                          limits = c(0,0.2))+
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))+
    ggtitle(y)
  # scale_color_gradient(low = "red",high = "blue")

  #ggsave(paste0("image/enrichR_out_",y,".pdf"),width = 7,height = 15)
})

wrap_plots(list_plot_DOWN)
ggsave("../../out/image/enrichR_DE_treatvsCSFctrl_cellIDwise_minpct5_DOWN.pdf",width = 40,height = 60,limitsize = FALSE)

# confirm some of the trends
list_enrichr_UP |>
  filter(comparison == "MG|CSF.MS.48h_CSF.ctrl.24h|UP") |>
  head()

# JUN seems to be an upregulated genes in MS.48h
df_test_subset |>
  filter(gene %in% c("JUN"),cell_id %in% c("MG"))

FeaturePlot(object = data.combined,features = "JUN",split.by = "treat_full",order = T)

# run GSEA on the same comparison -----------------------------------------
