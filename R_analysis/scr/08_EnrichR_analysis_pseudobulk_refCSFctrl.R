# libraries ---------------------------------------------------------------
library(tidyverse)
library(enrichR)
library(scales)
library(patchwork)

# run enrichr with the list of genes in the module
# DB selection ------------------------------------------------------------
dbs <- listEnrichrDbs()
# filter fo the db of interest
dbs %>%
  filter(str_detect(libraryName,pattern = "Atlas"))

dbs %>%
  filter(str_detect(libraryName,pattern = "Cell"))

dbs_db <- c("KEGG_2021_Human","MSigDB_Hallmark_2020","Reactome_2016","HDSigDB_Human_2021","GO_Biological_Process_2023")

# query -------------------------------------------------------------------
# select only the clusters with more than 10 genes as degs
# file <- dir("../../out/table/") %>%
#   str_subset(pattern = "res_") %>%
#   str_subset(pattern = ".txt") %>%
#   str_subset(pattern = "HUVEC") %>%
#   str_subset(pattern = "shr",negate = F)
# file 

# load the results 
# results <- lapply(paste0("../../out/table/",file),function(x){
#   read_tsv(x)  %>% 
#     filter(padj<0.05 & abs(log2FoldChange)>1&!is.na(symbol))
# }) %>%
#   setNames(str_remove_all(file,pattern = ".txt"))
list_genes_UP <- read_tsv("../../out/table/DE_treatvsCSFctrl_pseudobulk_MG_shr.tsv") %>%
  filter(padj<0.05 & log2FoldChange > 1&!is.na(symbol)) %>%
  split(f = .$conditionVsCSFctrl) %>%
  map(function(x){
    x %>%
      pull(symbol)
  })

list_genes_DOWN <- read_tsv("../../out/table/DE_treatvsCSFctrl_pseudobulk_MG_shr.tsv") %>%
  filter(padj<0.05 & log2FoldChange < (-1)&!is.na(symbol)) %>%
  split(f = .$conditionVsCSFctrl) %>%
  map(function(x){
    x %>%
      pull(symbol)
  })

# # pull the gene names dividing the up regulated from the downregulated
# list_genes <- list(list_UP = results$res_GMPvsHUVEC_shr %>% filter(log2FoldChange>0) %>% pull(symbol),
#                    list_DOWN = results$res_GMPvsHUVEC_shr %>% filter(log2FoldChange<0) %>% pull(symbol))

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

list_enrichr_UP %>%
  write_tsv("../../out/table/enrichR_DE_treatvsCSFctrl_pseudobulk_MG_UP_shr.tsv")

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

list_enrichr_DOWN %>%
  write_tsv("../../out/table/enrichR_DE_treatvsCSFctrl_pseudobulk_MG_DOWN_shr.tsv")

# see the top term from azimuth cell type
list_enrichr_UP %>%
  as.tibble() %>% 
  # group_by(annotation) %>% summarise()
  filter(annotation == "MSigDB_Hallmark_2020",
         comparison == "CSF.MS.48h_shr") %>% 
  # select(Term,Overlap,Genes,Adjusted.P.value,Odds.Ratio) %>% 
  data.frame()  %>% 
  head()

list_enrichr_UP %>%
  as.tibble() %>% 
  # group_by(annotation) %>% summarise()
  filter(annotation == "MSigDB_Hallmark_2020",
         comparison == "cytokine_shr") %>% 
  # select(Term,Overlap,Genes,Adjusted.P.value,Odds.Ratio) %>% 
  data.frame()  %>% 
  head()

list_enrichr_UP %>%
  as.tibble() %>% 
  # group_by(annotation) %>% summarise()
  filter(annotation == "MSigDB_Hallmark_2020",
         comparison == "Fe_shr") %>% 
  # select(Term,Overlap,Genes,Adjusted.P.value,Odds.Ratio) %>% 
  data.frame()  %>% 
  head()

# list  <- read_tsv("out/table/enrichR_out_all_filtered_clusters_CTRL_vs_CSF_harmonyMartina.tsv")

plot_list_UP <- list_enrichr_UP %>%
  split(f = .$comparison)

plot_list_DOWN <- list_enrichr_DOWN %>%
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

wrap_plots(list_plot_UP,nrow = 1)
ggsave("../../out/image/enrichR_DE_treatvsCSFctrl_pseudobulk_MG_UP_shr.pdf",width = 40,height = 12,limitsize = FALSE)

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

wrap_plots(list_plot_DOWN,nrow = 1)
ggsave("../../out/image/enrichR_DE_treatvsCSFctrl_pseudobulk_MG_DOWN_shr.pdf",width = 40,height = 12,limitsize = FALSE)

# -------------------------------------------------------------------------
read_tsv("../../out/table/DE_treatvsCSFctrl_pseudobulk_MG_shr.tsv") |> 
  filter(symbol %in% c("HEXB"))
