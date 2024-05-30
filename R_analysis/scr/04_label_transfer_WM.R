# # libraries ---------------------------------------------------------------
# library(tidyverse)
# library(Seurat)
# # library(SeuratData)
# # library(ggridges)
# # library(ComplexHeatmap)
# 
# # read in the data --------------------------------------------------------
# # read in the full annotations from the label transfer
# list_out <- readRDS("../../out/object/list_out_label_transfer_BS0102.rds")
# # read in the reference IMM WM
# ref_WM <- readRDS("../../data/all20_immune_model.rds")
# 
# DimPlot(ref_WM)
# Idents(ref_WM)
# 
# # save ccordinates fro reference of the UMAP and metadata
# ref_WM_meta <- left_join(ref_WM@reductions$umap@cell.embeddings %>%
#                            data.frame() %>%
#                            rownames_to_column("barcodes"),
#                          ref_WM@meta.data %>%
#                            data.frame() %>%
#                            rownames_to_column("barcodes"),by="barcodes")
# 
# write_tsv(ref_WM_meta,file = "../../out/table/ref_WM_meta.tsv")
# 
# # wrangling ---------------------------------------------------------------
# # focus on the 02_60 feature_filter. select only the cells from cluster 12 and 13
# df_all <- list_out %>% 
#   bind_rows(.id = "sample") %>% 
#   mutate(feature_filter = str_extract(sample,pattern = "02_60|06_60|10_60")) %>% 
#   dplyr::filter(feature_filter == "02_60",
#                 predicted.id %in% c("12","13"))
#   
# # read in the data and subset the immune only cells from the BS
# folder <- "../../out/object/"
# id_file <- dir(folder) %>%
#   str_subset(pattern = "datasc_fix_filter_norm_doublet_Singlet") %>% 
#   str_subset(pattern = "02_60") %>% 
#   str_remove_all(pattern = "datasc_fix_filter_norm_doublet_Singlet_|.rds")
# 
# id_file
# # x <- id_file[1]
# 
# # loop the filtering procedure
# list_out <- lapply(id_file,function(x){
#   print(x)
#   # read in the query dataset
#   query <- readRDS(paste0(folder,"datasc_fix_filter_norm_doublet_Singlet_",x,".rds"))
#   # fix the meta to add the general cell barcode
#   # load the full meta
#   meta_full <- list_out[[x]]
#   
#   meta_fix <- query@meta.data %>%
#     data.frame() %>%
#     rownames_to_column("barcodes") %>%
#     dplyr::select(barcodes) %>%
#     left_join(meta_full,by = "barcodes") %>%
#     column_to_rownames("barcodes")
#   
#   query@meta.data <- meta_fix
#   
#   # -------------------------------------------------------------------------
#   # subset the object
#   query_imm <- subset(query,subset = predicted.id%in%c("12","13"))
#   
#   # save the metadata with the new clsuter informations
#   return(query_imm)
# }) %>% 
#   setNames(id_file)
# 
# # save the list af all the immune cells
# saveRDS(list_out,"../../out/object/list_doublet_Singlet_02_60_IMM.rds")
# 
# # loop the label transfer procedure ---------------------------------------
# # load the list of the immune cells
# list_query <- readRDS("../../out/object/list_doublet_Singlet_02_60_IMM.rds")
# 
# x <- list_query[[1]]
# # loop the label transfer vs all the 
# list_label_transfer <- lapply(list_query,function(x){
#   print(x)
#   # read in the query dataset
#   query <- x
#   
#   # -------------------------------------------------------------------------
#   # to work the reference should be the integrated dataset
#   DefaultAssay(ref_WM)<-"integrated"
#   # the RNA slot do not work due to the missing variable features
#   # DefaultAssay(ref_WM)<-"RNA"
#   
#   # in the vignetted the query isn the RNA slot
#   # DefaultAssay(query) <- "integrated"
#   DefaultAssay(query) <- "RNA"
#   test.anchors_0 <- FindTransferAnchors(reference = ref_WM,
#                                         query = query,
#                                         dims = 1:30,
#                                         reference.reduction = "pca",k.filter = NA)
#   
#   predictions_0 <- TransferData(anchorset = test.anchors_0,
#                                 refdata = ref_WM$seurat_clusters,
#                                 dims = 1:30)
#   
#   # add the predictions ot the query dataset
#   query_transfer <- AddMetaData(query, metadata = predictions_0)
#   
#   # -------------------------------------------------------------------------
#   # add the meta to the coordinates
#   data_transfer <- left_join(
#     # get the coordinates of the query dataset
#     query_transfer@reductions$umap@cell.embeddings %>%
#       data.frame() %>%
#       rownames_to_column(var = "barcodes"),
#     
#     # get the metadata from the query dataset, after adding the prediciton information. use a threshold for the score
#     query_transfer@meta.data %>%
#       rownames_to_column(var = "barcodes") %>%
#       # this is based on martina's table, I need to ask her how she provided the accuracy of the imputation
#       mutate(robust_score = case_when(prediction.score.max>0.70~predicted.id,
#                                       T ~ "uncertain"))
#     ,"barcodes")
#   
#   # Unimodal UMAP Projection ------------------------------------------------
#   # ref_WM <- RunUMAP(object = ref_WM,reduction = "harmony",dims = 1:10, return.model = TRUE)
#   query_transfer <- MapQuery(anchorset = test.anchors_0,
#                              reference = ref_WM,
#                              query = query_transfer,
#                              refdata = list(celltype = "seurat_clusters"),
#                              reference.reduction = "pca", reduction.model = "umap")
#   
#   # add the meta to the coordinates
#   data2_transfer <- left_join(
#     # get the coordinates of the query dataset
#     query_transfer@reductions$ref.umap@cell.embeddings %>%
#       data.frame() %>%
#       rownames_to_column(var = "barcodes"),
#     
#     # get the metadata from the query dataset, after adding the prediciton information. use a threshold for the score
#     query_transfer@meta.data %>%
#       rownames_to_column(var = "barcodes") %>%
#       # this is based on martina's table, I need to ask her how she provided the accuracy of the imputation
#       mutate(robust_score = case_when(prediction.score.max>0.70~predicted.id,
#                                       T ~ "uncertain")),
#     suffix = c(".WM",".BS"),
#     "barcodes")
#   
#   # save the metadata with the new clsuter informations
#   return(data2_transfer)
# })
# 
# # ggplot(label= TRUE) +
# #   # reference layer
# #   geom_point(data = ref_WM_meta %>% dplyr::select(UMAP_1,UMAP_2),aes(x=UMAP_1,y=UMAP_2),col="gray",alpha=0.3,size=0.1)+
# #   # transfer layer
# #   geom_point(data = data2_transfer,aes(x = refUMAP_1.WM,y = refUMAP_2.WM, col = robust_score),size=0.3,alpha=0.8) +
# #   # labs(color= "Clusters") +
# #   theme_bw() +
# #   facet_wrap(~ID)+
# #   guides(colour = guide_legend(override.aes = list(size=5),ncol = 1))+theme(strip.background = element_blank())
# 
# saveRDS(list_label_transfer,"../../out/object/list_out_label_transfer_02_60_IMM_WM.rds")
