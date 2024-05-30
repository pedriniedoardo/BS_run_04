# # libraries ---------------------------------------------------------------
# library(tidyverse)
# library(Seurat)
# library(SeuratData)
# library(ggridges)
# library(ComplexHeatmap)
# 
# # read in the data --------------------------------------------------------
# # in this case I want to use the Leng 2021 dataset as reference to map the cells of martina's dataset
# # read in the reference dataset for the leng dataset
# ref_BS <- readRDS("/beegfs/scratch/ric.cosr/pedrini.edoardo/scRNAseq_Brainsphere_Absinta/BS_drop_RR25CTRL/data/data.combined_NOT_annotated_norm_fix_DoubletSinglet_harmonyMartina.rds")
# # add the model to the obejct if it is missing
# ref_BS <- RunUMAP(object = ref_BS,reduction = "harmony",dims = 1:10, return.model = TRUE)
# 
# # save ccordinates fro reference of the UMAP and metadata
# ref_BS_meta <- left_join(ref_BS@reductions$umap@cell.embeddings %>%
#             data.frame() %>%
#             rownames_to_column("barcodes"),
#           ref_BS@meta.data %>%
#             data.frame() %>%
#             rownames_to_column("barcodes"),by="barcodes")
# 
# write_tsv(ref_BS_meta,file = "../../out/table/ref_BS_meta.tsv")
# 
# # DimPlot(ref_BS,label = T)
# 
# # loop the label transfer procedure ---------------------------------------
# folder <- "../../out/object/"
# file <- dir(folder) %>%
#   str_subset(pattern = "datasc_fix_filter_norm_doublet_Singlet")
# 
# file
# # x <- file[1]
# 
# list_out <- lapply(file,function(x){
#   print(x)
#   # read in the query dataset
#   query <- readRDS(paste0(folder,x))
#   # fix the meta to add the general cell annotation
#   query$seurat_clusters_fix1 <- paste0("cluster_",query$seurat_clusters)
#   
#   # -------------------------------------------------------------------------
#   # to work the reference should be the integrated dataset
#   DefaultAssay(ref_BS)<-"integrated"
#   # the RNA slot do not work due to the missing variable features
#   # DefaultAssay(ref_BS)<-"RNA"
#   
#   # in the vignetted the query isn the RNA slot
#   # DefaultAssay(query) <- "integrated"
#   DefaultAssay(query) <- "RNA"
#   test.anchors_0 <- FindTransferAnchors(reference = ref_BS,
#                                         query = query,
#                                         dims = 1:30,
#                                         reference.reduction = "pca")
#   
#   predictions_0 <- TransferData(anchorset = test.anchors_0,
#                                 refdata = ref_BS$seurat_clusters,
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
#   # ref_BS <- RunUMAP(object = ref_BS,reduction = "harmony",dims = 1:10, return.model = TRUE)
#   query_transfer <- MapQuery(anchorset = test.anchors_0,
#                              reference = ref_BS,
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
#                                       T ~ "uncertain"))
#     ,"barcodes")
#   
#   # save the metadata with the new clsuter informations
#   return(data2_transfer)
# }) %>% 
#   setNames(file %>% 
#              str_remove_all(pattern = "datasc_fix_filter_norm_doublet_Singlet_|.rds"))
# 
# saveRDS(list_out,"../../out/object/list_out_label_transfer_BS0102.rds")
