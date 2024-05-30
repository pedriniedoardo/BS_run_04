# AIM ---------------------------------------------------------------------
# build a referece dataset out of a costum analysis of another rds objec.
# In this case I want to build a reference out the of the Run 03 dataset
# notice it requires a big amount of RAM to fully run.
# olso to avoid issue I had to install SeuratObject v 5.0.1

# libraries ---------------------------------------------------------------
# library(Matrix)
library(Seurat)
library(Azimuth)
library(presto)
library(dplyr)
library(glmGamPoi)

# read in the object ------------------------------------------------------
# args <- commandArgs(trailingOnly = TRUE)

# read in the object for which I want to build a reference
ref <- readRDS("/beegfs/scratch/ric.cosr/pedrini.edoardo/scRNAseq_Brainsphere_Absinta/BS_batch_2023/out/object/sobj_processed_donor.rds")
DimPlot(ref,raster = T,group.by = "expertAnno.l1",label = T)

# ref.dir <- "reference/"
# ob.dir <- "seurat_objects/"
# ref <- readRDS(file = args[1])
# annotations <- readRDS(file = args[2])
# Idents(object = ref) <- annotations
Idents(object = ref) <- "expertAnno.l1"

# if ("remove" %in% levels(x = ref)) {
#   ref <- subset(x = ref, idents = "remove", invert = TRUE)
#   ref <- RunPCA(object = ref, verbose = FALSE)
# }
# ref$annotation.l1 <- Idents(object = ref)
# needs to have the model
ref_SCT <- SCTransform(ref, method = "glmGamPoi", vars.to.regress = c("percent.mt","nCount_RNA"), verbose = T)
ref2 <- RunUMAP(object = ref_SCT,reduction = "harmony",dims = 1:10, return.model = TRUE)

DimPlot(ref2,raster = T)
# ref2@reductions$umap@misc$model
full.ref <- ref2

full.ref$annotation.l1 <- Idents(object = full.ref)
colormap <- list(annotation.l1 = CreateColorMap(object = ref2, seed = 2))
colormap[["annotation.l1"]] <- colormap[["annotation.l1"]][sort(x = names(x = colormap[["annotation.l1"]]))]

ref_final <- AzimuthReference(
  object = full.ref,
  refUMAP = "umap",
  refDR = "pca",
  refAssay = "SCT",
  metadata = c("annotation.l1"),
  dims = 1:30,
  k.param = 31,
  colormap = colormap,
  reference.version = "1.0.0"
)

DimPlot(ref_final,raster=T)

SaveAnnoyIndex(object = ref_final[["refdr.annoy.neighbors"]], file = file.path("../../data/ref_BS_run_03/","idx.annoy"))
saveRDS(object = ref_final, file = file.path("../../data/ref_BS_run_03/", "ref.Rds"))
saveRDS(object = full.ref, file = file.path("../../data/ref_BS_run_03/", "fullref.Rds"))

# test reading azimuth reference ------------------------------------------
# try to load the new reference object created
reference <- LoadReference(path = "../../data/ref_BS_run_03/")
DimPlot(reference$plot,group.by = "annotation.l1",label = T)
