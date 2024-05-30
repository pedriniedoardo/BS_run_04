test <- readRDS(file = "../../out/object/test.rds")

test2 <- test %>%
  Seurat::NormalizeData(verbose = T) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 4000) %>%
  # I can scale the missing features afterwards now focus on the highly variable one for speed purposes
  # ScaleData(vars.to.regress = c("percent.mt","nCount_RNA","S.Score","G2M.Score"), verbose = T) %>% 
  ScaleData(verbose = T) %>%
  # run this if you want to scale all the variables
  # ScaleData(vars.to.regress = c("percent.mt.harmony","nCount_RNA.harmony","S.Score.harmony","G2M.Score.harmony"), verbose = T,features = all.genes) %>% 
  RunPCA(npcs = 30, verbose = T) %>% 
  RunUMAP(reduction = "pca", dims = 1:30,return.model = TRUE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.5) %>%
  identity()
