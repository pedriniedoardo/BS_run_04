# AIM ---------------------------------------------------------------------
# add the donor deconvolution annotation to the barcodes

# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(scales)
library(ComplexHeatmap)
library(limma)
library(multtest)
library(metap)
library(ggbreak)
library(patchwork)
# library(lemon)
library(future)
library(ggnewscale)

#  read in the data -------------------------------------------------------
test <- readRDS("../../out/object/data.combined_harmonySkipIntegration_AllSoupX_01000_06000_15.rds")
DefaultAssay(test)
DimPlot(test,raster = T,label = T)

# # defirne the pattern in the barcodes
# pattern_remove <- table(test@meta.data$orig.ident) |> 
#   names() |>
#   paste0("_") |> 
#   paste0(collapse = "|")

# extract the full meta from the scobj
meta <- test@meta.data |> 
  rownames_to_column("full_barcode")

# define the table of the sample identity after donor devonvolution
# martina also asked to add another treat varaible that divide the CSF.MS 24 from the CSF.MS 48 treatmentnt
LUT_sample <- data.frame(id_deconvolution = c("W8_CSF_tolebrutinib","W8_untreated","W8_low_cytokines","W8_IgG_tolebrutinib","W8_CSF_dasatinib","W8_CSF_ibudilast","W8_CSF","W8_IgG","W8_CSF_alfa_lipoic_acid"),
                         id_sample_scRNAseq = c("CSF.MS_TOLE","BASELINE","CYTOlow","IgG_TOLE","CSF.MS_DAS","CSF.MS_IBU","CSF.MS","IgG","CSF.MS_ALA"),
                         id_sample_short = c("CSF.MS_TOLE","BASELINE","CYTOlow","IgG_TOLE","CSF.MS_DAS","CSF.MS_IBU","CSF.MS","IgG","CSF.MS_ALA"),
                         treat_full = c("CSF.MS_TOLE","BASELINE","CYTOlow","IgG_TOLE","CSF.MS_DAS","CSF.MS_IBU","CSF.MS","IgG","CSF.MS_ALA"))

# From the classification approach built the lut for the donor identity. Use a common donor ID
LUT_donor <- data.frame(id_donor_sample = c("donor2|CSF.MS_TOLE","donor2|BASELINE","donor0|CYTOlow","donor0|IgG_TOLE","donor2|CSF.MS_DAS","donor1|CSF.MS_IBU","donor1|CSF.MS","donor1|IgG","donor2|CSF.MS_ALA",
                                            "donor0|CSF.MS_TOLE","donor1|BASELINE","donor1|CYTOlow","donor1|IgG_TOLE","donor1|CSF.MS_DAS","donor2|CSF.MS_IBU","donor0|CSF.MS","donor2|IgG","donor0|CSF.MS_ALA",
                                            "donor1|CSF.MS_TOLE","donor0|BASELINE","donor2|CYTOlow","donor2|IgG_TOLE","donor0|CSF.MS_DAS","donor0|CSF.MS_IBU","donor2|CSF.MS","donor0|IgG","donor1|CSF.MS_ALA"),
                        harmonized_donor = c(rep("don01",9),
                                             rep("don02",9),
                                             rep("don03",9))) |> 
  separate(id_donor_sample,into = c("donor","id_sample_short"),sep = "\\|",remove = F)


dir("../../data/02_vireo_cellranger7/")

# pull the deconvolution results from vireo
# read in one run of genotyping
df_donorIds <- lapply(LUT_sample$id_deconvolution,function(x){
  file <- paste0("../../data/02_vireo_cellranger7/",x,"/donor_ids.tsv")
  read_tsv(file)
}) %>% 
  setNames(LUT_sample$id_sample_short) |> 
  bind_rows(.id = "id_sample_short") |> 
  dplyr::select(id_sample_short,cell,donor_id)

# merge the info
meta_full <- left_join(df_donorIds,LUT_sample,by="id_sample_short") |> 
  # build the join variable
  mutate(id_donor_sample = paste0(donor_id,"|",id_sample_short)) |> 
  # join the summary from the LUT_donor
  left_join(LUT_donor,by = c("id_donor_sample","id_sample_short")) |> 
  # build the barcode column
  mutate(full_barcode = paste0(id_deconvolution,"_",cell)) |> 
  # if the hamonyzed donor is missing pull the category from the donor_id imputation
  mutate(harmonized_donor2 = case_when(is.na(harmonized_donor)~donor_id,
                   T~harmonized_donor))
  # from the pure samples do not allow for more than one donor for the baseline
  # mutate(harmonized_donor2 = case_when(id_deconvolution == "hBS_CTR4_MG"~"donRR16",
  #                                      id_deconvolution == "hBS_RR16_MG"~"donRR16",
  #                                      id_deconvolution == "hBS_RR24_MG"~"donRR24",
  #                                      id_deconvolution == "hBS_RR25_MG"~"donRR25",
  #                                      T~harmonized_donor2))

# add the full meta to the sc meta table
meta_full2 <- meta |> 
  left_join(meta_full,by = c("full_barcode"))

# confirm the dimension of the metadata
dim(meta_full2)
dim(meta)

# swap the metadata in the sc object
test@meta.data <- meta_full2 |> 
  mutate(rowname = full_barcode) |> 
  column_to_rownames("rowname")

# count the number of cells per harmonized donors
meta_full2 |> 
  group_by(harmonized_donor2,orig.ident) |> 
  summarise(n = n()) |> 
  pivot_wider(names_from = harmonized_donor2,values_from = n)

# martina wanted to have an idea fo the relaitve proportion of doublet and unassigned per sample
meta_full2 |> 
  group_by(harmonized_donor2,orig.ident) |> 
  summarise(n = n()) |> 
  group_by(orig.ident) |> 
  mutate(tot_sample = sum(n)) |> 
  ungroup() |> 
  mutate(prop = n/tot_sample) |>
  ggplot(aes(x=orig.ident,y=prop,fill=harmonized_donor2))+
  geom_col()+
  theme_void()+
  theme(axis.text.x = element_text(hjust = 1,vjust = 1,angle = 45), plot.margin = margin(0.5, 0.5, 2, 2, "cm"),axis.text.y = element_text())
ggsave("../../out/plot/sobj_processed_donor_proportions.pdf",width = 8,height = 6)

# plot it as a scatter plot
meta_full2 |> 
  group_by(harmonized_donor2,orig.ident) |> 
  summarise(n = n()) |> 
  group_by(orig.ident) |> 
  mutate(tot_sample = sum(n)) |> 
  ungroup() |> 
  mutate(prop = n/tot_sample) |>
  ggplot(aes(x = orig.ident,
             y = prop,col = harmonized_donor2,
             group = harmonized_donor2))+
  geom_point() +
  geom_line() +
  theme_bw() +
  facet_grid(harmonized_donor2~"prop over sample") +
  theme(axis.text.x = element_text(hjust = 1,vjust = 1,angle = 45),
        plot.margin = margin(0.5, 0.5, 2, 2, "cm"),
        strip.background = element_blank())+
  # scale_y_sqrt(breaks = seq(from=0,to=1,by=0.1))
  scale_y_sqrt(breaks = c(0,0.001,0.02,0.1,0.2,0.4,0.6,0.8,1))
ggsave("../../out/plot/sobj_processed_donor_proportions_scatter.pdf",width = 8,height = 10)

# confirm the swap of the metadata
DimPlot(test,split.by = "id_sample_short",raster = T,group.by = "harmonized_donor2")+facet_wrap(~id_sample_short,nrow=3)
ggsave("../../out/plot/UMAPClusterSplit_DonorDeconvolution_id_sample_short.pdf",width = 13,height = 12)
DimPlot(test,split.by = "harmonized_donor2",raster = T,group.by = "id_sample_short")+facet_wrap(~harmonized_donor2,nrow=2)
ggsave("../../out/plot/UMAPClusterSplit_DonorDeconvolution_harmonized_donor2.pdf",width = 16,height = 9)
DimPlot(test,split.by = "harmonized_donor2",raster = T,group.by = "treat_full")+facet_wrap(~harmonized_donor2,nrow=2)

# add also the attempted cell ID annotation -------------------------------
# add the full meta to the sc meta table
meta_full2

# swap the metadata in the sc object
LUT_annotation <- data.frame(seurat_clusters = unique(meta_full2$seurat_clusters)) |> 
  mutate(expertAnno.l1 = case_when(seurat_clusters %in% c(1,3,4,5,6,7,12,13)~"NEU",
                                   seurat_clusters %in% c(2,17,19,21,22)~"MG",
                                   seurat_clusters %in% c(8)~"OLIGOLINEAGE",
                                   seurat_clusters %in% c(9,18)~"PROG",
                                   seurat_clusters %in% c(10,11,14)~"GLIA_IMM",
                                   seurat_clusters %in% c(0,10,15,16,20)~"ASTRO"))

test@meta.data <- meta_full2 %>% 
  mutate(rowname = full_barcode) %>% 
  left_join(LUT_annotation,by = "seurat_clusters") %>%
  column_to_rownames("rowname")

# check the groping
DimPlot(test,group.by = "expertAnno.l1",raster = T,label = T)

# -------------------------------------------------------------------------
# save the object with the updated matadata
saveRDS(test,"../../out/object/sobj_processed_donor.rds")
