# libraries ---------------------------------------------------------------
library(tidyverse)
library(ggrepel)

# read in the data --------------------------------------------------------
folder <-  "../../out/table/modules_SENESCENCE/"
file <- dir(folder) %>% 
  str_subset(pattern = "^Module_score_")

df_modules <- lapply(file, function(x){
  read_tsv(paste0(folder,x))
}) %>% 
  bind_rows()

# run it on cluster-wise BASELINE -----------------------------------------
id_signature <- unique(df_modules$signature)

lapply(id_signature,function(x){
  # keep track of the signatures
  print(x)
  # try to reference to top10% based on the control dataset
  df_90 <- df_modules %>% 
    filter(signature %in% x,
           # disease %in% "CTRL",
           treat == "BASELINE") %>% 
    group_by(seurat_clusters) %>% 
    summarise(thr = quantile(signature_score1,prob=0.90))
  
  # enumerate the cells 090 from the control only samples
  df_en_090_1 <- df_modules %>% 
    select(barcode,signature,treat,seurat_clusters,clone,doxy,exposure,ID,cell_id,signature_score1) %>% 
    filter(signature %in% x) %>% 
    left_join(df_90,by = "seurat_clusters") %>% 
    mutate(pass_thr = signature_score1>thr) %>% 
    group_by(seurat_clusters,treat) %>% 
    summarise(tot_cell = n(),
              tot_sen = sum(pass_thr)) %>%
    ungroup() %>% 
    mutate(prop_sen = tot_sen/tot_cell)
  
  # count the pro of the contols
  df_en_090_2 <- df_en_090_1 %>% 
    filter(treat %in% "BASELINE") %>% 
    select(seurat_clusters,prop_sen_ctrl=prop_sen)
  
  df_en_090_final <- df_en_090_1 %>% 
    left_join(df_en_090_2,by = "seurat_clusters") %>% 
    mutate(FC = prop_sen/prop_sen_ctrl,
           log2FC = log2(FC))
  
  df_en_090_final %>%
    write_tsv(paste0("../../out/table/modules_SENESCENCE/enumeration_",x,"_090_threshold_clusterBASELINE.tsv"))
  
  # read_tsv("../../out/table/ManualClean/modules_SENESCENCE/HVG4000/enumeration_senmayo_090_threshold_MSStatus_refWM.tsv")
  
  # plot using the 0.90 quantile
  df_modules %>% 
    filter(signature %in% x) %>% 
    # mutate(BraakStage=as.factor(BraakStage)) %>% 
    ggplot(aes(x=signature_score1,y=treat,fill=treat))+
    ggridges::geom_density_ridges(alpha=0.5)+
    facet_wrap(~seurat_clusters,scales = "free")+
    theme_bw()+
    theme(strip.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA))+
    # scale_fill_manual(values = c("green","red"))+
    geom_vline(data = df_90,aes(xintercept = thr),linetype="dashed",col="red")
  ggsave(paste0("../../out/image/modules_SENESCENCE/dist_score_ridges_",x,"_090_threshold_clusterBASELINE.pdf"),width = 12,height =9)
  
  # plot the enumeration
  df_en_090_final %>% 
    filter(tot_sen>5) %>% 
    ggplot(aes(y=treat,x=log2FC))+
    geom_col(aes(fill=tot_sen))+
    # geom_col(aes(width = tot_sen))+
    facet_wrap(~seurat_clusters)+
    theme_bw()+
    theme(strip.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA))+
    scale_fill_gradientn(colours = viridis::turbo(10),trans = "log",breaks = 4^seq(1,5))+geom_vline(xintercept = 0,linetype="dashed",col="gray")
  ggsave(paste0("../../out/image/modules_SENESCENCE/enumeration_log2Prop_",x,"_090_threshold_above5_clusterBASELINE.pdf"),width = 12,height =9)
})

# run it on cluster-wise CSF.ctrl -----------------------------------------
# x <- "senmayo"
id_signature <- unique(df_modules$signature)

lapply(id_signature,function(x){
  # keep track of the signatures
  print(x)
  # try to reference to top10% based on the control dataset
  df_90 <- df_modules %>% 
    filter(signature %in% x,
           # disease %in% "CTRL",
           treat == "CSF.ctrl") %>% 
    group_by(seurat_clusters) %>% 
    summarise(thr = quantile(signature_score1,prob=0.90))
  
  # enumerate the cells 090 from the control only samples
  df_en_090_1 <- df_modules %>% 
    select(barcode,signature,treat,seurat_clusters,clone,doxy,exposure,ID,cell_id,signature_score1) %>% 
    filter(signature %in% x) %>% 
    left_join(df_90,by = "seurat_clusters") %>% 
    mutate(pass_thr = signature_score1>thr) %>% 
    group_by(seurat_clusters,treat) %>% 
    summarise(tot_cell = n(),
              tot_sen = sum(pass_thr)) %>%
    ungroup() %>% 
    mutate(prop_sen = tot_sen/tot_cell)
  
  # count the pro of the contols
  df_en_090_2 <- df_en_090_1 %>% 
    filter(treat %in% "CSF.ctrl") %>% 
    select(seurat_clusters,prop_sen_ctrl=prop_sen)
  
  df_en_090_final <- df_en_090_1 %>% 
    left_join(df_en_090_2,by = "seurat_clusters") %>% 
    mutate(FC = prop_sen/prop_sen_ctrl,
           log2FC = log2(FC))
  
  df_en_090_final %>%
    write_tsv(paste0("../../out/table/modules_SENESCENCE/enumeration_",x,"_090_threshold_clusterCSFctrl.tsv"))
  
  # read_tsv("../../out/table/ManualClean/modules_SENESCENCE/HVG4000/enumeration_senmayo_090_threshold_MSStatus_refWM.tsv")
  
  # plot using the 0.90 quantile
  df_modules %>% 
    filter(signature %in% x) %>% 
    # mutate(BraakStage=as.factor(BraakStage)) %>% 
    ggplot(aes(x=signature_score1,y=treat,fill=treat))+
    ggridges::geom_density_ridges(alpha=0.5)+
    facet_wrap(~seurat_clusters,scales = "free")+
    theme_bw()+
    theme(strip.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA))+
    # scale_fill_manual(values = c("green","red"))+
    geom_vline(data = df_90,aes(xintercept = thr),linetype="dashed",col="red")
  ggsave(paste0("../../out/image/modules_SENESCENCE/dist_score_ridges_",x,"_090_threshold_clusterCSFctrl.pdf"),width = 12,height =9)
  
  # plot the enumeration
  df_en_090_final %>% 
    filter(tot_sen>5) %>% 
    ggplot(aes(y=treat,x=log2FC))+
    geom_col(aes(fill=tot_sen))+
    # geom_col(aes(width = tot_sen))+
    facet_wrap(~seurat_clusters)+
    theme_bw()+
    theme(strip.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA))+
    scale_fill_gradientn(colours = viridis::turbo(10),trans = "log",breaks = 4^seq(1,5))+geom_vline(xintercept = 0,linetype="dashed",col="gray")
  ggsave(paste0("../../out/image/modules_SENESCENCE/enumeration_log2Prop_",x,"_090_threshold_above5_clusterCSFctrl.pdf"),width = 12,height =9)
})

# run it on cell_id BASELINE ----------------------------------------------
id_signature <- unique(df_modules$signature)
# x <- "senmayo"

lapply(id_signature,function(x){
  # keep track of the signatures
  print(x)
  # try to reference to top10% based on the control dataset
  df_90 <- df_modules %>% 
    filter(signature %in% x,
           # disease %in% "CTRL",
           treat == "BASELINE") %>% 
    group_by(cell_id) %>% 
    summarise(thr = quantile(signature_score1,prob=0.90))
  
  # enumerate the cells 090 from the control only samples
  df_en_090_1 <- df_modules %>% 
    select(barcode,signature,treat,cell_id,clone,doxy,exposure,ID,cell_id,signature_score1) %>% 
    filter(signature %in% x) %>% 
    left_join(df_90,by = "cell_id") %>% 
    mutate(pass_thr = signature_score1>thr) %>% 
    group_by(cell_id,treat) %>% 
    summarise(tot_cell = n(),
              tot_sen = sum(pass_thr)) %>%
    ungroup() %>% 
    mutate(prop_sen = tot_sen/tot_cell)
  
  # count the pro of the contols
  df_en_090_2 <- df_en_090_1 %>% 
    filter(treat %in% "BASELINE") %>% 
    select(cell_id,prop_sen_ctrl=prop_sen)
  
  df_en_090_final <- df_en_090_1 %>% 
    left_join(df_en_090_2,by = "cell_id") %>% 
    mutate(FC = prop_sen/prop_sen_ctrl,
           log2FC = log2(FC))
  
  df_en_090_final %>%
    write_tsv(paste0("../../out/table/modules_SENESCENCE/enumeration_",x,"_090_threshold_cellidBASELINE.tsv"))
  
  # read_tsv("../../out/table/ManualClean/modules_SENESCENCE/HVG4000/enumeration_senmayo_090_threshold_MSStatus_refWM.tsv")
  
  # plot using the 0.90 quantile
  df_modules %>% 
    filter(signature %in% x) %>% 
    # mutate(BraakStage=as.factor(BraakStage)) %>% 
    ggplot(aes(x=signature_score1,y=treat,fill=treat))+
    ggridges::geom_density_ridges(alpha=0.5)+
    facet_wrap(~cell_id,scales = "free")+
    theme_bw()+
    theme(strip.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA))+
    # scale_fill_manual(values = c("green","red"))+
    geom_vline(data = df_90,aes(xintercept = thr),linetype="dashed",col="red")
  ggsave(paste0("../../out/image/modules_SENESCENCE/dist_score_ridges_",x,"_090_threshold_cellidBASELINE.pdf"),width = 9,height =6)
  
  # plot the enumeration
  df_en_090_final %>% 
    filter(tot_sen>5) %>% 
    ggplot(aes(y=treat,x=log2FC))+
    geom_col(aes(fill=tot_sen))+
    # geom_col(aes(width = tot_sen))+
    facet_wrap(~cell_id)+
    theme_bw()+
    theme(strip.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA))+
    scale_fill_gradientn(colours = viridis::turbo(10),trans = "log",breaks = 4^seq(1,5))+geom_vline(xintercept = 0,linetype="dashed",col="gray")
  ggsave(paste0("../../out/image/modules_SENESCENCE/enumeration_log2Prop_",x,"_090_threshold_above5_cellidBASELINE.pdf"),width = 9,height =6)
})

# run it on cell_id CSF.ctrl ----------------------------------------------
id_signature <- unique(df_modules$signature)
# x <- "senmayo"

lapply(id_signature,function(x){
  # keep track of the signatures
  print(x)
  # try to reference to top10% based on the control dataset
  df_90 <- df_modules %>% 
    filter(signature %in% x,
           # disease %in% "CTRL",
           treat == "CSF.ctrl") %>% 
    group_by(cell_id) %>% 
    summarise(thr = quantile(signature_score1,prob=0.90))
  
  # enumerate the cells 090 from the control only samples
  df_en_090_1 <- df_modules %>% 
    select(barcode,signature,treat,cell_id,clone,doxy,exposure,ID,cell_id,signature_score1) %>% 
    filter(signature %in% x) %>% 
    left_join(df_90,by = "cell_id") %>% 
    mutate(pass_thr = signature_score1>thr) %>% 
    group_by(cell_id,treat) %>% 
    summarise(tot_cell = n(),
              tot_sen = sum(pass_thr)) %>%
    ungroup() %>% 
    mutate(prop_sen = tot_sen/tot_cell)
  
  # count the pro of the contols
  df_en_090_2 <- df_en_090_1 %>% 
    filter(treat %in% "CSF.ctrl") %>% 
    select(cell_id,prop_sen_ctrl=prop_sen)
  
  df_en_090_final <- df_en_090_1 %>% 
    left_join(df_en_090_2,by = "cell_id") %>% 
    mutate(FC = prop_sen/prop_sen_ctrl,
           log2FC = log2(FC))
  
  df_en_090_final %>%
    write_tsv(paste0("../../out/table/modules_SENESCENCE/enumeration_",x,"_090_threshold_cellidCSFctrl.tsv"))
  
  # read_tsv("../../out/table/ManualClean/modules_SENESCENCE/HVG4000/enumeration_senmayo_090_threshold_MSStatus_refWM.tsv")
  
  # plot using the 0.90 quantile
  df_modules %>% 
    filter(signature %in% x) %>% 
    # mutate(BraakStage=as.factor(BraakStage)) %>% 
    ggplot(aes(x=signature_score1,y=treat,fill=treat))+
    ggridges::geom_density_ridges(alpha=0.5)+
    facet_wrap(~cell_id,scales = "free")+
    theme_bw()+
    theme(strip.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA))+
    # scale_fill_manual(values = c("green","red"))+
    geom_vline(data = df_90,aes(xintercept = thr),linetype="dashed",col="red")
  ggsave(paste0("../../out/image/modules_SENESCENCE/dist_score_ridges_",x,"_090_threshold_cellidCSFctrl.pdf"),width = 9,height =6)
  
  # plot the enumeration
  df_en_090_final %>% 
    filter(tot_sen>5) %>% 
    ggplot(aes(y=treat,x=log2FC))+
    geom_col(aes(fill=tot_sen))+
    # geom_col(aes(width = tot_sen))+
    facet_wrap(~cell_id)+
    theme_bw()+
    theme(strip.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA))+
    scale_fill_gradientn(colours = viridis::turbo(10),trans = "log",breaks = 4^seq(1,5))+geom_vline(xintercept = 0,linetype="dashed",col="gray")
  ggsave(paste0("../../out/image/modules_SENESCENCE/enumeration_log2Prop_",x,"_090_threshold_above5_cellidCSFcytrl.pdf"),width = 9,height =6)
})

