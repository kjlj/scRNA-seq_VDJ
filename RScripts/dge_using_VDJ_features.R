#set working directory
setwd("~/Documents/GitHub/scRNA-seq_VDJ/")

#packages
library(Seurat)
library(tidyverse)
library(ggsci)
library(ggrepel)


#themes for plots
theme_custom <- function (base_size = 12) { 
  theme_bw(base_size=base_size) %+replace% 
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour="black"),
          axis.text.x = element_text(colour="black", size=base_size),
          axis.text.y = element_text(colour="black", size = base_size, hjust = 1),
          axis.title = element_text(colour="black", face="bold", size = base_size),
          strip.background = element_blank(),
          strip.text = element_text(colour="black", face="bold", size=base_size),
          legend.text = element_text(colour="black", size=base_size)
    )
}

theme_x_rotated <- function (base_size = 12) { 
  theme_custom(base_size=base_size) %+replace% 
    theme(axis.text.x = element_text(colour="black", angle = 90, hjust = 1, vjust = 0.5, size=base_size))
}

#loading seurat object
tumour.seurat.obj <- readRDS("seurat_objects/tumour_seurat-without-VDJ-genes-azimuth.rds")

#loading the seurat exported data
seurat.data <- read_tsv("data/pbmc-tumour_seruat_data.tsv")

#add the metadata to seurat object
metadata <- as.data.frame(seurat.data %>% filter(dataset == "tumour"))
rownames(metadata) <- metadata$barcode
metadata <- metadata %>% select(-barcode)

tumour.seurat.obj <- AddMetaData(tumour.seurat.obj, metadata)

#umap of T subsets cells by Azimuth annotation
azimuth_levels <- c(unique(seurat.data$azimuth_l2)[grepl("^CD[4|8]", unique(seurat.data$azimuth_l2))], "MAIT", "Treg")

seurat.data <- seurat.data %>% 
  mutate(cell_type = ifelse(azimuth_l2 %in% azimuth_levels, paste(azimuth_l2), "other"))

p1 <- ggplot(seurat.data, 
       aes(x = umap_1, y = umap_2)) +
  geom_point(data = seurat.data %>% filter(vdj_cell_type != "T"), 
             size = 0.5, colour = "lightgray") +
  geom_point(data = seurat.data %>% filter(vdj_cell_type == "T"),
             size = 0.8,
             aes(colour = cell_type)) +
  scale_colour_d3(palette = "category20",
                  name = "cell type") +
  theme_custom() +
  guides(colour = guide_legend(override.aes = list(size = 5)))
print(p1)

#can use the FindMarkers for comparing gene expression between groups of cells

cd4naive.vs.cd4tcm.markers <- FindMarkers(tumour.seurat.obj,  
                                          group.by = "azimuth_l2", 
                                          ident.1 = "CD4 Naive", ident.2 = "CD4 TCM") %>% 
  rownames_to_column("gene") %>%
  arrange(desc(avg_log2FC))

cd4naive.vs.cd4tcm.markers %>% 
  filter(p_val_adj < 0.05) %>% 
  head()

#show as a volcano plot
p1 <- ggplot(cd4naive.vs.cd4tcm.markers, 
       aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(colour = ifelse(p_val_adj < 0.05, "sig", "ns"))) +
  geom_text_repel(aes(label = ifelse(p_val_adj < 0.05, gene, NA))) +
  theme_custom() +
  scale_colour_manual(values = c("firebrick", "darkgray"),
                      name = "",
                      limits = c("sig", "ns")) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  labs(x = "average log2 fold change",
       y = "-log10(adjusted p-value)") 
plot(p1)  


#can compare expanded versus singleton clones
texp.vs.tsingle.markers <- FindMarkers(tumour.seurat.obj,  
                                          group.by = "t_expanded_cln", 
                                          ident.1 = "expanded", ident.2 = "singleton") %>% 
  rownames_to_column("gene") %>%
  arrange(desc(avg_log2FC))

texp.vs.tsingle.markers %>% head()

texp.vs.tsingle.markers %>% 
  filter(p_val_adj < 0.05) %>% 
  head()

#show as a volcano plot
p1 <- ggplot(texp.vs.tsingle.markers, 
       aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(colour = ifelse(p_val_adj < 0.05, "sig", "ns"))) +
  geom_text_repel(aes(label = ifelse(p_val_adj < 0.05, gene, NA))) +
  theme_custom() +
  scale_colour_manual(values = c("firebrick", "darkgray"),
                      name = "",
                      limits = c("sig", "ns")) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  labs(x = "average log2 fold change",
       y = "-log10(adjusted p-value)")
plot(p1)


#the largest T cell clone in the tumour is only 7 cells, there are instances where you want to compare a single clone to other clone(s)
no1.vs.no2.markers <- FindMarkers(tumour.seurat.obj,  
                                       group.by = "t_cln_rnk", 
                                       ident.1 = "1", ident.2 = "2") %>% 
  rownames_to_column("gene") %>%
  arrange(desc(avg_log2FC))

no1.vs.no2.markers %>% head()

no1.vs.no2.markers %>% 
  filter(p_val_adj < 0.05) %>% 
  head()

#nothing sig, expanding to demonstrate that vectors can be used to group
## this will compare the clones ranked 1 and 2 with those ranked 3 - 5
grp1.vs.grp2.markers <- FindMarkers(tumour.seurat.obj,  
                                  group.by = "t_cln_rnk", 
                                  ident.1 = c("1", "2"), ident.2 = c("3", "5")) %>% 
  rownames_to_column("gene") %>%
  arrange(desc(avg_log2FC))

grp1.vs.grp2.markers %>% head()

grp1.vs.grp2.markers %>% 
  filter(p_val_adj < 0.05) %>% 
  head()

#still nothing significant in this case

#using VH genes 
vgene1.vs.vgene2.markers <- FindMarkers(tumour.seurat.obj,  
                                    group.by = "v_igh", 
                                    ident.1 = c("IGHV1-69"), ident.2 = c("IGHV3-23")) %>% 
  rownames_to_column("gene") %>%
  arrange(desc(avg_log2FC))

vgene1.vs.vgene2.markers %>% head()

vgene1.vs.vgene2.markers %>% 
  filter(p_val_adj < 0.05) %>% 
  head()

p1 <- ggplot(vgene1.vs.vgene2.markers, 
             aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(colour = ifelse(p_val_adj < 0.05, "sig", "ns"))) +
  geom_text_repel(aes(label = ifelse(p_val_adj < 0.05, gene, NA))) +
  theme_custom() +
  scale_colour_manual(values = c("firebrick", "darkgray"),
                      name = "",
                      limits = c("sig", "ns")) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  labs(x = "average log2 fold change",
       y = "-log10(adjusted p-value)") 
plot(p1)  

#IgM and FCER2 more common among IGHV1-69 utilising B cells