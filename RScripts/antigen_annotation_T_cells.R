#loading R packages
#tidyverse packages for data manipulation and plotting
library(tidyverse)
#in order to work with the seurat objects from the scRNA-seq
library(Seurat)
#multi-panel plots
library(ggpubr)
#colour scales for scientific journals
library(ggsci)
#stats package for ggplot2
library(rstatix)

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

#loading datasets

## seurat dump that includes the UMAP coords & meta with VDJ
seurat.data <- read_tsv("data/pbmc-tumour_seruat_data.tsv")

## annotated TCRs
trb.ags.data <- read_tsv("tcr_ag_refs/trb_antigens.tsv")


#we are focussing on TRB chain so will be joining the TRB clonotype from 10x with the TRB clonotype from the Ags
seurat.data <- seurat.data %>%
  mutate(trb_clonotype = ifelse(grepl("TRBV", paired_clonotype), str_replace(paired_clonotype, ":.+$", ""), NA_character_)) %>%
  left_join(trb.ags.data, by = c("trb_clonotype" = "clonotype")) %>%
  mutate(tcr_antigen_class = case_when(is.na(antigen) ~ "none",
                                       grepl(",", antigen) ~ "multiple ags",
                                       TRUE ~ paste0(antigen)))

#summarise the number of T cells in the dataset that we annotationed for Ag specificity
ag.summary <- seurat.data %>%
  filter(vdj_cell_type == "T") %>%
  group_by(dataset, tcr_antigen_class) %>%
  summarise(cell_cnt = length(barcode)) %>%
  ungroup() %>%
  group_by(dataset) %>%
  mutate(cell_prop = cell_cnt / sum(cell_cnt),
         cell_percent = cell_prop * 100) %>%
  ungroup()

ag.summary

#plotting onto the UMAP
ag.umap.p1 <- ggplot(seurat.data, aes(x = umap_1, y = umap_2)) +
  geom_point(data = seurat.data %>% filter(vdj_cell_type != "T"), 
             size = 0.5, 
             colour = "lightgray") +
  geom_point(data = seurat.data %>% filter(vdj_cell_type == "T", tcr_antigen_class == "none"), 
             olour = "darkgray", 
             size = 0.8) +
  geom_point(data = seurat.data %>% filter(vdj_cell_type == "T", tcr_antigen_class != "none"),
             aes(colour = tcr_antigen_class),
             size = 0.8) +
  scale_colour_futurama(name = "antigen") +
  theme_custom() +
  facet_grid(~ dataset) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  labs(x = "UMAP_1",
       y = "UMAP_2",
       caption = "unannotated T cells in dark gray\nnon-T cells in light gray")
plot(ag.umap.p1)
