#setting working directory
#setwd('~/Documents/GitHub/scRNA-seq_VDJ/')

#installing R packages
#check if the package is already installed, if not, install it
if (sum(grepl("tidyverse", rownames(installed.packages()))) > 0) {
  print("tidyverse is already installed")
} else {
  install.packages("tidyverse")
}

if (sum(grepl("Seurat", rownames(installed.packages()))) > 0) {
  print("Seurat is already installed")
} else {
  install.packages("Seurat")
}

if (sum(grepl("ggpubr", rownames(installed.packages()))) > 0) {
  print("ggpubr is already installed")
} else {
  install.packages("ggpubr")
}

if (sum(grepl("ggsci", rownames(installed.packages()))) > 0) {
  print("ggsci is already installed")
} else {
  install.packages("ggsci")
}

if (sum(grepl("rstatix", rownames(installed.packages()))) > 0) {
  print("rstatix is already installed")
} else {
  install.packages("rstatix")
}

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


#loading VDJ data
## loading form git hub
#b.data <- read_tsv("https://raw.githubusercontent.com/kjlj/scRNA-seq_VDJ/main/data/pbmc-tumour_Ig_cellranger_igblast_per-barcode_clns.tsv")
#t.data <- read_tsv("https://raw.githubusercontent.com/kjlj/scRNA-seq_VDJ/main/data/pbmc-tumour_TR_cellranger_igblast_per-barcode.tsv")
## loading from local
b.data <- read_tsv("data/pbmc-tumour_Ig_cellranger_igblast_per-barcode_clns.tsv")
t.data <- read_tsv("data/pbmc-tumour_TR_cellranger_igblast_per-barcode.tsv")


#loading seurat objects
#loading from dropbox as large file need to extend the timeout
#options(timeout=600)
#pbmc.seurat.obj <- readRDS(url("https://www.dropbox.com/scl/fi/4s610vt1mgtmgibdvfsar/pbmc_seurat-without-VDJ-genes-azimuth.rds?rlkey=ftdkxi9mnxezhbb42dqftbqel&dl=1"))
#tumour.seurat.obj <- readRDS(url("https://www.dropbox.com/scl/fi/scik47zay4x27t4wxmo70/tumour_seurat-without-VDJ-genes-azimuth.rds?rlkey=z8ghoeoboaneniv82e2xywji3&dl=1"))

#loading from local
pbmc.seurat.obj <- readRDS("data/pbmc_seurat-without-VDJ-genes-azimuth.rds")
tumour.seurat.obj <- readRDS("data/tumour_seurat-without-VDJ-genes-azimuth.rds")


#up to now we have treated the B and T VDJ data within the same dataset independently but we are now going to unite them with the scRNA-seq data
#shouldn't have cells with both B and T VDJ for the same cell barcode but doublets (which haven't been filtered from the VDJ) can result in this artifact
#check for any cell barcodes that have both B and T cell data, likely to represent doublets and should be filtered
dual.bt.summary <- bind_rows(b.data %>% 
                               select(barcode, dataset) %>% 
                               mutate(chain = "B"),
                             t.data %>%
                               select(barcode, dataset) %>%
                               mutate(chain = "T")) %>%
  pivot_wider(id_cols = c("barcode", "dataset"),
              names_from = "chain",
              values_from = "chain",
              values_fn = length) %>%
  replace(is.na(.), 0) %>%
  mutate(dual_cell = ifelse(B == 1 & T == 1, TRUE, FALSE))

#get count for the number of cells with both B and T cell data from each dataset
with(dual.bt.summary, table(dataset, dual_cell))

#add the dual_cell column to the b and t data and filter the cells that are TRUE
b.data <- b.data %>%
  left_join(dual.bt.summary %>% select(-c(T, B)), by = c("barcode", "dataset")) %>%
  filter(dual_cell == FALSE)

table(b.data$dataset)

t.data <- t.data %>%
  left_join(dual.bt.summary %>% select(-c(T, B)), by = c("barcode", "dataset")) %>%
  filter(dual_cell == FALSE)

table(t.data$dataset)

rm(dual.bt.summary)

##there may also be cell barcodes for which we have VDJ data but the cell barcode was removed from GEX during filtering
##need to harmonise the cellbarcodes between the VDJ and GEX, and drop VDJ cell barcodes that don't have data in the GEX
pbmc.barcodes <- tibble(barcode = rownames(pbmc.seurat.obj@meta.data)) %>%
  mutate(dataset = "pbmc",
         in_gex = TRUE)
tumour.barcodes <- tibble(barcode = rownames(tumour.seurat.obj@meta.data)) %>%
  mutate(dataset = "tumour",
         in_gex = TRUE)
barcodes <- bind_rows(pbmc.barcodes, tumour.barcodes)

b.data <- b.data %>%
  left_join(barcodes) %>%
  mutate(in_gex = ifelse(is.na(in_gex), FALSE, TRUE))

t.data <- t.data %>%
  left_join(barcodes) %>%
  mutate(in_gex = ifelse(is.na(in_gex), FALSE, TRUE))

with(t.data, table(dataset, in_gex))
with(b.data, table(dataset, in_gex))

rm(list = c("pbmc.barcodes", "tumour.barcodes", "barcodes"))

#there are some additional data fields that could be of interest for our investigations that we can add to the B and T cell datasets
#add extra data fields to the b cell data
b.data <- b.data %>%
  filter(in_gex) %>% #remove VDJ for cells that aren't in the GEX
  mutate(bcr_paired = ifelse(!(is.na(v_call_igblast_igh)) & !(is.na(v_call_igblast_igkl)), TRUE, FALSE)) %>% #add a column to indicate if the b cells include paired chains
  group_by(dataset, clone_id, unq_clone_id) %>% #group cells by their clone ids within each dataset in order to add columns for # cells in a clone and whether or not a clone is expanded
  mutate(cln_size = length(barcode)) %>% #add the size for each clone, using mutate rather than summarise
  ungroup() %>% #remove the grouping
  mutate(b_expanded_cln = ifelse(cln_size == 1, "singleton", "expanded")) %>% #add a column indicating if the clone is expanded based on the cell cnt for the clone
  mutate(vh_mut = ifelse(is.na(v_identity_igblast_igh), NA_integer_, 100 - v_identity_igblast_igh),
         vkl_mut = ifelse(is.na(v_identity_igblast_igkl), NA_integer_, 100 - v_identity_igblast_igkl),
         b_mut_class = ifelse(vh_mut > 2 | vkl_mut > 2, "mutated", "unmutated")) %>% #for the Ig data add the SHM & whether a cell's Ig is mutated (using 2% threshold)
  mutate(isotype = ifelse(is.na(c_call_igblast_igh), NA_character_, str_replace(c_call_igblast_igh, "\\*.+$", ""))) #for the Ig data add the isotype subclass

#number of CELLs in expanded clones
with(b.data, table(dataset, b_expanded_cln))
#number of cells that are mutated vs unmutated
with(b.data, table(dataset, b_mut_class))
#number of cells for each isotype subclass
with(b.data, table(dataset, isotype))

#add additional fields to the t cell data
t.data <- t.data %>%
  filter(in_gex) %>% #removing cells that aren't in the GEX
  mutate(tcr_paired = ifelse( !(is.na(clonotype_tra)) & !(is.na(clonotype_trb)), TRUE, FALSE)) %>% #add a column to indicate if the t cells include paired chains
  mutate(paired_clonotype = case_when(is.na(clonotype_tra) & !(is.na(clonotype_trb)) ~ paste(clonotype_trb),
                                      !(is.na(clonotype_tra)) & is.na(clonotype_trb) ~ paste(clonotype_tra),
                                      !(is.na(clonotype_tra)) & !(is.na(clonotype_trb)) ~ paste(clonotype_trb, clonotype_tra, sep = ":"),
                                      TRUE ~ NA_character_)) %>% #for the t cell add a clonotype label that joins both the TRA and TRB, the clone clustering for the Ig already considered the IGH and IGKL
  group_by(dataset, paired_clonotype) %>% #group cells by their paired clonotype labels to add cols for # cells in each clonotype and whether clonotype is expanded
  mutate(cln_size = length(barcode)) %>% #add the size for each clone, using mutate rather than summarise
  ungroup() %>% #remove the grouping
  mutate(t_expanded_cln = ifelse(cln_size == 1, "singleton", "expanded")) %>% #add a column indicating if the clone is expanded based on the cell cnt for the clone
  mutate(isotype = ifelse(is.na(c_call_igblast_trb), NA_character_, str_replace(c_call_igblast_trb, "\\*.+$", ""))) #adding the TRB constant region type, not actually isotype but for plotting same col name

#number of CELLs in expanded clones
with(t.data, table(dataset, t_expanded_cln))
with(t.data, table(dataset, isotype))

## adding a rank for the clones/clonotypes
b.cln.rnks <- b.data %>%
  filter(!(is.na(unq_clone_id))) %>%
  group_by(dataset, unq_clone_id) %>%
  summarise(cln_size = length(barcode)) %>%
  ungroup() %>%
  mutate(expanded = ifelse(cln_size > 1, TRUE, FALSE)) %>%
  filter(expanded) %>%
  group_by(dataset) %>%
  mutate(b_cln_rnk = min_rank(desc(cln_size))) %>% 
  ungroup() %>%
  arrange(b_cln_rnk)

b.data <- b.data %>%
  left_join(b.cln.rnks %>% select(dataset, unq_clone_id, b_cln_rnk),
            by = c("unq_clone_id", "dataset"))

rm(b.cln.rnks)

t.cln.rnks <- t.data %>%
  filter(!(is.na(paired_clonotype))) %>%
  group_by(dataset, paired_clonotype) %>%
  summarise(cln_size = length(barcode)) %>%
  ungroup() %>%
  mutate(expanded = ifelse(cln_size > 1, TRUE, FALSE)) %>%
  filter(expanded) %>%
  group_by(dataset) %>%
  mutate(t_cln_rnk = min_rank(desc(cln_size))) %>% 
  ungroup() %>%
  arrange(t_cln_rnk)

t.data <- t.data %>%
  left_join(t.cln.rnks %>% select(dataset, paired_clonotype, t_cln_rnk),
            by = c("paired_clonotype", "dataset"))

rm(t.cln.rnks)

##format the metadata to merge with the seurat objects for the scRNA-seq analysis from the GEX
## the join will happen via the barcode and here we are selecting the subset of columns that want to add into the seurat metadata
## this can be altered
b.metadata <- b.data %>%
  select(barcode, dataset, v_igh, j_igh, v_igkl, j_igkl, unq_clone_id, bcr_paired, b_cln_rnk, b_expanded_cln,
         vh_mut, vkl_mut, b_mut_class, isotype)

t.metadata <- t.data %>%
  select(barcode, dataset, v_tra, j_tra, v_trb, j_trb, paired_clonotype, tcr_paired, t_cln_rnk, t_expanded_cln, isotype)

##to add the metadata need to convert from a tibble to a data frame and move the barcode to the row names, also need to split the PBMC and tumour data

##PBMC
#create a data.frame from combining the B and T cell data for the PBMC dataset
#adding an extra column after bring the data together to indicate whether a cell has B VDJ, T VDJ or neither
pbmc.metadata <- as.data.frame(b.metadata %>% filter(dataset == "pbmc") %>%
                                 full_join(t.metadata %>% filter(dataset == "pbmc")) %>%
                                 mutate(vdj_cell_type = case_when(!(is.na(paired_clonotype)) ~ "T",
                                                                  !(is.na(unq_clone_id)) ~ "B",
                                                                  TRUE ~ NA_character_)))
#add the barcodes as the row names from the data.frame
rownames(pbmc.metadata) <- pbmc.metadata$barcode
#remove the barcode column
pbmc.metadata <- pbmc.metadata %>% select(-barcode)

#add the metadata to the seurat object, as the rownames are the cell barcodes this will ensure everything gets matched up!
pbmc.seurat.obj <- AddMetaData(pbmc.seurat.obj, pbmc.metadata)

#check the meta data cols now in the seurat object
colnames(pbmc.seurat.obj@meta.data)

rm(list = c("pbmc.metadata"))

##tumour
#create a data.frame from B and T data for the tumour dataset
tumour.metadata <- as.data.frame(b.metadata %>% filter(dataset == "tumour") %>%
                                   full_join(t.metadata %>% filter(dataset == "tumour")) %>%
                                   mutate(vdj_cell_type = case_when(!(is.na(paired_clonotype)) ~ "T",
                                                                    !(is.na(unq_clone_id)) ~ "B",
                                                                    TRUE ~ NA_character_)))
#add the barcodes as the row names from the data.frame
#copy the barcodes to the row names
rownames(tumour.metadata) <- tumour.metadata$barcode
#remove the barcode column
tumour.metadata <- tumour.metadata %>% select(-barcode)

#add the metadata to the seurat object, as the rownames are the cell barcodes this will ensure everything gets matched up!
tumour.seurat.obj <- AddMetaData(tumour.seurat.obj, tumour.metadata)

#check the meta data cols now in the seurat object
colnames(tumour.seurat.obj@meta.data)

rm(list = c("tumour.metadata", "b.metadata", "t.metadata"))

#plot the metadata onto the umap within Seurat
p1 <- DimPlot(pbmc.seurat.obj, group.by = "vdj_cell_type", cols = pal_d3()(2),  na.value = "lightgray")
plot(p1)

p1 <- DimPlot(pbmc.seurat.obj, group.by = "isotype", cols = pal_d3("category20")(10), na.value = "lightgray")
plot(p1)

p1 <- DimPlot(subset(pbmc.seurat.obj, vdj_cell_type == "B"), 
              group.by = "b_expanded_cln", 
              cols = pal_d3("category20")(10), 
              pt.size = 0.8,
              na.value = "lightgray")
plot(p1)

p1 <- DimPlot(subset(pbmc.seurat.obj, vdj_cell_type == "B"), 
              group.by = c("b_cln_rnk"), 
              cols = pal_d3("category20")(10), 
              pt.size = 0.8,
              na.value = "lightgray")
plot(p1)


p1 <- DimPlot(subset(pbmc.seurat.obj, vdj_cell_type == "T"), 
              group.by = "t_expanded_cln", 
              cols = pal_d3("category20")(10), 
              pt.size = 0.8,
              na.value = "lightgray")
plot(p1)

# p1 <- DimPlot(subset(pbmc.seurat.obj, vdj_cell_type == "T"), 
#                        group.by = "t_cln_rnk", 
#                        cols = pal_d3("category20")(10), 
#                        pt.size = 0.8,
#                        na.value = "lightgray")
# plot(p1)

#tumour
p1 <- DimPlot(tumour.seurat.obj, group.by = "vdj_cell_type", cols = pal_d3()(2),  na.value = "lightgray")
plot(p1)

p1 <- DimPlot(tumour.seurat.obj, group.by = "isotype", cols = pal_d3("category20")(20), na.value = "lightgray")
plot(p1)

# p1 <- DimPlot(subset(tumour.seurat.obj, vdj_cell_type == "B"), 
#                        group.by = "b_expanded_cln", 
#                        cols = pal_d3("category20")(10), 
#                        pt.size = 0.8,
#                        na.value = "lightgray")
# plot(p1)

# p1 <- DimPlot(subset(tumour.seurat.obj, vdj_cell_type == "B"), 
#                        group.by = c("b_cln_rnk"), 
#                        cols = pal_d3("category20")(10), 
#                        pt.size = 0.8,
#                        na.value = "lightgray")
# plot(p1)


# p1 <- DimPlot(subset(tumour.seurat.obj, vdj_cell_type == "T"), 
#                        group.by = "t_expanded_cln", 
#                        cols = pal_d3("category20")(10), 
#                        pt.size = 0.8,
#                        na.value = "lightgray")
# plot(p1)

# p1 <- DimPlot(subset(tumour.seurat.obj, vdj_cell_type == "T"), 
#                        group.by = "t_cln_rnk", 
#                        cols = pal_d3("category20")(10), 
#                        pt.size = 0.8,
#                        na.value = "lightgray")
# plot(p1)

#can use the VDJ data when exploring GEX too
p1 <- VlnPlot(tumour.seurat.obj, group.by = "seurat_clusters", split.by = "t_expanded_cln", features = c("CD3E"))
plot(p1)

p1 <- FeaturePlot(tumour.seurat.obj, split.by = "t_expanded_cln", features = c("CD3E"))
plot(p1)

rm(p1)
gc()

#plotting with seurat is not very flexible, so let's regain some control
## PBMCs
pbmc.seurat.meta <- as.data.frame(pbmc.seurat.obj@meta.data) %>% 
  rownames_to_column("barcode") 
pbmc.seurat.umap <- as.data.frame(pbmc.seurat.obj@reductions$umap@cell.embeddings) %>%
  rownames_to_column("barcode")
pbmc.seurat.data <- pbmc.seurat.meta %>%
  left_join(pbmc.seurat.umap, by = "barcode") %>%
  as_tibble()

pbmc.seurat.data %>% head()

rm(list = c("pbmc.seurat.meta", "pbmc.seurat.umap"))

## Tumour
tumour.seurat.meta <- as.data.frame(tumour.seurat.obj@meta.data) %>% 
  rownames_to_column("barcode") 
tumour.seurat.umap <- as.data.frame(tumour.seurat.obj@reductions$umap@cell.embeddings) %>%
  rownames_to_column("barcode")
tumour.seurat.data <- tumour.seurat.meta %>%
  left_join(tumour.seurat.umap, by = "barcode") %>%
  as_tibble()

tumour.seurat.data %>% head()

rm(list = c("tumour.seurat.meta", "tumour.seurat.umap",
            "tumour.seurat.obj", "pbmc.seurat.obj"))
gc()

#now have greater control over the plotting & the format of the plots

#an example of plotting B cell SHM levels on the UMAP
pbmc.shm.umap.p1 <- ggplot(pbmc.seurat.data, aes(x = umap_1, y = umap_2)) +
  geom_point(data = pbmc.seurat.data %>% filter(vdj_cell_type != "B"),
             size = 0.5,
             colour = "lightgray") + #plotting cells that are NOT classed as B cell from VDJ on the bottom layer of the plot with a smaller point size
  geom_point(data = pbmc.seurat.data %>% filter(vdj_cell_type == "B"),
             size = 1,
             aes(colour = vh_mut)) + #plotting cells that are classed as B cells on top layer and colouring by IGH SHM level
  scale_color_viridis_c(name = "VH SHM %") + #using the continous viridis for the colour
  theme_custom() +
  labs(x = "UMAP_1",
       y = "UMAP_2",
       title = "PBMC IGHV SHM",
       caption = "non B cells in gray") #labels for plot
plot(pbmc.shm.umap.p1)


tumour.shm.umap.p1 <- ggplot(tumour.seurat.data, aes(x = umap_1, y = umap_2)) +
  geom_point(data = tumour.seurat.data %>% filter(vdj_cell_type != "B"),
             size = 0.5,
             colour = "lightgray") + #plotting cells that are NOT classed as B cell from VDJ on the bottom layer of the plot with a smaller point size
  geom_point(data = tumour.seurat.data %>% filter(vdj_cell_type == "B"),
             size = 1,
             aes(colour = vh_mut)) + #plotting cells that are classed as B cells on top layer and colouring by IGH SHM level
  scale_color_viridis_c(name = "VH SHM %") + #using the continuous viridis for the colour
  theme_custom() +
  labs(x = "UMAP_1",
       y = "UMAP_2",
       title = "tumour IGHV SHM",
       caption = "non B cells in gray") #labels for plot
plot(tumour.shm.umap.p1)

shm.umap.panel <- ggarrange(pbmc.shm.umap.p1, tumour.shm.umap.p1)
shm.umap.panel

#only plotting the B cells, but splitting by isotype
pbmc.shm.umap.p2 <- ggplot(pbmc.seurat.data %>% filter(vdj_cell_type == "B"), aes(x = umap_1, y = umap_2)) +
  geom_point(size = 1,
             aes(colour = vh_mut)) + 
  scale_color_viridis_c(name = "VH SHM %") + #using the continuous viridis for the colour
  theme_custom() +
  labs(x = "UMAP_1",
       y = "UMAP_2",
       title = "PBMC IGHV SHM",
       caption = "only showing cells with Ig VDJ") + #labels for plot
  scale_y_continuous(limits = c(-20, 15)) +
  scale_x_continuous(limits = c(-15, 10)) +
  facet_wrap(~ isotype, nrow = 2, scales = "free") #splitting by the isotype subclass
plot(pbmc.shm.umap.p2)

tumour.shm.umap.p2 <- ggplot(tumour.seurat.data %>% filter(vdj_cell_type == "B"), aes(x = umap_1, y = umap_2)) +
  geom_point(size = 1,
             aes(colour = vh_mut)) + 
  scale_color_viridis_c(name = "VH SHM %") + #using the continuous viridis for the colour
  theme_custom() +
  labs(x = "UMAP_1",
       y = "UMAP_2",
       title = "tumour IGHV SHM",
       caption = "only showing cells with Ig VDJ") + #labels for plot
  scale_y_continuous(limits = c(-20, 15)) +
  scale_x_continuous(limits = c(-15, 25)) +
  facet_wrap(~ isotype, nrow = 2, scales = "free") #splitting by the isotype subclass
plot(tumour.shm.umap.p2)

#splitting by seurat cluster
tumour.shm.umap.p3 <- ggplot(tumour.seurat.data %>% filter(vdj_cell_type == "B"), aes(x = umap_1, y = umap_2)) +
  geom_point(size = 1,
             aes(colour = vh_mut)) + 
  scale_color_viridis_c(name = "VH SHM %") + #using the continuous viridis for the colour
  theme_custom() +
  labs(x = "UMAP_1",
       y = "UMAP_2",
       title = "tumour IGHV SHM",
       caption = "only showing cells with Ig VDJ") + #labels for plot
  scale_y_continuous(limits = c(-20, 15)) +
  scale_x_continuous(limits = c(-15, 25)) +
  facet_wrap(~ seurat_clusters, nrow = 2, scales = "free") #splitting by the cluster label
plot(tumour.shm.umap.p3)

##while UMAP is good for visualising the clusters, can pair the UMAPs up with other summary plots
#combine the datasets and create factor levels for the isotypes

## isotype order by xsome location
isotype_order <- c("IGHM", "IGHD", "IGHG3", "IGHG1", "IGHA1", "IGHG2", "IGHG4", "IGHE", "IGHA2")

##combining the datasets and adding the isotype levels
seurat.data <- bind_rows(pbmc.seurat.data %>% mutate(dataset = "PBMC"),
                         tumour.seurat.data %>% mutate(dataset = "tumour")) %>%
  mutate(isotype_fac = factor(isotype, isotype_order))


shm.summary.p1 <- ggplot(seurat.data %>% filter(vdj_cell_type == "B"), 
                         aes(x = isotype_fac, y = vh_mut)) +
  geom_boxplot(aes(group = interaction(dataset, isotype_fac)),
               position = position_dodge(width = 0.9)) +
  geom_point(aes(colour = dataset, group = dataset),
             position = position_dodge(width = 0.9)) +
  theme_x_rotated() +
  labs(x = "isotype subclass",
       y = "VH SHM %") +
  scale_colour_aaas() +
  guides(colour = guide_legend(override.aes = list(size = 5)))
plot(shm.summary.p1)

#what about some stats for this
shm.summary.stats <- seurat.data %>% 
  filter(vdj_cell_type == "B") %>% #only using data from B cells by VDJ
  filter(isotype_fac %in% c("IGHM", "IGHD", "IGHG3", "IGHG1", "IGHA1", "IGHG2", "IGHA2")) %>% #not all isotype have necessary data points, so restrict
  group_by(isotype_fac) %>% #want to perform the tests for each isotype
  wilcox_test(vh_mut ~ dataset) %>% #do th wilcoxon test for difference in means between the datasets
  ungroup() %>% #remove the grouping by isotype
  adjust_pvalue(method = "bonferroni") %>% # do the p-value adjustment
  add_xy_position(x = "isotype_fac", dodge = 0.9) %>% #add the X & y-positions that are needed for plotting
  mutate(group1 = isotype_fac,
         group2 = isotype_fac) #adjust groups to get values in correct position on panel

#add these to the panel using ggpubr function stat_pvalue_manual
shm.summary.p2 <- ggplot(seurat.data %>% filter(vdj_cell_type == "B"), 
                         aes(x = isotype_fac, y = vh_mut)) +
  geom_boxplot(aes(group = interaction(dataset, isotype_fac)),
               position = position_dodge(width = 0.9)) +
  geom_point(aes(colour = dataset, group = dataset),
             position = position_dodge(width = 0.9)) +
  stat_pvalue_manual(data = shm.summary.stats,
                     hide.ns = FALSE) +
  theme_custom() +
  labs(x = "isotype subclass",
       y = "VH SHM %",
       caption = "p-values: Wilcoxon, bonf adjusted") +
  scale_colour_aaas() +
  guides(colour = guide_legend(override.aes = list(size = 5)))
plot(shm.summary.p2)

plot(ggarrange(ggarrange(pbmc.shm.umap.p1, tumour.shm.umap.p1, nrow = 1),
               shm.summary.p2, nrow = 2, ncol = 1))


#how many clones/clonotypes are expanded?
exp.umap.p1 <- ggplot(seurat.data, aes(x = umap_1, y = umap_2)) +
  geom_point(data = seurat.data %>% filter(is.na(vdj_cell_type)),
             size = 0.5, colour = "lightgray") + #bottom layer are the cell that don't have any VDJ
  geom_point(data = seurat.data %>% filter(vdj_cell_type == "B", b_expanded_cln != "expanded"),
             size = 0.8,
             aes(colour = b_expanded_cln)) + #B cells, singleton
  geom_point(data = seurat.data %>% filter(vdj_cell_type == "B", b_expanded_cln == "expanded"),
             size = 0.8,
             aes(colour = b_expanded_cln)) + #B cells, expanded
  geom_point(data = seurat.data %>% filter(vdj_cell_type == "T", t_expanded_cln != "expanded"),
             aes(colour = t_expanded_cln),
             size = 0.8) + #T cells, singleton
  geom_point(data = seurat.data %>% filter(vdj_cell_type == "T", t_expanded_cln == "expanded"),
             aes(colour = t_expanded_cln),
             size = 0.8) + #T cells, expanded
  theme_custom() +
  facet_grid(dataset ~ vdj_cell_type) +
  scale_colour_npg(name = "expansions") +
  guides(colour = guide_legend(override.aes = list(size = 5)))
exp.umap.p1


b.expanded.clns.summary <- b.data %>%
  group_by(dataset, b_expanded_cln) %>%
  summarise(cell_cnt = length(barcode),
            cln_cnt = length(unique(unq_clone_id))) %>%
  ungroup() %>%
  group_by(dataset) %>%
  mutate(cell_prop = cell_cnt / sum(cell_cnt),
         cln_prop = cln_cnt / sum(cln_cnt)) %>%
  ungroup()

b.expanded.clns.summary

b.exp.p1 <- ggplot(b.expanded.clns.summary, aes(x = dataset, y = cln_prop * 100)) +
  geom_bar(stat = "identity", position = "stack",
           aes(fill = b_expanded_cln)) +
  scale_fill_npg(name = "expanded clone") +
  theme_custom() +
  theme(legend.position = "bottom") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(y = "% clones")
plot(b.exp.p1)

b.exp.p2 <- ggplot(b.expanded.clns.summary, aes(x = dataset, y = cell_prop * 100)) +
  geom_bar(stat = "identity", position = "stack",
           aes(fill = b_expanded_cln)) +
  scale_fill_npg(name = "expanded clone") +
  theme_custom() +
  theme(legend.position = "bottom") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(y = "% cells")
plot(b.exp.p2)


t.expanded.clns.summary <- t.data %>%
  group_by(dataset, t_expanded_cln) %>%
  summarise(cell_cnt = length(barcode),
            cln_cnt = length(unique(paired_clonotype))) %>%
  ungroup() %>%
  group_by(dataset) %>%
  mutate(cell_prop = cell_cnt / sum(cell_cnt),
         cln_prop = cln_cnt / sum(cln_cnt)) %>%
  ungroup()

t.expanded.clns.summary

t.exp.p1 <- ggplot(t.expanded.clns.summary, aes(x = dataset, y = cln_prop * 100)) +
  geom_bar(stat = "identity", position = "stack",
           aes(fill = t_expanded_cln)) +
  scale_fill_npg(name = "expanded clone") +
  theme_custom() +
  theme(legend.position = "bottom") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(y = "% clonotypes")
plot(t.exp.p1)

t.exp.p2 <- ggplot(t.expanded.clns.summary, aes(x = dataset, y = cell_prop * 100)) +
  geom_bar(stat = "identity", position = "stack",
           aes(fill = t_expanded_cln)) +
  scale_fill_npg(name = "expanded clone") +
  theme_custom() +
  theme(legend.position = "bottom") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(y = "% cells")
plot(t.exp.p2)

##stacked bar charts for the expanded clonotypes
b.cln.size.summary <- seurat.data %>%
  filter(vdj_cell_type == "B") %>%
  mutate(label = ifelse(b_expanded_cln == "expanded", paste0(unq_clone_id), "singletons")) %>%
  group_by(dataset, label) %>%
  summarise(cell_cnt = length(barcode)) %>%
  ungroup() %>%
  group_by(dataset) %>%
  mutate(cell_prop = cell_cnt / sum(cell_cnt)) %>%
  ungroup()

b.cln.size.summary

b.cln.sizes.p1 <- ggplot(b.cln.size.summary %>% filter(label != "singletons"), 
                         aes(x = dataset, y = cell_cnt, group = cell_cnt)) +
  geom_bar(stat = "identity", 
           position = position_stack(),
           aes(fill = factor(label)),
           linewidth = 0.05, colour = "gray") +
  scale_fill_manual(values = colorRampPalette(pal_d3("category20")(20))(40),
                    name = "clone") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_custom() +
  theme(legend.position = "none") +
  labs(y = "# cells")
b.cln.sizes.p1


t.cln.size.summary <- seurat.data %>%
  filter(vdj_cell_type == "T") %>%
  mutate(label = ifelse(t_expanded_cln == "expanded", paste0(paired_clonotype), "singletons")) %>%
  group_by(dataset, label) %>%
  summarise(cell_cnt = length(barcode)) %>%
  ungroup() %>%
  group_by(dataset) %>%
  mutate(cell_prop = cell_cnt / sum(cell_cnt)) %>%
  ungroup()

t.cln.size.summary


t.cln.sizes.p1 <- ggplot(t.cln.size.summary %>% filter(label != "singletons"), 
                         aes(x = dataset, y = cell_cnt, group = cell_cnt)) +
  geom_bar(stat = "identity", 
           position = position_stack(),
           aes(fill = factor(label)),
           linewidth = 0.05, colour = "gray") +
  scale_fill_manual(values = colorRampPalette(pal_d3("category20")(20))(70),
                    name = "clone") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_custom() +
  theme(legend.position = "none") +
  labs(y = "# cells")
t.cln.sizes.p1


#build into a single panel to summarise the expansions across the two datasets
## turning off the legends as will be provided by the UMAP but adding a title to each for cell type
b.exp.panel <- ggarrange(b.exp.p1, b.exp.p2, common.legend = TRUE, legend = "none")
b.exp.panel <- annotate_figure(b.exp.panel, top = text_grob("B cells", color = "black", face = "bold", size = 14))
#plot(b.exp.panel)

t.exp.panel <- ggarrange(t.exp.p1, t.exp.p2, common.legend = TRUE, legend = "none")
t.exp.panel <- annotate_figure(t.exp.panel, top = text_grob("T cells", color = "black", face = "bold", size = 14))
#plot(t.exp.panel)

b.cln.sizes.p1 <- annotate_figure(b.cln.sizes.p1, top = text_grob("B expanded", color = "black", face = "bold", size = 14))
t.cln.sizes.p1 <- annotate_figure(t.cln.sizes.p1, top = text_grob("T expanded", color = "black", face = "bold", size = 14))

#combining the bar plots and the UMAP with the stacked bar charts for clonotype distribution
plot(ggarrange(exp.umap.p1 + theme(legend.position = "top"), 
               ggarrange(b.exp.panel, t.exp.panel, nrow = 2), 
               ggarrange(b.cln.sizes.p1, t.cln.sizes.p1, nrow = 2),
               nrow = 1, widths = c(4, 2, 1)))

#save the tibble that include the seurat with VDJ
write_tsv(seurat.data, "data/pbmc-tumour_seruat_data.tsv")

##save the b and t cell data that has the additional columns and filtering
write_tsv(b.data, "data/pbmc-tumour_b_cells_vdj.tsv")
write_tsv(t.data, "data/pbmc-tumour_t_cells_vdj.tsv")


##what about adding cluster information to the umap when we are already colouring cells by another feature?
## using hull plots to outline clusters
if (sum(grepl("ggforce", rownames(installed.packages()))) > 0) {
  print("ggforce is already installed")
} else {
  install.packages("ggforce")
}

library(ggforce)

#for hull plot to work need to exclude cells that are well outside their cluster
clst.outliers <- seurat.data %>% 
  group_by(dataset, seurat_clusters) %>%
  mutate(median_umap_1 = median(umap_1),
         median_umap_2 = median(umap_2)) %>%
  ungroup() %>%
  mutate(cluster_outlier = case_when((umap_1 < (median_umap_1 - 2.2)) | (umap_1 > (median_umap_1 + 2.2)) ~ TRUE,
                                     (umap_2 < (median_umap_2 - 2.2)) | (`umap_2` > (median_umap_2 + 2.2)) ~ TRUE,
                                     TRUE ~ FALSE))

##points coloured by seurat cluster
umap.clst.p1 <- ggplot(seurat.data, aes(x = umap_1, y = umap_2, color = factor(seurat_clusters))) +
  geom_point(size = 0.2) +
  scale_color_d3(palette = "category20",
                 name = "") +
  theme_custom() +
  facet_wrap(~ dataset) +
  guides(colour = guide_legend(override.aes = list(size = 5)))
umap.clst.p1

##points + hull plot for cluster regions
umap.clst.p2 <- ggplot(seurat.data, aes(x = umap_1, y = umap_2, color = factor(seurat_clusters))) +
  geom_point(size = 0.4) +
  geom_mark_hull(data = clst.outliers %>% filter(!cluster_outlier),
                 aes(fill = factor(seurat_clusters)), concavity = 3, expand = 0.01, radius = 0.01) +
  scale_color_d3(palette = "category20",
                 name = "") +
  scale_fill_d3(palette = "category20",
                name = "") +
  theme_custom() +
  facet_wrap(~ dataset) +
  guides(colour = guide_legend(override.aes = list(size = 5)))
umap.clst.p2

plot(ggarrange(umap.clst.p1, umap.clst.p2, nrow = 2))

##keeping hull outline, but showing expanded clones
p1 <- ggplot(seurat.data, aes(x = umap_1, y = umap_2)) +
  geom_point(size = 0.8, shape = 21, 
             aes(fill = vh_mut),
             stroke = NA) +
  geom_mark_hull(data = clst.outliers %>% filter(!cluster_outlier),
                 aes(colour = factor(seurat_clusters)), 
                 concavity = 3, expand = 0.01, radius = 0.01,
                 linetype = "dashed", linewidth = 1) +
  scale_fill_viridis_c(name = "VH SHM%") +
  scale_colour_d3(palette = "category20",
                  name = "Cluster") +
  theme_custom() +
  facet_wrap(~ dataset) +
  guides(colour = guide_legend(override.aes = list(size = 5)))
p1
