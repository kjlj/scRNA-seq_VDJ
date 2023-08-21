#set working directory
setwd("~/Documents/GitHub/scRNA-seq_VDJ/")

#installing packages, only for first instance
#install.packages("scoper")
#install.packages("shazam")

#loading R packages
library(tidyverse)
library(scoper)
library(shazam)

##themes for plots
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


#loading the B and T cell post-processed data
bcr.data <- read_tsv("https://raw.githubusercontent.com/kjlj/scRNA-seq_VDJ/main/pbmc-tumour_Ig_cellranger_igblast_per-barcode.tsv")

#split into tumour and pbmc
pbmc.data <- bcr.data %>% filter(dataset == "pbmc")
tumour.data <- bcr.data %>% filter(dataset == "tumour")

#first need to find the clustering threshold
## this is done by comparing the hamming distance between sequences which should result in a bimodal distribution that define 'clones' and 'unrelated' sequences
## for this, we can use the IGH data
## the shazam package uses the IgBLAST columns, so can just select those and rename appropriately
pbmc.igh.data <- pbmc.data %>% 
  select(barcode, sequence_id = contig_id_igh, ends_with("_igblast_igh")) %>%
  rename_with(~ str_replace(.x, "_igblast_igh", ""), everything()) %>%
  filter(!(is.na(cdr3))) #need to remove any unpaired entries as they will cause errors

tumour.igh.data <- tumour.data %>%
  select(barcode, sequence_id = contig_id_igh, ends_with("_igblast_igh")) %>%
  rename_with(~ str_replace(.x, "_igblast_igh", ""), everything()) %>%
  filter(!(is.na(cdr3))) #need to remove any unpaired entries as they will cause errors

#add the distance to the nearest other heavy chain in the dataset to the tibbles
pbmc.igh.data <- distToNearest(pbmc.igh.data, nproc = 1)
tumour.igh.data <- distToNearest(tumour.igh.data, nproc = 1)

#plotting the distribution
pbmc.dist.p1 <- ggplot(pbmc.igh.data %>% filter(!is.na(dist_nearest)),
             aes(x = dist_nearest)) +
  theme_bw() +
  xlab("Hamming distance") + ylab("Count") +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  geom_histogram(color = "white", binwidth = 0.02) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_custom()
plot(pbmc.dist.p1)

#add a line for threshold
pbmc.dist.p1 <- pbmc.dist.p1 +
  geom_vline(xintercept = 0.1, linetype = "dashed", colour = "firebrick")
plot(pbmc.dist.p1)

#plotting distribution for the tumour sample
tumour.dist.p1 <- ggplot(tumour.igh.data %>% filter(!is.na(dist_nearest)),
                       aes(x = dist_nearest)) +
  theme_bw() +
  xlab("Hamming distance") + ylab("Count") +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  geom_histogram(color = "white", binwidth = 0.02) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_custom()
plot(tumour.dist.p1)

#add a line for threshold
tumour.dist.p1 <- tumour.dist.p1 +
  geom_vline(xintercept = 0.1, linetype = "dashed", colour = "firebrick")
plot(tumour.dist.p1)


#looks like 0.1 would be suitable threshold for both samples
threshold <- 0.1

#to run the clone building in 'single-cell' mode need to reformat the tibble to have igh and igkl as separate rows
pbmc.igkl.data <- pbmc.data %>% 
  select(barcode, sequence_id = contig_id_igkl, ends_with("_igblast_igkl")) %>%
  rename_with(~ str_replace(.x, "_igblast_igkl", ""), everything()) %>%
  filter(!(is.na(cdr3))) #need to remove any unpaired entries as they will cause errors

tumour.igkl.data <- tumour.data %>%
  select(barcode, sequence_id = contig_id_igkl, ends_with("_igblast_igkl")) %>%
  rename_with(~ str_replace(.x, "_igblast_igkl", ""), everything()) %>%
  filter(!(is.na(cdr3))) #need to remove any unpaired entries as they will cause errors

pbmc.input <- bind_rows(pbmc.igh.data, pbmc.igkl.data) 
tumour.input <- bind_rows(tumour.igh.data, tumour.igkl.data )

#add the clone id columns
pbmc.input <- hierarchicalClones(pbmc.input, 
                              threshold = threshold,
                              cell_id = "barcode",
                              method = "nt",
                              linkage = "average",
                              junction = "cdr3",
                              only_heavy = FALSE, 
                              split_light = TRUE,
                              summarize_clones = FALSE)

#add the clone id columns
tumour.input <- hierarchicalClones(tumour.input, 
                                   threshold = threshold,
                                   cell_id = "barcode",
                                   method = "nt",
                                   linkage = "average",
                                   junction = "cdr3",
                                   only_heavy = FALSE, 
                                   split_light = TRUE,
                                   summarize_clones = FALSE)

#extract the clone_id and add back to the main tibble
pbmc.clns <- pbmc.input %>%
  select(barcode, clone_id) %>%
  unique() %>% #as data is formatted as IGH and IGK/L as separate entries the barcode + clone_id is duplicated for IGH and IGK/l
  mutate(dataset = "pbmc")

tumour.clns <- tumour.input %>%
  select(barcode, clone_id) %>%
  unique() %>%
  mutate(dataset = "tumour")

#adding back
cln.data <- bind_rows(pbmc.clns, tumour.clns) 

#merge the BCR data that originally loaded with the clone labels just created
bcr.cln.data <- bcr.data %>%
  left_join(cln.data, by = c("barcode", "dataset"))

#check that the number of rows match
bcr.data %>% nrow() == bcr.cln.data %>% nrow()

#add a unique labels for the clns that is the clone_id + dataset
bcr.cln.data <- bcr.cln.data %>%
  mutate(unq_clone_id = ifelse(is.na(clone_id), NA_character_, paste(dataset, clone_id, sep = "_")))

#quick summary of clone sizes
cln.size.summary <- bcr.cln.data %>%
  filter(!(is.na(unq_clone_id))) %>% #remove cells that lack a clone call
  group_by(dataset, unq_clone_id) %>% #group by clone_id & dataset
  summarise(cln_size = length(barcode)) %>% #count the number of barcodes in each clones
  ungroup() %>% #remove the grouping
  arrange(desc(cln_size)) #reorder the data from largest to smallest clone size across the two datasets

cln.size.summary

#write the bcr data to a file
write_tsv(bcr.cln.data, "pbmc-tumour_Ig_cellranger_igblast_per-barcode_clns.tsv")
