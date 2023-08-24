#setting working directory
#setwd('~/Documents/GitHub/scRNA-seq_VDJ/')

#loading R packages
#tidyverse packages for data manipulation and plotting
library(tidyverse)

##two sources to use
## immuneCODE for COV2 specific https://clients.adaptivebiotech.com/pub/covid-2020, reference: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7418738/
## vdjdb for more diverse antigens tsv download https://github.com/antigenomics/vdjdb-db/releases/tag/2023-06-01 
## web interface https://vdjdb.cdr3.net/, reference: https://www.nature.com/articles/s41592-022-01578-0

#obtaining TCRs of known specificity

#if running on local computer change to the directory of your choice, on posit cloud start in the /cloud/project/ directory which is the default
#cd ~/Documents/GitHub/scRNA-seq_VDJ/

#make a new directory to store the references
# mkdir -p tcr_ag_refs/vdjdb/
# mkdir -p tcr_ag_refs/immuneCODE/

#downloading the file with cURL, need to include the -L option to follow the redirect
# cd tcr_ag_refs/vdjdb/
# curl -L -O https://github.com/antigenomics/vdjdb-db/releases/download/2023-06-01/vdjdb-2023-06-01.zip
# unzip vdjdb-2023-06-01.zip
# rm vdjdb-2023-06-01.zip

# cd ../immuneCODE/
# curl -O https://adaptivepublic.blob.core.windows.net/publishedproject-supplements/covid-2020/ImmuneCODE-MIRA-Release002.1.zip
# unzip ImmuneCODE-MIRA-Release002.1.zip
# rm  ImmuneCODE-MIRA-Release002.1.zip
# mv */* ./ 
# rmdir ImmuneCODE-MIRA-Release002.1

##read in the different files
immunecode.data1 <- read_csv("tcr_ag_refs/immuneCODE/minigene-detail.csv")
immunecode.data2 <- read_csv("tcr_ag_refs/immuneCODE/peptide-detail-ci.csv")
immunecode.data3 <- read_csv("tcr_ag_refs/immuneCODE/peptide-detail-cii.csv")

##combine into a single tibble and separate out the TCR clonotype information that is in the TCR BioIndenity columns
immunecode.data <- bind_rows(immunecode.data1 %>% select(tcr = `TCR BioIdentity`, epitope = ORF),
                             immunecode.data2 %>% select(tcr = `TCR BioIdentity`, epitope = `ORF Coverage`),
                             immunecode.data3 %>% select(tcr = `TCR BioIdentity`, epitope = `ORF Coverage`)) %>%
  separate(tcr, into = c("junction_aa", "v", "j"), sep = "\\+", remove = FALSE) %>%
  separate(v, into = c("v1", "v2"), sep = "\\/", remove = FALSE, fill = "right") %>%
  mutate(v2 = ifelse(!(is.na(v2)), paste0("TRBV", v2), NA_character_)) %>%
  pivot_longer(cols = c(v1, v2),
               values_to = "trbv",
               names_to = "v_alt") %>%
  filter(!(is.na(trbv))) %>%
  select(-c(v_alt, tcr, v)) %>%
  mutate(antigen = "SARS-COV-2",
         trbv = str_replace(trbv, "TCRBV0", "TRBV"),
         trbv = str_replace(trbv, "TCRBV", "TRBV"),
         trbv = str_replace(trbv, "-0", "-"),
         j = str_replace(j, "TCRBJ0", "TRBJ"),
         j = str_replace(j, "-0", ""),
         trbv = str_replace(trbv, "TRBV20-X", "TRBV20-1"),
         cdr3_aa = str_replace_all(junction_aa, "^.|.$", ""),
         clonotype = paste(trbv, j, cdr3_aa, sep = "_"),
         source = "immuneCODE") %>%
  filter(!grepl("nproductiv", clonotype)) 
immunecode.data

##reformat the data
immunecode.data <- immunecode.data %>% 
  select(clonotype, gene = epitope, antigen, source) %>%
  unique()


##loading the vdjdb data
vdjdb.data <- read_tsv("tcr_ag_refs/vdjdb/vdjdb_full.txt") %>%
  filter(!is.na(cdr3.beta)) %>%
  filter(!is.na(v.beta)) %>%
  filter(!is.na(j.beta)) %>%
  filter(species == "HomoSapiens") %>%
  mutate(vb = str_replace(v.beta, "\\*.+$", ""),
         jb = str_replace(j.beta, "\\*.+$", ""),
         cdr3 = str_replace_all(cdr3.beta, "^.|.$", ""),
         clonotype = paste(vb, jb, cdr3, sep = "_")) %>%
  select(clonotype, antigen = antigen.species, gene = antigen.gene) %>%
  unique() %>%
  mutate(source = "vdjdb")

#combining the databases
tcr.antigen.data <- bind_rows(immunecode.data %>% mutate(database = "immuneCODE"), 
                              vdjdb.data %>% mutate(database = "vdjdb")) %>%
  group_by(clonotype) %>%
  summarise(gene = str_c(unique(gene), collapse = ","),
            antigen = str_c(unique(antigen), collapse = ","),
            database = str_c(unique(database), collapse = ",")) %>%
  ungroup()

tcr.antigen.data

tcr.antigen.data %>% nrow()

tcr.ag.tbl <- tcr.antigen.data %>% 
  group_by(antigen) %>%
  summarise(cnt = length(clonotype)) %>%
  ungroup() %>%
  arrange(desc(cnt))
tcr.ag.tbl

#save the tcr antigen file for later use
write_tsv(tcr.antigen.data, "tcr_ag_refs/trb_antigens.tsv")

##go back to 'Terminal' tab and remove the files that are no longer required
#cd /cloud/project/tcr_ag_refs/
#rm -rf immuneCODE
#rm -rf vdjdb