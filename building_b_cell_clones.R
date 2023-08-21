#set working directory
setwd("~/Documents/GitHub/scRNA-seq_VDJ/")

#loading R packages
library(tidyverse)

#loading the B and T cell post-processed data
bcr.data <- read_tsv()
