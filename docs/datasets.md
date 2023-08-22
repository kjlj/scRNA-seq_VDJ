# Datasets

This workshop uses two datasets available for the [10x genomics website](https://www.10xgenomics.com/). 

The first is a 10x GEX + VDJ for [PBMC sample](https://www.10xgenomics.com/resources/datasets/human-pbmc-from-a-healthy-donor-10-k-cells-v-2-2-standard-5-0-0) and the second is 10x GEX + VDJ for a melanoma [tumour sample](https://www.10xgenomics.com/resources/datasets/melanoma-tumor-derived-cells-v-2-2-standard-4-0-0).

For post-processing of the VDJ need to obtain the B and T cell filtered_contig.fasta files from the 10x website.

Open up terminal and download the fasta files from the command line.
```
cd ~/data/
mkdir scRNA-seq_workshop/
cd scRNA-seq_workshop/
mkdir pbmcs/
mkdir tumour/

#PBMC dataset from [https://www.10xgenomics.com/resources/datasets/human-pbmc-from-a-healthy-donor-10-k-cells-v-2-2-standard-5-0-0]
cd ~/data/scRNA-seq_workshop/pbmcs/
## Ig
curl -O https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_PBMC_10k/sc5p_v2_hs_PBMC_10k_b_filtered_contig.fasta

## TR
curl -O https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_PBMC_10k/sc5p_v2_hs_PBMC_10k_t_filtered_contig.fasta

#tumor dataset from [https://www.10xgenomics.com/resources/datasets/melanoma-tumor-derived-cells-v-2-2-standard-4-0-0]
cd ~/data/scRNA-seq_workshop/tumour/
## Ig
curl -O https://cf.10xgenomics.com/samples/cell-vdj/4.0.0/sc5p_v2_hs_melanoma_10k/sc5p_v2_hs_melanoma_10k_b_filtered_contig.fasta

## TR
curl -O https://cf.10xgenomics.com/samples/cell-vdj/4.0.0/sc5p_v2_hs_melanoma_10k/sc5p_v2_hs_melanoma_10k_t_filtered_contig.fasta

```

Alternatively, can also obtain the files direct from the 10x website; [PBMCs](https://www.10xgenomics.com/resources/datasets/human-pbmc-from-a-healthy-donor-10-k-cells-v-2-2-standard-5-0-0) and [tumour](https://www.10xgenomics.com/resources/datasets/melanoma-tumor-derived-cells-v-2-2-standard-4-0-0). For each dataset need to grab the filtered_contig.fasta files for the T and B cell assays.

[Return to main page.](../README.md)