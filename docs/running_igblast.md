# running IgBLAST on contigs

The IgBLAST command for Ig and TR differs.

## Ig contigs

Running IgBLAST on Ig:
```
#for the PBMCs
## navigate to the directory where the 10x data was downloaded
cd ~/data/scRNA-seq_workshop/pbmcs/

#IgBLAST requires the IGDATA environmental variable to be set
export IGDATA=~/data/apps/ncbi-igblast-1.21.0/
#can check with what IGDATA is get to with
#echo $IGDATA

##running IgBLAST against the Ig references
##output to tab-delimited AIRR-C format
~/data/apps/ncbi-igblast-1.21.0/bin/igblastn \
    -germline_db_V ~/data/apps/ncbi-igblast-1.21.0/references/imgt_ig_v_human.fa \
    -germline_db_D ~/data/apps/ncbi-igblast-1.21.0/references/imgt_ig_d_human.fa \
    -germline_db_J ~/data/apps/ncbi-igblast-1.21.0/references/imgt_ig_j_human.fa \
    -c_region_db ~/data/apps/ncbi-igblast-1.21.0/references/ncbi_human_c_genes \
    -auxiliary_data ~/data/apps/ncbi-igblast-1.19.0/optional_file/human_gl.aux \
    -domain_system imgt -ig_seqtype Ig -organism human \
    -outfmt 19 \
    -query sc5p_v2_hs_PBMC_10k_b_filtered_contig.fasta \
    -out sc5p_v2_hs_PBMC_10k_b_filtered_contig.ig.igblast.tsv

#for the tumour sample
cd ~/data/scRNA-seq_workshop/tumour/
~/data/apps/ncbi-igblast-1.21.0/bin/igblastn \
    -germline_db_V ~/data/apps/ncbi-igblast-1.21.0/references/imgt_ig_v_human.fa \
    -germline_db_D ~/data/apps/ncbi-igblast-1.21.0/references/imgt_ig_d_human.fa \
    -germline_db_J ~/data/apps/ncbi-igblast-1.21.0/references/imgt_ig_j_human.fa \
    -c_region_db ~/data/apps/ncbi-igblast-1.21.0/references/ncbi_human_c_genes \
    -auxiliary_data ~/data/apps/ncbi-igblast-1.19.0/optional_file/human_gl.aux \
    -domain_system imgt -ig_seqtype Ig -organism human \
    -outfmt 19 \
    -query sc5p_v2_hs_melanoma_10k_b_filtered_contig.fasta \
    -out sc5p_v2_hs_melanoma_10k_b_filtered_contig.ig.igblast.tsv
```

## TR contigs

Running IgBLAST on TR:
```
#for the PBMCs
## navigate to directory for the PBMC data
cd ~/data/scRNA-seq_workshop/pbmcs/

##set IGDATA
export IGDATA=~/data/apps/ncbi-igblast-1.21.0/

##output to tab-delimited AIRR-C format
~/data/apps/ncbi-igblast-1.21.0/bin/igblastn \
    -germline_db_V ~/data/apps/ncbi-igblast-1.21.0/references/imgt_tr_v_human.fa \
    -germline_db_D ~/data/apps/ncbi-igblast-1.21.0/references/imgt_tr_d_human.fa \
    -germline_db_J ~/data/apps/ncbi-igblast-1.21.0/references/imgt_tr_j_human.fa \
    -c_region_db ~/data/apps/ncbi-igblast-1.21.0/references/imgt_tr_c_human.fa \
    -auxiliary_data ~/data/apps/ncbi-igblast-1.19.0/optional_file/human_gl.aux \
    -domain_system imgt -ig_seqtype TCR -organism human \
    -outfmt 19 \
    -query sc5p_v2_hs_PBMC_10k_t_filtered_contig.fasta \
    -out sc5p_v2_hs_PBMC_10k_t_filtered_contig.tr.igblast.tsv


#for the tumour sample
## nabivigate to the dorectory
cd ~/data/scRNA-seq_workshop/tumour/
## running IgBLAST with AIRR-C tab-delim output
~/data/apps/ncbi-igblast-1.21.0/bin/igblastn \
    -germline_db_V ~/data/apps/ncbi-igblast-1.21.0/references/imgt_tr_v_human.fa \
    -germline_db_D ~/data/apps/ncbi-igblast-1.21.0/references/imgt_tr_d_human.fa \
    -germline_db_J ~/data/apps/ncbi-igblast-1.21.0/references/imgt_tr_j_human.fa \
    -c_region_db ~/data/apps/ncbi-igblast-1.21.0/references/imgt_tr_c_human.fa \
    -auxiliary_data ~/data/apps/ncbi-igblast-1.19.0/optional_file/human_gl.aux \
    -domain_system imgt -ig_seqtype TCR -organism human \
    -outfmt 19 \
    -query sc5p_v2_hs_melanoma_10k_t_filtered_contig.fasta \
    -out sc5p_v2_hs_melanoma_10k_t_filtered_contig.tr.igblast.tsv

```