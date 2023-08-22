# scRNA-seq workshop

## IgBLAST post-processing of 10x contigs

Cell ranger output does not include allele level gene assignment and doesn't provide any information about identity to the germline that is needed for calculating somatic hypermutation (SHM) levels for the B cell datasets. To obtain this extra level of detail the contigs from 10x VDJ are re-aligned against the [IMGT reference directory](https://www.imgt.org/vquest/refseqh.html) using [IgBLAST](https://www.ncbi.nlm.nih.gov/igblast/igblast.cgi). To get tab-delimited AIRR-C format for the VDJ alignments we will use the stand alone version of IgBLAST rather than the [web-based version](https://ncbi.github.io/igblast/).

In the interest of time, we will not undertake the IgBLAST post-processing during this session. Rather we will start with datasets that have already been processed. 

All steps for generating the datasets that we will be working with are documented:

Step one - [setting up IgBLAST and human reference datasets](docs/igblast_setup.md)
Step one - [obtaining datasets from 10x genomics](docs/datasets.md)
Step two - [running IgBLAST](docs/running_igblast.md)
Step three - [combining IgBLAST & cellranger VDJ data summarise per cell barcode](https://kjlj.github.io/scRNA-seq_VDJ/docs/joining_cellranger_igblast.html)
Step four - [clone clustering for B cells](https://kjlj.github.io/scRNA-seq_VDJ/docs/building_b_cell_clones.html)

## adding VDJ data to UMAP

### within Seurat

### outside Seurat 

## antigen annotation for TR data

### obtaining TCRs of known specificity

### annotating TR clonotypes with Ag specificity 

## differential gene expression between clones/clonotypes

## misc

### alternative reference sets

Alternatvies to IMGT exists for sets of germline genes to use for references. NCBI, OGRdb.

### rhapsody V(D)J

For Rhapsody data can use the VDJ_Dominant_Contigs_AIRR tab-delimed file as a starting point for post-processing. The Rhapsody output already includes the allele level information for the nearest germline match and appears to use a recent version of the IMGT reference, however it lacks percent identity columns in the output that can be used for SHM levels. The % identity can be calculated from the available output as the sequence and germline strings are included, but can optionally re-run IgBLAST to get the output.

To run IgBLAST, first need to extract the sequences from the tab-delimited file and create FASTA file(s) for input to IgBLAST.

TODO: add example of script to extract
