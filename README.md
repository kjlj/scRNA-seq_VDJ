# scRNA-seq workshop

## IgBLAST post-processing of 10x contigs

Cell ranger output does not include allele level gene assignment and doesn't provide any information about identity to the germline that is needed for calculating somatic hypermutation (SHM) levels for the B cell datasets. To obtain this extra level of detail the contigs from 10x VDJ are re-aligned against the [IMGT reference directory](https://www.imgt.org/vquest/refseqh.html) using [IgBLAST](https://ncbi.github.io/igblast/). To get tab-delimited [AIRR-C format](https://docs.airr-community.org/en/stable/datarep/overview.html) for the VDJ alignments we will use the stand alone version of IgBLAST rather than the [web-based version](https://www.ncbi.nlm.nih.gov/igblast/igblast.cgi).

In the interest of time, we will not undertake the IgBLAST post-processing during this session. Rather we will start with datasets that have already been processed. 

All steps for generating the datasets that we will be working with are documented:
1. [Setting up IgBLAST and human reference datasets](docs/igblast_setup.md)
2. [Obtaining datasets from 10x genomics](docs/datasets.md)
3. [Running IgBLAST](docs/running_igblast.md)
4. [Combining IgBLAST & cellranger VDJ data summarise per cell barcode](https://kjlj.github.io/scRNA-seq_VDJ/docs/joining_cellranger_igblast.html)
5. [Clone clustering for B cells](https://kjlj.github.io/scRNA-seq_VDJ/docs/building_b_cell_clones.html)

The results of running the above steps are two tab-delimited files; one for the T cells and one for B cells:
- B cell VDJ data with B cell clones defined: [pbmc-tumour_Ig_cellranger_igblast_per-barcode_clns.tsv](https://raw.githubusercontent.com/kjlj/scRNA-seq_VDJ/main/data/pbmc-tumour_Ig_cellranger_igblast_per-barcode_clns.tsv)
- T cell VDJ data with clonotypes defined: [pbmc-tumour_TR_cellranger_igblast_per-barcode.tsv](https://raw.githubusercontent.com/kjlj/scRNA-seq_VDJ/main/data/pbmc-tumour_TR_cellranger_igblast_per-barcode.tsv)

These are the two files that we will use for integrating VDJ data with scRNA-seq data.

## Integrating scRNA-seq and VDJ

This session with focus on utilising post-processed VDJ data from 10x VDJ for integration with scRNA-seq. We will build off the earlier sessions that generated the Seurat objects with the Azimuth annotation and the cell clustering. The Seurat objects are available from the git at [pbmc](dummy_location) and [tumour](dummy_location).

Topics:
- [Adding VDJ to scRNA-seq]() including displaying VDJs on UMAP and summarising T and B cell repertoires.
- [Antigen annotation for T cells]()
- [DGE between clones/clonotypes]()

## Additional information

### alternative reference sets

Alternatvies to IMGT exists for sets of germline genes to use for references. NCBI, OGRdb.

### rhapsody V(D)J

For Rhapsody data can use the VDJ_Dominant_Contigs_AIRR tab-delimed file as a starting point for post-processing. The Rhapsody output already includes the allele level information for the nearest germline match and appears to use a recent version of the IMGT reference, however it lacks percent identity columns in the output that can be used for SHM levels. The % identity can be calculated from the available output as the sequence and germline strings are included, but can optionally re-run IgBLAST to get the output.

To run IgBLAST, first need to extract the sequences from the tab-delimited file and create FASTA file(s) for input to IgBLAST.

TODO: add example of script to extract
