# scRNA-seq workshop

Slides for the scVDJ-seq session of the scRNA-seq workshop are available at [https://kjlj.github.io/scRNA-seq_VDJ/docs/scVDJ_seq_session_31Aug23_KJLJ.pdf](https://kjlj.github.io/scRNA-seq_VDJ/docs/scVDJ_seq_session_31Aug23_KJLJ.pdf).

## IgBLAST post-processing of 10x contigs

`Cell Ranger` output does not include allele level gene assignment and doesn't provide information about identity to the germline gene that is needed for calculating somatic hypermutation (SHM) levels for the B cell datasets. To obtain this extra level of detail the contigs from 10x VDJ are re-aligned against the [IMGT reference directory](https://www.imgt.org/vquest/refseqh.html) using [IgBLAST](https://ncbi.github.io/igblast/). To get tab-delimited [AIRR-C format](https://docs.airr-community.org/en/stable/datarep/overview.html) for the VDJ alignments we will use the stand-alone version of `IgBLAST` rather than the [web-based version](https://www.ncbi.nlm.nih.gov/igblast/igblast.cgi).

In the interest of time, we will not undertake the `IgBLAST` post-processing during this session. 

Rather we will start with datasets that have already been processed. 

The steps for generating the datasets that we will be working with are documented:
1. [Setting up IgBLAST and human reference databases](docs/igblast_setup.md)
2. [Obtaining VDJ datasets from 10x genomics](docs/datasets.md)
3. [Running IgBLAST](docs/running_igblast.md)
4. [Combining IgBLAST and Cell Ranger VDJ data and summarising by each cell barcode](https://kjlj.github.io/scRNA-seq_VDJ/docs/joining_cellranger_igblast.html)
5. [Clone clustering for B cells](https://kjlj.github.io/scRNA-seq_VDJ/docs/building_b_cell_clones.html)

The results of running the above steps are two tab-delimited files; one for the T cells and one for B cells:
- B cell VDJ data with B cell clones defined: [pbmc-tumour_Ig_cellranger_igblast_per-barcode_clns.tsv](https://raw.githubusercontent.com/kjlj/scRNA-seq_VDJ/main/data/pbmc-tumour_Ig_cellranger_igblast_per-barcode_clns.tsv)
- T cell VDJ data with clonotypes defined: [pbmc-tumour_TR_cellranger_igblast_per-barcode.tsv](https://raw.githubusercontent.com/kjlj/scRNA-seq_VDJ/main/data/pbmc-tumour_TR_cellranger_igblast_per-barcode.tsv)

These are the two files that we will use for integrating VDJ data with scRNA-seq data.

### Do I really have to post-process the datasets?

If working with T cell VDJ data only, then probably not, but for B cell VDJ data it is probably worth the effort, but you don't *have* to!

## Integrating scRNA-seq and VDJ

This session with focus on utilising post-processed VDJ data from 10x VDJ for integration with scRNA-seq. 

We will build off the earlier sessions that generated the Seurat objects with the Azimuth annotation and the cell clustering. The Seurat objects are available from the Dropbox at [pbmc](https://www.dropbox.com/scl/fi/4s610vt1mgtmgibdvfsar/pbmc_seurat-without-VDJ-genes-azimuth.rds?rlkey=ftdkxi9mnxezhbb42dqftbqel&dl=0) and [tumour](https://www.dropbox.com/scl/fi/scik47zay4x27t4wxmo70/tumour_seurat-without-VDJ-genes-azimuth.rds?rlkey=z8ghoeoboaneniv82e2xywji3&dl=0), but have been pre-loaded into RStudio in Posit.cloud.

Topics:
- [Adding VDJ to scRNA-seq](https://kjlj.github.io/scRNA-seq_VDJ/docs/combining_GEX_and_VDJ.html) including displaying VDJ features on UMAPs.
- Antigen annotation for T cells
  - [generating reference datatset](https://kjlj.github.io/scRNA-seq_VDJ/docs/generating_Ag_reference_for_TRB.html)
  - [annotating single cell clontoypes](https://kjlj.github.io/scRNA-seq_VDJ/docs/antigen_annotation_T_cells.html)
- [DGE between clones/clonotypes](https://kjlj.github.io/scRNA-seq_VDJ/docs/dge_using_VDJ_features.html)

## Additional information

### Alternative reference sets

Alternatives to IMGT exists for sets of germline genes to use with `IgBLAST`.

One useful resource from the AIRR-C (Adaptive Immune Receptor Repertoire Community) is [OGRDB](https://ogrdb.airr-community.org/). OGRDB offers FASTA formatted germline reference sets for [direct or API download](https://wordpress.vdjbase.org/index.php/ogrdb_news/downloading-germline-sets-from-the-command-line-or-api/) that can be fed into the workflow for creating reference databases for use with `IgBLAST` as described [here](docs/igblast_setup.md).

### Rhapsody V(D)J

Rhapsody data is already in a very similar tab-delimited format to what is acquired from post-processing the 10x contigs. Some brief notes [here](docs/rhapsody_vdj.md).

