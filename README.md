# scRNA-seq workshop

## IgBLAST post-processing

In the interest of time, we will not undertake the IgBLAST post-processing during this sessions. Rather we will start with datasets that have already been processed. All steps for generating the datasets that we will be working with are documented.

Step one - [obtaining datasets from 10x genomics](datasets.md)

## post-processing

Cell ranger output does not include allele level gene assignment and doesn't provide any information about identity to the germline that is needed for calculating somatic hypermutation (SHM) levels for the B cell datasets. To obtain this extra level of detail the contigs from 10x VDJ are re-aligned against the [IMGT reference directory](https://www.imgt.org/vquest/refseqh.html) using [IgBLAST](https://www.ncbi.nlm.nih.gov/igblast/igblast.cgi). To get tab-delimited AIRR-C format for the VDJ alignments we will use the stand alone version of IgBLAST rather than the [web-based version](https://ncbi.github.io/igblast/).

### setting up IgBLAST 

Instructions for installation are available [here](https://ncbi.github.io/igblast/cook/How-to-set-up.html). At time of writing the current version is v1.21. Installation packages for various operating systems can be found at [https://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/](https://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/).

For OsX to obtain stand alone IgBLAST software:
```
cd ~/data/apps/
#downloading the file
wget https://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/ncbi-igblast-1.21.0-x64-macosx.tar.gz

#or, for curl
#curl -O https://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/ncbi-igblast-1.21.0-x64-macosx.tar.gz

#uncompress and untar
tar xzvf ncbi-igblast-1.21.0-x64-macosx.tar.gz

#remove the downloaded compressed tar
rm ncbi-igblast-1.21.0-x64-macosx.tar.gz
```

### VDJ reference sets for IgBLAST

To obtain the IMGT Reference Directory to use with IgBLAST. For each chain gene type (V, D, J) need to create a single file for the Ig and another for the TR that combines all the genes from the different loci. For example, for Ig V genes, need to collect all the V genes from IGH, IGK and IGL together into a single fasta file. Similarly, for the TR the V genes from TRA, TRB, TRG and TRD need to be compiled into a file. The description lines for the files obtained from IMGT need to be edited to have the gene name proximal to the '>'. IgBLAST provides a script for this. The example below is for human.

First, download the fasta files from IMGT:
```
#nagivate to the folder where IgBLAST is located
cd ~/data/apps/ncbi-igblast-1.21.0/

#create a new directory to store the references and navigate to it
mkdir references/
cd references/

#download the Ig gene files from IMGT
## V genes
curl -O https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGHV.fasta
curl -O https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGKV.fasta
curl -O https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGLV.fasta

## D genes, only for the IGH
curl -O https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGHD.fasta

## J genes
curl -O https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGHJ.fasta
curl -O https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGKJ.fasta
curl -O https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGLJ.fasta

#download the TR gene files from IMGT
## V genes
curl -O https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TRAV.fasta
curl -O https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TRBV.fasta
curl -O https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TRGV.fasta
curl -O https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TRDV.fasta

## D genes, only for TRB and TRG
curl -O https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TRBD.fasta
curl -O https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TRDD.fasta

## J genes
curl -O https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TRAJ.fasta
curl -O https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TRBJ.fasta
curl -O https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TRGJ.fasta
curl -O https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/TR/TRDJ.fasta
```

The TRAV and TRDV genes are interspered within the some locus. A subset of the genes are able to rearranged with both TRAJ and TRDJ genes are therefore considered part of both and have names such as TRAV14/DV4\*01 meaning that this gene segement is TRAV14 and TRDV4. These genes will be dupicated when the TR V genes are later collected into a single file. BLAST databases cannot included duplicated genes, therefore these need to be removed prior to this step. There are many options for doing this, but the simplest is to open the file in a suitable text editor and manually delete:
```
#navigate the the directory with the downloads
cd ~/data/apps/ncbi-igblast-1.21.0/references/

#remove TRAV/DV entries from the TRDV file, they will still be included in the final reference from the entries in the TRAV file
open -a TextEdit TRDV.fasta
# remove the entries for the TRAV14/DV4*01, TRAV14/DV4*02, TRAV14/DV4*03, TRAV14/DV4*04, TRAV23/DV6*01, TRAV23/DV6*02
# TRAV23/DV6*03, TRAV23/DV6*04, TRAV29/DV5*01, TRAV29/DV5*02, TRAV36/DV7*01, TRAV36/DV7*02, TRAV36/DV7*03, TRAV36/DV7*04, TRAV38-2/DV8*01
# save the file
```

There is also the option to exclude the TRD and TRG as these are not amplified by the 10x primers (only targets the TRB and TRA constant regions), but best to have a complete reference set that can be used for all analysis. 

Next, compile each gene into a single file for Ig and TR. This uses the cat command which prints the contents of a file. By using this command with a wildcard (the \* means 'any character' so IG\*V will apply the cat command to the IGKV, IGHV and IGLV) and directing the output to a new file (cat by default prints to the screen, but adding the > symbol redirects this to the filename provided) is telling we can print the contents of the files to a single new file.
```
#navigate to directory where the references were downloaded
cd ~/data/apps/ncbi-igblast-1.21.0/references/

#the Ig genes
cat IG*V.fasta > imgt_ig_v_human.fasta
## not strictly necessary as only a single file but good to keep the naming conventions
cat IG*D.fasta > imgt_ig_d_human.fasta
cat IG*J.fasta > imgt_ig_j_human.fasta

#the TR genes
cat TR*V.fasta > imgt_tr_v_human.fasta
cat TR*D.fasta > imgt_tr_d_human.fasta
cat TR*J.fasta > imgt_tr_j_human.fasta
```

Finally, to reformat the fasta description line to the format needed for IgBLAST:
```
cd ~/data/apps/ncbi-igblast-1.21.0/references/

#use the perl script provided by IgBLAST to do the reformatting
../bin/edit_imgt_file.pl imgt_ig_v_human.fasta > imgt_ig_v_human.fa
../bin/edit_imgt_file.pl imgt_ig_d_human.fasta > imgt_ig_d_human.fa
../bin/edit_imgt_file.pl imgt_ig_j_human.fasta > imgt_ig_j_human.fa

../bin/edit_imgt_file.pl imgt_tr_v_human.fasta > imgt_tr_v_human.fa
../bin/edit_imgt_file.pl imgt_tr_d_human.fasta > imgt_tr_d_human.fa
../bin/edit_imgt_file.pl imgt_tr_j_human.fasta > imgt_tr_j_human.fa

#or, one line version
#for f in *imgt*human.fasta; do echo $f; ../bin/edit_imgt_file.pl $f > ${f%.fasta}.fa; done;
```

These fasta files now need to be converted to blast database file for use with IgBLAST:
```
#navigate to the directory with the FASTA formatted gene sets
cd ~/data/apps/ncbi-igblast-1.21.0/references/

#for each gene set run the command to generate a BLAST database
../bin/makeblastdb -parse_seqids -dbtype nucl -in imgt_ig_v_human.fa
../bin/makeblastdb -parse_seqids -dbtype nucl -in imgt_ig_d_human.fa
../bin/makeblastdb -parse_seqids -dbtype nucl -in imgt_ig_j_human.fa

../bin/makeblastdb -parse_seqids -dbtype nucl -in imgt_tr_v_human.fa
../bin/makeblastdb -parse_seqids -dbtype nucl -in imgt_tr_d_human.fa
../bin/makeblastdb -parse_seqids -dbtype nucl -in imgt_tr_j_human.fa

#or, one line version
#for f in imgt_*_human.fa; do echo $f; ../bin/makeblastdb -parse_seqids -dbtype nucl -in $f; done;
```

Optionally, can clean up the directory to remove the intermediate files and keep only the BLAST databases that are need for IgBLAST:
```
cd ~/data/apps/ncbi-igblast-1.21.0/references/

#remove the files that end with .fa
rm *.fa
#remove the files that end with .fasta
rm *.fasta
```

### constant region references 

Newer versions of IgBLAST align the constant region in addition to the V(D)J. The primers for 10x V(D)J amplicons sit within the first exon of each of the constant regions, for Ig this is the CH1 region and for TR this is termed EX1. Unlike the V(D)J references, which are available as FASTA downloads from IMGT [https://www.imgt.org/vquest/refseqh.html], the constant regions are in HTML format and require cut&paste to collect the sequences of interest into a FASTA file. 

For the Ig can obtain blast database formatted human constant region gene sets from NCBI:
```
#navigate to the reference directory
cd ~/data/apps/ncbi-igblast-1.21.0/references/

#download the tar from NCBI
curl -O https://ftp.ncbi.nih.gov/blast/executables/igblast/release/database/ncbi_human_c_genes.tar

#or, using wget
#wget https://ftp.ncbi.nih.gov/blast/executables/igblast/release/database/ncbi_human_c_genes.tar

#untar
tar xvf ncbi_human_c_genes.tar

#remove the tar file
rm ncbi_human_c_genes.tar
```

Note, these are already in BLAST database format so ready to use with IgBLAST.

This database doesn't include the TR constant regions, so these need to be obtained from IMGT using the cut&paste method. 
```
#navigate to the references directory
cd ~/data/apps/ncbi-igblast-1.21.0/references/

#create an empty file
touch imgt_tr_c_human.fa
#open the empty file
open -a TextEdit imgt_tr_c_human.fa

#visit the following URLS to obtain EX1, and manually reformat the description lines to include only the gene name
#eg. change from 
#>X02883|TRAC*01|Homo sapiens|F|EX1|n,273..544|273 nt|1|+1|-1|| |273+0=273| | | 
#to
#>TRAC*01

#https://www.imgt.org/genedb/GENElect?query=7.2+TRAC&species=Homo+sapiens
#https://www.imgt.org/genedb/GENElect?query=7.2+TRBC&species=Homo+sapiens
#https://www.imgt.org/genedb/GENElect?query=7.2+TRGC&species=Homo+sapiens
#https://www.imgt.org/genedb/GENElect?query=7.2+TRDC&species=Homo+sapiens

#convert the FASTA to a blastdb
../bin/makeblastdb -parse_seqids -dbtype nucl -in imgt_tr_c_human.fa
```

### running IgBLAST on contigs

The IgBLAST command for Ig and TR differs.

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

## combining cell ranger and IgBLAST

For each dataset now have the V(D)J data from Cell Ranger (the filtered_contig_annotations.csv files) and IgBLAST (the igblast.tsv files). 

Combining these files in R to summarise the paired chains for each cell. Running R in RStudio on [Posit cloud](https://posit.cloud/). After logging in, can create a new project and within the project a new RScript (joining_cellranger_igblast.R).

Workflow for combining files is described [here](https://kjlj.github.io/scRNA-seq_VDJ/joining_cellranger_igblast.html). The RScript is available on [github](https://raw.githubusercontent.com/kjlj/scRNA-seq_VDJ/main/joining_cellranger_igblast.R) in the [workshop respository](https://github.com/kjlj/scRNA-seq_VDJ).
[

## defining expanded clones/clonotypes 

For the T cells the clonotype labels have already been created by joining the V, J and CDR3AA sequence. Note that this uses the IMGT CDR3 which excludes the 104-Cys and 118-P/W, many reported CDR3 incorporate these residues. To get an equivalent 'cdr3' can use the junction and junction_aa columns reported by IgBLAST. These will include the anchor residues.

As the B cell chains undergo somatic hypermutation (SHM) need to use a method for defining clonal groups that permits grouping of similar but not exact sequences.

### clustering for Ig chains 

Clustering of the Ig sequences is intended to identify groups of cells that carry sequences that are derived from the same original progentitor B cell. These cells will carry heavy and light chains with similar nucleotide sequences that have diverged due to the accumulation of point mutations. Clustering must be undertaken for each individual as a clonal lineage can only exist within a single individual. Similar Igs from different individuals are 'convergent' antibodies likely responding to the same antigen rather than 'divergent' members of the same clone. The same type of clustering can be used to explore convergence (aka public clones) in antibody responses but this is in addition to the clonal lineages.

In the example datasets the pbmc and tumour samples are from different individuals so clones will be built separately. If your 10x dataset includes multiplexed samples you must split the B cell V(D)J data before building clones. 

Generally, heirachical clustering is used on the CDR3 nucleotides sequences to group sequences at a specific threshold. First, the sequences are binned by V, J and CDR3 length as any clone will share these three features. The process within each bin then involves seeding the first cluster with the first sequence, then the next sequence is compared to whether or not it is within an identity threshold, if it is then it is added to the cluster, if not, it becomes the first member of a new cluster. This process is repeated until all sequences are assigned to clusters. Generally a centroid approach is taken, this means that adding a new sequence to a cluster doesn't increase average distance between cluster members beyond a threshold. This aims to prevent 'daisy-chaining', where a sequence may be within the threshold of a single member of a cluster but is increasing distance from other members.

An excellent resource for post-processing of 10x data including building clonal relationships can be found at [https://immcantation.readthedocs.io/en/stable/tutorials/10x_tutorial.html](https://immcantation.readthedocs.io/en/stable/tutorials/10x_tutorial.html). This approach can be undertaken in R, but can also export data and use alternative tools such as [cd-hit](https://github.com/weizhongli/cdhit) but this requires some additional coding to wrangle the cd-hit output into a format that can be joined to the cellranger + IgBLAST cell barcode level information. There is some information available [here](https://github.com/kjlj/10x-VDJ_post-processing).

The workflow for clone building in R is available [here][https://kjlj.github.io/scRNA-seq_VDJ/building_b_cell_clones.html]. This builds upon the Ig data generated from the combining IgBLAST and cellranger workflow. The RScript is available on [github](https://raw.githubusercontent.com/kjlj/scRNA-seq_VDJ/main/building_b_cell_clones.R) in the [workshop respository](https://github.com/kjlj/scRNA-seq_VDJ).


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

## additonal software

### wget

To install [wget](https://formulae.brew.sh/formula/wget0 on osX via [Homebrew](https://brew.sh/). Homebrew is a command line package manager for OsX that can be handy for installing a variety of software such as wget:
```
#installing homebrew
cd ~/
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

#using homebrew to install wget
brew install wget
```
