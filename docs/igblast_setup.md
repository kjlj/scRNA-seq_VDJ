# Setting up IgBLAST 

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

## VDJ reference sets for IgBLAST

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

## constant region references 

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

[Return to main page.](../README.md)

[Go to Running IgBLAST.](running_igblast.md)