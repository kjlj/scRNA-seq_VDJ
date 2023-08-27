# VDJ from Rhapsody

For Rhapsody data the `VDJ_Dominant_Contigs_AIRR` tab-delimited file(s) can be used as a starting point for post-processing. 

The Rhapsody output already includes the allele level information for the nearest germline match and appears to use a recent version of the IMGT reference, however it lacks percent identity columns in the output that can be used for SHM levels. The % SHM can be calculated from the available output as the sequence and germline strings are included, but can optionally re-run IgBLAST to get the v_identity column and derive the % SHM from there.

## adding SHM in R

To add SHM within R using the dominant contigs tab-delim file as input:
```
library(tidyverse)

dom.data <- read_tsv("VDJ_Dominant_Contigs_AIRR.tsv")

dom.data <- dom.data %>%
  mutate(v_seq_len = v_sequence_end - v_sequence_start,
         v_sequence = substr(sequence_alignment, 1, v_seq_len),
         v_sequence_germline = substr(germline_alignment, 1, v_seq_len)) %>%
  rowwise() %>%
  mutate(v_mut_cnt = StrDist(v_sequence, v_sequence_germline, method = "hamming")[1]) %>%
  ungroup() %>%
  mutate(v_mut_percent = v_mut_cnt / v_seq_len)

```


## running IgBLAST

To run IgBLAST, it is necessary to extract the sequences from the tab-delimited file and create FASTA file(s) for input to IgBLAST. An example of a [perl script](https://github.com/kjlj/scRNA-seq_VDJ/blob/main/scripts/extract_fasta_frm_contigs.pl) to do this conversion from tsv to FASTA is available on the [github repository](https://github.com/kjlj/scRNA-seq_VDJ). This script will create a separate FASTA file for each Ig/TR locus.

Here is a quick example workflow to get IgBLAST from Rhapsody:
```
cd ~/data/rhapsody/

perl extract_fasta_frm_contigs.pl VDJ_Dominant_Contigs_AIRR.tsv VDJ_Dominant_Contigs_AIRR

#combine to get one Ig file and one TR file for each set
cat VDJ_Dominant_Contigs_AIRR_IG*.fa > VDJ_Dominant_Contigs_AIRR_IG.fasta
cat VDJ_Dominant_Contigs_AIRR_TR*.fa > VDJ_Dominant_Contigs_AIRR_TR.fasta

#check the number of sequences in the joined FASTA files
grep -c ">" *.fasta
#VDJ_Dominant_Contigs_AIRR_IG.fasta:XXXX
#VDJ_Dominant_Contigs_AIRR_TR.fasta:XXXX
```

IgBLAST can then be run on the `VDJ_Dominant_Contigs_AIRR_IG.fasta` and `VDJ_Dominant_Contigs_AIRR_TR.fasta` as described [here](running_igblast.md).

[Return to main page](https://kjlj.github.io/scRNA-seq_VDJ/index.html)
