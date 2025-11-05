# E2G-jamboree-Oct2025
Author: Kilian Salomon (kilian.salomon@bih-charite.de)
date: 2025-10-23*
- Install required packages (bedtools)
    - or using "repeatmasker_fraction_overlap" environment
- Download RepeatMasker annotation for hg38
- I downloaded repeatmasker annotations from UCSC: 2025-10-23
    - downloading will be done automatically in the script (from UCSC)
### About the feature:
- it is an number between 0 and 1 indicating the fraction of the element region that overlaps with Repeatmasker annotated regions
- if no overlap is found, the value is 0
- if no region is found (e.g. element length = 0), the value is 0
### About Repeatmasker:
- Number of different repeat families:
  - all: 15600
  - without: (xyz)n: 1370
  - since no knowledge about the different families I decided to include all families and just compute the overlap with the given region universe and report it

- download RepeatMasker annotation for hg38
```bash
wget -O rmsk.hg38.tsv.gz "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz"
```
- generate a bed file
```bash
zcat rmsk.hg38.tsv.gz | awk 'BEGIN{OFS="\t"} {print $6, $7, $8, $11, $2, $10}' > rmsk.hg38.bed
```
- python command for testing
- download: https://www.synapse.org/Synapse:syn68648234
- change path of e2g data path
```bash
python repeatmasker_fraction_overlap.py   --e2g path/to/test/data   --rmsk rmsk.hg38.bed   --output element_gene_with_repeat_feature.tsv.gz
```
