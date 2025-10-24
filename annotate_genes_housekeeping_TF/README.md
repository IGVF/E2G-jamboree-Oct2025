# Annotating genes as housekeeping genes or TF-encoding genes


## Measure of a housekeeping genes using ENCODE RNA-Seq

`label_gene_ENCODE_stats.py`

Dependencies: Pandas

Calculates gene statistics from ENCODE RNA-seq data grouped by biosample type. Calculates `mean`, `std`, and `disp = mean / (std + 1)` mean normalized TPM for each gene. Genes not labeled in table (not many) are set to default values 1, 1, 0.5 respectively. Requires ENCODE RNA-seq data saved in `ENCODE_RNA.tsv.gz` and `ENCODE_RNA_info.tsv`. The idea is that housekeeping genes (low `std`, low `disp`) may not be effected by perturbations/variants of nearby regulatory elements. Thus, including this information may help avoid false positives.

Usage: 
```
python3 label_gene_ENCODE_stats.py fname.tsv.gz
```


## Labeling TF-encoding genes

`label_TFs.py`

Dependencies: Pandas

Labels genes as TF (1) or not TF (0). TFs are taken from list curated by Lambert et al. Cell 2018 https://pubmed.ncbi.nlm.nih.gov/29425488/. The idea is that genes encoding TFs may correlate with peaks because the TF binds to the peak rather than the peak regulating the TF-encoding gene. Thus, it may reasonable to remove TF-encoding genes when using correlation-based features or mask correlation-based features for TF-encoding genes.

Usage: 
```
python3 label_TFs.py fname.tsv.gz
```
