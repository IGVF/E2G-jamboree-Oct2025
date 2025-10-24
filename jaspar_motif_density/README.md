# Motif density features using PWMs from Jaspar2024 core non-redundant set

Required input files:
- E2G pair universe file
- [JASPAR hg38 annotations bigBed file](https://hgdownload.soe.ucsc.edu/gbdb/hg38/jaspar/JASPAR2024.bb)

Instructions for generating features:

```
conda env create -f generate_motif_densities_jaspar.yml

conda activate generate_motif_densities_jaspar

python generate_motif_densities_jaspar.py --e2g_pairs [PAIRS_FILE] --jaspar-bigbed JASPAR2024.bb --output-file [OUTPUT_FILE]
```

Feature columns:
| Feature name | Description |
| ------------ | ----------- |
| MotifDensityJaspar2024 | Number of total PWM hits with pvalue<10^-5 divided by element length |
| UniqueMotifDensityJaspar2024 | Number of unique PWM hits with pvalue<10^-5 divided by element length |
| CTCFMotifHitJaspar2024 | Binary presence/absence of at least one CTCF PWM hit with pvalue<10^-5 within element |