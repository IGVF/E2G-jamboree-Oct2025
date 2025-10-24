1. Use the .yaml file to make the conda environment.
2. The python script takes 3 arguments:
- File path to input E2G universe
- File path to gnomAD gene-level info (https://gnomad.broadinstitute.org/data#v2-constraint)
- File path to output file with pLI added

I am using gnomAD v2.1.1 scores. Citation: (https://www.nature.com/articles/s41586-020-2308-7). For genes with missing pLI, I am imputing the mean.

