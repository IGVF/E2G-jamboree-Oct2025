1. Use the .yaml file to make the conda environment.
2. The python script takes 3 arguments:
- File path to input E2G universe
- File path to GeneBayes supplementary tables (https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-024-01820-9/MediaObjects/41588_2024_1820_MOESM4_ESM.xlsx)
- File path to output file with Shet (should end with .tsv.gz)

Citation: (https://www.nature.com/articles/s41588-024-01820-9). For genes with missing Shet, I am imputing the mean.

The GeneBayes info should be an Excel table. I downloaded it from the paper directly.

