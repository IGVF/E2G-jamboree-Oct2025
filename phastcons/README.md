1. Use the .yaml file to make the conda environment.
2. The python script takes 3 arguments:
- File path to input E2G universe
- File path to UCSC phastCons bigWig (https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.phastCons100way.bw)
- File path to output file with mean phastCons score across enhancer added as features

Citation: https://genome.cshlp.org/content/15/8/1034. For some regions, it was not possible to compute a mean phastCons. For these regions, I imputed the mean of mean phastCons scores. These were a very small number ~400.

