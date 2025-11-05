1. Use the .yaml file to make the conda environment.
2. The python script takes 3 arguments:
- File path to input E2G universe
- File path to Zoonomia phyloP bigWig (https://cgl.gi.ucsc.edu/data/cactus/241-mammalian-2020v2-hub/Homo_sapiens/241-mammalian-2020v2.bigWig)
- File path to output file with max phyloP and min phyloP across enhancer added as features

Citation: https://www.science.org/doi/10.1126/science.abl8189. For some regions, it was not possible to compute a max or min phyloP. For these regions, I imputed the mean max or min phyloP, respectively. These were a very small number ~300.

