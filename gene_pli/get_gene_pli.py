# Gets pLI of genes using gnomAD v2.1.1. Fills in missing pLIs with mean pLI
# across genes from gnomAD.
#
# Author: Karthik Guruvayurappan

import numpy as np
import pandas as pd

import argparse

# set up argparse to read in command line arguments
parser = argparse.ArgumentParser(description="Input file and output file paths")
parser.add_argument("e2g_universe", help="Path to E2G universe")
parser.add_argument("gnomad_file", help="Path to gnomAD gene-level info")
parser.add_argument("output_path", help="Output path for E2G universe file with feature added")
args = parser.parse_args()

# read in E2G pair universe
e2g_universe = pd.read_csv(
    args.e2g_universe,
    sep = '\t'
)

# read in gene pLI info
gene_info = pd.read_table(
    args.gnomad_file
)

# filter for gene name and pLI info
gene_info = gene_info[['gene', 'pLI']]

# rename columns for merge
gene_info.columns = ['GeneSymbol', 'pLI']

# get mean pLI (for imputation)
mean_pli = gene_info['pLI'].mean()

# merge pLI scores with E2G universe
e2g_universe_pli = e2g_universe.merge(
    gene_info,
    how = 'left'
)

# impute missing pLIs with mean pLI across all genes
e2g_universe_pli['pLI'] = e2g_universe_pli['pLI'].fillna(mean_pli)

# write to output file
e2g_universe_pli.to_csv(
    args.output_path,
    index = False
)

