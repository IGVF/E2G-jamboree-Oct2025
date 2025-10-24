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

# read in gene constraint info from gnomAD
gene_info = pd.read_table(
    args.gnomad_file
)

# filter for gene name and pLI and LOEUF info
gene_info = gene_info[['gene', 'pLI', 'oe_lof_upper']]

# rename columns for merge
gene_info.columns = ['GeneSymbol', 'pLI', 'LOEUF']

# get mean pLI and LOEUF (for imputation)
mean_pli = gene_info['pLI'].mean()
mean_loeuf = gene_info['LOEUF'].mean()

# merge pLI and LOEUF scores with E2G universe
e2g_universe_merged = e2g_universe.merge(
    gene_info,
    how = 'left'
)

# impute missing pLIs and LOEUFs with mean pLI and LOEUFs across all genes
e2g_universe_merged['pLI'] = e2g_universe_merged['pLI'].fillna(mean_pli)
e2g_universe_merged['LOEUF'] = e2g_universe_merged['LOEUF'].fillna(mean_loeuf)

# write to output file
e2g_universe_merged.to_csv(
    args.output_path,
    index = False
)

