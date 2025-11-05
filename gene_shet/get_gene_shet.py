import numpy as np
import pandas as pd

import argparse

# set up argparse to read in command line arguments
parser = argparse.ArgumentParser(description="Input file and output file paths")
parser.add_argument("e2g_universe", help="Path to E2G universe")
parser.add_argument("genebayes_table", help="Path to supplementary tables from the GeneBayes paper")
parser.add_argument("output_path", help="Output path for E2G universe file with feature added")
args = parser.parse_args()

# read in E2G universe
e2g_universe = pd.read_csv(
    args.e2g_universe,
    sep = '\t'
)

# read in GeneBayes supplementary table
genebayes_results = pd.read_excel(
    args.genebayes_table,
    sheet_name = 'Supplementary Table 1'
)

# filter for gene ENSG and posterior mean (point estimate)
genebayes_results = genebayes_results[['ensg', 'post_mean']]

# rename columns before merge
genebayes_results.columns = ['GeneEnsemblID', 'genebayes_shet']

# merge
e2g_universe = e2g_universe.merge(
    genebayes_results,
    how = 'left'
)

# get mean Shet
mean_shet = genebayes_results['genebayes_shet'].mean()

# impute missing values using the mean Shet
e2g_universe['genebayes_shet'] = e2g_universe['genebayes_shet'].fillna(mean_shet)

# write to output file
e2g_universe.to_csv(
    args.output_path,
    sep = '\t',
    index = False
)



