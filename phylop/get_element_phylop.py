import pyBigWig, pandas as pd
from tqdm import tqdm

import argparse

# set up argparse to read in command line arguments
parser = argparse.ArgumentParser(description="Input file and output file paths")
parser.add_argument("e2g_universe", help="Path to E2G universe")
parser.add_argument("phylop_file", help="bigWig file from Zoonomia with phyloP information")
parser.add_argument("output_path", help="Output path for E2G universe file with features added")
args = parser.parse_args()

# read in E2G universe
e2g_universe = pd.read_csv(
    args.e2g_universe,
    sep = '\t'
)

# get unique enhancers
enhancers = e2g_universe[['ElementChr', 'ElementStart', 'ElementEnd']].drop_duplicates()

# read in bigWig file containing the Zoonomia phyloP scores
bw = pyBigWig.open(args.phylop_file)

# get minimum and maximum PhyloP across each unique enhancer
enhancers["min_phyloP"] = [bw.stats(c, s, e, type="min")[0]  for c,s,e in tqdm(enhancers[["ElementChr","ElementStart","ElementEnd"]].to_numpy())]
enhancers["max_phyloP"] = [bw.stats(c, s, e, type="max")[0]  for c,s,e in tqdm(enhancers[["ElementChr","ElementStart","ElementEnd"]].to_numpy())]

# close the bigWig connection
bw.close()

# get mean of min and max phyloP (for imputation)
mean_max_phyloP = enhancers['max_phyloP'].mean()
mean_min_phyloP = enhancers['min_phyloP'].mean()

# impute missing values with the mean
enhancers['max_phyloP'] = enhancers['max_phyloP'].fillna(mean_max_phyloP)
enhancers['min_phyloP'] = enhancers['min_phyloP'].fillna(mean_min_phyloP)

# merge with E2G universe
e2g_universe = e2g_universe.merge(
    enhancers,
    how = 'left'
)

# write E2G universe with added features to output file
e2g_universe.to_csv(
    args.output_path,
    sep = '\t',
    index = False
)

