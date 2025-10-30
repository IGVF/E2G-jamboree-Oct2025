import pyBigWig, pandas as pd
from tqdm import tqdm

import argparse

# set up argparse to read in command line arguments
parser = argparse.ArgumentParser(description="Input file and output file paths")
parser.add_argument("e2g_universe", help="Path to E2G universe")
parser.add_argument("phastcons_file", help="bigWig file from UCSC with phastCons information")
parser.add_argument("output_path", help="Output path for E2G universe file with features added")
args = parser.parse_args()

# read in E2G universe
e2g_universe = pd.read_csv(
    args.e2g_universe,
    sep = '\t'
)

# get unique enhancers
enhancers = e2g_universe[['ElementChr', 'ElementStart', 'ElementEnd']].drop_duplicates()

# read in bigWig file containing the phastCons scores
bw = pyBigWig.open(args.phastcons_file)

# get mean phastCons score across each enhancer
enhancers["mean_phastCons"] = [bw.stats(c, s, e, type="mean")[0]  for c,s,e in tqdm(enhancers[["ElementChr","ElementStart","ElementEnd"]].to_numpy())]

# close the bigWig connection
bw.close()

# get mean of the phastCons means for imputation
mean_mean_phastCons = enhancers['mean_phastCons'].mean()

# impute missing values with the mean
enhancers['mean_phastCons'] = enhancers['mean_phastCons'].fillna(mean_mean_phastCons)

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

