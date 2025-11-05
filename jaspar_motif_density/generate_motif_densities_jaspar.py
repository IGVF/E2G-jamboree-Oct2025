# Generate motif densities for GM12878 and K562 Multiome E-G pair universes using Jaspar 2024 set of motifs

import argparse
from pathlib import Path
import pandas as pd
import pyBigWig

def get_motif_features(motif_bigbed: pyBigWig, chrom, start, end, TF_name='CTCF'):
    hits = motif_bigbed.entries(chrom, start, end) or []
    filtered_hits = [item for item in hits if int(item[2].split('\t')[1]) >= 500]
    n_hits = len(filtered_hits)
    motif_density = n_hits / (end - start)

    uniq_motif_ids = list(set(item[2].split('\t')[0] for item in filtered_hits))
    uniq_n_hits = len(uniq_motif_ids)
    uniq_motif_density = uniq_n_hits / (end - start)

    ctcf_motif_hit = int(any(TF_name in item[2].split('\t')[3] for item in filtered_hits))

    return motif_density, n_hits, uniq_motif_density, uniq_n_hits, ctcf_motif_hit

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Generate motif densities for E-G pair universes.")
    parser.add_argument('--e2g-pairs', type=str, help="Path to E2G pair universe file.", required=True)
    parser.add_argument('--jaspar-bigbed', type=str, help="Path to JASPAR motif bigbed file.", required=True)
    parser.add_argument('--output-file', type=str, help="Output file path/name for processed E2G pair universe with motif features.", required=True)
    args = parser.parse_args()

    # Load E2G pair tables
    pairs_df = pd.read_table(args.e2g_pairs, compression='gzip')

    # Get unique elements
    uniq_elements_df = pairs_df[['ElementChr', 'ElementStart', 'ElementEnd']].drop_duplicates()

    # Load JASPAR2024 motif bigbed
    jaspar_bb = pyBigWig.open('data/raw/JASPAR2024.bb')

    motif_features_list = ['MotifDensityJaspar2024', 'MotifCountsJaspar2024', 'UniqueMotifDensityJaspar2024', 'UniqueMotifCountsJaspar2024', 'CTCFMotifHitJaspar2024']

    # Calculate motif densities
    uniq_elements_df[motif_features_list] = pd.DataFrame(
        uniq_elements_df.apply(
            lambda row: get_motif_features(jaspar_bb, row['ElementChr'], row['ElementStart'], row['ElementEnd'], 'CTCF'), axis=1).to_list(),
            index=uniq_elements_df.index)
    processed_pairs_df = pairs_df.merge(uniq_elements_df, on=['ElementChr', 'ElementStart', 'ElementEnd'], how='left')

    processed_pairs_df.to_csv(args.output_file, sep='\t', index=False, compression='gzip')