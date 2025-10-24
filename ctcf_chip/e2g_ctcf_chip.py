import pandas as pd
import argparse
import gzip

# Define your check_overlap function
def check_overlap(ctcf_row, elem_start, elem_end):
    return (elem_start < ctcf_row['CTCFEnd']) and (elem_end > ctcf_row['CTCFStart'])

def main(candidate_file, ctcf_file, output_file):
    # Read input files
    candidate_df = pd.read_csv(candidate_file, sep='\t')
    ctcf_df = pd.read_csv(ctcf_file, sep='\t')

    # Initialize a list to hold results for each chromosome
    results = []

    # Get unique candidate chromosomes
    unique_chromosomes = candidate_df['ElementChr'].unique()

    # Process each chromosome individually
    for chr in unique_chromosomes:
        print(f"Processing chromosome: {chr}")

        # Get candidates for this chromosome
        unique_candidates = candidate_df[candidate_df['ElementChr'] == chr][['ElementChr', 'ElementStart', 'ElementEnd', 'ElementName']].drop_duplicates()

        # Filter CTCF data for this chromosome
        ctcf_chr_df = ctcf_df[ctcf_df['chr'] == chr].rename(columns={'chr': 'ElementChr', 'start': 'CTCFStart', 'end': 'CTCFEnd'})

        # Check for overlaps with merged DataFrame
        merged_df = pd.merge(unique_candidates, ctcf_chr_df, on='ElementChr', how='left')

        # Check overlaps using vectorized operations. Store as 1 for True, 0 for False
        merged_df['Score'] = (
            (merged_df['ElementStart'] < merged_df['CTCFEnd']) &
            (merged_df['ElementEnd'] > merged_df['CTCFStart'])
        ).astype(int)  # Convert boolean to integer (1 for True, 0 for False)

        # Aggregate results to determine if any CTCF overlap exists for each candidate
        score_results = merged_df.groupby('ElementName')['Score'].sum()

        # Format the results DataFrame
        score_results_df = score_results.reset_index().rename(columns={'Score': 'Score'})

        # Merge with original unique_candidates to maintain the overall structure
        unique_candidates = unique_candidates.merge(score_results_df, on='ElementName', how='left')

        # Append results for this chromosome
        results.append(unique_candidates)

    # Concatenate results from all chromosomes
    final_result_df = pd.concat(results, ignore_index=True)

    # Merge back with original candidate_df
    result_df = candidate_df.merge(final_result_df[['ElementName', 'Score']], on='ElementName', how='left')

    # Fill NaN values with 0
    result_df['Score'] = result_df['Score'].fillna(0)

    # Save the final result to a gzipped file
    with gzip.open(output_file, 'wt') as gz_file:
        result_df.to_csv(gz_file, sep='\t', index=False)

    unique_elements = result_df[['ElementName', 'Score']].drop_duplicates()
    ctcf_enh = unique_candidates[unique_candidates['Score'] == 1]

    print(f'Among {len(unique_elements)} candidate enhancer elements, {len(ctcf_enh)} enhancer candidates overlap with CTCF ChIP-seq peaks.')

if __name__ == "__main__":
    # Set up command line argument parsing
    parser = argparse.ArgumentParser(description='Process candidate and CTCF files for overlaps.')
    parser.add_argument('-c', '--candidate_file', type=str, required=True, help='Name of the input candidate file (TSV)')
    parser.add_argument('-t', '--ctcf_file', type=str, required=True, help='Name of the input CTCF file (BED)')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Name of the output gzipped file')

    args = parser.parse_args()
    
    # Call the main function with input arguments
    main(args.candidate_file, args.ctcf_file, args.output_file)