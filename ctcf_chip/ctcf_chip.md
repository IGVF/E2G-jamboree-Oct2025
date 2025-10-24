# CTCF_chip

## Introduction
This script takes two file types: the E2G candidate enhancer file and CTCF Chip-seq file and annotate the candidate enhancers as overlapping with CTCF peak (1) or non-overlapping (0). The output is a zgipped tab-delimited file in the IGVF E2G format.

## Example usage
python e2g_ctcf_chip.py -c CharacterizationMcGinnis_Dataset1_K562_candidate_e2g_pairs.tsv -t ENCFF406IQL.bed -o K562_enhancers_CTCF.bed.gz

## Input type
- enhancer candidate file: tab-delimited file
- CTCF peak file: tab-delimited file (such as ENCODE bed files)

## Output
- One column "Score" is added to the original enhancer candidate file, showing each enhancer element as "overlapping with CTCF" (1) or 'Non-overlapping" (0)