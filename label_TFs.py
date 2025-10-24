# ANDREW ROJNUCKARIN 10.23.25 AROJNUC1@JH.EDU AROJNUCK@GMAIL.COM

"""
Lambert SA, Jolma A, Campitelli LF, Das PK, Yin Y, Albu M, Chen X, Taipale J, Hughes TR, Weirauch MT. The Human Transcription Factors. Cell. 2018 Feb 8;172(4):650-665. doi: 10.1016/j.cell.2018.01.029. Erratum in: Cell. 2018 Oct 4;175(2):598-599. doi: 10.1016/j.cell.2018.09.045. PMID: 29425488.
"""

import pandas as pd
import sys

print("reading Lambert Cell 2018 TF Ensembl IDs . . .")

LAMBERT_tfs = pd.read_table("LAMBERT_CELL_2018_hg38_TFs.tsv")

print("reading CRE-gene pairs . . .")

pair_info = pd.read_table(sys.argv[1])

tf_ensembl = [x.split(".")[0] for x in LAMBERT_tfs["gene_id"].to_numpy()]

print("labeling TFs . . .")

pair_info["is_tf"] = 0

pair_info.loc[pair_info["GeneEnsemblID"].isin(tf_ensembl), "is_tf"] = 1

print("writing output file to " + sys.argv[1].split(".")[0] + "_tf_genes_marked.tsv.gz . . .")
pair_info.to_csv(sys.argv[1].split(".")[0] + "_tf_genes_marked.tsv.gz", sep="\t", index=False)