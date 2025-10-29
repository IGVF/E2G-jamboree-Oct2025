# ANDREW ROJNUCKARIN 10.23.25 AROJNUC1@JH.EDU AROJNUCK@GMAIL.COM

import pandas as pd
import sys

print("reading ENCODE RNA-seq data from tables . . .")
ENCODE_rna = pd.read_table("ENCODE_RNA.tsv.gz")
ENCODE_info = pd.read_table("ENCODE_RNA_info.tsv")

pair_info = pd.read_table(sys.argv[1])

print("coalesce ENCODE RNA-seq data by biosample name . . .")
sample_dict = {x: [] for x in ENCODE_info["Biosample.term.name"].unique()}
for i in ENCODE_info.index:
    sample_dict[ENCODE_info.loc[i, "Biosample.term.name"]].append(ENCODE_info.loc[i, "tsv"])

ENCODE_rna_collapsed = ENCODE_rna.loc[:, ["gene_id"]].copy()

for i, sample in enumerate(list(sample_dict.keys())):
    ENCODE_rna_collapsed[sample] = ENCODE_rna.loc[:, sample_dict[sample]].mean(axis=1)
    if i % 50 == 0:
        ENCODE_rna_collapsed = ENCODE_rna_collapsed.copy()

ENCODE_rna_collapsed.iloc[:, 1:] = ENCODE_rna_collapsed.iloc[:, 1:] / ENCODE_rna_collapsed.iloc[:, 1:].mean(axis=0)

print("calculate ENCODE RNA-seq data statistics . . .")
ENCODE_rna_collapsed["std"] = ENCODE_rna_collapsed.iloc[:, 1:].std(axis=1)
ENCODE_rna_collapsed["mean"] = ENCODE_rna_collapsed.iloc[:, 1:].mean(axis=1)
ENCODE_rna_collapsed["disp"] = ENCODE_rna_collapsed["std"] / (ENCODE_rna_collapsed["mean"] + 1)

ENCODE_rna_collapsed["gene_id"] = [x.split(".")[0] for x in ENCODE_rna_collapsed["gene_id"]]

drop_id = []
for i, group in ENCODE_rna_collapsed.groupby("gene_id"):
    if len(group.index) > 1:
        drop_id.append(i)

ENCODE_rna_collapsed = ENCODE_rna_collapsed.loc[~ENCODE_rna_collapsed["gene_id"].isin(drop_id), :]
ENCODE_rna_collapsed = ENCODE_rna_collapsed.set_index(ENCODE_rna_collapsed["gene_id"].to_numpy())

pair_info["std"] = 1.0
pair_info["mean"] = 1.0
pair_info["disp"] = 0.5

pair_info.loc[pair_info["GeneEnsemblID"].isin(ENCODE_rna_collapsed.index), ["std", "mean", "disp"]] = ENCODE_rna_collapsed.loc[pair_info.loc[pair_info["GeneEnsemblID"].isin(ENCODE_rna_collapsed.index), "GeneEnsemblID"].to_numpy(), ["std", "mean", "disp"]].to_numpy()

print("writing output file to " + sys.argv[1].split(".")[0] + "_gene_ENCODE_stats.tsv.gz . . .")
pair_info.to_csv(sys.argv[1].split(".")[0] + "_gene_ENCODE_stats.tsv.gz", sep="\t", index=False)