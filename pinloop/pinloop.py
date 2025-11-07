# ANDREW ROJNUCKARIN 11.07.25 AROJNUC1@JH.EDU AROJNUCK@GMAIL.COM

"""
Currently Pinloop is cell-type agnostic (assuming most CTCF binding is constitutive), a representative set of ChIA-PET data was downloaded from the ENCODE portal (ENCFF377RDA.bedpe ENCFF519OAV.bedpe ENCFF743ZWY.bedpe) and aggregated

Luo, R., Yan, J., Oh, J.W. et al. Dynamic network-guided CRISPRi screen identifies CTCF-loop-constrained nonlinear enhancer gene regulatory activity during cell state transitions. Nat Genet 55, 1336–1346 (2023). https://doi.org/10.1038/s41588-023-01450-7

"""

import numpy as np
import pandas as pd
import sys
import os

print("reading ChIA-PET data . . .")
chiapet = pd.read_table("ChIA-PET_final.tsv")
print("reading CRE-gene pairs . . .\n")
pair_info = pd.read_table(sys.argv[1])
pair_info["Pinloop"] = 0.0


def alloc_pinloop(df_rna, df_atac, df_chiapet, loc_keys=["GeneTSS", "center"]):
    df_atac["center"] = (df_atac["ElementStart"] + df_atac["ElementEnd"]) / 2
    loops = df_chiapet.assign(**{"atac_start": 0, "atac_end": 0, "rnaseq_start": 0, "rnaseq_end": 0})
        
    if len(df_atac.index) != 0:
        loops.loc[df_chiapet.index, "atac_start"] = np.searchsorted(
            df_atac[loc_keys[1]], 
            df_chiapet["center1"]
        ) + df_atac.index[0]

        loops.loc[df_chiapet.index, "atac_end"] = np.searchsorted(
            df_atac[loc_keys[1]], 
            df_chiapet["center2"]
        ) + df_atac.index[0]
    
    if len(df_rna.index) != 0:
        loops.loc[df_chiapet.index, "rnaseq_start"] = np.searchsorted(
            df_rna[loc_keys[0]], 
            df_chiapet["center1"]
        ) + df_rna.index[0]

        loops.loc[df_chiapet.index, "rnaseq_end"] = np.searchsorted(
            df_rna[loc_keys[0]], 
            df_chiapet["center2"]
        ) + df_rna.index[0]
    return loops


def p_inloop_dense(prm_ids, enh_ids, df_chiapet):
    prm_id_dict = dict(zip(np.append(prm_ids, np.max(prm_ids)+1), np.arange(len(prm_ids)+1)))
    enh_id_dict = dict(zip(np.append(enh_ids, np.max(enh_ids)+1), np.arange(len(enh_ids)+1)))
    
    p_inloop = np.zeros(shape=(len(prm_ids), len(enh_ids)), dtype=float)
    prm_inloop = np.zeros(shape=len(prm_ids), dtype=float)
    for idx, loop in df_chiapet.iterrows():
        prm_start = prm_id_dict.get(loop["rnaseq_start"])
        prm_end = prm_id_dict.get(loop["rnaseq_end"])
        enh_start = enh_id_dict.get(loop["atac_start"])
        enh_end = enh_id_dict.get(loop["atac_end"])
        
        prm_inloop[prm_start:prm_end] += loop["count"]
        p_inloop[prm_start:prm_end, enh_start:enh_end] += loop["count"]

    return np.divide(p_inloop, prm_inloop[:, None], out=np.zeros_like(p_inloop), where=prm_inloop[:, None] != 0)


for chrm in pair_info["ElementChr"].unique():
    print("processing " + chrm + " . . .", flush=True)
    pairs_on_chr = pair_info.loc[pair_info["ElementChr"] == chrm, :]
    genes_on_chr = pairs_on_chr.loc[:, ["GeneEnsemblID", "GeneSymbol", "GeneTSS"]].drop_duplicates().reset_index(drop=True).copy()
    elements_on_chr = pairs_on_chr.loc[:, ["ElementName", "ElementStart", "ElementEnd"]].drop_duplicates().reset_index(drop=True).copy()
    chiapet_on_chr = chiapet.loc[chiapet["chr"] == chrm, :].reset_index(drop=True).copy()

    print("\tsorting loops, elements and genes . . .", flush=True)
    alloc_loops = alloc_pinloop(genes_on_chr, elements_on_chr, chiapet_on_chr.loc[chiapet_on_chr["loop_length"] < 1e6, :].reset_index(drop=True))
    print("\tcalculating Pinloop . . .", flush=True)
    pinloop = p_inloop_dense(genes_on_chr.index, elements_on_chr.index, alloc_loops).T
    lookup = pd.DataFrame(data=pinloop, index=elements_on_chr["ElementName"].to_numpy(), columns=genes_on_chr["GeneEnsemblID"].to_numpy())
    print("\tassigning Pinloop to CRE-gene pairs . . .\n", flush=True)
    pair_info.loc[pairs_on_chr.index, "Pinloop"] = [lookup.loc[r, c] for r, c in zip(pairs_on_chr["ElementName"].to_numpy(), pairs_on_chr["GeneEnsemblID"].to_numpy())]

print("writing to file " + sys.argv[1].split(".")[0] + "_pinloop.tsv.gz . . .")
pair_info.to_csv(sys.argv[1].split(".")[0] + "_pinloop.tsv.gz", sep="\t", index=False)