# Pinloop

`pinloop.py`

Dependencies: Pandas

Pinloop is a metric that approximates the probability that a CRE falls in the same CTCF-delineated TAD as a gene's TSS. Pinloop uses ChIA-PET data and is calculated for a TSS,CRE pair by taking the sum of loop counts enclosing both the TSS and CRE and dividing by the loop counts enclosing the TSS. ChIA-PET loops are filtered to only include loops spanning less than 1Mb. Currently Pinloop is cell-type agnostic (assumes most CTCF binding is constitutive). A representative set of ChIA-PET data was downloaded from the ENCODE portal (ENCFF377RDA.bedpe ENCFF519OAV.bedpe ENCFF743ZWY.bedpe) and aggregated in the provided file `ChIA-PET_final_filtered.tsv`. Pinloop was shown to be more informative than ABC and distance in predicting CRISPRi hits at CREs around genes encoding regulators of the embryonic stem cell to definitive endoderm differentiation in 2023.

Luo, R., Yan, J., Oh, J.W. et al. Dynamic network-guided CRISPRi screen identifies CTCF-loop-constrained nonlinear enhancer gene regulatory activity during cell state transitions. Nat Genet 55, 1336–1346 (2023). https://doi.org/10.1038/s41588-023-01450-7

Usage: 
```
python3 pinloop.py fname.tsv.gz
```

