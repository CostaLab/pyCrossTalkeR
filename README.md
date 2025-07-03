# pyCrossTalkeR

<img src="https://github.com/CostaLab/pyCrossTalkeR/blob/main/logo1.png" align="right" width="200" />

James S. Nagai<sup>1</sup>,
Vanessa Kloeker<sup>1</sup>,
Ruthvik Koppala<sup>1</sup>,
Nils B. Leimkühler<sup>2</sup>,
Michael T. Schaub <sup>3</sup>,
Rebekka K. Schneider<sup>4,5,6</sup>,
Ivan G. Costa<sup>1*</sup>

<sup>1</sup>Institute for Computational Genomics, Faculty of Medicine, RWTH Aachen University, Aachen, 52074 Germany

<sup>2</sup>Department of Hematology and Stem Cell Transplantation, University Hospital Essen, Germany

<sup>3</sup>Department of Computer Science, RWTH Aachen University, Germany

<sup>4</sup>Department of Cell Biology, Institute for Biomedical Engineering, Faculty of Medicine,RWTH Aachen University, Pauwelsstrasse 30, 52074 Aachen, NRW, Germany

<sup>5</sup>Oncode Institute, Erasmus Medical Center, Rotterdam, 3015GD, the Netherlands

<sup>6</sup>Department of Hematology, Erasmus Medical Center, Rotterdam, 3015GD, the Netherlands

![Alt text](image.png)

**Motivation:** Ligand-receptor (LR) analysis allows the characterization of cellular crosstalk from single cell RNA-seq data. However, current LR methods provide limited approaches for prioritization of cell types, ligands or receptors or characterizing changes in crosstalk between two biological conditions.

**Results:** pyCrossTalkeR is a framework for network analysis and visualisation of LR networks. pyCrossTalkeR identifies relevant ligands, receptors and cell types contributing to changes in cell communication when contrasting two biological states: disease vs. homeostasis. A case study on scRNA-seq of human myeloproliferative neoplasms reinforces the strengths of pyCrossTalkeR for characterisation of changes in cellular crosstalk in disease state.

## Install

You can install pyCrossTalkeR with the simple commands below:

```{python}
pip install pycrosstalker
```

*Note: Please avoid to use the following characters in celltype name: '$'*

## Possible system dependencies

```
libudunits2-dev
libgdal-dev
gdal-bin
libproj-dev
proj-data
proj-bin
libgeos-dev
```
  

## Features v2.0.0


- Single and Comparative Reports
   - Cell Cell Interaction visualization
   - Sending and Receiving Cells Ranking
   - CCI and GCI PCA ranking
      - All measures and PC table
      - PC1 and PC2 based barplot
   - LR pair visualization plot can be done


# References

[1] CrossTalkeR: Analysis and Visualisation of Ligand Receptor Networks [link](https://doi.org/10.1093/bioinformatics/btab370)

[2] Heterogeneous bone-marrow stromal progenitors drive myelofibrosis via a druggable alarmin axis. [link](https://www.cell.com/cell-stem-cell/fulltext/S1934-5909(20)30542-7#secsectitle0115)

[3] Comparison of Resources and Methods to infer Cell-Cell Communication from Single-cell RNA Data [link](https://www.biorxiv.org/content/10.1101/2021.05.21.445160v1.full)
