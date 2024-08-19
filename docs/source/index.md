# pyCrossTalkeR

<img src="/logo1.png" align="right" width="200" />

James S. Nagai<sup>1</sup>,
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

**Results:** CrossTalkeP is a framework for network analysis and visualisation of LR networks. CrossTalkeP identifies relevant ligands, receptors and cell types contributing to changes in cell communication when contrasting two biological states: disease vs. homeostasis. A case study on scRNA-seq of human myeloproliferative neoplasms reinforces the strengths of CrossTalkeP for characterisation of changes in cellular crosstalk in disease state.

## Install

You can install CrossTalkeP with the simple commands below:

``{python}
pip install git+https://github.com/CostaLab/pyCrossTalkeR/
```
