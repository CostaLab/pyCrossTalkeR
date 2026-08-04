# pyCrossTalkeR

<img src="https://raw.githubusercontent.com/CostaLab/pyCrossTalkeR/main/docs/source/logo1.png" align="right" width="200" />

James S. Nagai<sup>1</sup>,
Vanessa Kloeker<sup>1</sup>,
Ruthvik Koppala<sup>1</sup>,
Nils B. Leimkühler<sup>2</sup>,
Michael T. Schaub <sup>3</sup>,
Rebekka K. Schneider<sup>4,5,6</sup>,
Ivan G. Costa<sup>1*</sup>

<small>

<sup>1</sup>Institute for Computational Genomics, Faculty of Medicine, RWTH Aachen University, Aachen, 52074 Germany

<sup>2</sup>Department of Hematology and Stem Cell Transplantation, University Hospital Essen, Germany

<sup>3</sup>Department of Computer Science, RWTH Aachen University, Germany

<sup>4</sup>Department of Cell Biology, Institute for Biomedical Engineering, Faculty of Medicine, RWTH Aachen University, Pauwelsstrasse 30, 52074 Aachen, NRW, Germany

<sup>5</sup>Oncode Institute, Erasmus Medical Center, Rotterdam, 3015GD, the Netherlands

<sup>6</sup>Department of Hematology, Erasmus Medical Center, Rotterdam, 3015GD, the Netherlands

</small>

![Alt text](image.png)


**Motivation:** Ligand-receptor (LR) analysis allows the characterization of cellular crosstalk from single cell RNA-seq data. However, current LR methods provide limited approaches for prioritization of cell types, ligands or receptors or characterizing changes in crosstalk between two biological conditions.

**Results:** pyCrossTalkeR is a framework for network analysis and visualisation of LR networks. pyCrossTalkeR identifies relevant ligands, receptors and cell types contributing to changes in cell communication when contrasting two biological states: disease vs. homeostasis. A case study on scRNA-seq of human myeloproliferative neoplasms reinforces the strengths of pyCrossTalkeR for characterisation of changes in cellular crosstalk in disease state.

## Install

You can install pyCrossTalkeR with the simple commands below:

```
pip install pycrosstallker
```



***Note:** Please avoid to use the following characters in celltype name: '$'*

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


## Features v2.1.0


- Single and Comparative Reports
   - Cell Cell Interaction visualization
   - Sending and Receiving Cells Ranking
   - CCI and GCI PCA ranking
      - All measures and PC table
      - PC1 and PC2 based barplot
   - LR pair visualization plot can be done
   - Store analysis results directly in AnnData and export to `.h5ad`


## Citation

If you use `pyCrossTalkeR` in your research, please cite our paper:

> **CrossTalkeR: Analysis and Visualisation of Ligand Receptor Networks**  
> James S Nagai, Nils B Leimkühler, Michael T Schaub, Rebekka K Schneider, Ivan G Costa.  
> *Bioinformatics*, Volume 37, Issue 22, 2021, Pages 4263–4265.  
> [https://doi.org/10.1093/bioinformatics/btab370](https://doi.org/10.1093/bioinformatics/btab370)

```bibtex
@article{nagai_crosstalker_2021,
  title = {{CrossTalkeR}: {Analysis} and {Visualisation} of {Ligand} {Receptor} {Networks}},
  author = {Nagai, James S and Leimkühler, Nils B and Schaub, Michael T and Schneider, Rebekka K and Costa, Ivan G},
  journal = {Bioinformatics},
  volume = {37},
  number = {22},
  pages = {4263--4265},
  year = {2021},
  doi = {10.1093/bioinformatics/btab370},
  url = {https://doi.org/10.1093/bioinformatics/btab370}
}

```

## Tutorials

* [pyCrossTalkeR Example - Human Myelofibrosis](notebooks/Human_Myelofibrosis)
* [pyCrossTalkeR Example - Human Myocardial Infarction](notebooks/Human_Myocardial_Infarction)
* [pyCrossTalkeR Example - LIANA+ Integration](notebooks/Integration_Liana_pyCrossTalkeR)


```{toctree}
:maxdepth: 2
:hidden:

installation
tutorials
```

## API Reference
```{toctree}
---
maxdepth: 1
---
api
```

```{toctree}
:maxdepth: 2
:hidden:

references
contributors
```