# Overview

This repository contains the data and scripts used to reproduce the results presented in the accompanying manuscript.

The manuscript is available at [eLife](https://elifesciences.org/reviewed-preprints/105466) and on [bioRxiv](https://www.biorxiv.org/content/10.1101/2024.07.12.603090v3).

---
# Parallel HIV-1 fitness landscapes shape viral dynamics in humans and macaques that develop broadly neutralizing antibodies 
Kai S. Shimagaki<sup>1,2</sup>, Rebecca M. Lynch<sup>3</sup>, and John P. Barton<sup>1,2,+</sup>
<sup>1</sup> Department of Physics and Astronomy, University of California, Riverside  
<sup>2</sup> Department of Computational and Systems Biology, University of Pittsburgh School of Medicine  
<sup>3</sup> Department of Microbiology, Immunology and Tropical Medicine, School of Medicine and Health Sciences, George Washington University  
<sup>+</sup> correspondence to [jpbarton@pitt.edu](mailto:jpbarton@pitt.edu)  

# Contents

* **`figures.ipynb`**: Reproduces all the figures from the manuscript and saves them in the `figures/` directory.
* **`Inference_selection.ipynb`**: Generates bash scripts to execute inference simulations.
* **`Make_summary_CSV.ipynb`**: Processes inferred selection coefficients and creates annotated CSV files.
* **`note/` directory**: Contains notebooks that generate intermediate CSV files used for figure generation.
* **`data/` directory**: Contains all necessary input data for analysis.

Data preprocessing and the main analyses (including CSV generation) are implemented in **Julia**, compatible with versions **1.1 to 1.8 and higher**. 
Larger datasets, including sequence files, are available on [Zenodo](https://zenodo.org/records/15685461).

---

# License

This repository is dual licensed as [GPL-3.0](LICENSE-GPL) (source code) and [CC0 1.0](LICENSE-CC0) (figures, documentation, and our presentation of the data).

