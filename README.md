# Overview

This repository contains and data and scripts for reproducing the results accompanying the manuscript  

### Parallel HIV-1 evolutionary dynamics in humans and rhesus macaques who develop broadly neutralizing antibodies
Kai S. Shimagaki<sup>1,2</sup>, Rebecca M. Lynch<sup>3</sup>, and John P. Barton<sup>1,2,#</sup>

<sup>1</sup> Department of Physics and Astronomy, University of California, Riverside  
<sup>2</sup> Department of Computational and Systems Biology, University of Pittsburgh School of Medicine  
<sup>3</sup> Department of Microbiology, Immunology and Tropical Medicine, School of Medicine and Health Sciences, George Washington University  
<sup>#</sup> correspondence to [jpbarton@pitt.edu](mailto:jpbarton@pitt.edu)  

The preprint is available at [bioRxiv](https://doi.org/10.1101/2024.07.12.603090).


# Contents

All the figures used in the manuscript can be reproduced using the `figures.ipynb` notebook, which generates figures 
in the `figures/` directory. The required CSV files are generated by notebooks in the `note/` directory, and necessary 
files are located in the `data` directories. Data preprocessing and the main analysis, including CSV file generation, 
are implemented in Julia, compatible with versions 1.1 through 1.8 and higher. The `figures.ipynb` notebook is written in Python 3. 
Larger data files, including sequences, are located on Zenodo. 


# License

This repository is dual licensed as [GPL-3.0](LICENSE-GPL) (source code) and [CC0 1.0](LICENSE-CC0) (figures, documentation, and our presentation of the data).
