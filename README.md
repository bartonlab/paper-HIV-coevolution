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
# Inference of Selection Coefficients and Fitness Values

This document describes the computational pipeline used to infer selection coefficients and compute fitness values from sequence data.
> The core inference code was adapted from the prior work by Sohail et al. (2021) [MPL-inference](https://github.com/bartonlab/paper-MPL-inference).

## Overview of Steps

### 0. Sequence Preprocessing

* Process raw sequences and generate multiple sequence alignments using **HIVAlign**.
* Ensure that each sequence is associated with a collection time point.
* Data preprocessing and additional details are described in the prior study by Sohail et al. (2021) [MPL-inference](https://github.com/bartonlab/paper-MPL-inference).

---

### 1. Model of Evolution

We aim to infer the **selection coefficients** as the **maximum a posteriori** estimate from the trajectory of allele frequencies over time.

#### Notation

* Let $L$ be the number of sites, and $q$ be the number of allele states per site (e.g., $q = 5$ for DNA, $q = 21$ for protein sequences).
* We consider a genetic population consisting of many individuals, each carrying a genotype represented by a sequence:

$$
\mathbf{g} = (g_1, g_2, \ldots, g_L) = (g_i)_{i=1}^L
$$

* This population evolves under **selection**, **mutation**, and **recombination**.
  
* The selection coefficient vector is defined as:

$$
s = (s_{1,1}, \ldots, s_{1,q}, s_{2,1}, \ldots, s_{L,q}) = \left( s_{i,a} \right)_{i=1,\ldots,L;\; a=1,\ldots,q}
$$



* We model **allele frequency dynamics**, defined as the population average at each time point $t_k$:

$$
x_{i,a}(t_k) = \langle \delta_{g_i, a} \rangle_{t_k}
$$


We assume allele frequencies evolve over time according to a distribution conditioned on selection coefficients $s$. The joint distribution is modeled as a Markov chain:

$$
P(x(t_1),\ldots,x(t_{K+1}) | s; x(t_0) )P(x(t_0)) = \prod_{k=1}^{K} P( x(t_{k+1}) | x(t_k); s)P(x(t_0))
$$

This population evolves under **selection**, **mutation**, and **recombination**.


### 2. Inference of Selection Coefficients

#### Bayesian Inference of Selection

We place a prior over selection coefficients, denoted $P^{\text{prior}}(s)$. Using Bayes’ rule, we infer the posterior distribution:

$$
P( s \mid x(t_1), \ldots, x(t_{K+1}) ) \propto P( ( x(t_k) )_{k=0}^{K+1} \mid s ) \cdot P^{\text{prior}}(s)
$$


We estimate the most probable selection coefficients (MAP estimate):

$$
\hat{s} = \arg\max_{s}  \sum_{k=0}^{K} \log P( x(t_{k+1}) \mid x(t_k); s ) + \log P^{\text{prior}}(s)  
$$


#### Core Equation

Under the diffusion limit, the selection coefficients are obtained using the formula:

$$
\hat{s} = (C + \gamma I)^{-1} (\Delta x - \Delta u)
$$

Where:

* $\Delta x$ is the **net change in allele frequencies** over the observed time course:

$$
\Delta x = \sum_{k=0}^{K} ( x(t_{k+1}) - x(t_k) ) = x(t_{K+1}) - x(t_0)
$$

  with time points $t_0, t_1, \ldots, t_{K+1}$.

* $\Delta u$ is the **expected net change in allele frequencies due to mutation**.

* $C$ is the **integrated covariance matrix** over time:

$$
C = \sum_{k=0}^{K} \Delta t_k ~ C(t_k), \quad \Delta t_k = t_{k+1} - t_k
$$

  Each covariance matrix $C(t_k)$ is defined element-wise as:

  $$
  C_{ij,ab} =
  \begin{cases}
  x_{ij,ab} - x_{i,a} x_{j,b} & \text{if } i \neq j \\
  x_{i,a}(1 - x_{i,a}) & \text{if } i = j
  \end{cases}
  $$

  Where:

  * $x_{i,a}$ is the marginal frequency of allele $a$ at site $i$
  * $x_{ij,ab}$ is the joint frequency of allele $a$ at site $i$ and allele $b$ at site $j$

#### Covariance Matrix Estimation

Accurate estimation of $C$ is challenging due to limited sample size and the large number of site pairs. To improve robustness, we apply **linear interpolation** of covariance between timepoints:

$$
\int_0^1 ( x_{ij}^{\tau} - x_i^{\tau} x_j^{\tau}   ) d\tau
$$

Where:

* $x_{ij}^{\tau} = \tau x_{ij}(t_{k+1}) + (1 - \tau) x_{ij}(t_k)$
* $x_i^{\tau} = \tau x_i(t_{k+1}) + (1 - \tau) x_i(t_k)$

This interpolation is done for each interval $[t_k, t_{k+1}]$, and the results are summed over time to obtain the final integrated covariance matrix.

#### Regularization Term

* $\gamma$ is a regularization parameter derived from the prior assumption on the selection coefficients.
* We assume a **Gaussian prior** on $\boldsymbol{s}$:

$$
p(\boldsymbol{s}) \propto \exp\left(-\frac{\gamma}{2} \| \boldsymbol{s} \|^2\right), \quad \text{with } \| \boldsymbol{s} \|^2 = \sum_{i,a} s_{i,a}^2
$$

* $d$ is the dimension of $\boldsymbol{s}$, i.e., $d = L \times q$.

---

### 3. Computation of Fitness Values

Using the inferred selection coefficients $\boldsymbol{s}$, we compute the **fitness** of any genotype $\boldsymbol{g} = (g_1, \ldots, g_L)$, where $g_i = 1, \ldots, q$, as:

$$
F(\boldsymbol{g}) = 1 + \sum_{i=1}^{L} s_i(g_i)
$$

> Note: The absolute offset (e.g., "1") in fitness values is arbitrary, as only **relative fitness** is meaningful in our context.

---

# License

This repository is dual licensed as [GPL-3.0](LICENSE-GPL) (source code) and [CC0 1.0](LICENSE-CC0) (figures, documentation, and our presentation of the data).

