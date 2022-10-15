# CMOT: Cross Modality Optimal Transport for multimodal inference

## Abstract
<p align="justify">
Biological mechanisms are complex spanning multiple facets, each providing a unique view of the underlying mechanism. Recent single-cell technologies such as scRNA-seq and scATAC-seq have facilitated parallel probing into each of these facets, providing multimodal measurements of individual cells. An integrative analysis of these single-cell modalities can therefore provide a comprehensive understanding of specific cellular and molecular mechanisms. Despite the technological advances, simultaneous profiling of multiple modalities of single cells continues to be challenging as opposed to single modality measurements, and therefore, it may not always be feasible to conduct such profiling experiments (e.g., missing modalities). Furthermore, even with the available multimodalities, data integration remains elusive since modalities may not always have paired samples, leaving partial to no correspondence information. 
To address those challenges, we developed Cross-Modality Optimal Transport (CMOT), a computational approach to infer missing modalities of single cells based on optimal transport (OT). CMOT first aligns a group of cells (source) within available multi-modal data onto a common latent space. Then, it applies optimal transport to map the cells with a single modality (target) to the aligned cells in source from the same modality by minimizing their cost of transportation using Wasserstein distance. This distance can be regularized by prior knowledge (e.g., cell types) or induced cells clusters to improve mapping, and an entropy of transport to speed up OT computations. Once transported, CMOT uses K-Nearest-Neighbors to infer the missing modality for the cells in target from another modality of mapped cells in source. Furthermore, the alignment of CMOT works for partially corresponding information from multi-modalities, i.e., allowing a fraction of unmatched cells in source. We evaluated CMOT on several emerging multi-modal datasets, e.g., gene and protein expression and chromatin accessibility in developing brain, cancers, and immunology. We found that not only does CMOT outperform existing state-of-art methods, but its inferred gene expression is biologically interpretable such as classifying cell types and cancer types. 
</p>

## Flow chart
![alt text](https://github.com/sayali7/CMOT/blob/main/src/Fig1.png)

## System requirements
* Python 3.6 or above
* Pandas 1.3.5
* Numpy 1.21.4
* Scikit-Learn 1.0.2
* [Python Optimal Transport (POT)](https://pythonot.github.io/) [1]

## Download code
The code has been test on Python versions 3.6 and above on Linux.
```r
git clone https://github.com/daifengwanglab/CMOT
cd CMOT/src
```
## Usage
We demonstrate CMOT on Pan-Cancer data [2] for gene expression inference from chromatin accessiblity.

### Download data
Pan-cancer data can be downloaded [here](https://github.com/sayali7/CMOT/tree/main/src/data/Pan-Cancer).

1. Train CMOT:
```r
python3 run_cmot.py --sourceX ./src/data/Pan-Cancer/PanCancerX.csv --sourceY ./src/data/Pan-Cancer/PanCancerY.csv --targetY ./src/data/Pan-Cancer/PanCancerY_hat.csv --K 5 --d 10 --W ./src/data/Pan-Cancer/W.npy --hc 3 --reg_e 5e03 --reg_cl 1e00 --topFeat 150 --k 40
```
The command line arguments are:
* sourceX: .csv file of size (s<sub>X</sub>, f<sub>X</sub>) for source modality X
* sourceY: .csv file of size (s<sub>Y</sub>, f<sub>Y</sub>) for source modality Y
* labels: .csv file of size (s<sub>Y</sub>, f<sub>Y</sub>) specifying prior knowledge of cell labels in source modality Y
* targetY: .csv file of size (s<sub>$\widehat{Y}$</sub>, f<sub>$\widehat{Y}$</sub>) target modality $\widehat{Y}$
* K: integer specifying the nearest neighbors for Non-linear manidold alignment (NMA) (Step A)
* d: integer specifying the latent dimension for NMA (Step A)
* W: binary numpy matrix (.npy) of size (s<sub>X</sub>,s<sub>Y</sub>) specifying the correspondence between cells of X and Y (Step A)
* hc: integer specifying hierarchical clusters if no prior knowledge is available
* reg_e: entropy regularization for optimal transport (Step B)
* reg_cl: label regularization for optimal transport (Step B)
* topFeat: integer specifying number of top variable features to use for K-nearest neighbors (Step C)
* k: integer specifying the K-nearest neighbors for cross-modality inference (Step C)
* outdir: output directory

s<sub>X</sub>: number of cells in X, f<sub>X</sub>: number of features in X, s<sub>Y</sub>: number of cells in Y, f<sub>Y</sub>: number of cells in Y, s<sub>$\widehat{Y}$</sub>: number of cells in widehat{Y}, f<sub>$\widehat{Y}$</sub>: number of features in $\widehat{Y}$
 
Output:
The result is the normalized inferred modality $\widehat{X}$ of size (s<sub>$\widehat{X}$</sub>, f<sub>$\widehat{X}$</sub>); s<sub>$\widehat{X}$</sub>: number of cells in inferred modality $\widehat{X}$, f<sub>$\widehat{X}$</sub>: number of features in modality $\widehat{X}$, such that f<sub>X</sub> = f<sub>$\widehat{X}$</sub>.

2. Evaluate CMOT's inferred modality:
```r
python3 cmot_evaluations.py --targetX_hat ./src/data/Pan-Cancer/PanCancerX_hat.csv --predX_hat .src/results/Norm_ModalityXhat.csv
```
* targetX_hat: .csv file of size (s<sub>$\widehat{X}$</sub>, f<sub>$\widehat{X}$</sub>) held-out target modality $\widehat{X}$
* predX_hat: .csv file of inferred target modality $\widehat{X}$ by CMOT

s<sub>$\widehat{X}$</sub>: number of cells in $\widehat{X}$, f<sub>$\widehat{X}$</sub>: number of features in $\widehat{X}$

Output:
The above evaluations give cell-wise Pearson correlations for the inferred modality $\widehat{X}$ versus held out target modality $\widehat{X}$


## References
<a id="1">[1]</a>
Rémi Flamary, Nicolas Courty, Alexandre Gramfort, Mokhtar Z. Alaya, Aurélie Boisbunon, Stanislas Chambon, Laetitia Chapel, Adrien Corenflos, Kilian Fatras, Nemo Fournier, Léo Gautheron, Nathalie T.H. Gayraud, Hicham Janati, Alain Rakotomamonjy, Ievgen Redko, Antoine Rolet, Antony Schutz, Vivien Seguy, Danica J. Sutherland, Romain Tavenard, Alexander Tong, Titouan Vayer, “POT Python Optimal Transport library,” JMLR, [Online]. Available: Website: https://pythonot.github.io/

<a id="1">[2]</a> 
L. Liu et al., “Deconvolution of single-cell multi-omics layers reveals regulatory heterogeneity,” Nat Commun, vol. 10, no. 1, p. 470, Dec. 2019, doi: 10.1038/s41467-018-08205-7.
