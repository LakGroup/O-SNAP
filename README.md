# O-SNAP

**O-SNAP** is a data analysis pipeline to automate the generation of morphological features derived from SMLM images of nuclear targets. It is conceptually related to **ECLiPSE**, another analysis pipeline developed from our group. It features:
- Automated generation of 144 morphological features
- Volcano analyis to measure fold-change in pair-wise manner between O-SNAP featurs
- Feature set enrichment analysis for changes in general trends between phenotype pairs
- Create suite of classifiers to compare phenotypes based on O-SNAP features

The manuscript describing this work, "O-SNAP: A comprehensive pipeline for spatial profiling of chromatin architecture" by Kim et al. is available on bioRxiv at DOI: https://doi.org/10.1101/2025.07.18.665612.

## System Requirements
O-SNAP is implemented in MATLAB R2024b and is currently only available for Windows 64-bit systems. The code was implemented and tested on a Intel® Xeon® Silver 4214R CPU, 2.40GHz (48 GB RAM) desktop tower.

The application additionally requires the Bioinformatics, Deep Learning, Paralllel Computing, Signal Processing, and Statistics and Machine Learning Toolboxes in MATLAB. For speedups, multiple cores and a Mex setup is required (see tutorial).

To perform pseudotimeline analysis, R (https://www.r-project.org/) and the dynverse package (https://github.com/dynverse/dyno) are required. The respective links contain installation instructions.

## Example data
Full example data used for the tutorial is available on figshare, DOI:  https://doi.org/10.6084/m9.figshare.29533940. The STORM data is of H3K27me3 localizations originally collected from Martinez-Sarmiento et al, Cell Reports (2024) https://doi.org/10.1016/j.celrep.2024.114170.

Running the dataset on the example data set "Example_H3K27me3_ALL" on a Intel® Xeon® Silver 4214R CPU, 2.40GHz (48 GB RAM) system using all 12 processes took 45 hours and 26 minutes.
