# Analysis of repeated plot measurements

<p align="center">
  <img src="img/compressed_focal_stack.jpg" alt="Focal stack" width="45%" />
  <img src="img/Schematic_figure_sampling_data.jpg" alt="Sampling procedure" width="45%" />
</p>
This repository contains code to reproduce the statistical analysis, numeric results, figures, and tables reported in [todo: add link to published version].

## Author

> Jonas Anderegg  
> Plant Pathology Group  
> ETH Zürich  

## Installation
Download the R-code as an archive or using git. Executing the scripts will first check for missing R-libraries and their dependences and install them, if needed. 

## Data
Data is available via the ETH Zürich publications and research data repository:


## Reproducing the Analysis
To reproduce the analysis, execute the scripts in the order indicated by their numerical prefixes. Many scripts depend on outputs generated in earlier steps, so it is important to follow the sequence.

### Overview
1. `01_growth_analysis.R` includes all data pre-processing and filtering, and modelling of univariate relationships between features and lesion growth. Imports functions from `utils/analyze.R`.
2. `02_rfe.R` implements the described feature selection strategy. Imports functions from `utils/feature_selection.R`
3. `03_random_reg.R` implements the random regression approach. 
4. `04_random_reg.R` implements the reported 2-stage approach for heritability estimation.  
5. `05_reference_STB.R` implements the described analysis of reference data.  
