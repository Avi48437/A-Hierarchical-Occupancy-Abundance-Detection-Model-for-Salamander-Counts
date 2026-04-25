# Hierarchical Occupancy–Abundance–Detection Model for Salamander Counts

This repository contains code for modeling salamander count data using zero-inflated models and a hierarchical occupancy–abundance–detection framework.

## Repository Structure

Codes/
  - data_prep.R              # Data loading and preprocessing
  - exploratory_plots.R      # Exploratory analysis
  - zinb_comparison.R        # ZIP vs ZINB model comparison
  - hierarchical_model.R     # Proposed hierarchical model
  - dharma_diagnostics.R     # Residual diagnostics

## Workflow

Run the scripts in the following order:

1. data_prep.R  
2. exploratory_plots.R  
3. zinb_comparison.R  
4. hierarchical_model.R  
5. dharma_diagnostics.R  

## Models

Zero-inflated models:
- Poisson and Negative Binomial variants
- Account for excess zeros and overdispersion

Hierarchical model:
- Occupancy component models presence/absence
- Abundance component models latent counts
- Detection component accounts for imperfect detection

## Requirements

Required R packages:

glmmTMB  
DHARMa  
ggplot2  
dplyr  

Install using:

install.packages(c("glmmTMB", "DHARMa", "ggplot2", "dplyr"))

## Notes

- Data and generated outputs are not included in the repository
- The focus is on model specification, fitting, and comparison
