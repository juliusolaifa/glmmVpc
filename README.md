# VPC: Variance Partitioning in Generalized Linear Mixed Models

## Package Overview

This R package is designed to estimate and infer the Variance Partition Coefficient (VPC) from Generalized Linear Mixed Models (GLMMs). The package handles multiple GLMM families (e.g., Negative Binomial, Tweedie, Gaussian, comPoisson) and provides utilities for data generation, model fitting, and parallel computation.

## File Structure

The package is organized into the following key files and functions:

```plaintext
R/
├── design_matrices.R        # Functions related to design matrices (random and fixed effects)
    ├── generateRandomInterceptMatrix    # Generates design matrix for random intercepts
    └── generateRandomDesignMatrices     # Generates random design matrices based on inputs
├── fixed_random_effects.R   # Functions for handling random effects in GLMMs
    ├── computeFixedEffects               # Computes fixed effects for the model
    └── computeRandomEffects              # Computes random effects for the model
├── link_functions.R         # Functions for GLMM link functions (identity, logit, etc.)
    └── computeMuFromEta                  # Computes expected mean from linear predictor
├── glm_glmm_utils.R         # General utilities for GLMM/GLM fitting and computation
    ├── computeMuGLMM                      # Computes the expected mean for GLMM models
    └── compute_mu_sig                     # Computes both mean and variance for the linear predictor
├── vpc_calculation.R        # Functions related to VPC (Variance Partition Coefficient) estimation
    ├── calculate_vpc_for_family           # Main function for VPC calculation for various GLMM families
    └── vpc                               # Wrapper function to calculate VPC from fitted models
├── glmm_fitting.R           # Functions for fitting GLMM models (e.g., singleGLMMFit)
    ├── extractParametersByFamily          # Extracts family-specific parameters from a fitted GLMM
    └── singleGLMMFit                      # Fits a GLMM model and extracts key model parameters
├── data_generation.R        # Functions for generating data based on GLMM families
    ├── generateDataByFamily               # Generates data according to a specified GLMM family
    ├── singleGLMMData                     # Generates a single data set for GLMM
    └── batchGLMMData                      # Generates multiple GLMM data sets
├── parallel_glmm.R          # Parallelization-related functions for large-scale GLMMs
    └── parallelbatchGLMMData              # Runs batch GLMM data generation in parallel
├── utils.R                  # General utility functions (e.g., log moments, checks)
    ├── mvn                               # Multivariate normal distribution generator
    ├── groups                            # Groups structure helper
    ├── rgen01                            # Random number generator for generating 0-1 values
    ├── vec2mat                           # Converts vector into a matrix
    ├── mat2vec                           # Converts matrix into a vector
    └── logNormMoment                     # Computes moments for log-normal distribution
└── family_utils.R           # Functions related to handling GLMM family parameters (e.g., glmmTMBfamily)
    └── get_glmmTMBfamily                 # Converts a family string into a corresponding GLMM family object
```
