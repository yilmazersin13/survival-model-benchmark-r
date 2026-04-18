# survival-model-benchmark-r

A modular benchmark framework in R for comparing classical, penalized, Bayesian, tree-based, and deep survival models under different censoring regimes.

## Overview

This repository contains a research-oriented benchmark pipeline for survival analysis. The codebase is organized as modular R scripts and is designed to evaluate multiple survival models under fold-aware preprocessing, repeated cross-validation, and controlled censoring augmentation settings.

The framework includes:

- standardized data loaders for benchmark datasets
- fold-aware preprocessing without data leakage
- censoring augmentation for controlled censoring regimes
- repeated outer cross-validation and inner resampling scaffolding
- model fitters for classical, penalized, Bayesian, tree-based, and deep survival models
- predictive metrics and diagnostic summaries
- aggregation, statistical comparison, and figure generation

The repository is intended for reproducible benchmarking and methodological comparison rather than as a packaged software library.

## Benchmark scope

The codebase benchmarks several survival modeling strategies, including:

- Cox proportional hazards model
- lasso-penalized Cox model
- Bayesian elastic net Cox model
- random survival forest
- DeepSurv
- Cox-Time

The deep survival components are accessed from R through `reticulate` and require a Python environment with the relevant packages installed.

## Repository structure

```text
survival-model-benchmark-r/
├── README.md
├── .gitignore
├── LICENSE

├── data/
│   ├── README.md
│   └── .gitkeep
├── results/
│   └── .gitkeep
├── R/
│   ├── 01_data_loaders.R
│   ├── 02_preprocessing.R
│   ├── 03_censoring_augmentation.R
│   ├── 04_cv_scaffolding.R
│   ├── 05_model_fitters.R
│   ├── 05b_pycox_wrappers.R
│   ├── 06_evaluation_metrics.R
│   ├── 07_diagnostics.R
│   ├── 08_aggregation_and_stats.R
│   ├── 09_generate_figures.R
│   └── downstream.R
└── scripts/
    └── main.R
